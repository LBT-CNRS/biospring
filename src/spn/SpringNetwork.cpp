#include "SpringNetwork.h"

#include <fstream>
#include "logging.h"
#include "measure.hpp"

#include "forcefield/constants.hpp"
#include "forcefield/shared/torsion_shared.h"
#include "forcefield/ForceField.h"
#include "forcefield/energy/hydrogenbond.hpp"
#include "forcefield/energy/stacking.hpp"
#include "forcefield/ForceFieldElectrostaticCoulombAndStericLennardJones_12_6Amber.h"
#include "forcefield/ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Lewitt.h"
#include "forcefield/ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Zacharias.h"
#include "forcefield/ForceFieldElectrostaticCoulombAndStericLinear.h"

#include "IO/CSVSampleWriter.h"
#include "IO/OpenDXReader.h"
#include "IO/PDBTrajectoryWriter.h"
#include "IO/XTCTrajWriter.h"
#include "interactor/Interactor.h"
#ifdef MDDRIVER_SUPPORT
	#include "interactor/mddriver/InteractorMDDriver.h"
#endif
#ifdef FREESASA_SUPPORT
    #include "interactor/freesasa/InteractorFreeSASA.h"
#endif

#include <algorithm>
#include <sstream>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <limits>
#include <math.h>
#include <memory>
#include <stdlib.h>
#include <string.h>
#include <utility>
#include <chrono>
#include <thread>

#ifdef OPENMP_SUPPORT
#include <omp.h>
#endif

#include "rigidbody/RigidBodiesManager.h"

namespace biospring
{
namespace spn
{

unsigned SpringNetwork::_currentstructid = 0;

SpringNetwork::~SpringNetwork() {}

void SpringNetwork::_resetEnergies()
{
    _energies.reset();
}

void SpringNetwork::_updateInsertionVector()
{
    {
        _insertionVector->computeVector();
        _insertionVector->computeAngle();
        _insertionVector->computeRollAngle();
    }
}

// Calculates spring forces and applies them to the particles.
// Updates global `_energies.spring` variable.
void SpringNetwork::computeSpringForces()
{
    float springenergy = 0.0f;
    _springForceScratch.resize(_dynamicsprings.size());

#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < _dynamicsprings.size(); ++i)
    {
        Spring & spring = getSpring(_dynamicsprings[i]);
        _springForceScratch[i] = spring.computeForce(*_ff);
    }

    // Particle forces and energies are accumulated in a deterministic serial
    // pass. The expensive spring evaluation remains parallel, while concurrent
    // writes and thread-count-dependent floating-point reductions are avoided.
    for (size_t i = 0; i < _dynamicsprings.size(); ++i)
    {
        Spring & spring = getSpring(_dynamicsprings[i]);
        const Vector3f & force = _springForceScratch[i];
        spring.getParticle1().addForce(force);
        spring.getParticle2().addForce(-force);
        springenergy += spring.getEnergy();
    }

    _energies.spring = springenergy;
}

// The DIHEDRAL family switches gate the tabulated torsions, which is all this
// term is now: there is no spring collection per family any more, only the
// tables the families index into.
void SpringNetwork::computeDihedralForces()
{
    _energies.dihedral = computeTorsionForces();
}

int SpringNetwork::_anyHydrogenBondPartner(size_t index) const
{
    if (index + 1 < _hbDonorOffset.size())
        for (size_t s = _hbDonorOffset[index]; s < _hbDonorOffset[index + 1]; ++s)
            if (_hbDonorSlot[s] >= 0)
                return _hbDonorSlot[s];
    if (index + 1 < _hbAcceptorOffset.size())
        for (size_t s = _hbAcceptorOffset[index]; s < _hbAcceptorOffset[index + 1]; ++s)
            if (_hbAcceptorSlot[s] >= 0)
                return _hbAcceptorSlot[s];
    return -1;
}

size_t SpringNetwork::getHydrogenBondCount() const
{
    size_t n = 0;
    for (int slot : _hbDonorSlot)
        if (slot >= 0)
            ++n;
    return n;
}

std::array<size_t, 4> SpringNetwork::getHydrogenBondPairCensus() const
{
    std::map<std::pair<long, long>, size_t> per_pair;
    for (size_t i = 0; i + 1 < _hbDonorOffset.size(); ++i)
        for (size_t s = _hbDonorOffset[i]; s < _hbDonorOffset[i + 1]; ++s)
        {
            const int j = _hbDonorSlot[s];
            if (j < 0)
                continue;
            const Particle & p = _particles[i];
            const Particle & q = _particles[static_cast<size_t>(j)];
            const long dres = static_cast<long>(p.getResId()) - static_cast<long>(q.getResId());
            if (p.getChainName() == q.getChainName() && dres >= -1 && dres <= 1)
                continue;
            // Chain folded into the key by its first character: enough to
            // keep the two strands of a duplex apart.
            const long a = static_cast<long>(p.getResId()) * 256 +
                           (p.getChainName().empty() ? 0 : p.getChainName()[0]);
            const long b = static_cast<long>(q.getResId()) * 256 +
                           (q.getChainName().empty() ? 0 : q.getChainName()[0]);
            per_pair[{std::min(a, b), std::max(a, b)}]++;
        }
    std::array<size_t, 4> census = {0, 0, 0, 0};
    for (const auto & entry : per_pair)
        census[std::min<size_t>(entry.second, 4) - 1]++;
    return census;
}

// Every candidate the assignment could consider, and what became of it. The
// point is to answer "why was this pair not formed?" with the code's own
// arithmetic rather than a transcription of it.
//
// Columns: the two ends, the distance, the angular weight of each side, the
// planarity factor, the resulting strength in kJ/mol, how many slots each end
// has free, and whether the pair is currently held.
void SpringNetwork::dumpHydrogenBondCandidates(const std::string & path, int step) const
{
    std::ofstream out(path, std::ios::app);
    if (!out)
    {
        logging::warning("Cannot open hydrogen bond candidate log '%s'.", path.c_str());
        return;
    }
    const int exponent = getHydrogenBondAngularExponent();
    const float cutoff = getHydrogenBondCutoff();
    auto free_slots = [&](size_t i, bool donor) {
        const auto & off = donor ? _hbDonorOffset : _hbAcceptorOffset;
        const auto & slot = donor ? _hbDonorSlot : _hbAcceptorSlot;
        if (i + 1 >= off.size())
            return 0;
        int n = 0;
        for (size_t s = off[i]; s < off[i + 1]; ++s)
            if (slot[s] < 0)
                n++;
        return n;
    };
    out << "# step " << step << "\n";
    out << "# chainA residA resA atomA chainB residB resB atomB distance wA wB planarity "
           "strength freeDonorA freeAccA freeDonorB freeAccB held\n";
    const std::vector<size_t> candidates = _donorAcceptorParticleIndexes();
    for (const size_t i : candidates)
    {
        const Particle & p = getParticle(i);
        _nsearch.hbond->for_each_neighbor(i, [&](size_t j) {
            if (j <= i || isProbeParticle(j))
                return;
            const Particle & q = getParticle(j);
            if (p.getResId() == q.getResId() && p.getChainName() == q.getChainName())
                return;
            const float d = Particle::distance(p, q);
            if (d >= cutoff || d < 1e-6f)
                return;
            // Roles must be complementary in at least one direction, which is
            // what the assignment requires too.
            if (!((p.isDonor() && q.isAcceptor()) || (q.isDonor() && p.isAcceptor())))
                return;
            const Vector3f vhat = (q.getPosition() - p.getPosition()) / d;
            const Vector3f hd = donorDirection(p, vhat);
            const Vector3f ha = donorDirection(q, -vhat);
            const float wd = hd.norm() > 1e-6f
                                 ? forcefield::hydrogen_bond_angular_factor(hd.dot(vhat), exponent)
                                 : 1.0f;
            const float wa = ha.norm() > 1e-6f
                                 ? forcefield::hydrogen_bond_angular_factor(ha.dot(-vhat), exponent)
                                 : 1.0f;
            const float pl = hydrogenBondPlanarity(i, j).factor;
            const HydrogenBondParameters mp = hydrogenBondParameters(i, j);
            const float strength =
                _ff->computeHydrogenBondEnergy(d, mp.welldepth, mp.equilibrium, mp.width) * wd * wa * pl;
            out << p.getChainName() << ' ' << p.getResId() << ' ' << p.getResName() << ' ' << p.getName()
                << ' ' << q.getChainName() << ' ' << q.getResId() << ' ' << q.getResName() << ' '
                << q.getName() << ' ' << d << ' ' << wd << ' ' << wa << ' ' << pl << ' ' << strength << ' '
                << free_slots(i, true) << ' ' << free_slots(i, false) << ' ' << free_slots(j, true) << ' '
                << free_slots(j, false) << ' ' << (areHydrogenBonded(i, j) ? 1 : 0) << '\n';
        });
    }
}

void SpringNetwork::dumpHydrogenBonds(const std::string & path, int step) const
{
    std::ofstream out(path, std::ios::app);
    if (!out)
    {
        logging::warning("Cannot open hydrogen bond log '%s'.", path.c_str());
        return;
    }
    const int exponent = getHydrogenBondAngularExponent();
    out << "# step " << step << "\n";
    for (size_t i = 0; i + 1 < _hbDonorOffset.size(); ++i)
        for (size_t s = _hbDonorOffset[i]; s < _hbDonorOffset[i + 1]; ++s)
        {
            const int j = _hbDonorSlot[s];
            if (j < 0)
                continue;
            const Particle & d = _particles[i];
            const Particle & a = _particles[static_cast<size_t>(j)];
            const Vector3f v = a.getPosition() - d.getPosition();
            const float distance = v.norm();
            // BOTH angular factors, as the force applies them. Reporting only
            // the donor's would say the acceptor term does nothing, which is
            // how the previous one-sided version read.
            float w = 1.0f;
            if (distance > 1e-6f)
            {
                const Vector3f vhat = v / distance;
                const Vector3f hd = donorDirection(d, vhat);
                if (hd.norm() > 1e-6f)
                    w *= forcefield::hydrogen_bond_angular_factor(hd.dot(vhat), exponent);
                const Vector3f ha = donorDirection(a, -vhat);
                if (ha.norm() > 1e-6f)
                    w *= forcefield::hydrogen_bond_angular_factor(ha.dot(-vhat), exponent);
            }
            // The planarity factor folded into the reported weight, so the log
            // says what the force applies, and appended as its own last column
            // so it can be read on its own. Appended rather than inserted: the
            // first eleven columns keep their meaning and their position.
            const float planarity = hydrogenBondPlanarity(i, static_cast<size_t>(j)).factor;
            w *= planarity;
            out << d.getChainName() << ' ' << d.getResId() << ' ' << d.getResName() << ' ' << d.getName() << ' '
                << a.getChainName() << ' ' << a.getResId() << ' ' << a.getResName() << ' ' << a.getName() << ' '
                << distance << ' ' << w << ' ' << _ff->computeHydrogenBondEnergy(distance) * w << ' '
                << planarity << '\n';
        }
}

SpringNetwork::HydrogenBondParameters SpringNetwork::hydrogenBondParameters(size_t i, size_t j) const
{
    const int gi = _particles[i].hydrogenBondGroup();
    if (gi > 0 && gi == _particles[j].hydrogenBondGroup() && static_cast<size_t>(gi) <= _hbgroups.size())
        return _hbgroups[static_cast<size_t>(gi) - 1];
    return {_ff->getHydrogenBondWellDepth(), _ff->getHydrogenBondEquilibrium(), _ff->getHydrogenBondWidth()};
}

bool SpringNetwork::areHydrogenBonded(size_t a, size_t b) const
{
    if (a + 1 < _hbDonorOffset.size())
        for (size_t s = _hbDonorOffset[a]; s < _hbDonorOffset[a + 1]; ++s)
            if (_hbDonorSlot[s] == static_cast<int>(b))
                return true;
    if (b + 1 < _hbDonorOffset.size())
        for (size_t s = _hbDonorOffset[b]; s < _hbDonorOffset[b + 1]; ++s)
            if (_hbDonorSlot[s] == static_cast<int>(a))
                return true;
    return false;
}

// HiRE-RNA v3's planarity factor (arXiv:1404.0568, equation 5) and its
// gradient. The factor is
//
//     Pl = Fi Fj,   Fi = (1/3) sum_k g(s_k),   g(s) = e^(-(s/delta)^2)
//
// where s_k is the distance of the k-th plane atom of residue j out of the
// PLANE of residue i, and Fj the same with the roles swapped. It multiplies the
// well, so a pair of bases that meet edge to edge in a shared plane keeps its
// full depth and a skewed contact loses most of it.
//
// THE GRADIENT. Each F depends on six atoms: the three whose out-of-plane
// distances it measures, and the three that define the plane it measures them
// from -- through the centroid AND through the normal. Writing e1 = p1 - p0,
// e2 = p2 - p0, m = e1 x e2, n = m/|m|, o = (p0+p1+p2)/3 and
//
//     Si = sum_k g'(s_k),   Vi = (1/3) sum_k g'(s_k) (q_k - o_i),
//     Gi = (Vi - n_i (n_i.Vi)) / |m_i|
//
// the three plane atoms of i each take -(Si/9) n_i from the centroid, and from
// the normal p1 takes e2 x Gi, p2 takes Gi x e1 and p0 minus their sum -- so
// the normal's three forces cancel exactly, as they must. The partner's atom
// q_k takes (1/3) g'(s_k) n_i. Checked against the numerical gradient.
SpringNetwork::HydrogenBondPlanarity SpringNetwork::hydrogenBondPlanarity(size_t a, size_t b) const
{
    HydrogenBondPlanarity out;
    if (!isHydrogenBondPlanarityEnabled())
        return out;
    const int fa = getBaseFrameOfParticle(a), fb = getBaseFrameOfParticle(b);
    // No plane at one end is the ordinary case, not an error: it is every amino
    // acid, every sugar and every phosphate. Same residue would mean one base
    // paired with itself.
    if (fa < 0 || fb < 0 || fa == fb)
        return out;

    const float delta = static_cast<float>(_config.hbond.planaritywidth);
    const BaseFrame & Fa = _baseFrames[static_cast<size_t>(fa)];
    const BaseFrame & Fb = _baseFrames[static_cast<size_t>(fb)];

    Vector3f pa[3], pb[3];
    for (unsigned k = 0; k < 3; ++k)
    {
        pa[k] = getParticle(Fa.atoms[k]).getPosition();
        pb[k] = getParticle(Fb.atoms[k]).getPosition();
        out.atoms[k] = Fa.atoms[k];
        out.atoms[3 + k] = Fb.atoms[k];
    }

    const Vector3f e1a = pa[1] - pa[0], e2a = pa[2] - pa[0];
    const Vector3f e1b = pb[1] - pb[0], e2b = pb[2] - pb[0];
    const Vector3f ma = e1a ^ e2a, mb = e1b ^ e2b;
    const float lma = ma.norm(), lmb = mb.norm();
    if (lma < BASE_PLANE_MIN_AREA || lmb < BASE_PLANE_MIN_AREA)
    {
        out.degenerate = true;
        return out;
    }
    const Vector3f na = ma / lma, nb = mb / lmb;
    const Vector3f oa = (pa[0] + pa[1] + pa[2]) / 3.0f;
    const Vector3f ob = (pb[0] + pb[1] + pb[2]) / 3.0f;

    // s[k]: plane atom k of b, out of a's plane. t[k]: the other way round.
    float s[3], t[3], gs[3], gt[3], dgs[3], dgt[3];
    float Fi = 0.0f, Fj = 0.0f;
    for (unsigned k = 0; k < 3; ++k)
    {
        s[k] = (pb[k] - oa).dot(na);
        t[k] = (pa[k] - ob).dot(nb);
        gs[k] = forcefield::hydrogen_bond_planarity_term(s[k], delta);
        gt[k] = forcefield::hydrogen_bond_planarity_term(t[k], delta);
        dgs[k] = forcefield::hydrogen_bond_planarity_derivative(s[k], delta);
        dgt[k] = forcefield::hydrogen_bond_planarity_derivative(t[k], delta);
        Fi += gs[k] / 3.0f;
        Fj += gt[k] / 3.0f;
    }
    out.active = true;
    out.factor = Fi * Fj;

    float Si = 0.0f, Sj = 0.0f;
    Vector3f Vi, Vj;
    for (unsigned k = 0; k < 3; ++k)
    {
        Si += dgs[k];
        Sj += dgt[k];
        Vi = Vi + (pb[k] - oa) * (dgs[k] / 3.0f);
        Vj = Vj + (pa[k] - ob) * (dgt[k] / 3.0f);
    }
    const Vector3f Gi = (Vi - na * na.dot(Vi)) / lma;
    const Vector3f Gj = (Vj - nb * nb.dot(Vj)) / lmb;

    // dFi/d(atoms of a): a's own plane, through its centroid and its normal.
    Vector3f dFi_a[3], dFj_b[3];
    const Vector3f ai1 = e2a ^ Gi, ai2 = Gi ^ e1a;
    const Vector3f bj1 = e2b ^ Gj, bj2 = Gj ^ e1b;
    dFi_a[0] = na * (-Si / 9.0f) - ai1 - ai2;
    dFi_a[1] = na * (-Si / 9.0f) + ai1;
    dFi_a[2] = na * (-Si / 9.0f) + ai2;
    dFj_b[0] = nb * (-Sj / 9.0f) - bj1 - bj2;
    dFj_b[1] = nb * (-Sj / 9.0f) + bj1;
    dFj_b[2] = nb * (-Sj / 9.0f) + bj2;

    // d(Fi Fj)/dx, with each F's dependence on the OTHER residue's atoms being
    // the plain (1/3) g'(.) n of the out-of-plane distance.
    for (unsigned k = 0; k < 3; ++k)
    {
        out.gradient[k] = dFi_a[k] * Fj + nb * (dgt[k] / 3.0f) * Fi;
        out.gradient[3 + k] = dFj_b[k] * Fi + na * (dgs[k] / 3.0f) * Fj;
    }
    return out;
}

// Applies the Morse potential to every currently active hydrogen bond --
// see _assignHydrogenBondPairs, called separately beforehand to (re)build
// the slots. Walks the DONOR slots, so every bond is seen exactly once and
// with its direction known, which the angular weight needs.
//
// E = M(d) * w(cos t), M the Morse well in the distance between the two
// heavy atoms and w = max(cos t, 0)^2 the angular weight, t the angle at the
// donor between (donor - antecedent) and (acceptor - donor). Three bodies
// carry force, not two: the antecedent enters through cos t alone. With
// u = donor - antecedent and v = acceptor - donor,
//
//     dE/d(acceptor)   = M'(d) w vhat + M(d) w'(c) (uhat - c vhat)/|v|
//     dE/d(antecedent) = -M(d) w'(c) (vhat - c uhat)/|u|
//     dE/d(donor)      = -(the other two)      [translation invariance]
//
// A donor with no resolved antecedent keeps w = 1 and the old two-body
// behaviour exactly, so a .hbond table without the fifth column is unchanged
// by this. A particle can now appear in several bonds at once, so the
// apply pass stays serial -- it is a genuine write conflict now, not just
// summation order.
// 62 degrees is the DEFAULT lobe angle, and it is where a planar sp2 site's two
// hydrogens or its two lone pairs sit relative to its single heavy bond. Measured
// over 56 amine hydrogens of a B-DNA duplex the angle is 62.0 degrees on guanine N2,
// 62.1 on adenine N6 and 62.9 on cytosine N4, and over the carbonyls the acceptor
// side comes out at 52.7 to 57.6 from the axis. One value serves both: reconstructing
// the lobes at 62 degrees puts every one of the five Watson-Crick bond types within
// 8.2 to 15.3 degrees of its real partner, where naming the real hydrogen gives 8.5
// to 14.0. The reconstruction is as good as knowing where the hydrogen is.
//
// The angle is per site now, because the OTHER two-lobe geometry in this chemistry
// wants 105: the Burgi-Dunitz approach of a nucleophile to a carbonyl carbon, which
// leaves the sp2 plane instead of lying in it. The .hbond table carries both.
static const float HBOND_LOBE_DEFAULT_DEGREES = 62.0f;

SpringNetwork::HydrogenBondSite SpringNetwork::hydrogenBondSite(const Particle & self, int a1, int a2, int lobes,
                                                               float lcos, float lsin,
                                                               const Vector3f & towards) const
{
    HydrogenBondSite site;
    if (a1 < 0)
        return site;
    const Vector3f u1 = self.getPosition() - _particles[static_cast<size_t>(a1)].getPosition();
    site.l1 = u1.norm();
    if (site.l1 <= 1e-6f)
        return site;
    site.e1 = u1 / site.l1;

    if (lobes != Particle::HBOND_LOBES_NONE && a2 >= 0)
    {
        // The plane atom, taken from the SAME atom the bisector form would have
        // used, but read as a plane reference: r runs antecedent1 -> antecedent2, so
        // the three atoms span the group's plane. IN PLANE takes the Gram-Schmidt of
        // r against the axis; OUT OF PLANE takes the plane's normal. Either way shat
        // is unit and orthogonal to e1.
        site.r = _particles[static_cast<size_t>(a2)].getPosition() - _particles[static_cast<size_t>(a1)].getPosition();
        const Vector3f q = lobes == Particle::HBOND_LOBES_OUT_OF_PLANE
                               ? (site.e1 ^ site.r)
                               : site.r - site.e1 * site.e1.dot(site.r);
        site.m = q.norm();
        if (site.m <= 1e-6f)
        {
            // The three atoms are collinear, so no plane exists. Fall back to
            // the axis, which is what a single antecedent would have given.
            site.hhat = site.e1;
            site.hlen = 1.0f;
            site.valid = true;
            return site;
        }
        site.shat = q / site.m;
        site.lobes = lobes;
        site.lcos = lcos;
        site.lsin = lsin;
        // Both lobes, and the bond takes the one it fits. Not the slot index:
        // the slot assignment does not know about lobes, so binding them to
        // slots could seat an acceptor against the lobe pointing away from it.
        // The maximum is continuous, and the two are equal only on the axis,
        // where a lobe pair is at its weakest and the force smallest.
        const Vector3f plus = site.e1 * lcos + site.shat * lsin;
        const Vector3f minus = site.e1 * lcos - site.shat * lsin;
        site.sigma = plus.dot(towards) >= minus.dot(towards) ? 1.0f : -1.0f;
        site.hhat = site.e1 * lcos + site.shat * (site.sigma * lsin);
        site.hlen = 1.0f; // exactly, since e1 and shat are orthonormal
        site.valid = true;
        return site;
    }

    Vector3f h = site.e1;
    if (a2 >= 0)
    {
        const Vector3f u2 = self.getPosition() - _particles[static_cast<size_t>(a2)].getPosition();
        site.l2 = u2.norm();
        if (site.l2 > 1e-6f)
        {
            site.e2 = u2 / site.l2;
            h = h + site.e2;
        }
    }
    site.hlen = h.norm();
    if (site.hlen <= 1e-6f)
        return site;
    site.hhat = h / site.hlen;
    site.valid = true;
    return site;
}

// The unit vector along which a donor's hydrogen points, or a zero vector
// when the donor has no antecedent at all. `towards` picks the lobe on a site
// that has two; it is ignored otherwise.
//
// It exists because the direction was computed in three separate places (the
// force, the candidate ranking, the log) and fixing one of them left the
// other two reporting the old answer, which read as "the fix does nothing".
Vector3f SpringNetwork::donorDirection(const Particle & p, const Vector3f & towards) const
{
    const HydrogenBondSite s = hydrogenBondSite(p, p.antecedentIndex(), p.antecedentIndex2(), p.lobeMode(),
                                                p.lobeCos(), p.lobeSin(), towards);
    return s.valid ? s.hhat : Vector3f();
}

void SpringNetwork::computeHydrogenBondForces()
{
    struct Bond
    {
        size_t donor, acceptor;
        int antecedent, antecedent2;
        // The acceptor's own antecedents. A carbonyl oxygen's lone pairs are
        // no less directional than the donor's hydrogen: measured on
        // ubiquitin's native helix bonds, the angle between C->O and
        // O->N(donor) is 25 degrees median, never past 30, so cos^2 about
        // the C=O direction describes it as well as it describes the donor.
        int acceptorAntecedent, acceptorAntecedent2;
        // Whether each side's second antecedent is a plane reference rather
        // than a second bond: see ParticleProperties::hasLobes.
        int donorLobes, acceptorLobes;
    };
    const int exponent = getHydrogenBondAngularExponent();
    std::vector<Bond> bonds;
    bonds.reserve(_hbDonorSlot.size());
    for (size_t i = 0; i + 1 < _hbDonorOffset.size(); ++i)
        for (size_t s = _hbDonorOffset[i]; s < _hbDonorOffset[i + 1]; ++s)
            if (_hbDonorSlot[s] >= 0)
                bonds.push_back({i, static_cast<size_t>(_hbDonorSlot[s]), getParticle(i).antecedentIndex(),
                                 getParticle(i).antecedentIndex2(),
                                 getParticle(static_cast<size_t>(_hbDonorSlot[s])).antecedentIndex(),
                                 getParticle(static_cast<size_t>(_hbDonorSlot[s])).antecedentIndex2(),
                                 getParticle(i).lobeMode(),
                                 getParticle(static_cast<size_t>(_hbDonorSlot[s])).lobeMode()});

    // Twelve slots per bond: two antecedents on each side, the donor, the
    // acceptor -- the angular weight is a product of two factors, one per side,
    // so its gradient reaches every one of them -- and then the three base
    // plane atoms of each end's residue, which the planarity factor reaches.
    // The last six are unused, and left at zero, whenever that factor does not
    // apply, which is every hydrogen bond outside a base pair.
    _hydrogenBondForceScratch.resize(bonds.size() * 12);
    std::vector<float> pairEnergyScratch(bonds.size());
    std::vector<HydrogenBondPlanarity> planarityScratch(bonds.size());

#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
    for (size_t k = 0; k < bonds.size(); ++k)
    {
        const Particle & pd = getParticle(bonds[k].donor);
        const Particle & pa = getParticle(bonds[k].acceptor);
        Vector3f v = pa.getPosition() - pd.getPosition();
        const float distance = v.norm();
        if (distance < 1e-6f)
        {
            for (unsigned slot = 0; slot < 12; ++slot)
                _hydrogenBondForceScratch[12 * k + slot] = Vector3f();
            pairEnergyScratch[k] = 0.0f;
            continue;
        }
        const Vector3f vhat = v / distance;

        // Per PAIR, not per force field: a couple whose two ends name the same
        // group carries that group's well. See hydrogenBondParameters.
        const HydrogenBondParameters mp = hydrogenBondParameters(bonds[k].donor, bonds[k].acceptor);
        const float morse = _ff->computeHydrogenBondEnergy(distance, mp.welldepth, mp.equilibrium, mp.width);
        const float dmorse =
            _ff->computeHydrogenBondForceModule(distance, mp.welldepth, mp.equilibrium, mp.width);

        // A hydrogen bond is directional at BOTH ends. The donor's hydrogen
        // points somewhere, and so do the acceptor's lone pairs: measured on
        // ubiquitin's native helix bonds, the angle between C->O and the
        // direction of the donor is 25 degrees median and never past 30, so
        // the same cos^2 form describes it. Weighting only the donor side
        // overestimates every bond by the factor the acceptor term would have
        // applied -- about 0.82 on ideal geometry.
        //
        // Each side's direction comes from hydrogenBondSite, so the force, the
        // candidate ranking and the log cannot disagree about it.
        const HydrogenBondSite sd = hydrogenBondSite(pd, bonds[k].antecedent, bonds[k].antecedent2,
                                                     bonds[k].donorLobes, pd.lobeCos(), pd.lobeSin(), vhat);
        const HydrogenBondSite sa = hydrogenBondSite(pa, bonds[k].acceptorAntecedent,
                                                     bonds[k].acceptorAntecedent2, bonds[k].acceptorLobes,
                                                     pa.lobeCos(), pa.lobeSin(), -vhat);
        const bool hasD = sd.valid, hasA = sa.valid;
        const Vector3f dhat = sd.hhat, ahat = sa.hhat;

        const float cd = hasD ? dhat.dot(vhat) : 1.0f;
        const float ca = hasA ? ahat.dot(-vhat) : 1.0f;
        const float wd = hasD ? forcefield::hydrogen_bond_angular_factor(cd, exponent) : 1.0f;
        const float wa = hasA ? forcefield::hydrogen_bond_angular_factor(ca, exponent) : 1.0f;
        const float dwd = hasD ? forcefield::hydrogen_bond_angular_derivative(cd, exponent) : 0.0f;
        const float dwa = hasA ? forcefield::hydrogen_bond_angular_derivative(ca, exponent) : 0.0f;
        const float w = wd * wa;

        // Three sub-terms, each balanced on its OWN atoms rather than letting
        // one global "donor takes the rest" absorb everything. The acceptor
        // term acts on the DONOR through the same v that the donor term acts
        // on the acceptor through, and folding both into one balance puts that
        // reaction on the wrong atom with the wrong sign -- which is precisely
        // what the five-atom finite-difference test caught.
        //
        // computeHydrogenBondForceModule's own convention is "positive =
        // attractive, on the FIRST particle" -- here the donor, since v runs
        // donor -> acceptor.
        const float conv = static_cast<float>(forcefield::GLOBAL_SPRING_FORCE_CONVERT);
        Vector3f f_donor, f_acceptor, f_d1, f_d2, f_a1, f_a2;

        // The forces one side's angular factor puts on the ATOMS ITS DIRECTION IS
        // BUILT FROM -- the two antecedents. `towards` is the unit vector from
        // this side to its partner (vhat on the donor, -vhat on the acceptor),
        // `c` is hhat.towards and `k` the scalar the chain rule already carries.
        // The self and the partner are balanced by the caller.
        //
        // BISECTOR FORM. hhat = (e1 + e2)/|e1 + e2|, and with
        // t = towards - hhat c the force on antecedent i is the part of t
        // perpendicular to ei, scaled by k/(|e1+e2| li). e1.t is zero when there
        // is only one antecedent, which is why the same expression serves both.
        //
        // LOBE FORM. hhat = cos a e1 + sigma sin a shat, with shat either the
        // Gram-Schmidt of r = a2 - a1 against e1 (in plane) or e1 x r (out of
        // plane). For the in-plane form, writing
        //     ta = (towards - e1 (e1.towards)) / l1
        //     ts = (towards - shat (shat.towards)) / m
        //     z  = -[ (e1.r) ts + (ts.e1) r ]
        //     P(x) = x - e1 (e1.x)
        //     A = cos62 ta + sigma sin62 P(z)/l1
        //     B = sigma sin62 P(ts)
        // the force is k(A - B) on a1 and +kB on a2. With sin62 set to zero this
        // collapses to A = ta and B = 0, which is the single-antecedent bisector
        // case, so the two forms agree in that limit. Both the signs and the limit
        // are checked by finite differences, on two geometries of opposite sigma --
        // the algebra here went through a Gram-Schmidt projection and the sign of
        // the B term came out inverted on the first attempt.
        auto angular = [&](const HydrogenBondSite & st, const Vector3f & towards, float c, float k,
                           Vector3f & f1, Vector3f & f2) {
            f1 = Vector3f();
            f2 = Vector3f();
            if (!st.lobes)
            {
                const Vector3f t = towards - st.hhat * c;
                const float g = k / st.hlen;
                f1 = (t - st.e1 * st.e1.dot(t)) * (g / st.l1);
                if (st.l2 > 1e-6f)
                    f2 = (t - st.e2 * st.e2.dot(t)) * (g / st.l2);
                return;
            }
            const Vector3f ta = (towards - st.e1 * st.e1.dot(towards)) / st.l1;
            const Vector3f ts = (towards - st.shat * st.shat.dot(towards)) / st.m;
            const auto perp = [&](const Vector3f & x) { return x - st.e1 * st.e1.dot(x); };
            const float ss = st.sigma * st.lsin;
            Vector3f A, B;
            if (st.lobes == Particle::HBOND_LOBES_OUT_OF_PLANE)
            {
                // shat = normalise(e1 x r). d(e1 x r) = de1 x r + e1 x dr, and
                // ts.(de1 x r) = de1.(r x ts) while ts.(e1 x dr) = dr.(ts x e1).
                A = ta * st.lcos + perp(st.r ^ ts) * (ss / st.l1);
                B = (ts ^ st.e1) * ss;
            }
            else
            {
                // shat = normalise(r - e1 (e1.r)), the Gram-Schmidt.
                const Vector3f z = -(ts * st.e1.dot(st.r) + st.r * ts.dot(st.e1));
                A = ta * st.lcos + perp(z) * (ss / st.l1);
                B = perp(ts) * ss;
            }
            f1 = (A + B) * k;
            f2 = B * (-k);
        };

        // (A) radial, on the pair.
        const Vector3f radial = vhat * (-dmorse * w);
        f_acceptor += radial;
        f_donor -= radial;

        // (B) the donor's angular factor: explicit on its antecedents, its
        // v-dependence on the acceptor, the donor balancing the three.
        if (hasD && dwd != 0.0f)
        {
            angular(sd, vhat, cd, morse * wa * dwd * conv, f_d1, f_d2);
            const Vector3f on_acceptor = (dhat - vhat * cd) * (-morse * wa * dwd * conv / distance);
            f_acceptor += on_acceptor;
            f_donor -= (f_d1 + f_d2 + on_acceptor);
        }

        // (C) the acceptor's angular factor, the same expressions with the
        // roles swapped: its partner direction is -vhat and its v-dependence
        // lands on the DONOR, the acceptor balancing the three.
        if (hasA && dwa != 0.0f)
        {
            angular(sa, -vhat, ca, morse * wd * dwa * conv, f_a1, f_a2);
            const Vector3f on_donor = (ahat + vhat * ca) * (-morse * wd * dwa * conv / distance);
            f_donor += on_donor;
            f_acceptor -= (f_a1 + f_a2 + on_donor);
        }

        // (D) the planarity factor, which MULTIPLIES everything above: the six
        // forces already computed are scaled by it, and its own gradient lands
        // on the six base plane atoms. E = M(d) w Pl, so d(M w Pl) = Pl d(M w)
        // + M w dPl, which is exactly these two lines of bookkeeping.
        const HydrogenBondPlanarity pl = hydrogenBondPlanarity(bonds[k].donor, bonds[k].acceptor);
        if (pl.active)
        {
            f_d1 = f_d1 * pl.factor;
            f_d2 = f_d2 * pl.factor;
            f_a1 = f_a1 * pl.factor;
            f_a2 = f_a2 * pl.factor;
            f_donor = f_donor * pl.factor;
            f_acceptor = f_acceptor * pl.factor;
            for (unsigned slot = 0; slot < 6; ++slot)
                _hydrogenBondForceScratch[12 * k + 6 + slot] = pl.gradient[slot] * (-morse * w * conv);
        }
        else
            for (unsigned slot = 0; slot < 6; ++slot)
                _hydrogenBondForceScratch[12 * k + 6 + slot] = Vector3f();

        _hydrogenBondForceScratch[12 * k + 0] = f_d1;
        _hydrogenBondForceScratch[12 * k + 1] = f_d2;
        _hydrogenBondForceScratch[12 * k + 2] = f_a1;
        _hydrogenBondForceScratch[12 * k + 3] = f_a2;
        _hydrogenBondForceScratch[12 * k + 4] = f_donor;
        _hydrogenBondForceScratch[12 * k + 5] = f_acceptor;
        pairEnergyScratch[k] = morse * w * pl.factor;
        planarityScratch[k] = pl;
    }

    float hbond_energy = 0.0f;
    size_t degenerate_planes = 0;
    for (size_t k = 0; k < bonds.size(); ++k)
    {
        if (bonds[k].antecedent >= 0)
            getParticle(static_cast<size_t>(bonds[k].antecedent)).addForce(_hydrogenBondForceScratch[12 * k + 0]);
        if (bonds[k].antecedent2 >= 0)
            getParticle(static_cast<size_t>(bonds[k].antecedent2)).addForce(_hydrogenBondForceScratch[12 * k + 1]);
        if (bonds[k].acceptorAntecedent >= 0)
            getParticle(static_cast<size_t>(bonds[k].acceptorAntecedent)).addForce(_hydrogenBondForceScratch[12 * k + 2]);
        if (bonds[k].acceptorAntecedent2 >= 0)
            getParticle(static_cast<size_t>(bonds[k].acceptorAntecedent2)).addForce(_hydrogenBondForceScratch[12 * k + 3]);
        getParticle(bonds[k].donor).addForce(_hydrogenBondForceScratch[12 * k + 4]);
        getParticle(bonds[k].acceptor).addForce(_hydrogenBondForceScratch[12 * k + 5]);
        // A plane atom can also be the acceptor or an antecedent -- guanine's
        // N3 is both -- so these are separate additions, not a replacement.
        if (planarityScratch[k].active)
            for (unsigned slot = 0; slot < 6; ++slot)
                getParticle(planarityScratch[k].atoms[slot])
                    .addForce(_hydrogenBondForceScratch[12 * k + 6 + slot]);
        if (planarityScratch[k].degenerate)
            degenerate_planes++;
        hbond_energy += pairEnergyScratch[k];
    }

    if (degenerate_planes)
        logging::warning("Hydrogen bond planarity: %zu bond(s) took a factor of 1 because a base plane frame "
                         "has collapsed below %.2f A2 of area.",
                         degenerate_planes, static_cast<double>(BASE_PLANE_MIN_AREA));
    _energies.hbond = hbond_energy;}

// Base stacking, after HiRE-RNA v3 (Cragnolini, Laurin, Derreumaux & Pasquali,
// arXiv:1404.0568), their equation (4). Per pair of base planes,
//
//     E = -eps G(d) P A_i A_j
//
// with d the distance between the two frame CENTRES, G the Gaussian,
// P = (n_i.n_j)^2 the parallel factor and A = u^2(2-u^2), u = n.rhat, the
// face-to-face factor of each base. The forms and the measured parameters are
// in forcefield/shared/stacking_shared.h; what is here is the geometry and the
// gradient.
//
// SIX BODIES, NOT TWO. A frame's normal and its centre are both built from all
// three of its atoms, so every one of the six carries force. Writing
// b = p1 - p0, c = p2 - p0, m = b x c and n = m/|m|, the centre is the plain
// mean of the three and
//
//     dE/dm       = (N - n (n.N)) / |m|,     N = dE/dn
//     dE/dp1      = c x (dE/dm)
//     dE/dp2      = (dE/dm) x b
//     dE/dp0      = -(the other two)
//
// the last line making each frame's three normal-forces sum to zero exactly, by
// construction rather than by balancing. The centre's share is dE/dv split in
// three equal parts, with the sign of which frame it belongs to. Translation
// invariance therefore holds to the bit, and the finite-difference test checks
// the rest.
void SpringNetwork::computeStackingForces()
{
    _energies.stacking = 0.0f;
    if (_baseFrames.empty() || !_nsearch.stacking)
        return;

    const float cutoff = getStackingCutoff();
    // `scale` multiplies the well depth, which is what it means everywhere
    // else: one knob to weigh the whole family without touching its shape.
    const float eps = static_cast<float>(_config.stacking.epsilon * _config.stacking.scale);
    const float r_st = static_cast<float>(_config.stacking.equilibrium);
    const float sigma = static_cast<float>(_config.stacking.width);
    const float conv = static_cast<float>(forcefield::GLOBAL_SPRING_FORCE_CONVERT);

    // The pair list first, so the force loop below is a flat parallel loop over
    // pairs rather than a nested walk with a shared accumulator.
    struct Pair
    {
        size_t i, j;
    };
    std::vector<Pair> pairs;
    for (size_t fi = 0; fi < _baseFrames.size(); ++fi)
        _nsearch.stacking->for_each_neighbor(_baseFrameAnchors[fi], [&](size_t candidate) {
            const int fj = _frameOfAnchor[candidate];
            // Each pair once: the grid returns both directions.
            if (fj > static_cast<int>(fi))
                pairs.push_back({fi, static_cast<size_t>(fj)});
        });

    _stackingForceScratch.assign(pairs.size() * 6, Vector3f());
    std::vector<float> pairEnergy(pairs.size(), 0.0f);
    std::vector<unsigned char> degenerateScratch(pairs.size(), 0);

#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
    for (size_t k = 0; k < pairs.size(); ++k)
    {
        const BaseFrame & fi = _baseFrames[pairs[k].i];
        const BaseFrame & fj = _baseFrames[pairs[k].j];

        Vector3f pi[3], pj[3];
        for (unsigned a = 0; a < 3; ++a)
        {
            pi[a] = getParticle(fi.atoms[a]).getPosition();
            pj[a] = getParticle(fj.atoms[a]).getPosition();
        }

        const Vector3f bi = pi[1] - pi[0], ci = pi[2] - pi[0];
        const Vector3f bj = pj[1] - pj[0], cj = pj[2] - pj[0];
        const Vector3f mi = bi ^ ci, mj = bj ^ cj;
        const float lmi = mi.norm(), lmj = mj.norm();
        // A frame whose three atoms have nearly collapsed onto a line has no
        // usable normal, and the gradient below divides by this quantity -- so
        // the threshold is the PHYSICAL one, not an epsilon. See
        // BASE_PLANE_MIN_AREA; the tables pick triples whose smallest angle has
        // a sine above 0.7, so this does not fire on healthy geometry.
        if (lmi < BASE_PLANE_MIN_AREA || lmj < BASE_PLANE_MIN_AREA)
        {
            degenerateScratch[k] = 1;
            continue;
        }
        const Vector3f ni = mi / lmi, nj = mj / lmj;

        const Vector3f centre_i = (pi[0] + pi[1] + pi[2]) / 3.0f;
        const Vector3f centre_j = (pj[0] + pj[1] + pj[2]) / 3.0f;
        const Vector3f v = centre_j - centre_i;
        const float d = v.norm();
        if (d < 1e-6f || d >= cutoff)
            continue;
        const Vector3f rhat = v / d;

        const float w = ni.dot(nj);
        const float P = w * w;
        const float ui = ni.dot(rhat), uj = nj.dot(rhat);
        const float Ai = forcefield::stacking_alignment_factor(ui);
        const float Aj = forcefield::stacking_alignment_factor(uj);
        const float dAi = forcefield::stacking_alignment_derivative(ui);
        const float dAj = forcefield::stacking_alignment_derivative(uj);
        const float G = forcefield::stacking_distance_factor(d, r_st, sigma);
        const float dG = forcefield::stacking_distance_derivative(d, r_st, sigma);

        pairEnergy[k] = -eps * G * P * Ai * Aj;

        // dE/dv: the distance through G, and rhat through both alignments.
        const Vector3f dEdv = (rhat * (dG * P * Ai * Aj) + (ni - rhat * ui) * (G * P * dAi * Aj / d) +
                               (nj - rhat * uj) * (G * P * Ai * dAj / d)) *
                              (-eps);

        // dE/dn of each side: the parallel factor through the OTHER normal, and
        // its own alignment through rhat.
        const Vector3f Ni = (nj * (2.0f * w * Ai) + rhat * (P * dAi)) * (-eps * G * Aj);
        const Vector3f Nj = (ni * (2.0f * w * Aj) + rhat * (P * dAj)) * (-eps * G * Ai);
        const Vector3f gi = (Ni - ni * ni.dot(Ni)) / lmi;
        const Vector3f gj = (Nj - nj * nj.dot(Nj)) / lmj;

        // Forces are -dE/dx, converted from kJ.mol-1.A-1 to the integrator's
        // own unit the way every molar term in this code is.
        // v = centre_j - centre_i, so dE/dp is -dEdv/3 on each atom of i and
        // +dEdv/3 on each atom of j; the force is minus that.
        const Vector3f share_i = dEdv * (conv / 3.0f);
        const Vector3f share_j = dEdv * (-conv / 3.0f);

        const Vector3f fi1 = (ci ^ gi) * -conv;
        const Vector3f fi2 = (gi ^ bi) * -conv;
        const Vector3f fj1 = (cj ^ gj) * -conv;
        const Vector3f fj2 = (gj ^ bj) * -conv;

        _stackingForceScratch[6 * k + 0] = share_i - fi1 - fi2;
        _stackingForceScratch[6 * k + 1] = share_i + fi1;
        _stackingForceScratch[6 * k + 2] = share_i + fi2;
        _stackingForceScratch[6 * k + 3] = share_j - fj1 - fj2;
        _stackingForceScratch[6 * k + 4] = share_j + fj1;
        _stackingForceScratch[6 * k + 5] = share_j + fj2;
    }

    // Serial, like the hydrogen bond's: one particle can be in many stacking
    // pairs at once, so this is a genuine write conflict rather than a
    // question of summation order.
    float stacking_energy = 0.0f;
    for (size_t k = 0; k < pairs.size(); ++k)
    {
        for (unsigned a = 0; a < 3; ++a)
        {
            getParticle(_baseFrames[pairs[k].i].atoms[a]).addForce(_stackingForceScratch[6 * k + a]);
            getParticle(_baseFrames[pairs[k].j].atoms[a]).addForce(_stackingForceScratch[6 * k + 3 + a]);
        }
        stacking_energy += pairEnergy[k];
    }
    _energies.stacking = stacking_energy;
    const size_t degenerate = std::count(degenerateScratch.begin(), degenerateScratch.end(), 1);
    if (degenerate)
        logging::warning("Stacking: %zu base pair(s) skipped because a plane frame has collapsed below %.2f A2 "
                         "of area. The geometry has already gone wrong where that happens.",
                         degenerate, static_cast<double>(BASE_PLANE_MIN_AREA));
}

std::pair<size_t, size_t> SpringNetwork::getStackingPairCensus() const
{
    size_t within = 0, sprung = 0;
    if (_baseFrames.empty() || !_nsearch.stacking)
        return {0, 0};
    const float cutoff = getStackingCutoff();
    auto centre = [&](const BaseFrame & f) {
        return (getParticle(f.atoms[0]).getPosition() + getParticle(f.atoms[1]).getPosition() +
                getParticle(f.atoms[2]).getPosition()) /
               3.0f;
    };
    for (size_t fi = 0; fi < _baseFrames.size(); ++fi)
        _nsearch.stacking->for_each_neighbor(_baseFrameAnchors[fi], [&](size_t candidate) {
            const int fj = _frameOfAnchor[candidate];
            if (fj <= static_cast<int>(fi))
                return;
            if ((centre(_baseFrames[static_cast<size_t>(fj)]) - centre(_baseFrames[fi])).norm() >= cutoff)
                return;
            within++;
            for (unsigned a = 0; a < 3; ++a)
                for (unsigned b = 0; b < 3; ++b)
                    if (getParticle(_baseFrames[fi].atoms[a])
                            .isInSpringNeighbors(
                                static_cast<unsigned>(_baseFrames[static_cast<size_t>(fj)].atoms[b])))
                    {
                        sprung++;
                        return;
                    }
        });
    return {within, sprung};
}

void SpringNetwork::_setupStacking()
{
    if (!isStackingEnabled())
        return;
    if (getStackingCutoff() < 1e-6)
        throw std::runtime_error("Stacking cutoff must be > 0");
    if (_config.stacking.width <= 0.0)
        logging::die("stacking.width (sigma) must be > 0: the distance factor is e^(-(r-r_st)^2/sigma).");
    if (_baseFrames.empty())
    {
        logging::warning("stacking.enable is set but no base plane frame was installed: the term will do "
                         "nothing. Check stacking.path and that the topology carries nucleic residues.");
        return;
    }

    _baseFrameAnchors.clear();
    _frameOfAnchor.assign(_particles.size(), -1);
    for (size_t f = 0; f < _baseFrames.size(); ++f)
    {
        const size_t anchor = _baseFrames[f].atoms[0];
        if (_frameOfAnchor[anchor] >= 0)
            logging::die("Stacking: particle %zu is the first plane atom of two different base frames.", anchor);
        _frameOfAnchor[anchor] = static_cast<int>(f);
        _baseFrameAnchors.push_back(anchor);
    }

    // The grid is queried AT A PARTICLE -- atoms[0] -- while the term measures
    // centre to centre, so the search radius has to cover the cutoff plus the
    // offset at BOTH ends. Measured here rather than assumed: a nucleobase ring
    // is rigid, so this offset is a property of the frames the table chose and
    // cannot grow during a run. The margin is for the ring breathing under the
    // mesh.
    float widest = 0.0f;
    for (const BaseFrame & f : _baseFrames)
    {
        const Vector3f centre = (getParticle(f.atoms[0]).getPosition() + getParticle(f.atoms[1]).getPosition() +
                                 getParticle(f.atoms[2]).getPosition()) /
                                3.0f;
        widest = std::max(widest, (getParticle(f.atoms[0]).getPosition() - centre).norm());
    }
    const float margin = 0.5f;
    const float radius = getStackingCutoff() + 2.0f * (widest + margin);

    _nsearch.stacking = make_nsearch(_particles, radius, _baseFrameAnchors, getNeighborSkin(),
                                     getCellWidthFor(radius));
    _excludeProbeFromNeighborSearch(*_nsearch.stacking);

    logging::info("Stacking: %zu base frame(s), well %.3g kJ.mol-1 at %.2f A (sigma %.3g A2), cutoff %.2f A, "
                  "grid radius %.2f A (anchor offset up to %.2f A).",
                  _baseFrames.size(), _config.stacking.epsilon * _config.stacking.scale,
                  _config.stacking.equilibrium, _config.stacking.width, getStackingCutoff(), radius, widest);
}

// Calculate forces that apply on dynamic particles.
void SpringNetwork::computeParticleForces()
{
    float electrostatic_energy = 0.0f;
    float steric_energy = 0.0f;
    float imp_energy = 0.0f;
    float hydrophobic_energy = 0.0f;
    float hbond_core_repulsion_energy = 0.0f;

    _resizeNonbondedPairScratch();

#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < _dynamicparticules.size(); ++i)
    {
        Particle & p = getParticle(_dynamicparticules[i]);

        if (isAnyElectrostaticEnabled())
        {
            if (isElectrostaticCoulombEnabled() && p.isCharged() && _nsearch.electrostatic)
                p.addElectrostaticForce(_electrostaticPairScratch[i]);
            if (isElectrostaticFieldEnabled())
                p.addElectrostaticFieldForce();
        }

        if (isDensityGridEnabled())
            p.addDensityFieldForce();

        if (isStericEnabled())
            p.addStericForce(_stericPairScratch[i]);

        // Two paths on purpose. With the thermostat OFF this is the one this
        // code has always taken -- friction folded into the force, integrated
        // with everything else -- so every existing .msp keeps the trajectory
        // it had, to the bit. With the thermostat ON the friction moves to the
        // velocity update, where a bath has to act, and the two are no longer
        // the same arithmetic: the exact update damps the force it has just
        // added, this one does not. Both are legitimate splittings; only one
        // of them is what 45 validated example modes were measured against.
        if (isViscosityEnabled() && !isThermostatEnabled())
            p.applyViscosity(getViscosity());

        if (isIMPEnabled())
            p.addIMPForce();

        if (isHydrophobicityEnabled() && p.isHydrophobic() && _nsearch.hydrophobic)
            p.addHydrophobicityForce(_hydrophobicPairScratch[i]);

        if (isHydrogenBondEnabled() && (p.isDonor() || p.isAcceptor()) && _nsearch.hbond)
            p.addHydrogenBondCoreRepulsion(_hydrogenBondCoreRepulsionPairScratch[i]);

    }

    // Applies the deferred "other side" of each unique nonbonded pair
    // (Newton's third law) serially, since two threads may have deferred a
    // contribution to the same target particle. Must run before the force
    // read by rigid-body torque aggregation and setPreviousForce() below, and
    // before summing per-particle energies, since it feeds both.
    _applyNonbondedPairScratch(_stericPairScratch, steric_energy);
    _applyNonbondedPairScratch(_electrostaticPairScratch, electrostatic_energy);
    _applyNonbondedPairScratch(_hydrophobicPairScratch, hydrophobic_energy);
    _applyNonbondedPairScratch(_hydrogenBondCoreRepulsionPairScratch, hbond_core_repulsion_energy);

    // Sum per-particle energies in particle order to keep results reproducible
    // across OpenMP thread counts.
    for (const unsigned particle_id : _dynamicparticules)
    {
        const Particle & p = getParticle(particle_id);
        electrostatic_energy += p.getElectrostaticEnergy();
        steric_energy += p.getStericEnergy();
        imp_energy += p.getIMPEnergy();
        hydrophobic_energy += p.getHydrophobicityEnergy();
        hbond_core_repulsion_energy += p.getHydrogenBondCoreRepulsionEnergy();
    }

    // The probe is shared by every particle, therefore probe interactions must
    // not mutate it from an OpenMP loop. It is integrated exactly once per step.
    if (isProbeEnabled())
    {
        _probeparticule.resetForce();

        for (const unsigned particle_id : _dynamicparticules)
        {
            Particle & p = getParticle(particle_id);
            if (isProbeStericEnabled())
                steric_energy += p.addStericProbeForce(_probeparticule);
            if (isProbeElectrostaticEnabled())
                electrostatic_energy += p.addElectrostaticProbeForce(_probeparticule);
        }

        if (isProbeElectrostaticFieldEnabled())
        {
            const float previous_energy = _probeparticule.getElectrostaticEnergy();
            _probeparticule.addElectrostaticFieldForce();
            electrostatic_energy += _probeparticule.getElectrostaticEnergy() - previous_energy;
        }

        if (isViscosityEnabled())
            _probeparticule.applyViscosity(getViscosity());

        _probeparticule.IntegrateEuler(getTimeStep());
        _energies.kinetic += _probeparticule.getKineticEnergy();
        _syncProbeParticle();
    }

    // Rigid-body accumulators are shared between their particles. Aggregate
    // them serially after all per-particle forces are complete.
    for (const unsigned particle_id : _dynamicparticules)
    {
        Particle & p = getParticle(particle_id);
        if (isRigidBodyEnabled() && p.isRigid() && !isImpalaSamplingEnabled() && !isMonteCarloEnabled())
            rigidbody::RigidBody::computeParticleForceAndTorque(p);
        p.setPreviousForce();
    }

    _energies.electrostatic = electrostatic_energy;
    _energies.steric = steric_energy;
    _energies.imp = imp_energy;
    _energies.hydrophobic = hydrophobic_energy;
    // computeHydrogenBondForces() (called earlier in computeForces(), before
    // this method) already set _energies.hbond to the exclusive engaged-pair
    // energy -- add the always-on core-repulsion contribution on top rather
    // than overwrite it.
    _energies.hbond += hbond_core_repulsion_energy;
}

// Update the positions of the particles.
// Update global `_energies.kinetic` variable.
void SpringNetwork::updateParticlePositions()
{
    float kinetic_energy_particle = 0.0;
    // Read once rather than per particle: they are the same for everybody, and
    // the noise is keyed on the particle's own index so two threads never draw
    // the same number.
    const bool thermostatting = isThermostatEnabled();
    const float viscosity = isViscosityEnabled() ? getViscosity() : 0.0f;
    const float boltzmanntemperature = getBoltzmannTemperature();
    _thermostatstep++;
#ifdef OPENMP_SUPPORT
#pragma omp parallel default(shared)
#endif
    {
#ifdef OPENMP_SUPPORT
#pragma omp for reduction(+ : kinetic_energy_particle) schedule(static)
#endif
        // i stays a signed int: MSVC only supports OpenMP 2.0, which requires
        // a signed loop counter for #pragma omp parallel for.
        for (int i = 0; i < (int)_dynamicparticules.size(); i++)
        {
            Particle & p = getParticle(_dynamicparticules[static_cast<size_t>(i)]);
            if (p.isRigid())
                rigidbody::RigidBody::integrateParticleVelocity(p, i, getTimeStep());
            else if (thermostatting)
                p.IntegrateEulerLangevin(getTimeStep(), viscosity, boltzmanntemperature,
                                         _thermostatstep,
                                         _dynamicparticules[static_cast<size_t>(i)], THERMOSTAT_SEED);
            else
                p.IntegrateEuler(getTimeStep());

            // Check if position exploses, if one of float is NaN
            float x = p.getPosition().getX();
            float y = p.getPosition().getY();
            float z = p.getPosition().getZ();
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
            {
                logging::die("Found non-finite position for particle %d.", p.getId());
            }

            kinetic_energy_particle += p.getKineticEnergy();
            p.resetForce();
        } // omp for loop
    }     // omp parallel
    _energies.kinetic += kinetic_energy_particle;

    // The spatial grids must follow particle motion. Rebuild them once here,
    // after all particle positions have been integrated for the current step.
    _markNeighborSearchesDirty();
    _updateNeighborSearches();
}

/// @brief Compute particles force and, if activated, springs forces.
void SpringNetwork::computeForces()
{
    if (isSpringEnabled())
    {
        computeSpringForces();
        computeDihedralForces();
    }
    // Outside the isSpringEnabled() guard, unlike every family above: a
    // hydrogen bond is not a spring of the network. It is re-assigned each
    // step from the current geometry (see _assignHydrogenBondPairs), so it
    // has its own switch and its own lifetime.
    if (isHydrogenBondEnabled())
    {
        _assignHydrogenBondPairs();
        computeHydrogenBondForces();
    }
    // Outside isSpringEnabled() too, and with no assignment phase of its own:
    // unlike a hydrogen bond, stacking holds no slots and needs no partner to
    // be chosen. Every pair of base planes within the cutoff contributes what
    // its geometry is worth, the way the steric term does.
    if (isStackingEnabled())
        computeStackingForces();
    computeParticleForces();
}

// BAOAB (Leimkuhler & Matthews), the splitting that samples configurations
// correctly at a timestep where the naive one does not.
//
// Its shape is forced by the physics, not by taste: the closing half kick uses
// the force at the position the step ENDS at, so the force evaluation sits in
// the middle rather than at the top. It is still exactly ONE evaluation per
// step -- the force computed here serves the closing kick of this step and the
// opening kick of the next.
//
// idleRun() moves with it, which means an interactor and a written frame both
// see the position this step reports rather than the one it started from.
void SpringNetwork::computeStepBAOAB()
{
    _meanConstraintsDistances = 0.0;

    // B A O A, on the force left standing by the previous step. initRun() has
    // evaluated it for the first one.
    _baoabKickDriftBathDrift();

    idleRun();          // interactors, output, energies reset, step counter

    computeForces();    // at the position the drift just reached

    if (isConstraintEnabled())
        applyConstraints();
    if (isRigidBodyEnabled())
        rigidbody::RigidBodiesManager::SolveRigidBodiesDynamic();

    _baoabFinalKick();

    // Same place in the step as the ordinary path's, and it has to be repeated
    // here: computeStep() hands over to this function and RETURNS, so a call
    // added only there never runs under a thermostat. Which is the one case
    // that matters -- a bond forms because something is warm enough to bring
    // its two ends together, and a thermostat is how that is asked for.
    _peptidebond.update(*this, static_cast<unsigned>(_nbiter));

    if (isInsertionVectorEnabled())
        _updateInsertionVector();
}

void SpringNetwork::_baoabKickDriftBathDrift()
{
    const float viscosity = isViscosityEnabled() ? getViscosity() : 0.0f;
    const float boltzmanntemperature = getBoltzmannTemperature();
    _thermostatstep++;

#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < (int)_dynamicparticules.size(); i++)
    {
        Particle & p = getParticle(_dynamicparticules[static_cast<size_t>(i)]);
        if (p.isRigid())
            continue;   // the rigid-body solver owns these
        p.baoabKickDriftBathDrift(getTimeStep(), viscosity, boltzmanntemperature, _thermostatstep,
                                  _dynamicparticules[static_cast<size_t>(i)], THERMOSTAT_SEED);
    }

    // The grids have to follow the motion before the force is evaluated on it.
    _markNeighborSearchesDirty();
    _updateNeighborSearches();
}

void SpringNetwork::_baoabFinalKick()
{
    float kinetic_energy_particle = 0.0f;
#ifdef OPENMP_SUPPORT
#pragma omp parallel for reduction(+ : kinetic_energy_particle) schedule(static)
#endif
    for (int i = 0; i < (int)_dynamicparticules.size(); i++)
    {
        Particle & p = getParticle(_dynamicparticules[static_cast<size_t>(i)]);
        if (p.isRigid())
            continue;
        p.baoabFinalKick(getTimeStep());

        const float x = p.getPosition().getX();
        const float y = p.getPosition().getY();
        const float z = p.getPosition().getZ();
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
            logging::die("Found non-finite position for particle %d.", p.getId());

        kinetic_energy_particle += p.getKineticEnergy();
    }
    _energies.kinetic += kinetic_energy_particle;
}

void SpringNetwork::computeStep()
{
    if (isThermostatEnabled())
    {
        computeStepBAOAB();
        return;
    }

    idleRun();
    _meanConstraintsDistances = 0.0;

    computeForces();

    if (isConstraintEnabled())
        applyConstraints();

    if (isRigidBodyEnabled())
        rigidbody::RigidBodiesManager::SolveRigidBodiesDynamic();

    updateParticlePositions();

    // On the MOVED positions, and after the forces: a spring created here is
    // born at the length it already has, so it applies no force until the next
    // step evaluates it.
    _peptidebond.update(*this, static_cast<unsigned>(_nbiter));

    if (isInsertionVectorEnabled())
        _updateInsertionVector();
}

void SpringNetwork::run()
{
    initRun();

#ifdef MDDRIVER_SUPPORT
    // getInteractorInstance returns null when no interactor of that type was
    // registered, which is every caller that drives a SpringNetwork directly
    // rather than through biospring-cli -- the OpenCL parity tests, for one.
    // Dereferencing it unconditionally segfaulted all of them as soon as
    // MDDriver support was compiled in, and only then, which is why it went
    // unnoticed: the builds those tests were written against had it off.
    // The wait belongs with the announcement: it exists to let the server
    // finish opening its port, and there is no port without an interactor.
    if (interactor::InteractorMDDriver * interactormddriver =
            getInteractorInstance<interactor::InteractorMDDriver>())
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        logging::info("    MDDriver parameters:");
        logging::info("      port: %d is open for connection.", interactormddriver->getPort());
    }
#endif

    if (isRigidBodyEnabled())
    {
        // Set all structure to rigid.
        rigidbody::RigidBodiesManager::InitRigidBodies(this, getDynamicParticles());

        // If write csv and automatic sampling enabled, do not write relative 
        // to the option csvsampling.frequency (set to a big number)
        // The information will be written at the end of every protein rotation
        // around the insertion vector axis
        if (isImpalaSamplingEnabled() && _config.csvsample.enable)
            _config.csvsample.frequency = 1000000;
    }

    while (!isEnd())
    {
        computeStep();
        if (_isTimeToLogData())
        {
            _updateFrameRate();
            _displayFrameData();
        }
    }

    // Stop measuring time and calculate the elapsed time.
    _profiler["main"].stop();
    float elapsed = _profiler["main"].elapsed_seconds();
    logging::info("Total time measured: %5.2f seconds.", elapsed);
    logging::info("Average framerate: %5.2f sec-1.", getMaxIteration() / elapsed);

    endRun();
}

void SpringNetwork::initRun()
{
    // Starts measuring time.
    _profiler["main"].start();
    _profiler["samplerate"].start();

    _nbiter = 0;

    for (Interactor* interactor : getInteractors()) 
    {
        if (interactor != nullptr)
        {
            interactor->startInteractionThread();
        }
    }

    // BAOAB opens each step with a half kick on the force left standing by the
    // previous one. The first step has no previous one, so the force at the
    // starting structure is evaluated here. Without it the first half step
    // would kick on zero and the trajectory would start half a step behind.
    if (isThermostatEnabled())
    {
        _updateNeighborSearches();
        computeForces();
    }
}

void SpringNetwork::endRun()
{
    for (Interactor* interactor : getInteractors()) 
    {
        if (interactor != nullptr)
            interactor->stopInteractionThread();
    }
    for (Interactor* interactor : getInteractors())
    {
        if (interactor != nullptr)
            interactor->waitForInteractionThread();
    }
}

void SpringNetwork::idleRun()
{
    while (_pause)
        std::this_thread::sleep_for(std::chrono::milliseconds(10));

    for (Interactor* interactor : getInteractors()) 
    {
        if (interactor != nullptr)
        {
            interactor->syncSystemStateData();
        }
    }

    _writeNextStep();

    _resetEnergies();

    _nbiter++;

    if (_hasReachedEndOfRun())
        setEnd(true);
}

void SpringNetwork::_updateFrameRate()
{
    _framerate = getSampleRate() / _profiler["samplerate"].elapsed_seconds();
    _profiler["samplerate"].reset();
}

void SpringNetwork::_displayFrameData()
{
    logging::info("Step: %5d", _nbiter);
    logging::info("Framerate: %5.2f", _framerate);
    logging::info("Kinetic energy: %5.2f kJ.mol-1", _energies.kinetic);
    logging::info("Temperature: %5.1f K over %zu integrated particles", getInstantaneousTemperature(),
                  _dynamicparticules.size());
    if (isSpringEnabled())
    {
        logging::info("Spring energy: %5.2f kJ.mol-1", _energies.spring);
        logging::info("Dihedral energy: %5.2f kJ.mol-1", _energies.dihedral);
    }
    if (isAnyElectrostaticEnabled())
        logging::info("Electrostatic energy: %5.2f kJ.mol-1", _energies.electrostatic);
    if (isStericEnabled())
        logging::info("Steric energy: %5.2f kJ.mol-1", _energies.steric);
    if (isIMPEnabled())
        logging::info("IMP energy: %5.2f kJ.mol-1", _energies.imp);
    if (isHydrophobicityEnabled())
        logging::info("Hydrophobic energy: %5.2f kJ.mol-1", _energies.hydrophobic);
    if (isHydrogenBondEnabled())
    {
        const std::array<size_t, 4> census = getHydrogenBondPairCensus();
        logging::info("Hydrogen bond energy: %5.2f kJ.mol-1 over %zu bond(s); residue pairs holding "
                      "1/2/3/4+ bonds: %zu/%zu/%zu/%zu",
                      _energies.hbond, getHydrogenBondCount(), census[0], census[1], census[2], census[3]);
        if (!_config.hbond.log.empty())
            dumpHydrogenBonds(_config.hbond.log, _nbiter);
        if (!_config.hbond.candidatelog.empty())
            dumpHydrogenBondCandidates(_config.hbond.candidatelog, _nbiter);
    }
    if (isStackingEnabled())
    {
        // The sprung count beside the total: a cutoff-built elastic network
        // puts springs between stacked bases, and such a spring is a crude
        // stand-in for exactly this term, so a report that hid the overlap
        // would let the two double count unnoticed.
        const std::pair<size_t, size_t> census = getStackingPairCensus();
        logging::info("Stacking energy: %5.2f kJ.mol-1 over %zu base pair(s) within the cutoff, %zu of them "
                      "also joined by a spring",
                      _energies.stacking, census.first, census.second);
    }
    if (isPeptideBondEnabled())
    {
        // The two thresholds, beside what the structure is actually offering:
        // a run that forms nothing has to be able to say WHICH of the two it is
        // failing, and a bond count of zero says neither.
        logging::info("Peptide bonds formed: %u; closest candidate %5.2f A (forms below %5.2f), best angular "
                      "weight %5.3f (needs %5.3f)",
                      _peptidebond.getNumberOfBonds(), _peptidebond.getClosestApproach(), _config.peptidebond.distance,
                      _peptidebond.getBestWeight(), _config.peptidebond.weight);
    }
    if (isInsertionVectorEnabled())
    {
        logging::info("Insertion angle: %5.2lf °", _insertionVector->getAngle());
        logging::info("Roll angle: %5.2lf °", _insertionVector->getRollAngle());
        logging::info("Insertion depth: %5.2lf.", _insertionVector->getInsertionDepth());
    }
    if (isConstraintEnabled())
        logging::info("Constraints mean distance: %5.2f A", _meanConstraintsDistances);
#ifdef FREESASA_SUPPORT
    // Only when a FreeSASA interactor was actually registered. Printed
    // unconditionally it reported 0 on every run that has none, and on the
    // runs that do have one but are not dynamic it reports a figure computed
    // once at setup -- so say which, rather than letting a frozen number look
    // like a live measurement.
    if (getInteractorInstance<interactor::InteractorFreeSASA>())
        logging::info("Total SASA (%s): %5.2lf A2",
                      _freesasaState.isDynamic ? "updated" : "computed once at setup",
                      _freesasaState.sasaTotal);
#endif

}

// _nbiter stays a signed int: it is compared against getMaxIteration()
// (config's simulation.nbsteps, whose -1 means "infinite run"), and logged
// with %d in several places. Cast explicitly here (always >= 0 in practice,
// only incremented from 0).
void SpringNetwork::_writeNextStep() { _trajectories.write_step(static_cast<size_t>(_nbiter)); }

void SpringNetwork::writeNextStepNow() { _trajectories.write_step(); }

std::vector<Particle>::const_reference SpringNetwork::getParticleFromId(unsigned id) const
{
    for (unsigned i = 0; i < _particles.size(); i++)
    {
        if (_particles[i].getId() == static_cast<int>(id))
            return _particles[i];
    }
    throw std::out_of_range("SpringNetwork::getParticleFromId: Particle id not found.");
}

std::vector<Particle>::reference SpringNetwork::getParticleFromId(unsigned id)
{
    for (unsigned i = 0; i < _particles.size(); i++)
    {
        if (_particles[i].getId() == static_cast<int>(id))
            return _particles[i];
    }
    throw std::out_of_range("SpringNetwork::getParticleFromId: Particle id not found.");
}

void SpringNetwork::getParticlePosition(unsigned i, float position[3]) const
{
    const Particle & p = getParticle(i);
    position[0] = p.getX();
    position[1] = p.getY();
    position[2] = p.getZ();
}

void SpringNetwork::setForce(unsigned i, float force[3])
{
    Vector3f f(force[0], force[1], force[2]);
    getParticle(i).addForce(f);
}

// =====================================================================================
//
// Modification Methods.
//
// Should be used only to convert `Topology` to `SpringNetwork`.
//
// =====================================================================================

void SpringNetwork::addSpring(unsigned id1, unsigned id2, float equilibrium, float stiffness)
{
    if (id1 >= _particles.size() || id2 >= _particles.size())
        throw std::out_of_range("SpringNetwork::addSpring: particle index out of range");
    if (id1 == id2)
        throw std::invalid_argument("SpringNetwork::addSpring: a spring requires two distinct particles");

    // If the spring does not already exist.
    if (!_particles[id1].isInSpringNeighbors(id2))
    {
        // setup() may add one probe particle after springs have been created.
        // Reserve that slot before references to particles are stored in Spring.
        if (_springs.empty() && _particles.capacity() == _particles.size())
            _particles.reserve(_particles.size() + 1);

        Spring * old_storage = _springs.data();
        Particle & p1 = _particles[id1];
        Particle & p2 = _particles[id2];
        _springs.emplace_back(p1, p2, equilibrium, stiffness);
        _springs.back().setId(static_cast<unsigned>(_springs.size() - 1));

        if (old_storage != _springs.data())
            _rebuildSpringNeighbors();
        else
        {
            p1.addToSpringNeighbors(id2, &_springs.back());
            p2.addToSpringNeighbors(id1, &_springs.back());
        }

        if (p1.isStatic() && p2.isStatic())
            addStaticSpring(_springs.back().getId());
        else
            addDynamicSpring(_springs.back().getId());
    }
}

void SpringNetwork::releaseSpring(unsigned id)
{
    if (id >= _springs.size() || _releasedsprings.count(id))
        return;
    Spring & spring = _springs[id];
    spring.setStiffness(0.0f);
    const unsigned i = static_cast<unsigned>(spring.getParticle1().getId());
    const unsigned j = static_cast<unsigned>(spring.getParticle2().getId());
    // Only when no OTHER spring holds the same pair: two rigid-body groups
    // legitimately overlap on a pair (the _PHI/_PSI overlap is the documented
    // case), and dropping the exclusion because one of them broke would let a
    // pair that is still bonded start repelling itself.
    bool held = false;
    for (unsigned k = 0; k < _springs.size(); ++k)
    {
        if (k == id || _releasedsprings.count(k))
            continue;
        const unsigned a = static_cast<unsigned>(_springs[k].getParticle1().getId());
        const unsigned b = static_cast<unsigned>(_springs[k].getParticle2().getId());
        if ((a == i && b == j) || (a == j && b == i))
        {
            held = true;
            break;
        }
    }
    if (!held)
    {
        _particles[i].removeSpringNeighbor(j);
        _particles[j].removeSpringNeighbor(i);
    }
    _releasedsprings.insert(id);
    removeDynamicSpring(id);
    removeStaticSpring(id);
}

void SpringNetwork::updateSpringState(unsigned id, bool isStatic) {
    if (isStatic) {
        removeDynamicSpring(id);
        addStaticSpring(id);
    } else {
        removeStaticSpring(id);
        addDynamicSpring(id);
    }
}

void SpringNetwork::addParticle(const Particle & source)
{
    if (!_springs.empty())
        throw std::logic_error("SpringNetwork::addParticle: particles must be added before springs");

    Particle p = Particle(source);
    p.setSpringNetwork(this);
    // Particle::getId()/setId() stay signed (used as an "unassigned" sentinel
    // for the probe particle, see isProbeParticle()), so cast explicitly at
    // this array-index assignment and below where it's consumed as an
    // unsigned index (always >= 0 once assigned here).
    p.setId(static_cast<int>(_particles.size()));

    p.setInternalStructId(_structid);

    if (p.isStatic())
        addStaticParticle(static_cast<unsigned>(p.getId()));
    else
        addDynamicParticle(static_cast<unsigned>(p.getId()));

    if (p.isCharged())
        _chargedparticules.push_back(static_cast<unsigned>(p.getId()));

    if (p.isHydrophobic())
        _hydrophobicparticules.push_back(static_cast<unsigned>(p.getId()));

    _particles.push_back(p);
    _initparticles.push_back(p);
    _markNeighborSearchesDirty();
}

// AMBER's own V(phi), applied as a couple.
//
// V(phi) = sum_n V_n * (1 + cos(n*phi - gamma_n)), so the torque about the axis
// is -dV/dphi = sum_n n*V_n*sin(n*phi - gamma_n). A rigid side answers to that
// axial torque and to nothing else, which is why handing it to one atom per
// side is exact here rather than a distribution choice: the mesh has already
// removed every other freedom. The reaction on the two axis atoms is left to
// the same closed-form pass that serves the rings.
float SpringNetwork::computeTorsionForces()
{
    if (_torsions.empty())
        return 0.0f;

    const bool enabled[DIHEDRAL_FAMILY_COUNT] = {
        isDihedralPhiEnabled(),        isDihedralPsiEnabled(),           isDihedralOmegaEnabled(),
        isDihedralChiEnabled(),        isDihedralPlanarityEnabled(),     isDihedralNucleicBackboneEnabled(),
        isDihedralNucleicChiEnabled(), isDihedralNucleicSugarEnabled()};

    // Per family, and this is the only way to weigh the torsion term against
    // the rest of the model. `spring.scale` cannot do it: the torsion reads it,
    // but so does every rigid-body spring, so turning the torsions down turns
    // the mesh that holds bonds and angles down with them.
    //
    // The parameter already existed and was already accepted -- EnergySetting
    // carries a scale and every dihedral family is one, so `dihedralphi.scale`
    // parsed, printed in the configuration dump, and changed nothing at all.
    //
    // It is worth having because the two terms disagree about a helix. Measured
    // on a 34-residue beta-alpha-beta motif at 50 K, backbone hydrogen bonds
    // held out of 21: mesh alone 9, mesh + torsions 9, mesh + dynamic hbond 21,
    // all three together 15. The torsion term pulls phi/psi towards AMBER's own
    // minima, which are not the helical ones, and a hydrogen bond strong enough
    // to win that argument has to be worth -38 kJ/mol, which no backbone
    // hydrogen bond is.
    float familyscale[DIHEDRAL_FAMILY_COUNT];
    getDihedralFamilyScales(familyscale);

    const float PI = static_cast<float>(M_PI);
    const float unit = _ff->getSpringScale() * static_cast<float>(forcefield::GLOBAL_SPRING_FORCE_CONVERT);
    float energy = 0.0f;

    for (const Torsion & t : _torsions)
    {
        if (!enabled[t.family])
            continue;

        Particle & p1 = getParticle(t.atoms[0]);
        Particle & p2 = getParticle(t.atoms[1]);
        Particle & p3 = getParticle(t.atoms[2]);
        Particle & p4 = getParticle(t.atoms[3]);

        const Vector3f b1 = p2.getPosition() - p1.getPosition();
        const Vector3f b2 = p3.getPosition() - p2.getPosition();
        const Vector3f b3 = p4.getPosition() - p3.getPosition();

        const Vector3f n1 = b1 ^ b2;
        const Vector3f n2 = b2 ^ b3;
        const float n1sq = n1.dot(n1);
        const float n2sq = n2.dot(n2);
        const float b2len = b2.norm();
        if (n1sq < 1e-12f || n2sq < 1e-12f || b2len < 1e-6f)
            continue; // three atoms in line: the dihedral is not defined

        const float phi = std::atan2(b2len * b1.dot(n2), n1.dot(n2));

        // One atan2, one index, two lerps: the harmonics were summed once,
        // offline, and the table holds both what the energy is worth and how
        // hard it pulls.
        const TorsionTable & tab = _torsiontables[t.table];
        const int b = biospring_torsion_bin(phi, PI, static_cast<int>(tab.bins));
        const float f = biospring_torsion_fraction(phi, PI, static_cast<int>(tab.bins));
        const float fs = familyscale[t.family];
        energy += fs * (tab.energy[b] + f * (tab.energy[b + 1] - tab.energy[b]));
        const float torque = fs * (tab.torque[b] + f * (tab.torque[b + 1] - tab.torque[b]));
        if (torque == 0.0f)
            continue;

        // The exact gradient of phi, on all four atoms. A couple on the two
        // outer ones carries the same AXIAL torque and would be exact if the
        // two sides were rigid -- the mesh at 500 kJ.mol-1.A-2 is not, and the
        // difference measurably displaced the backbone by 3 degrees. This
        // costs a handful of cross products, needs no reaction on the axis and
        // conserves momentum by construction: the four forces sum to zero
        // identically, since F2 and F3 are built from F1 and F4.
        const float scale = unit * torque;
        const Vector3f F1 = n1 * biospring_torsion_k1(scale, b2len, n1sq);
        const Vector3f F4 = n2 * biospring_torsion_k4(scale, b2len, n2sq);
        const float b2lensq = b2len * b2len;
        const float c1 = biospring_torsion_c(b1.dot(b2), b2lensq);
        const float c3 = biospring_torsion_c(b3.dot(b2), b2lensq);
        // Signs matter here and are not guessable: the decomposition is written
        // in terms of r_ij = r_i - r_j and r_kl = r_k - r_l, which are the
        // NEGATIVES of b1 and b3 as spelled above. Checked against a finite
        // difference of phi itself -- 2.7e-10 -- after the other sign
        // convention passed the sum-to-zero test while being wrong by 5 rad/A.
        const Vector3f F2 = F1 * (-(c1 + 1.0f)) + F4 * c3;
        const Vector3f F3 = F1 * c1 - F4 * (c3 + 1.0f);

        p1.addForce(F1);
        p2.addForce(F2);
        p3.addForce(F3);
        p4.addForce(F4);
    }
    return energy;
}

void SpringNetwork::updateParticleState(unsigned id, bool isStatic) {
    if (isStatic) {
        removeDynamicParticle(id);
        addStaticParticle(id);
    } else {
        removeStaticParticle(id);
        addDynamicParticle(id);
    }
    _particles[id].setStatic(isStatic);
}

void SpringNetwork::clear()
{
    _initparticles.clear();
    _particles.clear();
    _staticparticules.clear();
    _dynamicparticules.clear();
    _chargedparticules.clear();
    _hydrophobicparticules.clear();
    _springs.clear();
    _staticsprings.clear();
    _dynamicsprings.clear();
    _springForceScratch.clear();
    _nsearch.steric.reset();
    _nsearch.electrostatic.reset();
    _nsearch.hydrophobic.reset();
    _nsearch.hbond.reset();
    _hbDonorSlot.clear();
    _hbDonorOffset.clear();
    _hbAcceptorSlot.clear();
    _hbAcceptorOffset.clear();
    _neighborSearchesDirty = false;
    _insertionVector.reset();
    _probeparticule = Particle();
}

// TODO: implement this function in `Topology`
// int SpringNetwork::removeSpringsBetweenSelections(std::vector<int> sel1, std::vector<int> sel2)
// {
//     int nbRemove = 0;
//     for (vector<int>::iterator it1 = sel1.begin(); it1 != sel1.end(); it1++)
//     {
//         for (vector<int>::iterator it2 = sel2.begin(); it2 != sel2.end(); it2++)
//         {
//             Spring * s = getSpringFromParticlesIds(*it1, *it2);
//             if (s)
//             {
//                 if (s->getParticle1().getId() == *it1)
//                     s->getParticle1().removeSpringNeighbor(*it2);
//                 else
//                     s->getParticle2().removeSpringNeighbor(*it2);

//                 if (s->getParticle1().getId() == *it2)
//                     s->getParticle1().removeSpringNeighbor(*it1);
//                 else
//                     s->getParticle2().removeSpringNeighbor(*it1);
//                 _springs.erase(std::find(_springs.begin(), _springs.end(), s));
//                 nbRemove++;
//             }
//         }
//     }

//     return nbRemove;
// }

// =====================================================================================

void SpringNetwork::setInsertionVector(unsigned aa1, unsigned aa2)
{
    try
    {
        getParticleFromId(aa1);
        getParticleFromId(aa2);
    }
    catch (const std::out_of_range & e)
    {
        logging::die("Insertion vector could not be initialized. One or both of the pair of particles "
                     "defining the insertion vector have an invalid identifier.");
    }

    //   Check that the two elements of the vector are identical
    if (aa1 == aa2)
    {
        logging::warning("The two values defining the insertion vector are identical: %zu and %zu",
                         aa1, aa2);
    }

    _insertionVector =
        std::make_unique<InsertionVector>(*this, getParticleFromId(aa1), getParticleFromId(aa2));
    
    const Particle p1 = _insertionVector->getParticle(0);
    const Particle p2 = _insertionVector->getParticle(1);
    logging::info("Insertion vector set between %s%d:%s and %s%d:%s",
        p1.getResName().c_str(), p1.getResId(), p1.getName().c_str(),
        p2.getResName().c_str(), p2.getResId(), p2.getName().c_str());
}

void SpringNetwork::applyConstraints()
{
    if (_constraints.empty())
    {
        _meanConstraintsDistances = 0.0f;
        return;
    }

    float sumDistances = 0.0;
    for (unsigned i = 0; i < _constraints.size(); i++)
    {
        _constraints[i]->apply();
        sumDistances += _constraints[i]->getDistance();
    }
    _meanConstraintsDistances = sumDistances / _constraints.size();
}

void SpringNetwork::clearParticles(void) { clear(); }

// =====================================================================================
//
// Setup methods
//
// =====================================================================================

void SpringNetwork::setup(const configuration::Configuration & conf)
{
    _config = conf;

    _setupForceField();
    _setupProbe();
    _setupSteric();
    _setupElectrostatic();
    _setupHydrophobic();
    _setupHydrogenBond();
    _setupStacking();
    // After the hydrogen bond term, whose out-of-plane lobe sites are what
    // brings a nucleophile into the attack conformation this then acts on.
    _peptidebond.setup(*this, conf.peptidebond);
    _setupDensityGrid();
    _setupInsertionVector();
    _setupTrajectories();
    // After everything else: it converts particles and adds springs, so it
    // needs the network whole and the configuration already stored.
    _neighborSearchesDirty = false;
    // _setupConstraints();
    // _setupSelections();
}

// Picks the cell width, in A, from the stencil radius the longest-reaching term
// should walk -- and from nothing else.
//
// WHY NOTHING ELSE, AND PARTICULARLY NOT A DENSITY. The obvious model is to
// price a width at "volume searched x particles per A3, plus what the cells
// cost to visit", and that is what this did. It does not work here, because
// BioSpring has no solvent: a structure is a molecule in vacuum, so there is no
// density to speak of. 76% of the capsid's bounding box is empty and so is 74%
// of the nucleosome's, and a fit of the cell-visit cost against measured times
// came out a factor of 4 apart on the CPU between the two, and a factor of 7 on
// the device. A number that moves by 7x between two examples is not a constant.
//
// What IS stable is the stencil radius. Measured on those same two examples --
// 25069 beads at 0.0175 A-3 with cutoffs of 9 and 16 A, and 37200 beads at
// 0.0015 A-3, twelve times sparser, with cutoffs of 5 and 16 A -- the best width
// is the one that puts the longest term at a radius of 2 cells on the CPU and 3
// on the device, in both cases, to within 1% of the measured optimum. It is a
// dimensionless quantity, which is why it survives the change of system where
// an energy per unit volume does not.
//
// Among the widths that give the longest term exactly that radius, the one that
// walks the fewest cells in total. That is what separates 9 A from 8 A on the
// nucleosome: both put the 16 A coulomb at two cells, but 9 A also puts the 9 A
// steric term at one, 27 cells against 125, and it measured 5% faster.

void SpringNetwork::_setupSteric()
{
    if (isStericEnabled())
    {
        if (getStericCutoff() < 1e-6)
            throw std::runtime_error("Steric cutoff must be > 0");
        _nsearch.steric = make_nsearch(_particles, getStericCutoff(), getNeighborSkin(), getCellWidthFor(getStericCutoff()));
        _excludeProbeFromNeighborSearch(*_nsearch.steric);
    }
}

void SpringNetwork::_setupHydrophobic()
{
    if (isHydrophobicityEnabled())
    {
        if (getHydrophobicCutoff() < 1e-6)
            throw std::runtime_error("Hydrophobic cutoff must be > 0");
        const std::vector<size_t> hydrophobic_particles = _hydrophobicParticleIndexes();
        if (!hydrophobic_particles.empty())
        {
            _nsearch.hydrophobic =
                make_nsearch(_particles, getHydrophobicCutoff(), hydrophobic_particles, getNeighborSkin(), getCellWidthFor(getHydrophobicCutoff()));
            _excludeProbeFromNeighborSearch(*_nsearch.hydrophobic);
        }
    }
}

void SpringNetwork::_setupHydrogenBond()
{
    // The thermostat's friction IS viscosity.value. Asking for one without the
    // other is asking a solvent to kick without dragging, which holds no
    // temperature at all -- the velocity decay would be 1, the kick would come
    // out of sqrt(1 - 1) = 0, and the run would silently be the unthermostatted
    // one. Better to say so than to let someone report a temperature they
    // never had.
    if (isThermostatEnabled() && !isViscosityEnabled())
        logging::die("thermostat.enable needs viscosity.enable: the thermostat's friction is viscosity.value, "
                     "and without it there is nothing for the random kicks to balance.");
    if (isThermostatEnabled() && getTemperature() <= 0.0f)
        logging::die("thermostat.temperature must be > 0 K when the thermostat is enabled "
                     "(use thermostat.enable = 0 for a pure viscous brake).");
    if (isThermostatEnabled())
        logging::info("Langevin thermostat: %.1f K, friction %.4g Da.fs-1.", getTemperature(), getViscosity());

    if (isHydrogenBondEnabled())
    {
        if (getHydrogenBondCutoff() < 1e-6)
            throw std::runtime_error("Hydrogen bond cutoff must be > 0");
        // One grid for the term, holding donors AND acceptors together; both
        // query it and filter candidates by role themselves (see
        // _assignHydrogenBondPairs), rather than pre-splitting them into two.
        // Its cells are the size of its own search radius, like every other
        // term's -- see getCellWidthFor.
        const std::vector<size_t> donor_acceptor_particles = _donorAcceptorParticleIndexes();
        if (!donor_acceptor_particles.empty())
        {
            _nsearch.hbond =
                make_nsearch(_particles, getHydrogenBondCutoff(), donor_acceptor_particles, getNeighborSkin(),
                             getCellWidthFor(getHydrogenBondCutoff()));
            _excludeProbeFromNeighborSearch(*_nsearch.hbond);
        }

        // Every particle starts unbonded; bonds are matched dynamically as
        // the simulation runs (see _assignHydrogenBondPairs). One slot per
        // donatable hydrogen and one per lone pair, laid out CSR so a
        // particle's slots are contiguous.
        const size_t n = _particles.size();
        _hbDonorOffset.assign(n + 1, 0);
        _hbAcceptorOffset.assign(n + 1, 0);
        for (size_t i = 0; i < n; ++i)
        {
            _hbDonorOffset[i + 1] = _hbDonorOffset[i] + _particles[i].donorCapacity();
            _hbAcceptorOffset[i + 1] = _hbAcceptorOffset[i] + _particles[i].acceptorCapacity();
        }
        _hbDonorSlot.assign(_hbDonorOffset[n], -1);
        _hbAcceptorSlot.assign(_hbAcceptorOffset[n], -1);

        logging::info("Hydrogen bond slots: %zu donor, %zu acceptor.", _hbDonorSlot.size(), _hbAcceptorSlot.size());

        // How many candidate pairs a SPRING has already pre-empted.
        //
        // The assignment refuses any pair its two atoms are sprung by, and that
        // is right -- a spring and a bond on the same pair would count one
        // interaction twice. But it was silent, and silence cost a great deal:
        // a 70S assembled with `mergespn --cutoff 7.0` carried 1780 springs
        // across its own base-edge atoms, at the merge's default stiffness of
        // 1.0 against 650 for the all-atom mesh. Removing all of them moved the
        // spring energy from 40.03 to 40.07 kJ/mol -- they held nothing -- and
        // took the real base pairs held from 16 of 49 to 48 of 49.
        //
        // Mechanically invisible, functionally decisive, and reported by
        // nothing: the spring energy looked normal and the bond count looked
        // plausible. One line here would have said it immediately.
        if (_nsearch.hbond)
        {
            size_t preempted = 0;
            for (const size_t i : donor_acceptor_particles)
            {
                const Particle & p = getParticle(i);
                _nsearch.hbond->for_each_neighbor(i, [&](size_t j) {
                    if (j <= i || isProbeParticle(j))
                        return;
                    const Particle & q = getParticle(j);
                    if (p.getResId() == q.getResId() && p.getChainName() == q.getChainName())
                        return;
                    if (!((p.isDonor() && q.isAcceptor()) || (q.isDonor() && p.isAcceptor())))
                        return;
                    if (Particle::distance(p, q) < getHydrogenBondCutoff() &&
                        p.isInSpringNeighbors(static_cast<unsigned>(q.getId())))
                        preempted++;
                });
            }
            if (preempted)
                logging::warning("%zu donor/acceptor pair(s) within the hydrogen bond cutoff are joined by a "
                                 "SPRING, so the term will decline them -- a spring and a bond on one pair "
                                 "would count it twice. If this network was assembled with a cutoff, that "
                                 "cutoff is welding the pairs the term is meant to form.",
                                 preempted);
        }
    }
}

void SpringNetwork::_setupForceField()
{
    const std::string steric = _config.steric.mode;

    if (steric == "lennard-jones-8-6Lewitt")
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Lewitt>();
    else if (steric == "lennard-jones-8-6Zacharias")
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLennardJones_8_6Zacharias>();
    else if (steric == "lennard-jones-12-6Amber")
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLennardJones_12_6Amber>();
    else
        _ff = std::make_unique<forcefield::ForceFieldElectrostaticCoulombAndStericLinear>();

    // "default" leaves the law's own rule, which each subclass sets in its
    // constructor. Anything else overrides it, which is what a .ff paired with a
    // law it did not come from needs -- the `radius` column does not hold the
    // same quantity in every file. Zacharias' law does not read it.
    const std::string radiusrule = _config.steric.radiusrule;
    if (radiusrule == "sum")
        _ff->setRadiusRule(BIOSPRING_RADIUS_SUM);
    else if (radiusrule == "geometric-mean")
        _ff->setRadiusRule(BIOSPRING_RADIUS_GEOMETRIC_MEAN);
    else if (radiusrule == "arithmetic-mean")
        _ff->setRadiusRule(BIOSPRING_RADIUS_ARITHMETIC_MEAN);

    _ff->setStericScale(_config.steric.gridscale);
    _ff->setCoulombScale(_config.electrostatic.scale);
    _ff->setDielectric(_config.electrostatic.dielectric);
    _ff->setDistanceDependentDielectric(_config.electrostatic.distancedependent);
    // electrostaticgrid, not potentialgrid: the grid was renamed for what it
    // holds (948a0dd), and the old name no longer exists.
    _ff->setForceFieldScale(_config.electrostaticgrid.scale);
    _ff->setSpringScale(_config.spring.scale);
    _ff->setIMPScale(_config.imp.scale);
    _ff->setHydrophobicityScale(_config.hydrophobicity.scale);
    _ff->setHydrophobicityDecayLength(_config.hydrophobicity.decaylength);
    _ff->setHydrogenBondScale(_config.hbond.scale);
}

void SpringNetwork::_setupElectrostatic()
{
    if (!isAnyElectrostaticEnabled())
        return;

    if (isElectrostaticFieldEnabled())
    {
        const std::string dxpath = _config.electrostaticgrid.path;
        logging::info("Reading electrostatic map from DX file '%s'", dxpath.c_str());
        _grids.electrostatic = opendx::readGrid(dxpath);
    }

    if (isElectrostaticCoulombEnabled())
    {
        if (getElectrostaticCutoff() < 1e-6)
            throw std::runtime_error("Electrostatic cutoff must be > 0");

        const std::vector<size_t> charged_particles = _chargedParticleIndexes();
        if (!charged_particles.empty())
        {
            _nsearch.electrostatic =
                make_nsearch(_particles, getElectrostaticCutoff(), charged_particles, getNeighborSkin(), getCellWidthFor(getElectrostaticCutoff()));
            _excludeProbeFromNeighborSearch(*_nsearch.electrostatic);
        }
    }
}

void SpringNetwork::_setupDensityGrid()
{
    if (isDensityGridEnabled())
    {
        const std::string dxpath = _config.densitygrid.path;
        logging::info("Reading density grid from DX file '%s'", dxpath.c_str());
        _grids.density = opendx::readGrid(dxpath);
    }
}

void SpringNetwork::_setupProbe()
{
    if (isProbeEnabled())
    {
        _probeparticule.setSpringNetwork(this);
        _probeparticule.setName("PRB");
        _probeparticule.setResName("PRB");
        _probeparticule.setCharge(_config.probe.charge);
        _probeparticule.setRadius(_config.probe.radius);
        _probeparticule.setEpsilon(_config.probe.epsilon);
        _probeparticule.setMass(_config.probe.mass);
        _probeparticule.setX(_config.probe.x);
        _probeparticule.setY(_config.probe.y);
        _probeparticule.setZ(_config.probe.z);
        _probeparticule.setId(static_cast<int>(_particles.size()));
        _particles.push_back(_probeparticule);
    }
}

void SpringNetwork::_setupTrajectories()
{
    if (_config.pdbtraj.enable)
        _trajectories.add_writer(
            std::make_unique<io::modern::PDBTrajectoryWriter>(_config.pdbtraj.path, *this, _config.pdbtraj.frequency));
    if (_config.xtctraj.enable)
        _trajectories.add_writer(
            std::make_unique<io::modern::XTCTrajectoryWriter>(_config.xtctraj.path, *this, _config.xtctraj.frequency));
    if (_config.csvsample.enable)
        _trajectories.add_writer(std::make_unique<io::modern::CSVTrajectoryWriter>(_config.csvsample.path, *this,
                                                                                   _config.csvsample.frequency));
}

void SpringNetwork::_setupInsertionVector()
{
    if (_config.ivector.enable)
        setInsertionVector(_config.ivector.vector[0], _config.ivector.vector[1]);
}


std::vector<size_t> SpringNetwork::_chargedParticleIndexes() const
{
    std::vector<size_t> indexes;
    indexes.reserve(_particles.size());

    for (size_t i = 0; i < _particles.size(); ++i)
    {
        if (_particles[i].isCharged())
            indexes.push_back(i);
    }

    return indexes;
}

std::vector<size_t> SpringNetwork::_hydrophobicParticleIndexes() const
{
    std::vector<size_t> indexes;
    indexes.reserve(_particles.size());

    for (size_t i = 0; i < _particles.size(); ++i)
    {
        if (_particles[i].isHydrophobic())
            indexes.push_back(i);
    }

    return indexes;
}

std::vector<size_t> SpringNetwork::_donorAcceptorParticleIndexes() const
{
    std::vector<size_t> indexes;
    indexes.reserve(_particles.size());

    for (size_t i = 0; i < _particles.size(); ++i)
    {
        if (_particles[i].isDonor() || _particles[i].isAcceptor())
            indexes.push_back(i);
    }

    return indexes;
}

void SpringNetwork::_assignHydrogenBondPairs()
{
    const float cutoff = getHydrogenBondCutoff();
    const int exponent = getHydrogenBondAngularExponent();
    const size_t n = _particles.size();
    if (_hbDonorOffset.size() != n + 1)
        return;

    // Step 1: free any slot whose bond has drifted beyond the cutoff, or whose
    // two ends have become SPRUNG. Only the donor side is walked -- the
    // acceptor side is cleared with it, so the two never disagree.
    //
    // The sprung test is not symmetry with the proposal step for its own sake.
    // Nothing used to change the springs during a run, so a pair could only be
    // sprung before it was ever proposed, and refusing it there was enough.
    // PeptideBondFormation does change them, and a bond held from before the
    // change stays held: at the peptide C-N length of 1.33 A the Morse well is
    // deep into its repulsive wall, so the hydrogen bond term ends up pushing
    // apart exactly the pair the new springs are holding together -- measured
    // at +4.31 kJ/mol on two alanines, a positive "hydrogen bond energy".
    for (size_t i = 0; i < n; ++i)
        for (size_t sd = _hbDonorOffset[i]; sd < _hbDonorOffset[i + 1]; ++sd)
        {
            const int j = _hbDonorSlot[sd];
            if (j < 0)
                continue;
            // Beyond the cutoff, newly sprung, or no longer WORTH its slot.
            // The third test matters because a slot is held until released: a
            // contact that was worth taking when it formed, and that a factor
            // has since reduced to nothing, would otherwise keep the atom for
            // the rest of the run.
            const bool far =
                Particle::distance(getParticle(i), getParticle(static_cast<size_t>(j))) > cutoff;
            const bool sprung = getParticle(i).isInSpringNeighbors(static_cast<unsigned>(j));
            bool worthless = false;
            if (!far && !sprung && getHydrogenBondMinimumStrength() > 0.0f)
            {
                const HydrogenBondParameters mp = hydrogenBondParameters(i, static_cast<size_t>(j));
                const float d = Particle::distance(getParticle(i), getParticle(static_cast<size_t>(j)));
                float w = 1.0f;
                if (d > 1e-6f)
                {
                    const Vector3f vhat =
                        (getParticle(static_cast<size_t>(j)).getPosition() - getParticle(i).getPosition()) / d;
                    const Vector3f hd = donorDirection(getParticle(i), vhat);
                    if (hd.norm() > 1e-6f)
                        w *= forcefield::hydrogen_bond_angular_factor(hd.dot(vhat), exponent);
                    const Vector3f ha = donorDirection(getParticle(static_cast<size_t>(j)), -vhat);
                    if (ha.norm() > 1e-6f)
                        w *= forcefield::hydrogen_bond_angular_factor(ha.dot(-vhat), exponent);
                }
                const float strength =
                    _ff->computeHydrogenBondEnergy(d, mp.welldepth, mp.equilibrium, mp.width) * w *
                    hydrogenBondPlanarity(i, static_cast<size_t>(j)).factor;
                worthless = strength > -getHydrogenBondMinimumStrength();
            }
            if (!far && !sprung && !worthless)
                continue;
            _hbDonorSlot[sd] = -1;
            for (size_t sa = _hbAcceptorOffset[j]; sa < _hbAcceptorOffset[j + 1]; ++sa)
                if (_hbAcceptorSlot[sa] == static_cast<int>(i))
                {
                    _hbAcceptorSlot[sa] = -1;
                    break;
                }
        }

    auto free_donor = [&](size_t i) {
        for (size_t s = _hbDonorOffset[i]; s < _hbDonorOffset[i + 1]; ++s)
            if (_hbDonorSlot[s] < 0)
                return true;
        return false;
    };
    auto free_acceptor = [&](size_t i) {
        for (size_t s = _hbAcceptorOffset[i]; s < _hbAcceptorOffset[i + 1]; ++s)
            if (_hbAcceptorSlot[s] < 0)
                return true;
        return false;
    };

    const std::vector<size_t> candidates = _donorAcceptorParticleIndexes();

    // Steps 2 and 3, repeated: propose, then confirm mutual proposals. One
    // round fills one slot per particle, so as many rounds as the largest
    // capacity are enough -- 2 or 3 in practice. The loop stops early as
    // soon as a round confirms nothing, so the common single-slot case
    // costs exactly what it did before.
    //
    // AN OPEN QUESTION, measured but not explained. On 8PKL's all-atom region
    // this holds 16 of the 49 base pairs the crystal has. The 33 it misses are
    // not weak: asked through dumpHydrogenBondCandidates -- which uses these
    // same functions, so it cannot disagree -- their best atom couple scores
    // -15.7 kJ/mol at the median against -10.8 for the pairs that DO form, with
    // better distances, better angular weights and better planarity.
    //
    // Three explanations were tried and none survived:
    //   25 of the 33 have no free DONOR facing a free ACCEPTOR, so slots are
    //     part of it -- but freeing 320 parasite slots with
    //     hbond.minimumstrength recovered none of them;
    //   the remaining 8 are fully pairable and still unformed, which looks like
    //     the reciprocal-best rule stranding a pair whose ends prefer others --
    //     but adding a Gale-Shapley fallback (a refused particle tries its next
    //     choice) gained nothing and cost one real pair, 16 -> 15;
    //   raising the round count from 4 to 24 changed the result not at all,
    //     bit for bit.
    // Whatever is missing, it is none of those. Worth settling on a small
    // controlled case where the assignment can be watched, not on a ribosome.
    for (unsigned round = 0; round < 4; ++round)
    {
        std::vector<int> nearest(n, -1);
        // Holds the best bond STRENGTH seen so far, in kJ/mol. It starts at
        // MINUS THE FLOOR rather than at zero: a candidate whose angular or
        // planarity factor kills it is worth exactly nothing and must not take
        // a slot from a real partner, which is what the previous version of
        // this comment said -- but a threshold of zero exactly admits a bond
        // worth -0.07 kJ/mol, and that bond then holds the slot while a real
        // pair worth -15 finds the atom taken. Measured on the ribosome: 330
        // slots held by stacked sequence neighbours worth 0.000 kJ/mol at the
        // median. See hbond.minimumstrength, which is zero by default.
        const float floor = -getHydrogenBondMinimumStrength();
        std::vector<float> nearest_distance(n, floor);

#ifdef OPENMP_SUPPORT
#pragma omp parallel for schedule(static)
#endif
        for (size_t k = 0; k < candidates.size(); ++k)
        {
            const size_t i = candidates[k];
            const bool i_donor = free_donor(i);
            const bool i_acceptor = free_acceptor(i);
            if (!i_donor && !i_acceptor)
                continue;

            const Particle & p = getParticle(i);
            _nsearch.hbond->for_each_neighbor(p, [&](size_t j) {
                if (isProbeParticle(j))
                    return;

                // A free donor slot needs a free acceptor slot facing it. A
                // particle that is both (a hydroxyl) can pair either way,
                // but never twice with the SAME partner -- that would be one
                // bond counted as two.
                const bool roles_match = (i_donor && free_acceptor(j)) || (i_acceptor && free_donor(j));
                if (!roles_match || areHydrogenBonded(i, j))
                    return;

                const Particle & q = getParticle(j);

                // A residue's own backbone N and O (or a side chain's own
                // donor and acceptor atom, e.g. Asn's ND2/OD1) sit at a
                // fixed, short covalent-geometry distance -- not a real
                // hydrogen bond, and not necessarily an explicit spring
                // neighbour either. Without this, it is invariably the
                // closest candidate and starves the real inter-residue bond.
                if (p.getResId() == q.getResId() && p.getChainName() == q.getChainName())
                    return;

                if (isSpringEnabled() && p.isInSpringNeighbors(static_cast<unsigned>(q.getId())))
                    return;

                const float distance = Particle::distance(p, q);
                if (distance >= cutoff)
                    return;

                // Ranked by what the bond is actually WORTH -- the Morse well
                // times the angular factor -- and not by distance alone.
                //
                // Distance alone hands the slot to whichever candidate is
                // nearest even when the angle says the bond cannot exist, and
                // the angular factor then reduces it to zero energy and zero
                // force. It still occupies the slot. Measured on 071 before
                // this: 47 bonds held, 11 of them with an angular weight of
                // exactly zero, 43 between residues one or two apart, and NOT
                // ONE at i,i+-4 -- the neighbours a chain always has nearby
                // had saturated every donor, so no helical partner could ever
                // win one. Scoring by the real strength drops them without a
                // residue-separation rule, which would have thrown away the
                // legitimate side-chain bonds between adjacent residues too.
                float weight = 1.0f;
                if (distance > 1e-6f)
                {
                    const Vector3f vhat = (q.getPosition() - p.getPosition()) / distance;
                    const Vector3f hd = donorDirection(p, vhat);
                    if (hd.norm() > 1e-6f)
                        weight *= forcefield::hydrogen_bond_angular_factor(hd.dot(vhat), exponent);
                    // The acceptor's own geometry decides too: a lone pair
                    // pointing away forbids the bond as surely as a hydrogen
                    // pointing away does.
                    const Vector3f ha = donorDirection(q, -vhat);
                    if (ha.norm() > 1e-6f)
                        weight *= forcefield::hydrogen_bond_angular_factor(ha.dot(-vhat), exponent);
                }
                // Ranked on the well this COUPLE would use, not the force
                // field's: a pair that carries a group is worth what its own
                // group says, and ranking it as an ordinary hydrogen bond
                // would let a shallower neighbour take the slot.
                const HydrogenBondParameters mp = hydrogenBondParameters(i, j);
                // And on the planarity factor, so the ranking sees the same
                // quantity the energy and the force will. Ranking without it
                // would hand a slot to a skewed contact that the force term
                // then values at almost nothing -- which is the failure the
                // angular factor was added to this ranking to fix, repeated.
                const float strength =
                    _ff->computeHydrogenBondEnergy(distance, mp.welldepth, mp.equilibrium, mp.width) * weight *
                    hydrogenBondPlanarity(i, j).factor;
                if (strength < nearest_distance[i])
                {
                    nearest_distance[i] = strength;
                    nearest[i] = static_cast<int>(j);
                }
            });
        }

        // Reciprocal best hit, the same criterion as before -- `nearest` is
        // fixed by this point, so the outcome does not depend on the
        // (serial, deterministic) order pairs are confirmed in.
        unsigned confirmed = 0;
        for (size_t k = 0; k < candidates.size(); ++k)
        {
            const size_t i = candidates[k];
            const int j = nearest[i];
            if (j < 0 || nearest[static_cast<size_t>(j)] != static_cast<int>(i))
                continue;
            const size_t jj = static_cast<size_t>(j);
            if (areHydrogenBonded(i, jj))
                continue;

            // Whichever way round still has both slots free. i first, so a
            // hydroxyl pair resolves deterministically.
            size_t donor = i;
            size_t acceptor = jj;
            if (!(free_donor(i) && free_acceptor(jj)))
            {
                if (!(free_donor(jj) && free_acceptor(i)))
                    continue;
                donor = jj;
                acceptor = i;
            }
            for (size_t sd = _hbDonorOffset[donor]; sd < _hbDonorOffset[donor + 1]; ++sd)
                if (_hbDonorSlot[sd] < 0)
                {
                    _hbDonorSlot[sd] = static_cast<int>(acceptor);
                    break;
                }
            for (size_t sa = _hbAcceptorOffset[acceptor]; sa < _hbAcceptorOffset[acceptor + 1]; ++sa)
                if (_hbAcceptorSlot[sa] < 0)
                {
                    _hbAcceptorSlot[sa] = static_cast<int>(donor);
                    break;
                }
            ++confirmed;
        }
        if (confirmed == 0)
            break;
    }
}

void SpringNetwork::_excludeProbeFromNeighborSearch(NeighborSearch::Searcher & searcher)
{
    if (!isProbeEnabled())
        return;

    const int probe_id = _probeparticule.getId();
    if (probe_id >= 0)
        searcher.exclude_index(static_cast<size_t>(probe_id));
}

void SpringNetwork::_markNeighborSearchesDirty()
{
    if (_nsearch.steric || _nsearch.electrostatic || _nsearch.hydrophobic || _nsearch.hbond)
        _neighborSearchesDirty = true;
}

void SpringNetwork::_updateNeighborSearches()
{
    if (!_neighborSearchesDirty)
        return;

    if (_nsearch.steric)
        _nsearch.steric->update();
    if (_nsearch.electrostatic)
        _nsearch.electrostatic->update();
    if (_nsearch.hydrophobic)
        _nsearch.hydrophobic->update();
    if (_nsearch.hbond)
        _nsearch.hbond->update();

    _neighborSearchesDirty = false;
}

void SpringNetwork::_resizeNonbondedPairScratch()
{
    const size_t n = _dynamicparticules.size();

    _stericPairScratch.resize(n);
    _electrostaticPairScratch.resize(n);
    _hydrophobicPairScratch.resize(n);
    _hydrogenBondCoreRepulsionPairScratch.resize(n);

    // Clears logical contents but keeps each bucket's capacity, so the
    // simulation loop does not reallocate every step.
    for (auto & bucket : _stericPairScratch)
        bucket.clear();
    for (auto & bucket : _electrostaticPairScratch)
        bucket.clear();
    for (auto & bucket : _hydrophobicPairScratch)
        bucket.clear();
    for (auto & bucket : _hydrogenBondCoreRepulsionPairScratch)
        bucket.clear();
}

void SpringNetwork::_applyNonbondedPairScratch(
    const std::vector<std::vector<spn::DeferredNonbondedContribution>> & scratch, float & energy)
{
    for (const auto & bucket : scratch)
    {
        for (const auto & contribution : bucket)
        {
            getParticle(contribution.target).addForce(contribution.force);
            energy += contribution.energy;
        }
    }
}

void SpringNetwork::_syncProbeParticle()
{
    if (!isProbeEnabled())
        return;

    const int probe_id = _probeparticule.getId();
    if (probe_id < 0 || static_cast<size_t>(probe_id) >= _particles.size())
        return;

    _particles[static_cast<size_t>(probe_id)] = _probeparticule;
    _particles[static_cast<size_t>(probe_id)].setSpringNetwork(this);
}

void SpringNetwork::_rebuildSpringNeighbors()
{
    for (Particle & particle : _particles)
        particle.clearSpringNeighbors();

    for (unsigned id = 0; id < _springs.size(); ++id)
    {
        if (_releasedsprings.count(id))
            continue; // broken: it holds nothing and excludes nothing
        Spring & spring = _springs[id];
        Particle & p1 = spring.getParticle1();
        Particle & p2 = spring.getParticle2();
        p1.addToSpringNeighbors(static_cast<unsigned>(p2.getId()), &spring);
        p2.addToSpringNeighbors(static_cast<unsigned>(p1.getId()), &spring);
    }
}

void SpringNetwork::_setupSelections() { throw "Not Implemented Error"; }
void SpringNetwork::_setupConstraints() { throw "Not Implemented Error"; }

} // namespace spn
} // namespace biospring
