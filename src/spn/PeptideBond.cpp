#include "spn/PeptideBond.h"

#include <algorithm>
#include <cmath>

#include "IO/ReduceRuleReader.h"
#include "IO/RigidBodyRuleReader.h"
#include "forcefield/energy/hydrogenbond.hpp"
#include "logging.h"
#include "spn/SpringNetwork.h"
#include "utils.hpp"

namespace biospring
{
namespace spn
{

namespace
{

// The peptide bond's OWN internal coordinates, Engh & Huber (1991). Only the
// three that describe the BOND: everything else about the two residues is
// already in the structure, and taking it from a table instead is exactly what
// made the product strained -- see referenceGeometry below.
const float PEPTIDE_CN = 1.329f;         // A
const float PEPTIDE_CA_C_N = 116.2f;     // deg, at the carbonyl carbon
const float PEPTIDE_C_N_CA = 121.7f;     // deg, at the amine nitrogen
// omega (CA-C-N-CA) is 180 by construction below, and the unit is trans.

// Where the fifteen springs of a new peptide plane should END UP, as POSITIONS
// for the group's atoms.
//
// WHY NOT A TABLE OF FIFTEEN DISTANCES. Six points in space have fifteen
// pairwise distances, which is exactly enough to fix them up to a reflection --
// so a table determines omega completely, and the first version of this used
// one, taken from a fully idealised trans plane. It came out 17.3 degrees from
// trans, with the six atoms 0.195 A off coplanar.
//
// The reason is that six of those fifteen pairs are INTERNAL to one residue or
// the other, are already sprung by the _CA group, and therefore keep the rest
// length the structure gave them -- they are never ramped. The idealised table
// set the other nine from a geometry those six do not have. Measured on two
// alanines: the residue's own amide hydrogen sits at 108.6 degrees of N-CA
// where the table assumed 119.1, so H...CA was 2.023 A against a table that had
// built its cross terms around 2.104. The fifteen distances were then not
// realisable by any arrangement of six points, and the mesh settled in the
// least-bad strained one, which is what twisted omega. Making that ONE distance
// consistent took omega from 162.7 to -174.2 degrees and the out-of-plane error
// from 0.195 to 0.049 A.
//
// So the reference is BUILT, from what is actually there: each residue keeps
// its own geometry exactly, and only their relative placement comes from the
// table. The fifteen distances are then realisable by construction.
//
// Returns false when the group does not name the four atoms the construction
// needs, which is how a rule that is not a peptide plane opts out of the ramp
// rather than being drawn towards a geometry that means nothing for it.
// normalize() is in place and returns void, so this keeps each direction one
// expression instead of three statements.
Vector3f unitOf(const Vector3f & v)
{
    Vector3f u = v;
    u.normalize();
    return u;
}

bool referenceGeometry(const std::vector<Vector3f> & position, const std::vector<std::string> & refname,
                       const std::string & electrophile, const std::string & nucleophile,
                       std::vector<Vector3f> & reference)
{
    const auto find = [&](const std::string & want) {
        for (size_t k = 0; k < refname.size(); ++k)
            if (refname[k] == want)
                return static_cast<int>(k);
        return -1;
    };
    const int iCA1 = find("CA"), iC = find(electrophile), iO = find("O");
    const int iN = find("+" + nucleophile), iCA2 = find("+CA");
    if (iCA1 < 0 || iC < 0 || iN < 0 || iCA2 < 0)
        return false;

    const Vector3f C = position[iC];
    const Vector3f uCA = unitOf(position[iCA1] - C);
    // The carbonyl's own sp2 plane, which the new bond lies in. Without an
    // oxygen in the group there is no plane to inherit, so any perpendicular
    // will do: the result is still planar, just not tied to the carbonyl.
    Vector3f normal = iO >= 0 ? (uCA ^ (position[iO] - C)) : Vector3f(0.0f, 0.0f, 0.0f);
    if (normal.norm() < 1e-4f)
    {
        const Vector3f trial = std::abs(uCA.getX()) < 0.9f ? Vector3f(1.0f, 0.0f, 0.0f) : Vector3f(0.0f, 1.0f, 0.0f);
        normal = uCA ^ trial;
    }
    normal = unitOf(normal);
    const Vector3f inplane = unitOf(normal ^ uCA);

    const float a1 = PEPTIDE_CA_C_N * static_cast<float>(M_PI) / 180.0f;
    // Two directions make that angle with C->CA; the nitrogen takes the one
    // AWAY from the oxygen, which is the only thing that distinguishes them.
    Vector3f dirN = uCA * std::cos(a1) + inplane * std::sin(a1);
    if (iO >= 0 && dirN.dot(unitOf(position[iO] - C)) > 0.0f)
        dirN = uCA * std::cos(a1) - inplane * std::sin(a1);
    const Vector3f N = C + dirN * PEPTIDE_CN;

    // The second alpha carbon, at the residue's OWN N-CA length, trans across
    // the new bond: of the two in-plane choices, the one that puts it on the
    // far side from the first alpha carbon.
    const Vector3f uNC = unitOf(C - N);
    const Vector3f inplane2 = unitOf(normal ^ uNC);
    const float a2 = PEPTIDE_C_N_CA * static_cast<float>(M_PI) / 180.0f;
    const float lNCA = (position[iCA2] - position[iN]).norm();
    const Vector3f plus = N + (uNC * std::cos(a2) + inplane2 * std::sin(a2)) * lNCA;
    const Vector3f minus = N + (uNC * std::cos(a2) - inplane2 * std::sin(a2)) * lNCA;
    const Vector3f CA2 = (plus - position[iCA1]).norm() > (minus - position[iCA1]).norm() ? plus : minus;

    // Everything else on the amine side rides along: that residue is a rigid
    // body here, so its remaining atoms are placed by the same rotation that
    // takes its own (N, CA) onto the two positions just built. The spin about
    // that axis is fixed by putting the amide hydrogen IN the peptide plane,
    // on the side away from the second alpha carbon -- which is the trans amide.
    const Vector3f e1 = unitOf(position[iCA2] - position[iN]);
    const Vector3f E1 = unitOf(CA2 - N);
    Vector3f e2, E2;
    int spin = -1;
    for (size_t k = 0; k < refname.size(); ++k)
        if (refname[k].size() > 1 && refname[k][0] == '+' && static_cast<int>(k) != iN &&
            static_cast<int>(k) != iCA2)
        {
            spin = static_cast<int>(k);
            break;
        }
    if (spin >= 0)
    {
        const Vector3f v = position[spin] - position[iN];
        e2 = v - e1 * e1.dot(v);
        if (e2.norm() < 1e-4f)
            spin = -1;
        else
            e2 = unitOf(e2);
    }
    if (spin < 0)
    {
        // Nothing to orient by: any perpendicular keeps the construction valid.
        const Vector3f trial = std::abs(e1.getX()) < 0.9f ? Vector3f(1.0f, 0.0f, 0.0f) : Vector3f(0.0f, 1.0f, 0.0f);
        e2 = unitOf(e1 ^ trial);
        E2 = unitOf(E1 ^ trial);
    }
    else
    {
        // Which of the two in-plane perpendiculars, and the answer is NOT
        // "the one pointing at the first alpha carbon". E2 is perpendicular to
        // N->CA, while the side an amide hydrogen sits on is a rotation about
        // C->N -- a different axis. Testing against the wrong one put the
        // hydrogen trans to the first alpha carbon instead of cis, i.e. on top
        // of the second one, and the mesh then had to compromise between that
        // and omega: it settled 15 degrees off trans.
        //
        // So both signs are built and the one that makes the spin atom CIS to
        // the first alpha carbon across the new bond wins -- omega(CA, C, N, H)
        // near 0, which is what a trans amide is.
        const Vector3f perp = unitOf(normal ^ E1);
        float best = 0.0f;
        for (int sign = 0; sign < 2; ++sign)
        {
            const Vector3f candidate = sign == 0 ? perp : perp * -1.0f;
            const Vector3f third = E1 ^ candidate;
            const Vector3f v = position[spin] - position[iN];
            const Vector3f placed =
                N + E1 * v.dot(e1) + candidate * v.dot(e2) + third * v.dot(e1 ^ e2);
            // cos of the CA-C-N-spin dihedral, which is all the sign needs.
            const Vector3f axis = unitOf(N - C);
            const Vector3f a = position[iCA1] - C;
            const Vector3f b = placed - N;
            const Vector3f pa = a - axis * a.dot(axis);
            const Vector3f pb = b - axis * b.dot(axis);
            const float cosine = pa.norm() > 1e-6f && pb.norm() > 1e-6f ? pa.dot(pb) / (pa.norm() * pb.norm()) : 0.0f;
            if (sign == 0 || cosine > best)
            {
                best = cosine;
                E2 = candidate;
            }
        }
    }
    const Vector3f e3 = e1 ^ e2;
    const Vector3f E3 = E1 ^ E2;

    reference.assign(position.begin(), position.end());
    reference[iN] = N;
    reference[iCA2] = CA2;
    for (size_t k = 0; k < refname.size(); ++k)
    {
        if (refname[k].empty() || refname[k][0] != '+' || static_cast<int>(k) == iN || static_cast<int>(k) == iCA2)
            continue;
        const Vector3f v = position[k] - position[iN];
        reference[k] = N + E1 * v.dot(e1) + E2 * v.dot(e2) + E3 * v.dot(e3);
    }
    return true;
}

unsigned long long pairKey(unsigned a, unsigned b)
{
    return (static_cast<unsigned long long>(a) << 32) | static_cast<unsigned long long>(b);
}

} // namespace

void PeptideBondFormation::setup(SpringNetwork & network, const configuration::PeptideBondSetting & settings)
{
    _enabled = false;
    if (!settings.enable)
        return;

    _settings = settings;

    if (settings.path.empty())
        logging::die("peptidebond.enable is on but peptidebond.path names no .rbody file: the new springs are made "
                     "by the model's own rules, so there is nothing to make them with.");

    const auto fields = utils::string::split(settings.bond, ":");
    if (fields.size() != 2 || fields[0].empty() || fields[1].empty())
        logging::die("peptidebond.bond = '%s': expected '<electrophile>:<nucleophile>', e.g. 'C:N'.",
                     settings.bond.c_str());
    _electrophilename = fields[0];
    _nucleophilename = fields[1];
    _suffix = "_" + settings.group;

    _hastranslation = false;
    if (!settings.naming.empty())
    {
        reduce::ReduceRuleReader naming(settings.naming);
        naming.read();
        _translation = naming.rules();
        _hastranslation = true;
    }

    rigidbodygroup::RigidBodyRuleReader reader(settings.path);
    reader.read();
    _rules = reader.rules();

    // Group the particles by residue, in topology order. The same grouping
    // RigidBodyBuilder makes at build time, and for the same reason: a rule
    // names atoms within a residue, so a residue has to be a set of indices
    // before any name can be resolved.
    _residues.clear();
    _residueof.assign(network.getNumberOfParticles(), 0);
    for (unsigned i = 0; i < network.getNumberOfParticles(); ++i)
    {
        const Particle & p = network.getParticle(i);
        const bool same = !_residues.empty() &&
                          network.getParticle(_residues.back().back()).getResId() == p.getResId() &&
                          network.getParticle(_residues.back().back()).getChainName() == p.getChainName();
        if (!same)
            _residues.emplace_back();
        _residues.back().push_back(i);
        _residueof[i] = static_cast<unsigned>(_residues.size() - 1);
    }

    // The candidates: an atom of the right name, in a residue the rules cover,
    // whose terminus is not already used.
    _electrophiles.clear();
    _nucleophiles.clear();
    for (unsigned i = 0; i < network.getNumberOfParticles(); ++i)
    {
        if (_isNamed(network, i, _electrophilename) && _isFree(network, i, _nucleophilename))
            _electrophiles.push_back(i);
        if (_isNamed(network, i, _nucleophilename) && _isFree(network, i, _electrophilename))
            _nucleophiles.push_back(i);
    }

    logging::info("Peptide bond formation: %zu free %s and %zu free %s, rules '%s' group '%s'.", _electrophiles.size(),
                  _electrophilename.c_str(), _nucleophiles.size(), _nucleophilename.c_str(), settings.path.c_str(),
                  settings.group.c_str());
    if (_electrophiles.empty() || _nucleophiles.empty())
        logging::warning("Peptide bond formation: one of the two candidate sets is empty, so no bond can ever form. "
                      "Check peptidebond.bond against the atom names this network actually carries.");

    // Without a site the angular factor is 1 and the only condition left is the
    // distance -- which would bond a nucleophile arriving edge-on, through the
    // sp2 plane, where it cannot reach the pi* at all. Said out loud because a
    // run configured that way does not fail: it bonds MORE readily.
    unsigned sited = 0;
    for (unsigned i : _electrophiles)
        if (network.getParticle(i).antecedentIndex() >= 0 && network.getParticle(i).hasLobes())
            ++sited;
    if (sited < _electrophiles.size())
        logging::warning("Peptide bond formation: %zu of %zu %s atoms carry no out-of-plane lobe site, so for those "
                         "the Burgi-Dunitz condition cannot be evaluated and peptidebond.weight is inert. Declare "
                         "them in the .hbond table as '<res> %s 0 1 O ^CA:75'.",
                         _electrophiles.size() - sited, _electrophiles.size(), _electrophilename.c_str(),
                         _electrophilename.c_str());

    _leaving.clear();
    for (const std::string & name : utils::string::split(settings.leaving, " "))
        if (!name.empty())
            _leaving.push_back(name);

    _learnTorsionPattern(network);

    _enabled = true;
}

// WHY OMEGA NEEDS ITS OWN TERM AT ALL, since the fifteen springs already fix
// the plane: they do not fix omega, and the arithmetic says so rather than the
// measurement. Rotating the amine side about the new bond changes the fifteen
// distances only at SECOND order -- CA...CA, the most sensitive of them, moves
// 0.016 A over fifteen degrees -- so the whole mesh resists a fifteen-degree
// twist with 0.237 kJ/mol, while one residue deforming by the 0.015 A measured
// under an interactive pull costs 0.073 kJ/mol on its own. Omega is free, and
// the product settled 14 degrees off trans however good the rest lengths were.
// Same evenness as a ring's out-of-plane coordinate; see the ring planarity
// note. A torsion is the only term with any stiffness there.
//
// AMBER puts FOUR torsions on a peptide bond, not one -- {CA, O} x {+CA, +H} --
// with two different tables, and the table INDEX is this network's own. So the
// pattern is LEARNED from a peptide bond the structure already has rather than
// written down here: whatever pdb2spn built, this builds the same.
void PeptideBondFormation::_learnTorsionPattern(const SpringNetwork & network)
{
    _torsionpattern.clear();
    for (const SpringNetwork::Torsion & t : network.getTorsions())
    {
        // The bond is the middle pair, by the .bi.ff's own convention
        // (X - C - +N - Y). Only a torsion straddling a real electrophile /
        // nucleophile couple describes the bond being made.
        if (!_isNamed(network, t.atoms[1], _electrophilename) ||
            !_isNamed(network, t.atoms[2], _nucleophilename))
            continue;
        if (_residueof[t.atoms[1]] == _residueof[t.atoms[2]])
            continue;
        TorsionPattern p;
        // Stored as the reduction renamed them is wrong: the next bond may be
        // between two other residue types, whose names differ. Recovered to the
        // ORIGINAL name, which _isNamed then translates again per residue.
        p.before = _original(network, t.atoms[0]);
        p.after = _original(network, t.atoms[3]);
        p.family = t.family;
        p.table = t.table;
        const bool known = std::any_of(_torsionpattern.begin(), _torsionpattern.end(),
                                       [&](const TorsionPattern & q) {
                                           return q.before == p.before && q.after == p.after && q.table == p.table;
                                       });
        if (!known)
            _torsionpattern.push_back(p);
    }
    if (_torsionpattern.empty())
        logging::warning("Peptide bond formation: this network carries no torsion across an existing %s-%s bond, so "
                         "there is no pattern to copy and a new bond's omega will be FREE -- the fifteen springs "
                         "resist a 15 degree twist with 0.24 kJ/mol. Build the network with -bondedinteraction and "
                         "-dihedralbackbone, on a structure that already contains one such bond.",
                         _electrophilename.c_str(), _nucleophilename.c_str());
    else
        logging::info("Peptide bond formation: %zu torsion(s) will be created per bond, copied from this network's "
                      "own %s-%s torsions.",
                      _torsionpattern.size(), _electrophilename.c_str(), _nucleophilename.c_str());
}

void PeptideBondFormation::_addTorsions(SpringNetwork & network, unsigned electrophile, unsigned nucleophile)
{
    const unsigned own = _residueof[electrophile];
    const unsigned next = _residueof[nucleophile];
    for (const TorsionPattern & p : _torsionpattern)
    {
        const int a = _atomInResidue(network, own, p.before);
        const int d = _atomInResidue(network, next, p.after);
        if (a < 0 || d < 0)
            continue;
        network.addTorsion(p.family, static_cast<unsigned>(a), electrophile, nucleophile,
                           static_cast<unsigned>(d), p.table);
    }
}

std::string PeptideBondFormation::_translate(const std::string & resname, const std::string & name) const
{
    if (!_hastranslation)
        return "";
    for (const auto & rule : _translation.get_rules_for_residue(resname))
        if (rule.hasAtomNamed(name))
            return rule.getName();
    return "";
}

std::string PeptideBondFormation::_original(const SpringNetwork & network, unsigned particle) const
{
    const Particle & p = network.getParticle(particle);
    if (!_hastranslation)
        return p.getName();
    for (const auto & rule : _translation.get_rules_for_residue(p.getResName()))
        if (rule.getName() == p.getName() && rule.getNumberOfAtoms() == 1)
            return *rule.getAtomNames().begin();
    return p.getName();
}

bool PeptideBondFormation::_isNamed(const SpringNetwork & network, unsigned particle, const std::string & name) const
{
    const Particle & p = network.getParticle(particle);
    if (p.getName() == name)
        return true;
    const std::string renamed = _translate(p.getResName(), name);
    return !renamed.empty() && p.getName() == renamed;
}

bool PeptideBondFormation::_isFree(const SpringNetwork & network, unsigned particle, const std::string & other) const
{
    const Particle & p = network.getParticle(particle);
    for (const auto & entry : p.getSpringNeighbors())
    {
        const Particle & q = network.getParticle(entry.first);
        if (!_isNamed(network, entry.first, other))
            continue;
        if (q.getResId() == p.getResId() && q.getChainName() == p.getChainName())
            continue;
        return false;
    }
    return true;
}

int PeptideBondFormation::_atomInResidue(const SpringNetwork & network, unsigned residue,
                                         const std::string & name) const
{
    for (unsigned i : _residues[residue])
        if (_isNamed(network, i, name))
            return static_cast<int>(i);
    return -1;
}

unsigned PeptideBondFormation::update(SpringNetwork & network, unsigned iteration)
{
    if (!_enabled)
        return 0;

    _advanceRamps(network, iteration);

    _closest = -1.0f;
    _bestweight = 0.0f;
    unsigned made = 0;

    // Pairs that reach the threshold on this step, collected before any of them
    // is acted on: forming a bond invalidates both candidate lists, and a loop
    // that mutates what it is walking is the kind of thing that works until two
    // bonds happen to form on the same step.
    std::vector<std::pair<unsigned, unsigned>> ready;

    for (unsigned e : _electrophiles)
    {
        const Particle & pe = network.getParticle(e);
        for (unsigned n : _nucleophiles)
        {
            if (_residueof[e] == _residueof[n] || _bonded.count(pairKey(e, n)))
                continue;
            const Particle & pn = network.getParticle(n);

            const Vector3f axis = pn.getPosition() - pe.getPosition();
            const float distance = axis.norm();
            if (_closest < 0.0f || distance < _closest)
                _closest = distance;
            if (distance < 1e-6f)
            {
                _dwell.erase(pairKey(e, n));
                continue;
            }
            const Vector3f vhat = axis / distance;

            // The SAME angular factor the hydrogen bond term scores with, from
            // the same site builder, at both ends. "The configuration is right"
            // therefore means exactly "the term that brought them here likes
            // this pair", rather than a second opinion that could disagree with
            // the first.
            //
            // Computed BEFORE the distance is judged, and for every pair, so
            // that the two reported numbers answer different questions. Gated
            // on the distance first, they both collapsed together: an approach
            // still six A out reported an angular weight of 0, which reads as
            // "the angle is wrong" when it means "no pair was close enough for
            // anyone to have looked at its angle".
            const SpringNetwork::HydrogenBondSite se =
                network.hydrogenBondSite(pe, pe.antecedentIndex(), pe.antecedentIndex2(), pe.lobeMode(), pe.lobeCos(),
                                         pe.lobeSin(), vhat);
            const SpringNetwork::HydrogenBondSite sn =
                network.hydrogenBondSite(pn, pn.antecedentIndex(), pn.antecedentIndex2(), pn.lobeMode(), pn.lobeCos(),
                                         pn.lobeSin(), -vhat);
            float weight = 1.0f;
            if (se.valid)
                weight *= forcefield::hydrogen_bond_angular_factor(se.hhat.dot(vhat));
            if (sn.valid)
                weight *= forcefield::hydrogen_bond_angular_factor(sn.hhat.dot(-vhat));
            _bestweight = std::max(_bestweight, weight);

            if (distance > _settings.distance || weight < _settings.weight)
            {
                _dwell.erase(pairKey(e, n));
                continue;
            }

            const unsigned held = ++_dwell[pairKey(e, n)];
            if (held >= _settings.dwell)
                ready.emplace_back(e, n);
        }
    }

    for (const auto & pair : ready)
    {
        if (_bonded.count(pairKey(pair.first, pair.second)))
            continue;
        // A second bond on the same terminus would be a pentavalent nitrogen.
        // Both ends are re-checked here because `ready` was collected before
        // any of it was acted on.
        if (!_isFree(network, pair.first, _nucleophilename) || !_isFree(network, pair.second, _electrophilename))
            continue;
        _form(network, pair.first, pair.second, iteration);
        ++made;
    }

    if (made > 0)
    {
        const auto used = [&](unsigned i, const std::string & other)
        { return !_isFree(network, i, other); };
        _electrophiles.erase(std::remove_if(_electrophiles.begin(), _electrophiles.end(),
                                            [&](unsigned i) { return used(i, _nucleophilename); }),
                             _electrophiles.end());
        _nucleophiles.erase(std::remove_if(_nucleophiles.begin(), _nucleophiles.end(),
                                           [&](unsigned i) { return used(i, _electrophilename); }),
                            _nucleophiles.end());
        _dwell.clear();
    }

    return made;
}

void PeptideBondFormation::_form(SpringNetwork & network, unsigned electrophile, unsigned nucleophile,
                                 unsigned iteration)
{
    const std::string & resname = network.getParticle(electrophile).getResName();

    // BY VALUE. get_rules_for_residue builds and returns a fresh container, so a
    // pointer into it dangles the moment the loop that produced it ends -- which
    // read as an empty atom set and made a bond with no springs at all.
    rigidbodygroup::RigidBodyRule chosen;
    bool found = false;
    for (const rigidbodygroup::RigidBodyRule & rule : _rules.get_rules_for_residue(resname))
    {
        const std::string & name = rule.getName();
        if (name.size() >= _suffix.size() && name.compare(name.size() - _suffix.size(), _suffix.size(), _suffix) == 0)
        {
            chosen = rule;
            found = true;
            break;
        }
    }
    if (!found)
    {
        logging::warning("Peptide bond formation: residue '%s' has no '*%s' rule, so the bond it just made carries no "
                      "springs. Nothing is bonded.",
                      resname.c_str(), _suffix.c_str());
        return;
    }

    // Which reference atom each resolved particle is, so the ramp knows where to
    // draw it to. Built alongside the resolution rather than from the particle
    // names, because the "+" prefix is what distinguishes the two CA.
    const unsigned own = _residueof[electrophile];
    const unsigned next = _residueof[nucleophile];
    std::vector<unsigned> resolved;
    std::vector<std::string> refname;
    for (const std::string & atomname : chosen.getAtomNames())
    {
        if (!atomname.empty() && atomname[0] == '-')
            continue;
        const bool plus = !atomname.empty() && atomname[0] == '+';
        const int found = _atomInResidue(network, plus ? next : own, plus ? atomname.substr(1) : atomname);
        if (found < 0)
            continue;
        resolved.push_back(static_cast<unsigned>(found));
        refname.push_back(atomname);
    }

    Bond bond;
    bond.electrophile = electrophile;
    bond.nucleophile = nucleophile;
    bond.step = iteration;

    // Where the group should END UP, built from these very residues rather
    // than read off an idealised plane -- see referenceGeometry for why that
    // distinction is the difference between a trans peptide and a strained one.
    std::vector<Vector3f> here;
    here.reserve(resolved.size());
    for (unsigned index : resolved)
        here.push_back(network.getParticle(index).getPosition());
    std::vector<Vector3f> reference;
    const bool ramped = referenceGeometry(here, refname, _electrophilename, _nucleophilename, reference);
    if (!ramped)
        logging::warning("Peptide bond formation: rule '%s' does not name the atoms a peptide plane is built from "
                         "(CA, %s, +%s, +CA), so its springs keep the length they were born at.",
                         chosen.getName().c_str(), _electrophilename.c_str(), _nucleophilename.c_str());

    for (size_t i = 0; i < resolved.size(); ++i)
    {
        for (size_t j = i + 1; j < resolved.size(); ++j)
        {
            if (network.getParticle(resolved[i]).isInSpringNeighbors(resolved[j]))
                continue;
            // Born at the length it already has, so the bond costs nothing at
            // the instant it is made. The chemistry is in where it is then
            // DRAWN to, not in a step change of the energy -- which at 650
            // kJ/mol/A^2 over the 1.6 A this has to travel would be 830 kJ/mol
            // appearing inside one timestep.
            const float length =
                (network.getParticle(resolved[i]).getPosition() - network.getParticle(resolved[j]).getPosition())
                    .norm();
            const unsigned before = network.getNumberOfSprings();
            network.addSpring(resolved[i], resolved[j], length, static_cast<float>(_settings.stiffness));
            if (network.getNumberOfSprings() == before)
                continue; // addSpring declined it
            bond.springs.push_back(network.getNumberOfSprings() - 1);
            bond.born.push_back(length);
            bond.target.push_back(ramped ? (reference[i] - reference[j]).norm() : length);
        }
    }

    logging::info("Peptide bond formed at step %u between %s %u%s and %s %u%s: %zu springs.", iteration,
                  network.getParticle(electrophile).getResName().c_str(),
                  network.getParticle(electrophile).getResId(),
                  network.getParticle(electrophile).getChainName().c_str(),
                  network.getParticle(nucleophile).getResName().c_str(),
                  network.getParticle(nucleophile).getResId(),
                  network.getParticle(nucleophile).getChainName().c_str(), bond.springs.size());

    if (bond.springs.empty())
        logging::warning("Peptide bond formation: rule '%s' created no spring for this pair -- every atom it names "
                         "was already sprung, or none of them was found. The bond is recorded and will not be "
                         "attempted again, but nothing holds it.",
                         chosen.getName().c_str());

    _addTorsions(network, electrophile, nucleophile);
    const unsigned going = _markLeavingGroup(network, electrophile, bond);
    if (going > 0)
        logging::info("Peptide bond formation: %u spring(s) of the leaving group will be drawn out and broken.",
                      going);

    _bonded.insert(pairKey(electrophile, nucleophile));
    _bonds.push_back(bond);
}

// A real carboxyl does not simply gain a fourth partner: it loses its hydroxyl,
// and a tRNA-esterified carbon loses the ester oxygen. Modelled as the only
// thing that distinguishes a bonded atom from a free one here -- its springs.
// Every spring between the leaving atom and its OWN residue is broken, which
// both stops it being held and stops it being excluded from the non-bonded
// terms, so what pushes it away afterwards is ordinary sterics.
//
// It does not become water: that would need a hydrogen this model has no
// mechanism to add. It leaves as the atom it was, with its own mass and radius.
unsigned PeptideBondFormation::_markLeavingGroup(SpringNetwork & network, unsigned electrophile, Bond & bond)
{
    if (_leaving.empty())
        return 0;
    const unsigned residue = _residueof[electrophile];
    for (const std::string & name : _leaving)
    {
        const int atom = _atomInResidue(network, residue, name);
        if (atom < 0)
            continue;
        for (const auto & entry : network.getParticle(static_cast<unsigned>(atom)).getSpringNeighbors())
        {
            if (entry.second == nullptr)
                continue;
            const unsigned id = entry.second->getId();
            const Particle & a = network.getParticle(static_cast<unsigned>(atom));
            const Particle & b = network.getParticle(entry.first);
            const float now = (a.getPosition() - b.getPosition()).norm();
            // Far enough out that the pair is no longer in each other's way:
            // the steric contact distance, which for AMBER is the sum of the
            // two R*. Breaking a spring already at that length costs nothing.
            const float apart = a.getRadius() + b.getRadius();
            bond.leaving.push_back(id);
            bond.leavingborn.push_back(now);
            bond.leavingtarget.push_back(std::max(apart, now * 1.1f));
        }
    }
    return static_cast<unsigned>(bond.leaving.size());
}

// A LEAVING GROUP CANNOT TELEPORT, and that is what broke the first version.
// Released the instant the bond formed, it sat 1.23 A from the carbon it had
// just left with nothing holding it -- and a broken spring is also a dropped
// exclusion, so it met the full steric wall of its own carbon in one step. On a
// concentrated soup the temperature went from 325 K to 3.8e+18 K on the step a
// bond formed.
//
// So it is drawn OUT over the same ramp that draws the new bond in, to the
// steric contact distance where the two no longer overlap, and only then let
// go. Bond made and bond broken take the same 2000 steps, which is what makes
// either of them integrable.
void PeptideBondFormation::_advanceRamps(SpringNetwork & network, unsigned iteration)
{
    for (Bond & bond : _bonds)
    {
        if (iteration <= bond.step)
            continue;
        if (iteration > bond.step + _settings.ramp)
        {
            if (!bond.released)
            {
                for (unsigned id : bond.leaving)
                    network.releaseSpring(id);
                bond.released = true;
                if (!bond.leaving.empty())
                    logging::info("Peptide bond formation: leaving group let go at step %u, %zu spring(s) broken.",
                                  iteration, bond.leaving.size());
            }
            continue;
        }
        if (_settings.ramp == 0)
            continue;
        const float t = static_cast<float>(iteration - bond.step) / static_cast<float>(_settings.ramp);
        for (size_t k = 0; k < bond.springs.size(); ++k)
            network.getSpring(bond.springs[k]).setEquilibrium(bond.born[k] + t * (bond.target[k] - bond.born[k]));
        for (size_t k = 0; k < bond.leaving.size(); ++k)
            network.getSpring(bond.leaving[k])
                .setEquilibrium(bond.leavingborn[k] + t * (bond.leavingtarget[k] - bond.leavingborn[k]));
    }
}

} // namespace spn
} // namespace biospring
