#include "spn/PeptideBond.h"

#include <algorithm>
#include <cmath>

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

// An ideal TRANS peptide plane, as the six atoms of an X_PSI group, in A. The
// two ends of a bond that has just formed sit ~2.9 A apart -- an attack
// distance, not a bond -- so every spring the rule creates is born far from
// where it belongs, and something has to say where that is.
//
// Placed from Engh & Huber's (1991) internal coordinates: C-N 1.329, C=O 1.231,
// N-H 1.010, CA-C 1.525, N-CA 1.458 A; CA-C-N 116.2, O-C-N 123.0, C-N-CA 121.7,
// C-N-H 119.2 deg; omega 180, all six coplanar. C is at the origin and N on +x;
// the two CA are on opposite sides of that axis, which is what trans means.
// The fifteen distances are computed from these rather than tabulated, so the
// geometry can be checked against its source in one place.
struct IdealAtom
{
    const char * name;
    float x, y;
};
const IdealAtom IDEAL_PEPTIDE_PLANE[] = {
    {"CA", -0.67250f, 1.36872f},  // the carbonyl side's alpha carbon
    {"C", 0.0f, 0.0f},            // the carbonyl carbon: the electrophile
    {"O", -0.67045f, -1.03240f},  //
    {"+N", 1.32900f, 0.0f},       // the amine nitrogen: the nucleophile
    {"+H", 1.82122f, 0.88193f},   //
    {"+CA", 2.09442f, -1.24093f}, // the amine side's alpha carbon
};

// The ideal distance between two of the six, or -1 when either name is not one
// of them -- which is how a group that is not a peptide plane opts out of the
// ramp instead of being drawn towards a geometry that means nothing for it.
float idealDistance(const std::string & a, const std::string & b)
{
    const IdealAtom * pa = nullptr;
    const IdealAtom * pb = nullptr;
    for (const IdealAtom & atom : IDEAL_PEPTIDE_PLANE)
    {
        if (a == atom.name)
            pa = &atom;
        if (b == atom.name)
            pb = &atom;
    }
    if (pa == nullptr || pb == nullptr)
        return -1.0f;
    const float dx = pa->x - pb->x;
    const float dy = pa->y - pb->y;
    return std::sqrt(dx * dx + dy * dy);
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
        const Particle & p = network.getParticle(i);
        if (p.getName() == _electrophilename && _isFree(network, i, _nucleophilename))
            _electrophiles.push_back(i);
        if (p.getName() == _nucleophilename && _isFree(network, i, _electrophilename))
            _nucleophiles.push_back(i);
    }

    logging::info("Peptide bond formation: %zu free %s and %zu free %s, rules '%s' group '%s'.", _electrophiles.size(),
                  _electrophilename.c_str(), _nucleophiles.size(), _nucleophilename.c_str(), settings.path.c_str(),
                  settings.group.c_str());
    if (_electrophiles.empty() || _nucleophiles.empty())
        logging::warning("Peptide bond formation: one of the two candidate sets is empty, so no bond can ever form. "
                      "Check peptidebond.bond against the atom names this network actually carries.");

    _enabled = true;
}

bool PeptideBondFormation::_isFree(const SpringNetwork & network, unsigned particle, const std::string & other) const
{
    const Particle & p = network.getParticle(particle);
    for (const auto & entry : p.getSpringNeighbors())
    {
        const Particle & q = network.getParticle(entry.first);
        if (q.getName() != other)
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
        if (network.getParticle(i).getName() == name)
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
            if (distance > _settings.distance || distance < 1e-6f)
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

            if (weight < _settings.weight)
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
            const float ideal = idealDistance(refname[i], refname[j]);
            bond.target.push_back(ideal > 0.0f ? ideal : length);
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

    _bonded.insert(pairKey(electrophile, nucleophile));
    _bonds.push_back(bond);
}

void PeptideBondFormation::_advanceRamps(SpringNetwork & network, unsigned iteration)
{
    if (_settings.ramp == 0)
        return;
    for (const Bond & bond : _bonds)
    {
        if (iteration <= bond.step || iteration > bond.step + _settings.ramp)
            continue;
        const float t = static_cast<float>(iteration - bond.step) / static_cast<float>(_settings.ramp);
        for (size_t k = 0; k < bond.springs.size(); ++k)
            network.getSpring(bond.springs[k]).setEquilibrium(bond.born[k] + t * (bond.target[k] - bond.born[k]));
    }
}

} // namespace spn
} // namespace biospring
