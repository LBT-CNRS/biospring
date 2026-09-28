#ifndef __SPN_PEPTIDEBOND_H__
#define __SPN_PEPTIDEBOND_H__

#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "configuration/Configuration.hpp"
#include "rigidbodygroup/RigidBodyRuleContainer.hpp"

namespace biospring
{
namespace spn
{

class SpringNetwork;

// Makes a covalent bond during a run, by rewriting the spring topology.
//
// WHY THIS IS NOT A POTENTIAL. The hydrogen bond term can already bring a
// nucleophile onto a carbonyl carbon at the Burgi-Dunitz angle -- that is what
// an out-of-plane lobe site is for, and it is the whole of stage one. What no
// well can express is the rest: after the bond is made the two residues are ONE
// molecule, with a peptide plane that resists rotation, and they do not come
// apart when the well is climbed. A deeper well is still a well. So stage two
// does the only thing that actually changes the molecule -- it changes the
// springs -- and does it once, for good.
//
// WHERE THE NEW SPRINGS COME FROM. Not from this file. They come from the same
// .rbody rule that would have created them had the bond been there from the
// start: X_PSI is {CA, C, O, +N, +H, +CA}, sprung pairwise, which IS the peptide
// plane in this model. The product is therefore the model's own peptide and not
// a special case that has to be maintained beside it.
//
// The one thing the rule cannot supply is WHICH residue is the "+" one.
// RigidBodyBuilder answers that from the order residues appear in the file,
// which is exactly right at build time and useless here: two amino acids that
// meet in solution are not neighbours in anybody's file. So the pair is chosen
// geometrically and the rule is resolved against it -- unprefixed names in the
// residue that brings the carbonyl, "+" names in the one that brings the amine.
class PeptideBondFormation
{
  public:
    // One bond that has been made. Kept for the ramp, for the log, and because
    // the pair must never be considered again.
    struct Bond
    {
        unsigned electrophile = 0; // particle index, the carbonyl carbon
        unsigned nucleophile = 0;  // particle index, the attacking nitrogen
        unsigned step = 0;         // the iteration it formed on
        // The springs this bond created, with where each one was born and where
        // it is being drawn to. Springs that already existed are not listed:
        // they belong to whatever made them.
        std::vector<unsigned> springs;
        std::vector<float> born;   // A
        std::vector<float> target; // A
    };

    // Reads the rules and resolves the two atom names. Does nothing, and stays
    // disabled, when the setting is off.
    void setup(SpringNetwork & network, const configuration::PeptideBondSetting & settings);

    bool isEnabled() const { return _enabled; }

    // One step: advances any ramp in progress, then looks for a pair that has
    // been in the attack conformation long enough and bonds it. Returns how
    // many bonds were made on this step.
    unsigned update(SpringNetwork & network, unsigned iteration);

    unsigned getNumberOfBonds() const { return static_cast<unsigned>(_bonds.size()); }
    const std::vector<Bond> & getBonds() const { return _bonds; }

    // How close the closest candidate pair currently is, in A, and the angular
    // factor it carries. Both are what the .msp thresholds are compared against,
    // so a run that never bonds can be asked which of the two conditions it is
    // failing rather than only that it failed.
    float getClosestApproach() const { return _closest; }
    float getBestWeight() const { return _bestweight; }

  protected:
    // The particle of `residue` named `name`, or -1.
    int _atomInResidue(const SpringNetwork & network, unsigned residue, const std::string & name) const;

    // Whether `particle` is already covalently bound to an atom named `other`
    // in a DIFFERENT residue -- which is what "this terminus is already used"
    // means, and what makes the formation irreversible without a flag to keep.
    bool _isFree(const SpringNetwork & network, unsigned particle, const std::string & other) const;

    // The torsions AMBER puts on a peptide bond, learned from one the network
    // already has rather than written down here -- see _learnTorsionPattern.
    struct TorsionPattern
    {
        std::string before; // the atom name on the carbonyl side
        std::string after;  // the atom name on the amine side
        unsigned family = 0;
        unsigned table = 0;
    };
    void _learnTorsionPattern(const SpringNetwork & network);
    void _addTorsions(SpringNetwork & network, unsigned electrophile, unsigned nucleophile);
    // Breaks the leaving group off the electrophile's residue. Returns how many
    // springs it broke.
    unsigned _releaseLeavingGroup(SpringNetwork & network, unsigned electrophile);

    void _form(SpringNetwork & network, unsigned electrophile, unsigned nucleophile, unsigned iteration);
    void _advanceRamps(SpringNetwork & network, unsigned iteration);

    bool _enabled = false;
    configuration::PeptideBondSetting _settings{"peptidebond"};
    rigidbodygroup::RigidBodyRuleContainer _rules;
    std::string _electrophilename; // "C"
    std::string _nucleophilename;  // "N"
    std::string _suffix;           // "_PSI"

    // Particles grouped by residue, in topology order, and the residue each
    // particle belongs to. Both are fixed for the run.
    std::vector<std::vector<unsigned>> _residues;
    std::vector<unsigned> _residueof;

    // The two candidate sets, pruned as bonds form.
    std::vector<unsigned> _electrophiles;
    std::vector<unsigned> _nucleophiles;

    // How many consecutive steps each candidate pair has satisfied both
    // conditions. A pair that fails either one is erased rather than decremented:
    // an approach that is interrupted is not a slower approach.
    std::unordered_map<unsigned long long, unsigned> _dwell;

    // Every pair that has bonded, whatever came of it. The candidate lists are
    // pruned by _isFree, which reads the springs -- so a formation that created
    // none would be retried for ever. Irreversibility is recorded, not inferred.
    std::set<unsigned long long> _bonded;

    std::vector<TorsionPattern> _torsionpattern;
    std::vector<std::string> _leaving;

    std::vector<Bond> _bonds;
    float _closest = 0.0f;
    float _bestweight = 0.0f;
};

} // namespace spn
} // namespace biospring

#endif
