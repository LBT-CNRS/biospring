#ifndef __IO_DONORACCEPTORRULEREADER_H__
#define __IO_DONORACCEPTORRULEREADER_H__

#include <map>
#include <vector>
#include <string>
#include <utility>

#include "IO/ReaderBase.h"
#include "reduce/ReduceRuleContainer.hpp"
#include "SpringNetwork.h"

namespace biospring
{
namespace io
{

// Hydrogen-bond role of one named atom in one residue type: how many bonds
// it can hold at once in each direction, and which heavy atom it hangs off.
struct DonorAcceptorRole
{
    unsigned donorCapacity = 0;
    unsigned acceptorCapacity = 0;
    // Name of the neighbouring heavy atom, in the same residue and in the
    // same naming convention as the second column. Empty when the file names
    // none, in which case the bond keeps its old, purely distance-based
    // behaviour for this atom. See ParticleProperties::antecedentIndex.
    std::string antecedent;

    // A SECOND antecedent, optional. One is not enough for a planar sp2
    // donor: a backbone amide nitrogen sits between CA and the previous
    // residue's C, and its hydrogen points opposite their bisector, 58
    // degrees off the CA->N direction a single antecedent gives. Naming both
    // recovers the true N-H direction to under a degree. Same for a guanine
    // N1 or a thymine N3 between two ring carbons.
    std::string antecedent2;

    // Is antecedent2 a PLANE reference rather than a second bond, and if so which
    // way do the lobes go?
    //
    //   ~atom        IN the plane of the three atoms, at 62 degrees by default.
    //                Where an sp2 centre's two hydrogens or two lone pairs sit: an
    //                exocyclic amine, a carbonyl oxygen.
    //   ^atom        OUT of that plane, along its normal. The Burgi-Dunitz approach
    //                to a carbonyl CARBON: a nucleophile attacks one of its two
    //                faces, at 105 degrees of the C=O axis. There is no sensible
    //                default angle for this, so it must be written.
    //   ~atom:ANGLE  either form with the angle in degrees.
    //
    // A bare '~' keeps 62 degrees so every table written before the angle existed
    // reads the same. See ParticleProperties::setLobes.
    int lobeMode = 0;        // ParticleProperties::HBOND_LOBES_*
    float lobeAngle = 62.0f; // degrees

    // The named Morse parameters this site belongs to, written '@name' anywhere
    // after the capacities. Empty means the force field's own -- which is what
    // a hydrogen bond is. See HydrogenBondGroup.
    std::string group;
};

// A named set of Morse parameters, declared in the table as
//
//     GROUP <name> <welldepth kJ.mol-1> <equilibrium A> <width A-1>
//
// and used by a PAIR only when BOTH its ends name the SAME group. That
// condition is the whole point and not a convenience: a peptide bond's well is
// two orders of magnitude deeper than a hydrogen bond's and sits at 1.33 A
// rather than 2.90, so a carbonyl carbon carrying it alone would drag every
// amine and every hydroxyl in reach down to a covalent distance. It is the
// COUPLE that reacts, so it is the couple that carries the parameters.
struct HydrogenBondGroup
{
    std::string name;
    float welldepth = 0.0f;   // kJ.mol-1
    float equilibrium = 0.0f; // A
    float width = 0.0f;       // A-1
};

// <resname> may be '*', which matches ANY residue. The protein backbone is the
// same in every amino acid, so its two lines are said once instead of forty
// times -- and a residue that differs says so explicitly and wins. Two do:
// proline, whose ring nitrogen has no amide hydrogen and is therefore NOT a
// donor, and the N-methyl cap, whose nitrogen hangs off its methyl rather than
// an alpha carbon. The acetyl cap needs no exception: it has no nitrogen for
// the wildcard to find.
//
// Parses a .hbond file: lines of
// "<resname> <atomname> <donor> <acceptor> [<antecedent> [[~|^]<antecedent2>[:angle]]]", donor/acceptor
// being CAPACITIES (0, 1, 2...) rather than flags -- 0/1 keeps its old
// meaning exactly, so a table written before capacities existed still reads
// the same. The fifth column is optional and names the heavy atom that gives
// the bond its direction -- see data/reducerules/ProteinDonorAcceptor.hbond
// (plain PDB atom names) and data/reducerules/amber.hbond (same
// classification, keyed by amber.grp's renamed types instead, for use with
// a --grp amber.grp --ff amber.ff -reduced topology) for the format and the
// classification itself.
//
// Unlike ReduceRuleReader/RigidBodyRuleReader, there is no grouping of
// several atoms under one named rule: every line stands on its own, keyed by
// (resname, atomname), so a plain lookup table is enough.
class DonorAcceptorRuleReader : public ReaderBase
{
  public:
    DonorAcceptorRuleReader() : ReaderBase() {}
    DonorAcceptorRuleReader(const std::string & path) : ReaderBase(path) {}
    DonorAcceptorRuleReader(const char * const path) : ReaderBase(path) {}

    void read();

    // The .grp the topology was REDUCED with, from simulation.naming. With it,
    // a table written in plain PDB names resolves against a renamed topology,
    // so ONE table serves both -- instead of a second copy of the same 60
    // entries keyed by amber.grp's types, which is what this file used to need.
    //
    // It also replaces a heuristic. A cross-residue antecedent like '-C' was
    // resolved by dropping the first character of the type (AC -> C, ACA ->
    // CA), which holds for amber.grp's one-letter prefixes and silently fails
    // elsewhere: the N-formyl cap's carbon, FRC, becomes RC instead of C. With
    // a real table the translation is exact; without one the heuristic is kept,
    // so nothing that worked before stops working.
    void setNaming(const std::string & path);

    // Sets the donor/acceptor capacities and resolves the antecedent on
    // every particle of `spn` whose (resname, name) matches an entry read
    // from the file. Particles with no matching entry (waters, ligands,
    // unknown residues, ...) are silently left untagged rather than treated
    // as an error. An antecedent that names an atom absent from the residue
    // is reported once and leaves that particle undirected, which is a
    // degraded but valid state, not a failure.
    void tagParticles(spn::SpringNetwork & spn) const;

    const std::vector<HydrogenBondGroup> & groups() const { return _groups; }

  protected:
    std::map<std::pair<std::string, std::string>, DonorAcceptorRole> _roles;
    std::vector<HydrogenBondGroup> _groups;
    reduce::ReduceRuleContainer _translation;
    bool _hastranslation = false;

    void _parse_line(const std::string & line, size_t line_id);

    // The particle's name in the convention the TABLE is written in: its own
    // name translated back through the .grp when there is one, and the
    // first-character heuristic when there is not.
    std::string _plainName(const spn::Particle & p) const;

    // The role that applies to a particle: its residue's own entry first, then
    // the '*' wildcard's. A residue that says something specific must win, or
    // the two lines that describe every backbone could not be overridden for
    // proline -- whose nitrogen has no hydrogen to give and must NOT be tagged
    // a donor.
    const DonorAcceptorRole * _roleFor(const spn::Particle & p) const;
};

} // namespace io
} // namespace biospring

#endif // __IO_DONORACCEPTORRULEREADER_H__
