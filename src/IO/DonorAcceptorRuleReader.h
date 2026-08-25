#ifndef __IO_DONORACCEPTORRULEREADER_H__
#define __IO_DONORACCEPTORRULEREADER_H__

#include <map>
#include <string>
#include <utility>

#include "IO/ReaderBase.h"
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
};

// Parses a .hbond file: lines of
// "<resname> <atomname> <donor> <acceptor> [<antecedent> [<antecedent2>]]", donor/acceptor
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

    // Sets the donor/acceptor capacities and resolves the antecedent on
    // every particle of `spn` whose (resname, name) matches an entry read
    // from the file. Particles with no matching entry (waters, ligands,
    // unknown residues, ...) are silently left untagged rather than treated
    // as an error. An antecedent that names an atom absent from the residue
    // is reported once and leaves that particle undirected, which is a
    // degraded but valid state, not a failure.
    void tagParticles(spn::SpringNetwork & spn) const;

  protected:
    std::map<std::pair<std::string, std::string>, DonorAcceptorRole> _roles;

    void _parse_line(const std::string & line, size_t line_id);
};

} // namespace io
} // namespace biospring

#endif // __IO_DONORACCEPTORRULEREADER_H__
