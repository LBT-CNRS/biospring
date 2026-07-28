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

// Hydrogen-bond donor/acceptor role of one named atom in one residue type.
struct DonorAcceptorRole
{
    bool isDonor = false;
    bool isAcceptor = false;
};

// Parses a .hbond file: lines of "<resname> <atomname> <donor> <acceptor>",
// donor/acceptor each 0 or 1 -- see data/reducerules/ProteinDonorAcceptor.hbond
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

    // Sets isDonor()/isAcceptor() on every particle of `spn` whose
    // (resname, name) matches an entry read from the file. Particles with no
    // matching entry (waters, ligands, unknown residues, ...) are silently
    // left untagged rather than treated as an error.
    void tagParticles(spn::SpringNetwork & spn) const;

  protected:
    std::map<std::pair<std::string, std::string>, DonorAcceptorRole> _roles;

    void _parse_line(const std::string & line, size_t line_id);
};

} // namespace io
} // namespace biospring

#endif // __IO_DONORACCEPTORRULEREADER_H__
