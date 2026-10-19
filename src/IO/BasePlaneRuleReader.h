#ifndef __IO_BASEPLANERULEREADER_H__
#define __IO_BASEPLANERULEREADER_H__

#include <array>
#include <map>
#include <string>

#include "IO/ReaderBase.h"
#include "SpringNetwork.h"

namespace biospring
{
namespace io
{

// Parses a .stack file: lines of
//
//     <resname> <atom0> <atom1> <atom2>
//
// naming three atoms of one nucleobase that span its plane. The stacking term
// takes that base's centre as their centroid and its normal as
// (atom1-atom0) x (atom2-atom0); which three atoms is a data decision, measured
// rather than assumed -- see data/reducerules/*.stack for the measurement.
//
// WHY THIS IS NOT THE .hbond TABLE, for the same reason .psite is not. A .hbond
// line classifies ONE ATOM as a donor or an acceptor and gives it capacities. A
// .stack line is a property of a RESIDUE: three atoms taken together, meaningless
// one at a time. Filing a residue-level declaration in an atom-level table is
// what put the peptide-bond attack sites in the .hbond table with capacities
// '0 0' and let either file silently overwrite the other's direction on an atom
// both named. Each term reads its own file into its own storage.
//
// Like PeptideSiteRuleReader, every line stands on its own, keyed by resname: a
// plain lookup table is enough.
class BasePlaneRuleReader : public ReaderBase
{
  public:
    BasePlaneRuleReader() : ReaderBase() {}
    BasePlaneRuleReader(const std::string & path) : ReaderBase(path) {}
    BasePlaneRuleReader(const char * const path) : ReaderBase(path) {}

    void read();

    const std::map<std::string, std::array<std::string, 3>> & rules() const { return _planes; }

    /// Resolves the table against a network and installs one frame per residue
    /// whose three atoms are all present. A residue the table does not name, or
    /// one missing an atom, simply gets no frame and so takes no part in the
    /// stacking term -- which is how a protein stays untouched by it.
    void buildFrames(spn::SpringNetwork & spn) const;

  protected:
    std::map<std::string, std::array<std::string, 3>> _planes;

    void _parse_line(const std::string & line, size_t line_id);
};

} // namespace io
} // namespace biospring

#endif // __IO_BASEPLANERULEREADER_H__
