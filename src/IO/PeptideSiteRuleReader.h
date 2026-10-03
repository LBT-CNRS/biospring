#ifndef __IO_PEPTIDESITERULEREADER_H__
#define __IO_PEPTIDESITERULEREADER_H__

#include <map>
#include <string>
#include <utility>

#include "IO/ReaderBase.h"

namespace biospring
{
namespace io
{

// The geometry of one reactive site: an atom of one residue type, the heavy
// atom that gives it a direction, and the plane that splits that direction into
// two faces.
//
// For a peptide bond's ELECTROPHILE this is the Burgi-Dunitz approach to a
// carbonyl carbon: the antecedent is the carbonyl oxygen, so the base direction
// is opposite C=O, and the plane reference is CA, so the two lobes leave the
// (C, O, CA) plane at +/- the angle. 75 degrees puts them 105 degrees off the
// C=O axis, which is where the pi* lobe is -- 1VY4 measures 100.6 degrees.
struct PeptideSite
{
    // The bonded heavy atom the direction hangs off, in the same naming
    // convention as the second column. Required: a site with no antecedent has
    // no direction, and the angular criterion could not be evaluated at all.
    std::string antecedent;
    // The second reference. Either another bonded neighbour, which makes the
    // direction the opposite of the bisector, or a PLANE marked '~' or '^'.
    std::string antecedent2;
    int lobeMode = 0;        // ParticleProperties::HBOND_LOBES_*
    float lobeAngle = 62.0f; // degrees
};

// Parses a .psite file: lines of
//
//     <resname> <atomname> <antecedent> [[~|^]<plane>[:angle]]
//
// WHY THIS IS NOT THE .hbond TABLE. The two files describe different
// chemistry and are consumed by different terms. A .hbond line classifies an
// atom as a hydrogen-bond donor or acceptor and gives it capacities; a .psite
// line says where a REACTION may be attacked from. They happened to share a
// column layout, and the attack sites were written into the .hbond table with
// capacities '0 0' so that a carbonyl carbon would not take a donor's slot --
// which worked, but filed a peptide-bond parameter under the hydrogen bond's
// name and let either file silently overwrite the other's direction on an atom
// both named. Each term now reads its own file into its own storage.
//
// The lobe column is parsed by the same code as the .hbond table's
// (io::parse_lobe_spec), so the two files cannot drift apart on the one thing
// they do share.
//
// Like DonorAcceptorRuleReader and unlike ReduceRuleReader, every line stands
// on its own, keyed by (resname, atomname): a plain lookup table is enough.
class PeptideSiteRuleReader : public ReaderBase
{
  public:
    PeptideSiteRuleReader() : ReaderBase() {}
    PeptideSiteRuleReader(const std::string & path) : ReaderBase(path) {}
    PeptideSiteRuleReader(const char * const path) : ReaderBase(path) {}

    void read();

    const std::map<std::pair<std::string, std::string>, PeptideSite> & rules() const { return _sites; }

  protected:
    std::map<std::pair<std::string, std::string>, PeptideSite> _sites;

    void _parse_line(const std::string & line, size_t line_id);
};

} // namespace io
} // namespace biospring

#endif // __IO_PEPTIDESITERULEREADER_H__
