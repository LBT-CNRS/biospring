#ifndef __IO_LOBESPEC_H__
#define __IO_LOBESPEC_H__

#include <cstddef>
#include <string>

namespace biospring
{
namespace io
{

// Parses the last column of a site declaration: the atom that gives the site
// its SECOND reference, and whether that atom is a second bond or a plane.
//
//   atom         a second bonded neighbour. The direction becomes the opposite
//                of the two antecedents' bisector, which is exact for a planar
//                sp2 site with two heavy neighbours.
//   ~atom        atom is a PLANE reference, and the two directions lie IN that
//                plane at +/- ANGLE. Where an sp2 centre's two hydrogens or two
//                lone pairs sit: an exocyclic amine, a carbonyl oxygen.
//                Defaults to 62 degrees, so a table written before the angle
//                existed reads the same.
//   ^atom        atom is a plane reference and the two directions leave that
//                plane along its normal, at +/- ANGLE. The Burgi-Dunitz
//                approach to a carbonyl CARBON. There is no canonical angle for
//                a face approach, so this form REQUIRES one.
//   ~atom:ANGLE  either form with the angle in degrees.
//
// Returns the bare atom name and writes `mode` (ParticleProperties::HBOND_LOBES_*)
// and `angle` in degrees. `reader` and `line_id` only appear in the diagnostics,
// so the same spelling and the same errors serve every file that carries this
// column.
std::string parse_lobe_spec(const std::string & token, int & mode, float & angle, const char * reader,
                            size_t line_id);

} // namespace io
} // namespace biospring

#endif // __IO_LOBESPEC_H__
