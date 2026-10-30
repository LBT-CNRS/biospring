#ifndef __BIOSPRING_HYBRID36_HPP__
#define __BIOSPRING_HYBRID36_HPP__

// PDB's atom serial is five columns wide and its residue number four, and
// BioSpring writes both with %d. A value that does not fit does not fail: it
// takes the column after it and shifts every field to its right by one
// character, coordinates included. Measured on 079.Ribosome70S_AllAtom, whose
// trajectory holds 279430 atoms and residue numbers from -9 to 91392:
//
//   field               width  atoms over  what breaks
//   serial   %5d            5      179431  coordinates shift, unreadable
//   resId    %4d            4       68261  eats the chain column, chain LOST
//
// And the read side does not notice. PDBReader takes line.substr(38, 8) and
// hands it to std::stof, which parses the leading digits and STOPS rather than
// throwing, so on atom 100000 "6 336.58" becomes 6 and the atom lands 560 A
// from where it belongs with no warning. Writer and reader are both ours, so
// BioSpring could not even read back its own file; every column-based reader
// outside -- VMD, PyMOL, MDAnalysis, mdtraj -- corrupts it the same way.
//
// Hybrid-36 is Grosse-Kunstleve, Sauter, Moriarty and Adams (2006), the
// encoding cctbx and PHENIX use and that recent VMD and PyMOL read. It keeps
// the field WIDTH and spends the characters PDB never used:
//
//   width 5   -9999 .. 99999      decimal, right-justified, UNCHANGED
//             100000 .. 43770015  base 36 on 0-9A-Z, so 100000 -> "A0000"
//             43770016 .. ...     base 36 on 0-9a-z
//   width 4    -999 .. 9999       decimal, UNCHANGED
//             10000 .. 1223055    base 36 on 0-9A-Z, so 10000 -> "A000"
//
// The first range is why this is safe to turn on everywhere: every file that
// fits today is written byte for byte as it was, and 079 is the FIRST network
// over 99999 particles -- the next biggest is 90849 -- so no existing output
// changes at all. 079 needs 279430 and 91392, both inside the upper range.
//
// The decimal branch keeps negative numbers, which the base-36 ranges cannot
// express and which real structures do use: 079 has residue -9.

#include "utils/string.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace biospring
{
namespace io
{
namespace hybrid36
{

static const char * const DIGITS_UPPER = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static const char * const DIGITS_LOWER = "0123456789abcdefghijklmnopqrstuvwxyz";

inline long ipow(const long base, unsigned exponent)
{
    long result = 1;
    while (exponent-- > 0)
        result *= base;
    return result;
}

/// The value written in base 36 over exactly `width` characters.
/// Callers only reach this with a value in [10*36^(width-1), 36^width), which
/// always occupies width digits, so the padding below never actually runs --
/// it is there so a future range change cannot silently produce a short field.
inline std::string _toBase36(const char * digits, long value, const unsigned width)
{
    std::string out;
    while (value > 0)
    {
        out.push_back(digits[value % 36]);
        value /= 36;
    }
    while (out.size() < width)
        out.push_back(digits[0]);
    std::reverse(out.begin(), out.end());
    return out;
}

/// @param width Field width in characters: 5 for an atom serial, 4 for a
///              residue number.
/// @param value The number to write.
/// @return Exactly `width` characters, right-justified when decimal.
/// @throws std::runtime_error when the value is outside everything the field
///         can express. Throwing is the point: the old code emitted a line it
///         could not read back.
inline std::string encode(const unsigned width, const long value)
{
    const long decimal_min = 1 - ipow(10, width - 1);
    const long decimal_max = ipow(10, width);
    if (value >= decimal_min && value < decimal_max)
        return utils::string::format("%*ld", static_cast<int>(width), value);

    const long block = 26 * ipow(36, width - 1);
    const long offset = 10 * ipow(36, width - 1);

    long i = value - decimal_max;
    if (i >= 0 && i < block)
        return _toBase36(DIGITS_UPPER, i + offset, width);
    i -= block;
    if (i >= 0 && i < block)
        return _toBase36(DIGITS_LOWER, i + offset, width);

    throw std::runtime_error(utils::string::format(
        "hybrid36: %ld does not fit a %u-character PDB field", value, width));
}

/// Reads back what encode() wrote, and also every plain decimal field written
/// before this existed -- the two agree below the limit, so an old file and a
/// new one parse identically.
/// @throws std::runtime_error on a field that is neither decimal nor base 36.
inline long decode(const unsigned width, const std::string & field)
{
    const std::string s = utils::string::trim(field);
    if (s.empty())
        return 0;

    bool upper = false;
    bool lower = false;
    for (const char c : s)
    {
        if (c >= 'A' && c <= 'Z')
            upper = true;
        else if (c >= 'a' && c <= 'z')
            lower = true;
    }

    // A decimal field never holds a letter, so this is the whole test. It also
    // means the sign and the leading spaces are handled by the decimal branch
    // rather than by a special case here.
    if (!upper && !lower)
    {
        // The whole field has to be consumed. std::stol parses as far as it
        // can and returns rather than throwing, which is precisely the
        // behaviour that let a shifted column pass for a number for so long.
        try
        {
            size_t used = 0;
            const long value = std::stol(s, &used);
            if (used == s.size())
                return value;
        }
        catch (const std::exception &)
        {
        }
        throw std::runtime_error(
            utils::string::format("hybrid36: \"%s\" is not a number", s.c_str()));
    }
    if (upper && lower)
        throw std::runtime_error(utils::string::format(
            "hybrid36: \"%s\" mixes the upper and lower base-36 ranges", s.c_str()));

    const char * const digits = upper ? DIGITS_UPPER : DIGITS_LOWER;
    long value = 0;
    for (const char c : s)
    {
        const char * const found = std::find(digits, digits + 36, c);
        if (found == digits + 36)
            throw std::runtime_error(utils::string::format(
                "hybrid36: \"%s\" is not a base-36 field", s.c_str()));
        value = value * 36 + (found - digits);
    }

    value += ipow(10, width) - 10 * ipow(36, width - 1);
    if (lower)
        value += 26 * ipow(36, width - 1);
    return value;
}

} // namespace hybrid36
} // namespace io
} // namespace biospring

#endif // __BIOSPRING_HYBRID36_HPP__
