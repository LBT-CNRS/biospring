#include "IO/hybrid36.hpp"

#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

using biospring::io::hybrid36::decode;
using biospring::io::hybrid36::encode;

// The decimal range is the whole reason this encoding is safe to switch on:
// every file that fit the old %5d and %4d is written byte for byte as before.
TEST(Hybrid36, DecimalRangeIsUnchanged)
{
    EXPECT_EQ(encode(5, 1), "    1");
    EXPECT_EQ(encode(5, 99999), "99999");
    EXPECT_EQ(encode(4, 1), "   1");
    EXPECT_EQ(encode(4, 9999), "9999");
}

// 079.Ribosome70S_AllAtom has residue -9. Negative numbers have no base-36
// form, so they have to stay in the decimal branch rather than raise.
TEST(Hybrid36, NegativeResidueNumbers)
{
    EXPECT_EQ(encode(4, -9), "  -9");
    EXPECT_EQ(decode(4, "  -9"), -9);
    EXPECT_EQ(encode(5, -999), " -999");
    EXPECT_EQ(decode(5, " -999"), -999);
}

// The published boundaries of Grosse-Kunstleve et al. (2006).
TEST(Hybrid36, PublishedBoundaries)
{
    EXPECT_EQ(encode(5, 100000), "A0000");
    EXPECT_EQ(encode(5, 100001), "A0001");
    EXPECT_EQ(encode(5, 43770015), "ZZZZZ");
    EXPECT_EQ(encode(5, 43770016), "a0000");
    EXPECT_EQ(encode(5, 87440031), "zzzzz");

    EXPECT_EQ(encode(4, 10000), "A000");
    EXPECT_EQ(encode(4, 1223055), "ZZZZ");
    EXPECT_EQ(encode(4, 1223056), "a000");
    EXPECT_EQ(encode(4, 2436111), "zzzz");
}

// The two values this was written for, measured on 079's own trajectory.
TEST(Hybrid36, TheValuesThatBroke079)
{
    EXPECT_EQ(encode(5, 279430).size(), 5u);
    EXPECT_EQ(decode(5, encode(5, 279430)), 279430);
    EXPECT_EQ(encode(4, 91392).size(), 4u);
    EXPECT_EQ(decode(4, encode(4, 91392)), 91392);
}

TEST(Hybrid36, RoundTripOverBothRanges)
{
    for (long v = -999; v < 12000; ++v)
        EXPECT_EQ(decode(4, encode(4, v)), v) << "width 4, value " << v;
    for (long v : {1L, 99998L, 99999L, 100000L, 250000L, 279430L, 1000000L, 43770015L, 43770016L, 50000000L})
        EXPECT_EQ(decode(5, encode(5, v)), v) << "width 5, value " << v;
}

// A field is always exactly its width, or the columns after it shift and the
// coordinates become unreadable -- which is the defect this replaces.
TEST(Hybrid36, FieldIsAlwaysItsFullWidth)
{
    for (long v : {-999L, 0L, 1L, 9999L, 10000L, 91392L, 1223055L})
        EXPECT_EQ(encode(4, v).size(), 4u) << "value " << v;
    for (long v : {-9999L, 0L, 1L, 99999L, 100000L, 279430L, 43770015L})
        EXPECT_EQ(encode(5, v).size(), 5u) << "value " << v;
}

// Out of range must throw. The old code emitted a line it could not read back;
// refusing is the improvement, not a wider range.
TEST(Hybrid36, OutOfRangeThrows)
{
    EXPECT_THROW(encode(5, 87440032), std::runtime_error);
    EXPECT_THROW(encode(4, 2436112), std::runtime_error);
    EXPECT_THROW(encode(5, -10000), std::runtime_error);
}

TEST(Hybrid36, MalformedFieldsThrow)
{
    EXPECT_THROW(decode(5, "A0z00"), std::runtime_error);  // mixed ranges
    EXPECT_THROW(decode(5, "12*45"), std::runtime_error);  // not a number
    EXPECT_THROW(decode(4, "A$00"), std::runtime_error);   // not base 36
}

// A blank field is how real PDB files spell "absent", and the old std::stoi
// threw on it. Reading it as 0 keeps those files loading.
TEST(Hybrid36, BlankFieldIsZero)
{
    EXPECT_EQ(decode(5, "     "), 0);
    EXPECT_EQ(decode(4, "    "), 0);
}
