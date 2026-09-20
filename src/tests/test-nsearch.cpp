#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "Particle.h"
#include "measure.hpp"
#include "nsearch.hpp"

using Particle = biospring::spn::Particle;

bool is_position_duplicate(const Particle & candidate, const std::vector<Particle> & system);
bool has_position_duplicate(const std::vector<Particle> & system);
std::array<double, 3> random_position();
std::vector<Particle> generate_random_particles(size_t n);
std::vector<Particle> generate_particles(size_t n);
std::vector<Particle> generate_particle_group(size_t n, const std::array<double, 3> & position);
std::vector<Particle> generate_particle_groups(size_t n);

// =====================================================================================
//
// Test for `NeighborSearch` class.
//
// =====================================================================================

TEST(TestNeighborSearch, NeighbourSearchingOneNeighbor)
{
    // 10 groups of 10 particles spaced 10 units apart.
    const auto particles = generate_particle_groups(10);
    double cutoff = 1.0;
    biospring::nsearch::NeighborSearch ns(particles, cutoff);

    for (size_t i = 0; i < particles.size() - 1; i++)
    {
        const auto & neighbors = ns.get_neighbors(particles[i]);
        // Groups have 10 particles so each particle should have 9 neighbors.
        EXPECT_EQ(neighbors.size(), 9);
    }
}

TEST(TestNeighborSearch, NeighbourSearchingNoNeighbor)
{
    const auto particles = generate_particles(10);
    biospring::nsearch::NeighborSearch ns(particles, 0.5);

    for (size_t i = 0; i < particles.size(); i++)
    {
        const auto & neighbors = ns.get_neighbors(particles[i]);
        EXPECT_TRUE(neighbors.empty());
    }
}

TEST(TestNeighborSearch, NeighbourSearchingRandom)
{
    const auto particles = generate_random_particles(1000);

    // makes sure there are no duplicates in the system.
    ASSERT_FALSE(has_position_duplicate(particles));

    double cutoff = 10.0;
    biospring::nsearch::NeighborSearch ns(particles, cutoff);

    for (size_t i = 0; i < particles.size() - 1; i++)
    {
        const auto & neighbors = ns.get_neighbors(particles[i]);
        for (size_t j = i + 1; j < particles.size(); j++)
        {
            const auto & p1 = particles[i];
            const auto & p2 = particles[j];
            double distance = biospring::measure::distance(p1, p2);

            bool neighbors_contains_j = std::find(neighbors.begin(), neighbors.end(), j) != neighbors.end();
            if (distance < cutoff)
                EXPECT_TRUE(neighbors_contains_j);
            else
                EXPECT_FALSE(neighbors_contains_j);
        }
    }
}

// =====================================================================================
//
// Tests for the cell width, which is not the cutoff.
//
// =====================================================================================

// The answer must not depend on how the space was cut up. Held against brute
// force -- which is what NeighborSearchO2 is for -- at widths from the whole
// cutoff down to a fifth of it, so every stencil radius from 1 to 5 is walked.
//
// This is the test the change needs: a stencil that is too small silently drops
// pairs, and nothing else in the suite would notice.
TEST(TestNeighborSearchCellWidth, EveryWidthFindsTheSameNeighbors)
{
    const auto particles = generate_random_particles(500);
    ASSERT_FALSE(has_position_duplicate(particles));

    const double cutoff = 10.0;
    biospring::nsearch::NeighborSearchO2 reference(particles, cutoff);

    for (int divisor = 1; divisor <= 5; divisor++)
    {
        const float width = static_cast<float>(cutoff) / static_cast<float>(divisor);
        biospring::nsearch::NeighborSearch ns(particles, cutoff, 0.0f, width);

        EXPECT_EQ(ns.stencil_radius(), divisor) << "at a width of " << width;

        for (size_t i = 0; i < particles.size(); i++)
        {
            auto expected = reference.get_neighbors(particles[i]);
            auto found = ns.get_neighbors(particles[i]);
            std::sort(expected.begin(), expected.end());
            std::sort(found.begin(), found.end());
            ASSERT_EQ(found, expected) << "particle " << i << " at a cell width of " << width;
        }
    }
}

// The skin widens the sphere the stencil has to cover, and the stencil radius
// has to grow with it. Getting this wrong loses exactly the pairs the skin was
// there to keep.
TEST(TestNeighborSearchCellWidth, TheStencilCoversTheSkinToo)
{
    const auto particles = generate_random_particles(500);
    const double cutoff = 8.0;
    const float skin = 4.0f;

    // Cells of 2 A: the cutoff alone reaches 4 cells, cutoff plus skin reaches 6.
    biospring::nsearch::NeighborSearch ns(particles, cutoff, skin, 2.0f);
    EXPECT_EQ(ns.stencil_radius(), 6);

    biospring::nsearch::NeighborSearchO2 reference(particles, cutoff);
    for (size_t i = 0; i < particles.size(); i++)
    {
        auto expected = reference.get_neighbors(particles[i]);
        auto found = ns.get_neighbors(particles[i]);
        std::sort(expected.begin(), expected.end());
        std::sort(found.begin(), found.end());
        ASSERT_EQ(found, expected) << "particle " << i;
    }
}

// A cell that the cutoff cannot reach into must be skipped, and one it can must
// not be. Checked against the definition directly, because it is the one piece
// of arithmetic the OpenCL kernels also compile.
TEST(TestNeighborSearchCellWidth, TheStencilDropsItsDeadCorners)
{
    const float width = 2.0f;
    const float cutoff = 6.0f;

    // Touching cells: no gap at all, whatever the cutoff.
    EXPECT_TRUE(biospring_cell_in_range(1, 1, 1, width, 0.001f));
    // Straight out along one axis: 3 cells away leaves a gap of 2 widths = 4 A.
    EXPECT_TRUE(biospring_cell_in_range(3, 0, 0, width, cutoff * cutoff));
    // The corner at the same radius is 4 * sqrt(3) = 6.93 A away, and is not.
    EXPECT_FALSE(biospring_cell_in_range(3, 3, 3, width, cutoff * cutoff));
    // Which is the whole point: the stencil is a cube and the cutoff is a ball.
    EXPECT_EQ(biospring_stencil_radius(cutoff, width), 3);
}

// =====================================================================================
//
// Tests for the cached pair list.
//
// =====================================================================================

// The list must answer exactly what a fresh search answers -- that is the whole
// claim, and it is what makes the list an optimisation rather than an
// approximation. Held against brute force after every move, while the particles
// drift far enough to force several rebuilds.
TEST(TestNeighborList, MatchesBruteForceWhileTheParticlesDrift)
{
    auto particles = generate_random_particles(400);
    const double cutoff = 10.0;
    // A narrow skin against a wide jitter, so that every round moves the
    // particles past half of it: a list that failed to rebuild would be stale
    // from the first round on, which is what this has to catch.
    const float skin = 1.0f;

    biospring::nsearch::NeighborSearch ns(particles, cutoff, skin);
    ASSERT_TRUE(ns.has_list()) << "a searcher with a skin should hold a list";

    unsigned state = 11u;
    const auto jitter = [&state]() {
        state = state * 1103515245u + 12345u;
        return (static_cast<float>((state >> 16) & 0x7fffu) / static_cast<float>(0x7fff) - 0.5f) * 1.6f;
    };

    for (int round = 0; round < 10; round++)
    {
        for (auto & p : particles)
            p.setPosition(p.getPosition() + Vector3f(jitter(), jitter(), jitter()));
        ns.update();

        biospring::nsearch::NeighborSearchO2 reference(particles, cutoff);
        for (size_t i = 0; i < particles.size(); i++)
        {
            auto expected = reference.get_neighbors(particles[i]);
            auto found = ns.get_neighbors(i);
            std::sort(expected.begin(), expected.end());
            std::sort(found.begin(), found.end());
            ASSERT_EQ(found, expected) << "round " << round << ", particle " << i;
        }
    }
}

// What the skin is FOR: the list has to stay right while nothing has moved far
// enough to rebuild it. A pair that was outside the cutoff at build time can
// cross into it before the next rebuild, and the list has to have been holding
// it all along -- which is why it is built at cutoff + skin and not at cutoff.
//
// The drift here stays under half the skin, so `update()` deliberately does
// NOT rebuild: what is tested is the list it kept.
TEST(TestNeighborList, HoldsPairsThatCrossInBeforeAnyRebuild)
{
    auto particles = generate_random_particles(400);
    const double cutoff = 25.0;
    const float skin = 4.0f;   // rebuilt only past a drift of 2 A

    biospring::nsearch::NeighborSearch ns(particles, cutoff, skin);
    ASSERT_TRUE(ns.has_list());
    const size_t built = ns.list_size();

    unsigned state = 29u;
    const auto jitter = [&state]() {
        state = state * 1103515245u + 12345u;
        return (static_cast<float>((state >> 16) & 0x7fffu) / static_cast<float>(0x7fff) - 0.5f) * 1.8f;
    };
    for (auto & p : particles)
        p.setPosition(p.getPosition() + Vector3f(jitter(), jitter(), jitter()));

    ns.update();
    ASSERT_EQ(ns.list_size(), built) << "nothing drifted past half the skin, so nothing should have been rebuilt";

    biospring::nsearch::NeighborSearchO2 reference(particles, cutoff);
    size_t pairs = 0;
    for (size_t i = 0; i < particles.size(); i++)
    {
        auto expected = reference.get_neighbors(particles[i]);
        auto found = ns.get_neighbors(i);
        std::sort(expected.begin(), expected.end());
        std::sort(found.begin(), found.end());
        ASSERT_EQ(found, expected) << "particle " << i << " after drifting inside the skin";
        pairs += expected.size();
    }
    EXPECT_GT(pairs, 100u) << "too sparse for the comparison to mean anything";
}

// A list is only worth having if it holds more than this step's neighbours:
// what it holds beyond the cutoff is exactly what lets it survive the next
// step. Without a skin there is nothing to survive on, so there is no list.
TEST(TestNeighborList, NoSkinMeansNoList)
{
    const auto particles = generate_random_particles(200);
    biospring::nsearch::NeighborSearch bare(particles, 10.0);
    EXPECT_FALSE(bare.has_list());

    biospring::nsearch::NeighborSearch skinned(particles, 10.0, 3.0f);
    ASSERT_TRUE(skinned.has_list());

    // And the answers are the same either way.
    for (size_t i = 0; i < particles.size(); i++)
    {
        auto a = bare.get_neighbors(i), b = skinned.get_neighbors(i);
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        ASSERT_EQ(a, b) << "particle " << i;
    }

    // The list carries the pairs the cutoff does not, which is the margin it
    // lives on.
    EXPECT_GT(skinned.list_size(), 0u);
}

// A subset searcher -- electrostatic holds only the charged particles -- must
// give every unlisted particle an empty slice, not a broken one, and must not
// let the offsets go backwards whatever order the indices arrive in.
TEST(TestNeighborList, ASubsetLeavesTheOthersEmpty)
{
    const auto particles = generate_random_particles(300);
    // Deliberately out of order, which is what a caller building an index list
    // from a filter may well produce.
    std::vector<size_t> included;
    for (size_t i = particles.size(); i-- > 0;)
        if (i % 3 == 0)
            included.push_back(i);

    biospring::nsearch::NeighborSearch ns(particles, 12.0, included, 3.0f);
    ASSERT_TRUE(ns.has_list());

    biospring::nsearch::NeighborSearchO2 reference(particles, 12.0);
    for (size_t i = 0; i < particles.size(); i++)
    {
        auto found = ns.get_neighbors(i);
        std::sort(found.begin(), found.end());

        // Only the included particles are in the grid, so only they can be
        // neighbours -- of anyone.
        auto expected = reference.get_neighbors(particles[i]);
        expected.erase(std::remove_if(expected.begin(), expected.end(),
                                      [&](size_t j) { return j % 3 != 0; }),
                       expected.end());
        if (i % 3 != 0)
            expected.clear(); // not a member: the grid was never asked about it
        std::sort(expected.begin(), expected.end());
        ASSERT_EQ(found, expected) << "particle " << i;
    }
}

// =====================================================================================
//
// Tests for the skin margin of `NeighborSearch` (deferred grid rebuilds).
//
// =====================================================================================

// A particle drifting by less than the skin margin must still be found (or correctly
// excluded) by distance, even though the grid was never rebuilt to reflect the move.
TEST(TestNeighborSearchSkin, NeighborsStayCorrectWithinSkinMargin)
{
    auto particles = generate_random_particles(200);
    ASSERT_FALSE(has_position_duplicate(particles));

    double cutoff = 5.0;
    double skin = 3.0;
    biospring::nsearch::NeighborSearch ns(particles, cutoff, skin);

    // Moves every particle by less than the skin margin, without ever calling update().
    for (auto & p : particles)
    {
        double dx = 2.0 * rand() / RAND_MAX - 1.0;
        double dy = 2.0 * rand() / RAND_MAX - 1.0;
        double dz = 2.0 * rand() / RAND_MAX - 1.0;
        double norm = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (norm < 1e-9)
        {
            dx = 1.0;
            norm = 1.0;
        }
        double magnitude = 0.9 * skin * (double)rand() / RAND_MAX;
        const auto & pos = p.getPosition();
        p.setPosition(Vector3f(pos.getX() + dx / norm * magnitude, pos.getY() + dy / norm * magnitude,
                                pos.getZ() + dz / norm * magnitude));
    }

    // Neighbor lists must still be exact against the (moved) positions.
    for (size_t i = 0; i < particles.size() - 1; i++)
    {
        const auto & neighbors = ns.get_neighbors(particles[i]);
        for (size_t j = i + 1; j < particles.size(); j++)
        {
            double distance = biospring::measure::distance(particles[i], particles[j]);
            bool neighbors_contains_j = std::find(neighbors.begin(), neighbors.end(), j) != neighbors.end();
            if (distance < cutoff)
                EXPECT_TRUE(neighbors_contains_j);
            else
                EXPECT_FALSE(neighbors_contains_j);
        }
    }
}

// A particle drifting past the skin margin can be missed by a stale grid; `update()` must
// detect the excess drift and rebuild to restore correctness.
TEST(TestNeighborSearchSkin, UpdateRebuildsAfterExceedingSkinMargin)
{
    // 10 groups of 10 particles spaced 10 units apart.
    auto particles = generate_particle_groups(10);
    double cutoff = 1.0;
    double skin = 2.0;
    biospring::nsearch::NeighborSearch ns(particles, cutoff, skin);

    ASSERT_EQ(ns.get_neighbors(particles[0]).size(), 9);

    // Moves the first particle of the second group right next to particle 0, well beyond
    // the skin margin, without calling update().
    particles[10].setPosition(Vector3f(0.5, 0.0, 0.0));
    ASSERT_LT(biospring::measure::distance(particles[0], particles[10]), cutoff);

    // The grid still places particle 10 in its old, far-away cell: it is missed.
    EXPECT_EQ(ns.get_neighbors(particles[0]).size(), 9);

    // update() must detect that particle 10 drifted past the skin and rebuild.
    ns.update();
    const auto & neighbors = ns.get_neighbors(particles[0]);
    EXPECT_EQ(neighbors.size(), 10);
    EXPECT_NE(std::find(neighbors.begin(), neighbors.end(), 10u), neighbors.end());
}

// =====================================================================================
//
// Test for `NeighborSearchBase` class.
//
// =====================================================================================

// Excepts a failure when the particle list is empty.
TEST(TestNeighborSearch, InitFailsEmptyParticleList)
{
    const auto particles = generate_particles(0);
    EXPECT_THROW(biospring::nsearch::NeighborSearch(particles, 1.0), std::invalid_argument);
}

// Excepts a failure when the cutoff distance is negative.
TEST(TestNeighborSearch, InitFailsNegativeCutoff)
{
    const auto particles = generate_particles(10);
    EXPECT_THROW(biospring::nsearch::NeighborSearch(particles, -1.0);, std::invalid_argument);
}

// Excepts a failure when the cutoff distance is too close to zero.
TEST(TestNeighborSearch, InitFailsCutoffIsZero)
{
    const auto particles = generate_particles(10);
    EXPECT_THROW(biospring::nsearch::NeighborSearch(particles, 1e-8), std::invalid_argument);
}

// =====================================================================================
//
// Helper functions.
//
// =====================================================================================

bool is_position_duplicate(const Particle & candidate, const std::vector<Particle> & system)
{
    for (size_t i = 0; i < system.size(); i++)
    {
        const auto & existing = system.at(i);
        const auto & candidate_position = candidate.getPosition();
        const auto & existing_position = existing.getPosition();

        // Checks that candidate and existing are not the same particle.
        if (&candidate == &existing)
            continue;

        if (std::abs(existing_position.getX() - candidate_position.getX()) < 1e-6 &&
            std::abs(existing_position.getY() - candidate_position.getY()) < 1e-6 &&
            std::abs(existing_position.getZ() - candidate_position.getZ()) < 1e-6)
        {
            return true; // The position is a duplicate
        }
    }
    return false; // The position is not a duplicate
}

bool has_position_duplicate(const std::vector<Particle> & system)
{
    for (size_t i = 0; i < system.size(); i++)
    {
        if (is_position_duplicate(system.at(i), system))
            return true;
    }
    return false;
}

// Returns a random 3D-position.
// Each coordinate is drawn from a uniform distribution on the interval [-100, 100].
std::array<double, 3> random_position()
{
    std::array<double, 3> pos;
    for (double & x : pos)
        x = 200.0 * (double)rand() / RAND_MAX - 100.0;
    return pos;
}

// Generates a collection of particles with random positions.
std::vector<Particle> generate_random_particles(size_t n)
{
    std::vector<Particle> particles(n);
    for (Particle & p : particles)
        p.setPosition(random_position());
    return particles;
}

// Generates a list of particles of size `n`.
// Particle coordinates are set to (0, 0, 0), (1, 1, 1), ..., (n - 1, n - 1, n - 1).
std::vector<Particle> generate_particles(size_t n)
{
    std::vector<Particle> particles(n);
    for (size_t i = 0; i < n; i++)
        particles[i].setPosition(Vector3f(i, i, i));
    return particles;
}

// Generate a list of particles of size `n` with the same coordinates.
std::vector<Particle> generate_particle_group(size_t n, const std::array<double, 3> & position)
{
    std::vector<Particle> particles(n);
    for (size_t i = 0; i < n; i++)
        particles[i].setPosition(position);
    return particles;
}

// Generates a list of n particle groups, spaced 10 units apart.
std::vector<Particle> generate_particle_groups(size_t n)
{
    std::vector<Particle> particles;
    for (size_t i = 0; i < n; i++)
    {
        std::array<double, 3> position = {10 * double(i), 0.0, 0.0};
        const auto group = generate_particle_group(10, position);
        particles.insert(particles.end(), group.begin(), group.end());
    }
    return particles;
}

// -- Main function  ----------------------------------------------------------
int main(int argc, char * argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

