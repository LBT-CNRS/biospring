
#include <gtest/gtest.h>

#include <filesystem>
#include <string>
#include <vector>

#include "IO/modern/XTCTrajectoryWriter.hpp"
#include "Particle.h"
#include "SpringNetwork.h"
#include "configuration/Configuration.hpp"

#include "xdrfile.h"
#include "xdrfile_xtc.h"

#include "TempPath.h"

using namespace biospring;

namespace
{

// Angstroms. The writer divides by ten, so a frame read back must carry these
// over ten -- which is also what tells a working writer from one that happens
// to produce a plausible-looking file.
constexpr float X0 = 1.0f;
constexpr float X1 = 24.5f;
constexpr float X2 = -7.25f;

struct TestXTCTrajectoryWriter : public ::testing::Test
{
    configuration::Configuration config;
    spn::SpringNetwork spn;
    std::string path;

    void SetUp() override
    {
        ::testing::Test::SetUp();
        path = biospring::tests::temp_path(".xtc");

        const float xs[3] = {X0, X1, X2};
        for (int i = 0; i < 3; ++i)
        {
            spn::Particle p;
            p.setPosition(Vector3f(xs[i], 0.0f, 0.0f));
            p.setMass(12.0f);
            p.setRadius(1.5f);
            spn.addParticle(p);
        }
        spn.setup(config);
    }

    void TearDown() override
    {
        std::error_code ignored;
        std::filesystem::remove(path, ignored);
        ::testing::Test::TearDown();
    }
};

// This writer used to inherit TrajectoryWriterBase's constructor. That
// constructor calls safe_open(), safe_open() is NOT virtual, so the call
// resolved to WriterBase's -- which opens an std::ofstream -- while
// XTCTrajectoryWriter::safe_open(), the one that opens the xdr handle, merely
// HID the base's and was never called from anywhere. Two consequences, both
// seen on every run that asked for an XTC: the file was created and truncated
// as a text stream nobody wrote to, so it stayed at 0 bytes, and _xdr was left
// uninitialised, so the first frame handed write_xtc whatever was on the stack
// and the process died of SIGSEGV.
//
// Reproduced on 077.Ribosome (5207 particles) and 079's 70S (279430), with and
// without the peptide bond term, which is why this test asks for nothing more
// than three particles: the defect was in construction, not in size.
TEST_F(TestXTCTrajectoryWriter, writes_a_readable_file_rather_than_crashing)
{
    {
        io::modern::XTCTrajectoryWriter writer(path, spn, 1);
        writer.write_step();
        writer.write_step();
    } // closed here, by the destructor

    ASSERT_TRUE(std::filesystem::exists(path));
    EXPECT_GT(std::filesystem::file_size(path), 0u) << "a 0-byte .xtc is the old failure, not a new one";

    int natoms = 0;
    ASSERT_EQ(read_xtc_natoms(const_cast<char *>(path.c_str()), &natoms), exdrOK);
    EXPECT_EQ(natoms, 3);

    XDRFILE * xdr = xdrfile_open(path.c_str(), "r");
    ASSERT_NE(xdr, nullptr);

    // rvec is float[3], so it cannot be a vector element: a flat buffer, cast
    // the way the writer itself casts.
    std::vector<float> frame(3 * static_cast<size_t>(natoms));
    matrix box;
    int step = 0;
    float time = 0.0f, precision = 0.0f;
    int frames = 0;
    while (read_xtc(xdr, natoms, &step, &time, box, reinterpret_cast<rvec *>(frame.data()), &precision) == exdrOK)
    {
        ++frames;
        // Angstroms in, nanometres out.
        EXPECT_NEAR(frame[0], X0 / 10.0f, 1e-3);
        EXPECT_NEAR(frame[3], X1 / 10.0f, 1e-3);
        EXPECT_NEAR(frame[6], X2 / 10.0f, 1e-3);
    }
    xdrfile_close(xdr);

    EXPECT_EQ(frames, 2) << "both write_step() calls must reach the file";
}

} // namespace
