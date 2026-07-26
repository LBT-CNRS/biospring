
#include "../logging.h"

#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace
{

// logging::warning writes to std::cerr; swap its buffer to capture the output.
class CerrCapture
{
  public:
    CerrCapture() : _saved(std::cerr.rdbuf(_buffer.rdbuf())) {}
    ~CerrCapture() { std::cerr.rdbuf(_saved); }

    std::string str() const { return _buffer.str(); }

    size_t count(const std::string & needle) const
    {
        const std::string haystack = _buffer.str();
        size_t n = 0;
        for (size_t pos = haystack.find(needle); pos != std::string::npos; pos = haystack.find(needle, pos + 1))
            n++;
        return n;
    }

  private:
    std::ostringstream _buffer;
    std::streambuf * _saved;
};

} // namespace

// The point of BIOSPRING_WARN_ONCE: a degraded path inside the simulation loop
// must announce itself, but exactly once, no matter how many times it is hit.
TEST(WarnOnce, FiresOnlyOnceForRepeatedCalls)
{
    CerrCapture capture;

    for (int i = 0; i < 100; ++i)
        BIOSPRING_WARN_ONCE("repeated degraded path");

    EXPECT_EQ(capture.count("repeated degraded path"), 1u);
}

// The flag is function-local, so two distinct sites must not silence each other.
TEST(WarnOnce, EachCallSiteThrottlesIndependently)
{
    CerrCapture capture;

    for (int i = 0; i < 10; ++i)
    {
        BIOSPRING_WARN_ONCE("first site");
        BIOSPRING_WARN_ONCE("second site");
    }

    EXPECT_EQ(capture.count("first site"), 1u);
    EXPECT_EQ(capture.count("second site"), 1u);
}

// Several of these sites straddle an interactor's worker thread and the main
// thread, so "once" must hold under concurrency, not just in a single thread.
TEST(WarnOnce, FiresOnlyOnceAcrossThreads)
{
    CerrCapture capture;

    std::vector<std::thread> threads;
    for (int t = 0; t < 8; ++t)
        threads.emplace_back([]() {
            for (int i = 0; i < 100; ++i)
                BIOSPRING_WARN_ONCE("concurrent degraded path");
        });

    for (std::thread & thread : threads)
        thread.join();

    EXPECT_EQ(capture.count("concurrent degraded path"), 1u);
}

// The message must still go through logging::warning, prefix included.
TEST(WarnOnce, UsesTheWarningChannel)
{
    CerrCapture capture;

    BIOSPRING_WARN_ONCE("formatted %s and %d", "string", 42);

    EXPECT_NE(capture.str().find(biospring::logging::WARNING_PREFIX), std::string::npos);
    EXPECT_NE(capture.str().find("formatted string and 42"), std::string::npos);
}
