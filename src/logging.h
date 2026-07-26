#ifndef __BIOSPRING_LOGGING_H__
#define __BIOSPRING_LOGGING_H__

#include <atomic>
#include <cstdlib>
#include <string>

#define asize(a) (sizeof(a) / sizeof((a)[0]))

namespace biospring
{

namespace logging
{

static const std::string DEBUG_COLOR = "\033[1;36m";
static const std::string ERROR_COLOR = "\033[1;31m";
static const std::string WARNING_COLOR = "\033[93m";
static const std::string STATUS_COLOR = "\033[95m";
static const std::string INFO_COLOR = "\033[95m";
static const std::string RESET_COLOR = "\033[0m";

static const std::string DEBUG_PREFIX = "!! DEBUG: ";
static const std::string ERROR_PREFIX = "!! ERROR: ";
static const std::string WARNING_PREFIX = "!! WARNING: ";
static const std::string INFO_PREFIX = "-- ";
static const std::string STATUS_PREFIX = "";

void debug(const char * fmt, ...);
void die(const char * fmt, ...);
void error(const char * fmt, ...);
void info(const char * fmt, ...);
void status(const char * fmt, ...);
void warning(const char * fmt, ...);

} // namespace logging

// Emits a warning the first time this call site is reached, then stays quiet.
// For degraded-but-recoverable paths inside the simulation loop, where an
// unconditional warning() would print once per particle per step and drown the
// output it is meant to draw attention to. The flag is function-local, so each
// call site throttles independently, and atomic because several of these sites
// straddle an interactor's worker thread and the main thread.
#define BIOSPRING_WARN_ONCE(...)                                                                                       \
    do                                                                                                                 \
    {                                                                                                                  \
        static std::atomic<bool> _biospring_warned{false};                                                             \
        if (!_biospring_warned.exchange(true, std::memory_order_relaxed))                                              \
            biospring::logging::warning(__VA_ARGS__);                                                                  \
    } while (false)

} // namespace biospring

#endif // __BIOSPRING_LOGGING_H__
