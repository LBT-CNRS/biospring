#ifndef __BIOSPRING_LOGGING_H__
#define __BIOSPRING_LOGGING_H__

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

} // namespace biospring

#endif // __BIOSPRING_LOGGING_H__
