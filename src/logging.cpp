

#include "logging.h"

#include <iostream>
#include <stdarg.h>
#include <stdio.h>
#include <vector>

#define STRLEN 1024

using std::cerr;
using std::cout;
using std::endl;

void biospring::logging::debug(const char * fmt, ...)
{
    va_list ap;
    char msg[STRLEN];

    va_start(ap, fmt);
    vsnprintf(msg, STRLEN, fmt, ap);
    va_end(ap);

    cerr << DEBUG_COLOR << DEBUG_PREFIX << msg << RESET_COLOR << endl;
}

void biospring::logging::die(const char * fmt, ...)
{
    va_list ap;
    char msg[STRLEN];

    va_start(ap, fmt);
    vsnprintf(msg, STRLEN, fmt, ap);
    va_end(ap);

    cerr << ERROR_COLOR << ERROR_PREFIX << msg << RESET_COLOR << endl;
    exit(EXIT_FAILURE);
}

void biospring::logging::status(const char * fmt, ...)
{
    va_list ap;
    char msg[STRLEN];

    va_start(ap, fmt);
    vsnprintf(msg, STRLEN, fmt, ap);
    va_end(ap);

    cerr << STATUS_COLOR << STATUS_PREFIX << msg << RESET_COLOR << endl;
}

void biospring::logging::warning(const char * fmt, ...)
{
    va_list ap;
    char msg[STRLEN];

    va_start(ap, fmt);
    vsnprintf(msg, STRLEN, fmt, ap);
    va_end(ap);

    cerr << WARNING_COLOR << WARNING_PREFIX << msg << RESET_COLOR << endl;
}

void biospring::logging::info(const char * fmt, ...)
{
    va_list ap;
    char msg[STRLEN];

    va_start(ap, fmt);
    vsnprintf(msg, STRLEN, fmt, ap);
    va_end(ap);

    cerr << INFO_COLOR << INFO_PREFIX << msg << RESET_COLOR << endl;
}

void biospring::logging::error(const char * fmt, ...)
{
    va_list ap;
    char msg[STRLEN];

    va_start(ap, fmt);
    vsnprintf(msg, STRLEN, fmt, ap);
    va_end(ap);

    cerr << ERROR_COLOR << ERROR_PREFIX << msg << RESET_COLOR << endl;
}
