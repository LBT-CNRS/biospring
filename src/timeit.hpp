#ifndef __TIMEIT_HPP__
#define __TIMEIT_HPP__

#include <chrono>
#include <iomanip> // std::setprecision
#include <iostream>
#include <map>
#include <ostream>
#include <string>

namespace biospring
{
namespace timeit
{

struct Timer
{
    using high_resolution_clock = std::chrono::high_resolution_clock;

    high_resolution_clock::time_point time_start;
    high_resolution_clock::time_point time_stop;

    Timer() { start(); }

    inline void start();
    inline void stop();
    inline void reset();

    // Returns the elapsed time between start and stop in milliseconds.
    inline unsigned long elapsed() const;
    inline double elapsed_seconds() const;
    inline static unsigned long elapsed(const high_resolution_clock::time_point & start,
                                        const high_resolution_clock::time_point & stop);

    inline friend std::ostream & operator<<(std::ostream & os, const Timer & t)
    {
        os << t.elapsed() << " ms";
        return os;
    }

    inline static auto now() { return high_resolution_clock::now(); }
};

// A Profiler holds several timers which allow to profile several portion of the code
// at the time.
class Profiler
{
  protected:
    std::map<std::string, Timer> _counters;

  public:
    Profiler() = default;
    virtual ~Profiler() {}

    inline Timer & operator[](const std::string & name);
    inline const Timer & operator[](const std::string & name) const;

    inline void create_timer(const std::string & name);
    inline void erase_timer(const std::string & name);

    inline void print_elapsed(const std::string & name, std::ostream & os = std::cerr) const;

  protected:
    inline void _assert_timer_exists(const std::string & name) const;
    inline void _assert_timer_does_not_exists(const std::string & name) const;
    inline bool _timer_exists(const std::string & name) const;

  private:
    Profiler(const Profiler &) = delete;
    Profiler & operator=(const Profiler &) = delete;
};

// ======================================================================================
//
// Timer
//
// ======================================================================================

void Timer::start() { time_start = now(); }
void Timer::stop() { time_stop = now(); }
void Timer::reset()
{
    stop();
    start();
}

// Returns the elapsed time between start and stop in milliseconds.
unsigned long Timer::elapsed() const { return elapsed(time_start, now()); }

// Returns the elapsed time between start and stop in seconds.
double Timer::elapsed_seconds() const { return elapsed() * 1e-3; }

unsigned long Timer::elapsed(const high_resolution_clock::time_point & start,
                             const high_resolution_clock::time_point & stop)
{
    if (stop == start)
        throw std::runtime_error("timer not started");
    if (stop < start)
        throw std::runtime_error("timer not stopped");
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    return ms.count();
}

// ======================================================================================
//
// Profiler
//
// ======================================================================================

void Profiler::print_elapsed(const std::string & name, std::ostream & os) const
{
    _assert_timer_exists(name);
    os << "Timer: " << name << _counters.at(name).elapsed() << " ms" << std::endl;
}

void Profiler::create_timer(const std::string & name)
{
    _assert_timer_does_not_exists(name);
    _counters[name] = Timer();
}

void Profiler::erase_timer(const std::string & name)
{
    _assert_timer_exists(name);
    _counters.erase(name);
}

Timer & Profiler::operator[](const std::string & name)
{
    _assert_timer_exists(name);
    return _counters[name];
}

const Timer & Profiler::operator[](const std::string & name) const
{
    _assert_timer_does_not_exists(name);
    return _counters.at(name);
}

// Throws std::invalid_argument if timer named `name` is not registered in Profiler instance.
void Profiler::_assert_timer_exists(const std::string & name) const
{
    if (not _timer_exists(name))
        throw std::invalid_argument("timer '" + name + "' does not exists");
}

// Throws std::invalid_argument if timer named `name` is registered in Profiler instance.
void Profiler::_assert_timer_does_not_exists(const std::string & name) const
{
    if (_timer_exists(name))
        throw std::invalid_argument("timer '" + name + "' already exists");
}

bool Profiler::_timer_exists(const std::string & name) const { return _counters.count(name); }

void timeit(const auto fn, size_t N = 10000)
{
    Timer timer;
    for (size_t i = 0; i < N; ++i)
        fn();
    std::cout << std::fixed << std::setprecision(3) << timer << "\n";
}

} // namespace timeit
} // namespace biospring

#endif // __TIMEIT_HPP__