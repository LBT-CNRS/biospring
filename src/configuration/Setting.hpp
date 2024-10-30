#ifndef __SETTING_H__
#define __SETTING_H__

#include <array>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "utils.hpp"

namespace biospring
{
namespace configuration
{

class ChoiceType
{
  public:
    std::string value;

    ChoiceType(const std::string value, const std::set<std::string> & valid) : value(value), _valid_choices(valid) {}

    ChoiceType & operator=(const std::string & v)
    {
        if (not _valid_choices.count(v))
            throw std::invalid_argument(utils::string::format("cannot assign value '%s': invalid choice", v.c_str()));
        value = v;
        return *this;
    }

    friend std::ostream & operator<<(std::ostream & os, const ChoiceType & choice)
    {
        os << choice.value;
        return os;
    }

    operator const std::string &() const { return value; }

    bool operator==(const ChoiceType & rhs) const { return value == rhs.value; }

  protected:
    std::set<std::string> _valid_choices;
};

struct _SettingFormatter
{
    std::string prefix;
    std::string separator;
    std::string assignment_operator;
    std::string eol;

    template <typename T> std::string format(const std::string & name, const T & value) const
    {
        std::stringstream ss;
        ss << prefix << separator << name << assignment_operator << value << eol;
        return ss.str();
    }

    template <typename T, size_t N> std::string format(const std::string & name, const std::array<T, N> & value) const
    {
        std::stringstream ss;
        ss << prefix << separator << name << assignment_operator;
        for (size_t i = 0; i < N - 1; ++i)
            ss << value[i] << " ";
        ss << value[N - 1] << eol;
        return ss.str();
    }

    template <typename T> void print(const std::string & name, const T & value, std::ostream & os = std::cout) const
    {
        os << format(name, value);
    }
};

class SettingBase
{
  public:
    std::string name;

    SettingBase() = delete;
    SettingBase(const std::string & name) : name(name), _parameterNames(), _mspFormatter({name, ".", " = ", "\n"}) {}

    virtual void setFromString(const std::string & param, const std::string & s) = 0;
    virtual void print(std::ostream & os = std::cout) const = 0;

    const std::vector<std::string> & getAllParameters() const { return _parameterNames; }

  protected:
    std::vector<std::string> _parameterNames;
    _SettingFormatter _mspFormatter;

    void _parse_bool(bool & b, const std::string & s, const std::string param)
    {
        bool success = utils::string::from_string<bool>(b, s);
        if (not success)
            logging::die("Configuration: %s.%s: cannot assign value '%s': cannot interpret boolean", name.c_str(),
                         param.c_str(), s.c_str());
    }
};

class EnergySetting : public SettingBase
{
  public:
    bool enable;
    double scale;
    double cutoff;

    EnergySetting(const std::string & name) : SettingBase(name), enable(false), scale(1.0), cutoff(0.0)
    {
        _parameterNames = {"enable", "scale", "cutoff"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "scale")
            utils::string::from_string<decltype(scale)>(scale, s);
        else if (param == "cutoff")
            utils::string::from_string<decltype(cutoff)>(cutoff, s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("scale", scale, os);
        _mspFormatter.print("cutoff", cutoff, os);
    }
};

class ElectrostaticSetting : public EnergySetting
{
  public:
    double dielectric;

    ElectrostaticSetting(const std::string & name) : EnergySetting(name), dielectric(1.0)
    {
        _parameterNames = {"enable", "scale", "cutoff", "dielectric"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "dielectric")
            utils::string::from_string<decltype(dielectric)>(dielectric, s);
        else
            EnergySetting::setFromString(param, s);
    }

    void print(std::ostream & os = std::cout) const override
    {
        EnergySetting::print(os);
        _mspFormatter.print("dielectric", dielectric, os);
    }
};

class TrajectorySetting : public SettingBase
{
  public:
    bool enable;
    std::string path;
    size_t frequency;

    TrajectorySetting(const std::string & name) : SettingBase(name), enable(false), path(), frequency(100)
    {
        _parameterNames = {"enable", "path", "frequency"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "frequency")
            utils::string::from_string<decltype(frequency)>(frequency, s);
        else if (param == "path")
            path = s;
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("path", path, os);
        _mspFormatter.print("frequency", frequency, os);
    }
};

class GridSetting : public SettingBase
{
  public:
    bool enable;
    std::string path;
    double scale;

    GridSetting(const std::string & name) : SettingBase(name), enable(false), path(), scale(1.0)
    {
        _parameterNames = {"enable", "path", "scale"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "scale")
            utils::string::from_string<decltype(scale)>(scale, s);
        else if (param == "path")
            path = s;
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("path", path, os);
        _mspFormatter.print("scale", scale, os);
    }
};

class DensityGridSetting : public SettingBase
{
  public:
    bool enable;
    std::string path;

    DensityGridSetting(const std::string & name) : SettingBase(name), enable(false), path()
    {
        _parameterNames = {"enable", "path"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "path")
            path = s;
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("path", path, os);
    }
};

class SimulationSetting : public SettingBase
{
  public:
    int nbsteps;
    double timestep;
    size_t samplerate;

    SimulationSetting(const std::string & name) : SettingBase(name), nbsteps(0), timestep(0.0), samplerate(1)
    {
        _parameterNames = {"nbsteps", "timestep", "samplerate"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "nbsteps")
            utils::string::from_string<decltype(nbsteps)>(nbsteps, s);
        else if (param == "timestep")
            utils::string::from_string<decltype(timestep)>(timestep, s);
        else if (param == "samplerate")
            utils::string::from_string<decltype(samplerate)>(samplerate, s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("nbsteps", nbsteps, os);
        _mspFormatter.print("timestep", timestep, os);
        _mspFormatter.print("samplerate", samplerate, os);
    }
};

class StericSetting : public SettingBase
{
  public:
    bool enable;
    double gridscale;
    double cutoff;
    ChoiceType mode;

    StericSetting(const std::string & name)
        : SettingBase(name), enable(false), gridscale(1.0), cutoff(0.0),
          mode("linear", {"linear", "lennard-jones-8-6Lewitt", "lennard-jones-8-6Zacharias", "lennard-jones-8-6Amber"})
    {
        _parameterNames = {"enable", "gridscale", "cutoff", "mode"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "gridscale")
            utils::string::from_string<decltype(gridscale)>(gridscale, s);
        else if (param == "cutoff")
            utils::string::from_string<decltype(cutoff)>(cutoff, s);
        else if (param == "mode")
            _parse_mode(s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("gridscale", gridscale, os);
        _mspFormatter.print("cutoff", cutoff, os);
        _mspFormatter.print("mode", mode, os);
    }

  protected:
    void _parse_mode(const std::string & value)
    {
        try
        {
            mode = value;
        }
        catch (const std::invalid_argument &)
        {
            logging::die("Configuration: %s: cannot assign value '%s': invalid choice", name.c_str(), value.c_str());
        }
    }
};

class ImpalaSetting : public SettingBase
{
  public:
    bool enable;
    double scale;

    ImpalaSetting(const std::string & name) : SettingBase(name), enable(false), scale(1.0)
    {
        _parameterNames = {"enable", "scale"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "scale")
            utils::string::from_string<decltype(scale)>(scale, s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("scale", scale, os);
    }
};

class ViscositySetting : public SettingBase
{
  public:
    bool enable;
    double value;

    ViscositySetting(const std::string & name) : SettingBase(name), enable(false), value(1.0)
    {
        _parameterNames = {"enable", "value"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "value")
            utils::string::from_string<decltype(value)>(value, s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("value", value, os);
    }
};

class InsertionVectorSetting : public SettingBase
{
  public:
    bool enable;
    std::array<size_t, 2> vector;

    InsertionVectorSetting(const std::string & name) : SettingBase(name), enable(false), vector({0, 0})
    {
        _parameterNames = {"enable", "vector"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "vector")
            _parseInsertionVector(s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("vector", vector, os);
    }

  protected:
    void _parseInsertionVector(const std::string & s)
    {
        const std::vector<std::string> tokens = utils::string::split(s);

        if (tokens.size() != 2)
            logging::die("Configuration: %s: Invalid value for argument insertionvector: \"%s\" (expected: <particle1 "
                         "external id> <particle2 external id> ",
                         name.c_str(), s.c_str());

        if (not utils::string::from_string<size_t>(vector[0], tokens[0]))
            logging::die("Configuraton: %s: Cannot parse first particle id \"%s\"", name.c_str(), tokens[1].c_str());

        if (not utils::string::from_string<size_t>(vector[1], tokens[1]))
            logging::die("Configuraton: %s: Cannot parse second particle id \"%s\"", name.c_str(), tokens[1].c_str());
    }
};

class ProbeSetting : public SettingBase
{
  public:
    bool enable;
    bool enableelectrostatic;
    bool enablesteric;
    double x;
    double y;
    double z;
    double mass;
    double epsilon;
    double radius;
    double charge;

    ProbeSetting(const std::string & name)
        : SettingBase(name), enable(false), enableelectrostatic(false), enablesteric(false), x(0.0), y(0.0), z(0.0),
          mass(1.0), epsilon(1.0), radius(1.0), charge(0.0)
    {
        _parameterNames = {"enable", "enableelectrostatic", "enablesteric", "x", "y", "z", "mass", "epsilon", "radius",
                           "charge"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "enableelectrostatic")
            _parse_bool(enableelectrostatic, s, param);
        else if (param == "enablesteric")
            _parse_bool(enablesteric, s, param);
        else if (param == "x")
            utils::string::from_string<decltype(x)>(x, s);
        else if (param == "y")
            utils::string::from_string<decltype(y)>(y, s);
        else if (param == "z")
            utils::string::from_string<decltype(z)>(z, s);
        else if (param == "mass")
            utils::string::from_string<decltype(mass)>(mass, s);
        else if (param == "epsilon")
            utils::string::from_string<decltype(epsilon)>(epsilon, s);
        else if (param == "radius")
            utils::string::from_string<decltype(radius)>(radius, s);
        else if (param == "charge")
            utils::string::from_string<decltype(charge)>(charge, s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("enableelectrostatic", enableelectrostatic, os);
        _mspFormatter.print("enablesteric", enablesteric, os);
        _mspFormatter.print("x", x, os);
        _mspFormatter.print("y", y, os);
        _mspFormatter.print("z", z, os);
        _mspFormatter.print("mass", mass, os);
        _mspFormatter.print("epsilon", epsilon, os);
        _mspFormatter.print("radius", radius, os);
        _mspFormatter.print("charge", charge, os);
    }
};

class RigidBodySetting : public SettingBase
{
  public:
    bool enable, enablesampling, enablemontecarlo;
    double montecarlo_translation_norm;  // random translation to apply each step in Å
    double montecarlo_rotation_norm;     // random rotation to apply each step in °
    double montecarlo_temperature;

    RigidBodySetting(const std::string & name) : SettingBase(name), 
        enable(false), enablesampling(false), enablemontecarlo(false),
        montecarlo_translation_norm(0.1), montecarlo_rotation_norm(0.1),
        montecarlo_temperature(298.1)
    {
        _parameterNames = {"enable", "enablesampling", "enablemontecarlo",
            "montecarlo_translation_norm", "montecarlo_rotation_norm",
            "montecarlo_temperature"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "enablesampling")
            _parse_bool(enablesampling, s, param);
        else if (param == "enablemontecarlo")
            _parse_bool(enablemontecarlo, s, param);
        else if (param == "montecarlo_translation_norm")
            utils::string::from_string<decltype(montecarlo_translation_norm)>(montecarlo_translation_norm, s);
        else if (param == "montecarlo_rotation_norm")
            utils::string::from_string<decltype(montecarlo_rotation_norm)>(montecarlo_rotation_norm, s);
        else if (param == "montecarlo_temperature")
            utils::string::from_string<decltype(montecarlo_temperature)>(montecarlo_temperature, s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("enablesampling", enablesampling, os);
        _mspFormatter.print("enablemontecarlo", enablemontecarlo, os);
        _mspFormatter.print("montecarlo_translation_norm", montecarlo_translation_norm, os);
        _mspFormatter.print("montecarlo_rotation_norm", montecarlo_rotation_norm, os);
        _mspFormatter.print("montecarlo_temperature", montecarlo_temperature, os);
    }
};

} // namespace configuration
} // namespace biospring

#endif // __SETTING_H__