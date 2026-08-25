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

    // Make the dielectric grow with distance -- epsilon(r) = dielectric * r,
    // r in Angstrom (Warshel & Levitt 1976). A constant 78 screens a pair at
    // contact exactly as hard as one 30 A away, and that is wrong in the
    // direction that matters: no water fits between two atoms 3 A apart, so
    // the short-range interaction is barely screened in reality. Measured on
    // a native alpha-helix hydrogen bond of ubiquitin, the amide/carbonyl
    // attraction is -13.61 kJ/mol unscreened and -0.17 at 78.
    //
    // OFF by default: it changes every electrostatic energy in the model, so
    // it is a modelling choice and not a correction.
    bool distancedependent;

    ElectrostaticSetting(const std::string & name)
        : EnergySetting(name), dielectric(1.0), distancedependent(false)
    {
        _parameterNames = {"enable", "scale", "cutoff", "dielectric", "distancedependent"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "dielectric")
            utils::string::from_string<decltype(dielectric)>(dielectric, s);
        else if (param == "distancedependent")
            _parse_bool(distancedependent, s, param);
        else
            EnergySetting::setFromString(param, s);
    }

    void print(std::ostream & os = std::cout) const override
    {
        EnergySetting::print(os);
        _mspFormatter.print("dielectric", dielectric, os);
        _mspFormatter.print("distancedependent", distancedependent, os);
    }
};

// How the dihedral ghost rings' forces reach the real atoms. Not an energy
// term of its own -- the families above already switch those -- so a small
// setting of its own rather than another EnergySetting.
class DihedralSetting : public SettingBase
{
  public:
    // Project each ring's reaction onto the tangential direction about its
    // axis, so a torsion pushes a substituent only AROUND that axis, as
    // AMBER's does. What it drops carries no torque about the axis at all,
    // and measurably nothing else: on ubiquitin it is 96 % of the applied
    // force, and it deforms the rigid-body mesh rather than twisting it --
    // spring energy after 20000 steps falls from 442.90 to 0.80 kJ/mol.
    // Per-atom agreement with AMBER's own dihedral forces goes from a
    // correlation of 0.000 to 0.650, with the median magnitude ratio moving
    // from 21x to 1.00.
    //
    // ON by default, because there is no case where the unprojected force is
    // the better answer -- it is not a trade-off between two models, it is a
    // defect and its fix. The switch stays only so a run can be reproduced
    // against results predating it, and so the two can be measured side by
    // side; it is not a modelling choice.
    bool tangentialonly;

    // Spread each axis's torsional torque over EVERY substituent bonded to
    // its two axis atoms, instead of leaving it on the one or two reference
    // atoms the ghost rings happen to be built from.
    //
    // The ring's total torque is already right (correlation 0.993 against
    // AMBER over 468 axes on ubiquitin, mean error 0.69 kJ/mol/rad, zero sign
    // flips) but it arrives CONCENTRATED: on 936 of 937 axes the ring pushes
    // fewer atoms than AMBER does -- 1.35 against 4.92 on average -- so each
    // atom it does push receives about 2.38x too much. AMBER splits a torsion
    // into one term per (substituent, substituent) pair, each carrying an
    // equal share, so every substituent on a side ends up with 1/N of that
    // side's torque. This reproduces that split directly: |dphi/dr_i| =
    // 1/(|r_i| sin theta), the inverse perpendicular radius, so a share of
    // torque becomes a force by dividing by that atom's own lever arm. Total
    // torque per side is unchanged by construction -- only its delivery is.
    //
    // The substituent lists are NOT derived at runtime. They come from
    // DIHEDRALAXIS records in the .bi.ff, written by the same generator pass
    // that emits the dihedral terms, and travel through the .nc -- the
    // force field knows which atoms are bonded, and no runtime rule
    // recovers that safely (a distance cutoff cannot separate 1-2 bonds
    // from 1-3 pairs: a methyl's H...H is 1.78 A, shorter than a C-S bond at
    // 1.81 A). A .nc written before those records existed simply carries no
    // lists, and this setting then has nothing to act on -- it warns and
    // leaves the ring's own reference atoms alone rather than guessing.
    //
    // OFF by default, unlike tangentialonly. The projection fixes a defect;
    // this refines how an already-correct total is shared out, and it is the
    // one place in the model where a force reaches an atom no spring
    // connects it to.
    bool distributetorque;

    DihedralSetting(const std::string & name) : SettingBase(name), tangentialonly(true), distributetorque(false)
    {
        _parameterNames = {"tangentialonly", "distributetorque"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "tangentialonly")
            _parse_bool(tangentialonly, s, param);
        else if (param == "distributetorque")
            _parse_bool(distributetorque, s, param);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("tangentialonly", tangentialonly, os);
        _mspFormatter.print("distributetorque", distributetorque, os);
    }
};

// Same enable/scale/cutoff triplet as EnergySetting, plus a `path` to the
// donor/acceptor table (see IO/DonorAcceptorRuleReader.h and
// data/reducerules/*.hbond) -- mirrors how GridSetting/TrajectorySetting
// carry their own input/output file path directly in the .msp file, rather
// than as a separate command-line-only option.
class HydrogenBondSetting : public EnergySetting
{
  public:
    std::string path;
    // Where to write the list of bonds actually held, one line each, at every
    // sample step. Empty (the default) writes nothing. A total energy and a
    // bond count cannot say WHICH bonds are held, and that is the only thing
    // that settles a disagreement about pairing: a fully extended chain
    // reports a perfectly healthy count made entirely of contacts between
    // neighbouring residues.
    std::string log;

    HydrogenBondSetting(const std::string & name) : EnergySetting(name), path(), log()
    {
        _parameterNames = {"enable", "scale", "cutoff", "path", "log"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "path")
            path = s;
        else if (param == "log")
            log = s;
        else
            EnergySetting::setFromString(param, s);
    }

    void print(std::ostream & os = std::cout) const override
    {
        EnergySetting::print(os);
        _mspFormatter.print("path", path, os);
        _mspFormatter.print("log", log, os);
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

class SimulationSetting : public SettingBase
{
  public:
    int nbsteps;
    double timestep;
    size_t samplerate;
    double neighborskin;

    SimulationSetting(const std::string & name)
        : SettingBase(name), nbsteps(0), timestep(0.0), samplerate(1), neighborskin(0.0)
    {
        _parameterNames = {"nbsteps", "timestep", "samplerate", "neighborskin"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "nbsteps")
            utils::string::from_string<decltype(nbsteps)>(nbsteps, s);
        else if (param == "timestep")
            utils::string::from_string<decltype(timestep)>(timestep, s);
        else if (param == "samplerate")
            utils::string::from_string<decltype(samplerate)>(samplerate, s);
        else if (param == "neighborskin")
            utils::string::from_string<decltype(neighborskin)>(neighborskin, s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("nbsteps", nbsteps, os);
        _mspFormatter.print("timestep", timestep, os);
        _mspFormatter.print("samplerate", samplerate, os);
        _mspFormatter.print("neighborskin", neighborskin, os);
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
          mode("linear", {"linear", "lennard-jones-8-6Lewitt", "lennard-jones-8-6Zacharias", "lennard-jones-12-6Amber"})
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
    // Former spelling of lennard-jones-12-6Amber, still carried by .msp files
    // in the wild. Translated to the current name rather than merely added to
    // the choice set: SpringNetwork::_setupForceField dispatches on the
    // string and falls back to "linear" for anything it does not recognise,
    // without complaining, so accepting the old spelling without mapping it
    // would silently run a different force field than the file asks for.
    static constexpr const char * _DEPRECATED_AMBER_MODE = "lennard-jones-8-6Amber";
    static constexpr const char * _AMBER_MODE = "lennard-jones-12-6Amber";

    void _parse_mode(const std::string & value)
    {
        if (value == _DEPRECATED_AMBER_MODE)
        {
            logging::warning("Configuration: %s.mode: '%s' was renamed '%s'; using the latter. "
                             "Update the .msp -- the old name will stop being accepted.",
                             name.c_str(), _DEPRECATED_AMBER_MODE, _AMBER_MODE);
            mode = _AMBER_MODE;
            return;
        }
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