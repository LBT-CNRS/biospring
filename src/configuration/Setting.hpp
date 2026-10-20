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

// The pairwise hydrophobic term, which needs one parameter the other energies
// do not: the distance over which it decays.
//
// It sits here rather than being written into the law for the same reason
// `dielectric` does: it is a property of the solvent that the modeller picks,
// not a universal constant. The default, 10 A, is the decay length
// Israelachvili & Pashley measured between two MACROSCOPIC hydrophobic
// surfaces (Nature 300:341, 1982, published as 22 exp(-D/10) mJ.m-2 with D in
// A), which is where the law in hydrophobic_shared.h comes from. Later work
// puts the range anywhere between 3 and 10 A depending on the system, and two
// coarse-grained beads are not two macroscopic plates -- so 10 is the value
// the published law calls for, not a value validated between beads.
class HydrophobicitySetting : public EnergySetting
{
  public:
    double decaylength;

    HydrophobicitySetting(const std::string & name) : EnergySetting(name), decaylength(10.0)
    {
        _parameterNames = {"enable", "scale", "cutoff", "decaylength"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "decaylength")
            utils::string::from_string<decltype(decaylength)>(decaylength, s);
        else
            EnergySetting::setFromString(param, s);
    }

    void print(std::ostream & os = std::cout) const override
    {
        EnergySetting::print(os);
        _mspFormatter.print("decaylength", decaylength, os);
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

// The Langevin thermostat: a switch of its own, so it can be turned off
// without touching the temperature and on without having to remember that a
// temperature of zero means "off".
//
// Its FRICTION is viscosity.value -- not a second coefficient. Friction and
// the random kicks are one physical effect, the solvent this model does not
// simulate, and the fluctuation-dissipation relation fixes the size of one
// from the other. Enabling the thermostat with the viscosity off would ask for
// kicks from a solvent that does not drag, which is not a thermostat; setup()
// refuses it rather than running a model nobody asked for.
class ThermostatSetting : public SettingBase
{
  public:
    bool enable;
    double temperature;   // K

    ThermostatSetting(const std::string & name) : SettingBase(name), enable(false), temperature(300.0)
    {
        _parameterNames = {"enable", "temperature"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "temperature")
            utils::string::from_string<decltype(temperature)>(temperature, s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("temperature", temperature, os);
    }
};

class SimulationSetting : public SettingBase
{
  public:
    int nbsteps;
    double timestep;
    size_t samplerate;
    double neighborskin;

    // Width of a neighbour-search cell, in A. Zero asks for the automatic
    // choice, which is what every run should want; see
    // SpringNetwork::chooseCellWidth. It is NOT a cutoff, and making it one
    // is the pessimal setting rather than the neutral one.
    double cellsize;

    SimulationSetting(const std::string & name)
        : SettingBase(name), nbsteps(0), timestep(0.0), samplerate(1), neighborskin(-1.0), cellsize(0.0)
    {
        _parameterNames = {"nbsteps", "timestep", "samplerate", "neighborskin", "cellsize"};
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
        else if (param == "cellsize")
            utils::string::from_string<decltype(cellsize)>(cellsize, s);
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("nbsteps", nbsteps, os);
        _mspFormatter.print("timestep", timestep, os);
        _mspFormatter.print("samplerate", samplerate, os);
        _mspFormatter.print("neighborskin", neighborskin, os);
        _mspFormatter.print("cellsize", cellsize, os);
    }
};

class StericSetting : public SettingBase
{
  public:
    bool enable;
    double gridscale;
    double cutoff;
    ChoiceType mode;
    // How two radii combine. "default" means the law's own rule, which is the
    // rule of the parameter set that law was published with -- see
    // STERIC_RADIUS_RULE_* in forcefield/energy/steric.hpp. Set it explicitly
    // only when pairing a .ff with a law it did not come from, because the
    // `radius` column does not hold the same quantity in every file: amber*.ff
    // holds R* = rmin/2 and its radii SUM, while CAonlyLewitt.ff holds a CA-CA
    // contact distance and its radii must be MEANED. The full table is in
    // forcefield/shared/steric_shared.h. Zacharias' law ignores this: its
    // product on the radius is inseparable from its product on epsilon.
    ChoiceType radiusrule;

    StericSetting(const std::string & name)
        : SettingBase(name), enable(false), gridscale(1.0), cutoff(0.0),
          mode("linear", {"linear", "lennard-jones-8-6Lewitt", "lennard-jones-8-6Zacharias", "lennard-jones-12-6Amber"}),
          radiusrule("default", {"default", "sum", "geometric-mean", "arithmetic-mean"})
    {
        _parameterNames = {"enable", "gridscale", "cutoff", "mode", "radiusrule"};
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
        else if (param == "radiusrule")
            radiusrule = s;
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("gridscale", gridscale, os);
        _mspFormatter.print("cutoff", cutoff, os);
        _mspFormatter.print("mode", mode, os);
        _mspFormatter.print("radiusrule", radiusrule, os);
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

    // HiRE-RNA v3's planarity factor (arXiv:1404.0568, equation 5), which
    // MULTIPLIES the well: a pair of bases scores 1 when every one of their
    // plane atoms lies in the other's plane and falls off as they skew. It
    // needs a base plane at BOTH ends (baseplane.path), so a residue without
    // one -- every amino acid -- keeps a factor of exactly 1 and behaves as
    // before, to the bit.
    //
    // This is deliberately NOT applied to protein backbone hydrogen bonds even
    // though their amide plane is available. In an alpha helix the donor's and
    // the acceptor's peptide planes are NOT coplanar, so requiring coplanarity
    // would unmake the helix. Coplanarity is a property of two aromatic rings
    // meeting edge to edge, not of an amide bond.
    bool planarity;
    // delta in equation (5), in A. 1.5 A, chosen by MEASUREMENT OF BEHAVIOUR on
    // three simple systems, after a first value of 2.5 taken from a structural
    // histogram alone turned out to be too wide.
    //
    // The structural measurement still holds and is worth keeping: over 8PKL's
    // 1893 real base pairs (two or more base-edge N/O contacts under 3.3 A,
    // excluding sequence neighbours, which are STACKED rather than paired) the
    // plane atoms sit 0.55 A out of the partner's plane at the median, against
    // 3.32 A for a stack and 2.80 A for a single skewed contact. What it does
    // NOT determine is where to put delta: its p90 tail at 2.95 A is mostly
    // non-canonical pairs the 2-contact criterion admitted, and reading that
    // tail as real pairs needing tolerance is what made 2.5 too wide.
    //
    // Scanned instead on what the term DOES, over 0.75 to 3.5 A:
    //
    //   21-mer DNA duplex, starting native -- all 20 native pairs survive at
    //     every width (-772 to -822 kJ/mol, a 6 % spread), while the parasitic
    //     pairing the hydrogen bond term's own comment documents falls from
    //     -7.4 kJ/mol at 3.5 A to exactly -0.0 at 1.5 A and below;
    //   RNA hairpin, starting unfolded, ten paired starts -- the native-register
    //     share of the pairing rises monotonically as delta narrows, 0.575 +/-
    //     0.085 at 3.5 A to 0.824 +/- 0.054 at 0.75 A;
    //   anticodon trimer, starting in register, three paired starts -- the
    //     OPPOSITE trend: in-register binding deepens as delta widens, -44.2
    //     +/- 11.1 kJ/mol at 1.0 A to -72.3 +/- 6.8 at 2.5 A. This is the
    //     stringent system, a 3-nt trimer with no neighbours, no continuous
    //     backbone and nothing else holding it.
    //
    // The three do not agree on an optimum, they bracket one. 1.5 A is the only
    // width satisfying all three: the narrowest that costs nothing measurable in
    // binding (18.8 kJ/mol below 2.5 A against a 10.8 combined standard error,
    // 1.7 sigma, where 1.0 A is 2.2 sigma worse and so genuinely worse), and the
    // widest that leaves the duplex parasite at exactly zero.
    //
    // A RIBOSOME WANTS 2.0 A, and the reason is physical rather than a fudge.
    // Narrowing delta on a 21-mer B-DNA duplex costs nothing -- all 20 native
    // pairs survive at every width -- because a B-form pair is ideal and
    // coplanar. Many ribosomal RNA pairs are NOT: non-canonical, Hoogsteen,
    // sheared, part of a base triple. They therefore pay a factor the ideal pair
    // does not, and on 8PKL's all-atom region the native pairing weight falls
    // 34 % between 2.5 and 1.5 A:
    //
    //   delta   native  correct share  native kept
    //   off       22.1          0.280           --
    //   2.5       18.5          0.614         84 %
    //   2.0       16.9          0.679         76 %
    //   1.5       14.5          0.748         66 %
    //
    // 2.0 A recovers 17 % more native pairing than 1.5 for 0.069 of
    // cleanliness. Set it per configuration; this default is the simple-nucleic
    // one.
    //
    // Two things that breakdown settled, worth not re-deriving: pairing weight
    // involving a COARSE bead is exactly zero (a one-bead residue has no
    // base-edge N or O and cannot pair at all, so the multi-resolution boundary
    // is not a source of error), and the crystal holds 111 pairs with a SINGLE
    // edge contact -- base triples, A-minor motifs -- which a two-contact ground
    // truth miscounts as parasitic. A ground truth demanding two contacts
    // understates a ribosome.
    double planaritywidth;
    // The exponent in the angular weight max(cos theta, 0)^n, applied at BOTH
    // ends. Two is what this code has always used; HiRE-RNA v3 uses SIX
    // (nu(alpha) = cos^6), a far narrower acceptance. Unlike the planarity
    // factor this says nothing about bases, so it is the one part of the
    // HiRE-RNA treatment that can change a PROTEIN hydrogen bond.
    //
    // The OpenCL kernel hard-codes two, so a value other than two is refused on
    // that backend rather than silently disagreeing with the CPU.
    int angularexponent;

    HydrogenBondSetting(const std::string & name)
        : EnergySetting(name), path(), log(), planarity(false), planaritywidth(1.5), angularexponent(2)
    {
        _parameterNames = {"enable",    "scale",          "cutoff", "path", "log", "planarity",
                           "planaritywidth", "angularexponent"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "path")
            path = s;
        else if (param == "log")
            log = s;
        else if (param == "planarity")
            _parse_bool(planarity, s, param);
        else if (param == "planaritywidth")
            utils::string::from_string<decltype(planaritywidth)>(planaritywidth, s);
        else if (param == "angularexponent")
            utils::string::from_string<decltype(angularexponent)>(angularexponent, s);
        else
            EnergySetting::setFromString(param, s);
    }

    void print(std::ostream & os = std::cout) const override
    {
        EnergySetting::print(os);
        _mspFormatter.print("path", path, os);
        _mspFormatter.print("log", log, os);
        _mspFormatter.print("planarity", planarity, os);
        _mspFormatter.print("planaritywidth", planaritywidth, os);
        _mspFormatter.print("angularexponent", angularexponent, os);
    }
};

// Where the base plane frames come from. Its own group because TWO terms read
// them: the stacking term, whose whole geometry they are, and the hydrogen
// bond's planarity factor. Filing the path inside either one would make the
// other depend on a setting named after something it is not.
class BasePlaneSetting : public SettingBase
{
  public:
    std::string path;

    BasePlaneSetting(const std::string & name) : SettingBase(name), path() { _parameterNames = {"path"}; }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "path")
            path = s;
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override { _mspFormatter.print("path", path, os); }
};

// Base stacking, after HiRE-RNA v3 (arXiv:1404.0568). Same enable/scale/cutoff
// triplet as EnergySetting; the frames come from baseplane.path.
//
// The well's three numbers live HERE rather than in the force field, unlike the
// hydrogen bond's, because there is one global parameter set: stacking has no
// per-couple named groups to resolve. Their defaults are measured -- see
// forcefield/energy/stacking.hpp for the measurement and its sample.
class StackingSetting : public EnergySetting
{
  public:
    double epsilon;     // well depth, kJ.mol-1
    double equilibrium; // r_st, A
    double width;       // sigma, A^2 -- the factor is e^(-(r-r_st)^2/sigma)

    StackingSetting(const std::string & name)
        : EnergySetting(name), epsilon(8.0), equilibrium(3.92), width(0.85)
    {
        _parameterNames = {"enable", "scale", "cutoff", "epsilon", "equilibrium", "width"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "epsilon")
            utils::string::from_string<decltype(epsilon)>(epsilon, s);
        else if (param == "equilibrium")
            utils::string::from_string<decltype(equilibrium)>(equilibrium, s);
        else if (param == "width")
            utils::string::from_string<decltype(width)>(width, s);
        else
            EnergySetting::setFromString(param, s);
    }

    void print(std::ostream & os = std::cout) const override
    {
        EnergySetting::print(os);
        _mspFormatter.print("epsilon", epsilon, os);
        _mspFormatter.print("equilibrium", equilibrium, os);
        _mspFormatter.print("width", width, os);
    }
};

// The SECOND stage of making a covalent bond at run time: once two residues
// have been held in the attack conformation long enough, the spring topology is
// REWRITTEN so the new bond exists and cannot come undone.
//
// The first stage -- bringing the nucleophile onto the electrophile at the
// Burgi-Dunitz angle -- is the hydrogen bond term with an out-of-plane lobe
// site (see spn/ParticleProperty.h), and nothing here duplicates it. What this
// adds is the irreversible part, which no potential well can express: a bond
// that has formed is a different molecule, not a deeper minimum of the same one.
class PeptideBondSetting : public SettingBase
{
  public:
    bool enable;
    // The .rbody the network was built from. The new springs are created by the
    // SAME rule that would have created them had the bond been there from the
    // start, so the product is the model's own peptide and not a special case.
    std::string path;
    // Which rule describes the new bond: the suffix after the residue code, so
    // "PSI" picks up A_PSI, R_PSI, ... -- {CA, C, O, +N, +H, +CA}, the peptide
    // plane ahead of a residue.
    std::string group;
    // Which PAIR inside that group is the bond being made, as
    // "<electrophile>:<nucleophile>". The rule itself cannot say: it is six
    // atoms sprung pairwise, with no distinguished pair among the fifteen.
    std::string bond;
    double stiffness; // kJ.mol-1.A-2, what the new springs get
    double distance;  // A, below which the nucleophile counts as attacking
    double weight;    // 0..1, the angular factor both ends must reach together
    unsigned dwell;   // consecutive steps the two conditions must hold
    // Over how many steps the new springs are drawn from the length they were
    // BORN at -- the attack distance, so no energy appears out of nothing -- to
    // the length an ideal trans peptide plane gives them. 0 leaves them where
    // they were born, which is a 2.9 A "peptide bond".
    unsigned ramp;
    // The names an atom may have to be the LEAVING GROUP, space-separated (e.g.
    // "OXT O3'"), and empty by default. It is looked for among the atoms BONDED
    // to the electrophile, wherever they live, and what breaks is every spring
    // between it and the electrophile's residue.
    //
    // That one rule covers both chemistries. A free acid's OXT is inside that
    // residue, so all its springs go and it departs whole -- the water that
    // leaves. An aminoacyl-tRNA's O3' belongs to the tRNA, so only the ester
    // C-O goes and the oxygen stays on its ribose: the peptide transfers to the
    // other tRNA and this one is left deacylated, which is peptidyl transfer.
    //
    // Breaking is both ways: a spring is also what keeps a pair out of the
    // steric and Coulomb terms, so what pushes the group away afterwards is
    // ordinary sterics.
    std::string leaving;
    // The .psite file that says WHERE a carbonyl carbon may be attacked from:
    // its antecedent and the plane that splits its two faces. Its own file and
    // its own reader on purpose -- a reactive site is not a hydrogen bond, and
    // filing it in the .hbond table meant either file could silently overwrite
    // the other's direction on an atom both named. Empty falls back to whatever
    // the .hbond table tagged, which is what every configuration written before
    // this setting existed relies on.
    std::string sites;
    // The .grp the network was REDUCED with, when it was. amber.grp renames
    // every atom per residue -- C becomes AC for alanine, RC for arginine --
    // so on a reduced network the names in `bond` and in the .rbody rules find
    // nothing without it. Same table, same purpose, as pdb2spn --rigidbody's
    // own naming translation.
    std::string naming;
    std::string log;

    PeptideBondSetting(const std::string & name)
        : SettingBase(name), enable(false), path(), group("PSI"), bond("C:N"), stiffness(650.0), distance(3.2),
          weight(0.5), dwell(200), ramp(2000), leaving(), sites(), naming(), log()
    {
        _parameterNames = {"enable",  "path",  "group",   "bond",   "stiffness", "distance",
                           "weight",  "dwell", "ramp",    "leaving", "sites",    "naming",
                           "log"};
    }

    void setFromString(const std::string & param, const std::string & s) override
    {
        if (param == "enable")
            _parse_bool(enable, s, param);
        else if (param == "path")
            path = s;
        else if (param == "group")
            group = s;
        else if (param == "bond")
            bond = s;
        else if (param == "stiffness")
            utils::string::from_string<decltype(stiffness)>(stiffness, s);
        else if (param == "distance")
            utils::string::from_string<decltype(distance)>(distance, s);
        else if (param == "weight")
            utils::string::from_string<decltype(weight)>(weight, s);
        else if (param == "dwell")
            utils::string::from_string<decltype(dwell)>(dwell, s);
        else if (param == "ramp")
            utils::string::from_string<decltype(ramp)>(ramp, s);
        else if (param == "leaving")
            leaving = s;
        else if (param == "sites")
            sites = s;
        else if (param == "naming")
            naming = s;
        else if (param == "log")
            log = s;
        else
            logging::die("%s: unknown parameter '%s'", name.c_str(), param.c_str());
    }

    void print(std::ostream & os = std::cout) const override
    {
        _mspFormatter.print("enable", enable, os);
        _mspFormatter.print("path", path, os);
        _mspFormatter.print("group", group, os);
        _mspFormatter.print("bond", bond, os);
        _mspFormatter.print("stiffness", stiffness, os);
        _mspFormatter.print("distance", distance, os);
        _mspFormatter.print("weight", weight, os);
        _mspFormatter.print("dwell", dwell, os);
        _mspFormatter.print("ramp", ramp, os);
        _mspFormatter.print("leaving", leaving, os);
        _mspFormatter.print("sites", sites, os);
        _mspFormatter.print("naming", naming, os);
        _mspFormatter.print("log", log, os);
    }
};

} // namespace configuration
} // namespace biospring

#endif // __SETTING_H__
