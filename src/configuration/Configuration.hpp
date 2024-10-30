
#ifndef __CONFIGURATION_H__
#define __CONFIGURATION_H__

#include "Setting.hpp"
#include "utils.hpp"

#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace biospring
{
namespace configuration
{

class Configuration
{
  public:
    SimulationSetting sim;
    StericSetting steric;
    EnergySetting spring;
    EnergySetting hydrophobicity;
    ElectrostaticSetting electrostatic;
    ImpalaSetting imp;
    InsertionVectorSetting ivector;
    ViscositySetting viscosity;
    TrajectorySetting pdbtraj;
    TrajectorySetting xtctraj;
    TrajectorySetting csvsample;
    GridSetting potentialgrid;
    DensityGridSetting densitygrid;
    ProbeSetting probe;
    RigidBodySetting rigidbody;

    Configuration()
        : sim("simulation"), steric("steric"), spring("spring"), hydrophobicity("hydrophobicity"),
          electrostatic("coulomb"), imp("impala"), ivector("insertionvector"), viscosity("viscosity"),
          pdbtraj("pdbtrajectory"), xtctraj("xtctrajectory"), csvsample("csvsampling"), potentialgrid("potentialgrid"),
          densitygrid("densitygrid"), probe("probe"), rigidbody("rigidbody")
    {
        _register(sim);
        _register(steric);
        _register(spring);
        _register(hydrophobicity);
        _register(electrostatic);
        _register(imp);
        _register(ivector);
        _register(viscosity);
        _register(pdbtraj);
        _register(xtctraj);
        _register(csvsample);
        _register(potentialgrid);
        _register(densitygrid);
        _register(probe);
        _register(rigidbody);
    }

    void print(std::ostream & os = std::cout) const
    {
        sim.print();
        os << "\n";
        steric.print();
        os << "\n";
        spring.print();
        os << "\n";
        hydrophobicity.print();
        os << "\n";
        electrostatic.print();
        os << "\n";
        imp.print();
        os << "\n";
        ivector.print();
        os << "\n";
        viscosity.print();
        os << "\n";
        pdbtraj.print();
        os << "\n";
        xtctraj.print();
        os << "\n";
        csvsample.print();
        os << "\n";
        potentialgrid.print();
        os << "\n";
        densitygrid.print();
        os << "\n";
        probe.print();
        os << "\n";
        rigidbody.print();
    }

    bool exists(const std::string & name) { return _allSettingNames.count(name); }

    // Sets a parameter from its full name.
    void setFromString(const std::string & param, const std::string & value)
    {
        if (not exists(param))
            throw std::invalid_argument(utils::string::format("invalid parameter '%s'", param.c_str()));

        const auto tokens = utils::string::split(param, ".");
        if (tokens.size() != 2)
            throw std::invalid_argument(
                utils::string::format("invalid parameter '%s': expected format <group>.<name>", param.c_str()));

        const std::string & group = tokens[0];
        const std::string & name = tokens[1];

        if (group == sim.name)
            sim.setFromString(name, value);
        else if (group == steric.name)
            steric.setFromString(name, value);
        else if (group == spring.name)
            spring.setFromString(name, value);
        else if (group == hydrophobicity.name)
            hydrophobicity.setFromString(name, value);
        else if (group == electrostatic.name)
            electrostatic.setFromString(name, value);
        else if (group == imp.name)
            imp.setFromString(name, value);
        else if (group == ivector.name)
            ivector.setFromString(name, value);
        else if (group == viscosity.name)
            viscosity.setFromString(name, value);
        else if (group == pdbtraj.name)
            pdbtraj.setFromString(name, value);
        else if (group == xtctraj.name)
            xtctraj.setFromString(name, value);
        else if (group == csvsample.name)
            csvsample.setFromString(name, value);
        else if (group == potentialgrid.name)
            potentialgrid.setFromString(name, value);
        else if (group == densitygrid.name)
            densitygrid.setFromString(name, value);
        else if (group == probe.name)
            probe.setFromString(name, value);
        else if (group == rigidbody.name)
            rigidbody.setFromString(name, value);
    }

  protected:
    std::set<std::string> _allSettingNames; // stores settings full names

    // Stores all settings sub-parameters in _allSettingNames.
    // This allows to assess quickly if a parameter exists or not.
    void _register(SettingBase & setting)
    {
        for (const std::string & name : setting.getAllParameters())
            _allSettingNames.insert(setting.name + "." + name);
    }
};

inline Configuration defaultConfiguration()
{
    Configuration config;

    config.sim.nbsteps = -1;
    config.sim.timestep = 0.01;
    config.sim.samplerate = 100;

    config.steric.enable = false;
    config.steric.gridscale = 1.0;
    config.steric.cutoff = 1.0;
    config.steric.mode = "linear";

    config.spring.enable = false;
    config.spring.cutoff = 15.0;
    config.spring.scale = 1.0;

    config.hydrophobicity.enable = false;
    config.hydrophobicity.cutoff = 15.0;
    config.hydrophobicity.scale = 1.0;

    config.electrostatic.enable = false;
    config.electrostatic.cutoff = 16.0;
    config.electrostatic.scale = 1.0;
    config.electrostatic.dielectric = 1.0;

    config.ivector.enable = false;
    config.ivector.vector = {0, 0};

    config.viscosity.enable = false;
    config.viscosity.value = 1.0;

    config.probe.enable = false;
    config.probe.enableelectrostatic = false;
    config.probe.enablesteric = false;
    config.probe.x = 1.0;
    config.probe.y = 1.0;
    config.probe.z = 1.0;
    config.probe.mass = 1.0;
    config.probe.epsilon = 1.0;
    config.probe.radius = 1.0;
    config.probe.charge = 0.0;

    config.rigidbody.enable = false;
    config.rigidbody.enablesampling = false;
    config.rigidbody.enablemontecarlo = false;
    config.rigidbody.montecarlo_translation_norm = 0.1;
    config.rigidbody.montecarlo_rotation_norm = 0.1;
    config.rigidbody.montecarlo_temperature = 298.1;


    config.imp.enable = false;
    config.imp.scale = 1.0;

    return config;
}

} // namespace configuration
} // namespace biospring

#endif // __CONFIGURATION_H__