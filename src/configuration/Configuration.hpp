
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
    GridSetting densitygrid;
    ProbeSetting probe;
    RigidBodySetting rigidbody;

    // Runtime (.msp) enable/disable for the bonded-force-field families
    // that BondedForceFieldReader::buildSprings already decided (at
    // -dihedral* build time, see pdb2spn-cli.cpp) to
    // actually create springs for. Unlike every EnergySetting above,
    // these default to enable=true (see defaultConfiguration()) -- built
    // springs are meant to be active by default, this is an opt-OUT
    // debugging knob (e.g. "dihedralomega.enable = 0" to isolate a
    // family's contribution), not an opt-in feature switch. dihedralchi
    // gates the same springs as the (SIDECHAIN family/chi1-4) collection;
    // "chi" is the name exposed here since that's the physically
    // meaningful term.
    EnergySetting dihedralphi;
    EnergySetting dihedralpsi;
    EnergySetting dihedralomega;
    EnergySetting dihedralchi;
    // The improper (aromatic-ring/His hub planarity) family. Gated like
    // the four proper ones above -- it was the only built dihedral family
    // with no toggle, which made "turn every dihedral off" quietly leave
    // ~40 kJ/mol of improper energy behind on GKinase.
    EnergySetting dihedralplanarity;
    // The nucleic-acid families (see
    // spn::SpringNetwork::DihedralFamilyIndex for why they are separate
    // families rather than SIDECHAIN). dihedralnucleicsugar gates the four
    // furanose ring bonds on their own: the sugar pucker is what selects
    // the A or B helical form, so being able to read -- and remove -- just
    // that contribution is the point of giving it its own family.
    EnergySetting dihedralnucleicbackbone;
    EnergySetting dihedralnucleicchi;
    EnergySetting dihedralnucleicsugar;

    Configuration()
        : sim("simulation"), steric("steric"), spring("spring"), hydrophobicity("hydrophobicity"),
          electrostatic("coulomb"), imp("impala"), ivector("insertionvector"), viscosity("viscosity"),
          pdbtraj("pdbtrajectory"), xtctraj("xtctrajectory"), csvsample("csvsampling"), potentialgrid("potentialgrid"),
          densitygrid("densitygrid"), probe("probe"), rigidbody("rigidbody"), dihedralphi("dihedralphi"),
          dihedralpsi("dihedralpsi"), dihedralomega("dihedralomega"), dihedralchi("dihedralchi"),
          dihedralplanarity("dihedralplanarity"), dihedralnucleicbackbone("dihedralnucleicbackbone"),
          dihedralnucleicchi("dihedralnucleicchi"), dihedralnucleicsugar("dihedralnucleicsugar")
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
        _register(dihedralphi);
        _register(dihedralpsi);
        _register(dihedralomega);
        _register(dihedralchi);
        _register(dihedralplanarity);
        _register(dihedralnucleicbackbone);
        _register(dihedralnucleicchi);
        _register(dihedralnucleicsugar);
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
        os << "\n";
        dihedralphi.print();
        os << "\n";
        dihedralpsi.print();
        os << "\n";
        dihedralomega.print();
        os << "\n";
        dihedralchi.print();
        dihedralplanarity.print();
        dihedralnucleicbackbone.print();
        dihedralnucleicchi.print();
        dihedralnucleicsugar.print();
        os << "\n";
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
        else if (group == dihedralphi.name)
            dihedralphi.setFromString(name, value);
        else if (group == dihedralpsi.name)
            dihedralpsi.setFromString(name, value);
        else if (group == dihedralomega.name)
            dihedralomega.setFromString(name, value);
        else if (group == dihedralchi.name)
            dihedralchi.setFromString(name, value);
        else if (group == dihedralplanarity.name)
            dihedralplanarity.setFromString(name, value);
        else if (group == dihedralnucleicbackbone.name)
            dihedralnucleicbackbone.setFromString(name, value);
        else if (group == dihedralnucleicchi.name)
            dihedralnucleicchi.setFromString(name, value);
        else if (group == dihedralnucleicsugar.name)
            dihedralnucleicsugar.setFromString(name, value);
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
    config.sim.neighborskin = 0.0;

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

    // Deliberately true, unlike every setting above: these only gate
    // springs that -dihedral* already decided to build (see
    // Configuration.hpp's own comment on these members) --
    // an opt-OUT debugging knob, not an opt-in feature switch, so an .msp
    // written before these existed keeps exactly the same behaviour.
    config.dihedralphi.enable = true;
    config.dihedralpsi.enable = true;
    config.dihedralomega.enable = true;
    config.dihedralchi.enable = true;
    config.dihedralplanarity.enable = true;
    config.dihedralnucleicbackbone.enable = true;
    config.dihedralnucleicchi.enable = true;
    config.dihedralnucleicsugar.enable = true;

    return config;
}

} // namespace configuration
} // namespace biospring

#endif // __CONFIGURATION_H__