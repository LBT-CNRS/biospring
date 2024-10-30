#include "NetCDFReader.h"
#include "Vector3f.h"
#include "logging.h"
#include "topology.hpp"

using namespace biospring;

void NetCDFReader::addSpringsToSpn()
{
    for (size_t i = 0; i < _sbuffer.number_of_springs; ++i)
    {
        _topology.add_spring(_sbuffer.springs[i][0], _sbuffer.springs[i][1], _sbuffer.springsequilibriums[i],
                             _sbuffer.springsstiffnesses[i]);
    }
}

void NetCDFReader::addParticlesToSpn()
{
    char buf[5] = "";
    for (size_t i = 0; i < _pbuffer.number_of_particles; ++i)
    {
        topology::Particle p;
        p.properties().set_position(
            Vector3f(_pbuffer.coordinates[i][0], _pbuffer.coordinates[i][1], _pbuffer.coordinates[i][2]));
        p.properties().set_atom_id(_pbuffer.ids[i]);
        p.properties().set_residue_id(_pbuffer.resids[i]);
        p.properties().set_charge(_pbuffer.charges[i]);
        p.properties().set_radius(_pbuffer.radii[i]);
        p.properties().set_epsilon(_pbuffer.epsilons[i]);
        p.properties().set_dynamic((_pbuffer.dynamic_states[i] == 0) ? false : true);
        p.properties().set_mass(_pbuffer.masses[i]);

        if (std::abs(p.properties().mass()) < 1e-6)
        {
            logging::warning("Mass of particle %d is 0.0. Changing it to 1.0", i);
            p.properties().set_mass(1.0);
        }

        memset(buf, 0, PARTICLE_NAME_LENGTH);
        strncpy(buf, _pbuffer.particlenames[i], PARTICLE_NAME_LENGTH);
        p.properties().set_name(buf);

        memset(buf, 0, RESIDUE_NAME_LENGTH);
        strncpy(buf, _pbuffer.resnames[i], RESIDUE_NAME_LENGTH);
        p.properties().set_residue_name(buf);

        memset(buf, 0, CHAIN_NAME_LENGTH);
        strncpy(buf, _pbuffer.chainnames[i], CHAIN_NAME_LENGTH);
        p.properties().set_chain_name(buf);

        /* we test if these arrays are NULL since they are not mandatory
         * in the Nc. */
        // topology::IMPProperties impProperties = topology::IMPProperties::build();
        if (_pbuffer.surface_accessibilities)
            p.properties().set_imp(topology::IMPProperties::build().solvent_accessible_surface(_pbuffer.surface_accessibilities[i]));
            // impProperties.set_solvent_accessible_surface(_pbuffer.surface_accessibilities[i]);
            // p.properties().imp().set_solvent_accessible_surface(_pbuffer.surface_accessibilities[i]);

        if (_pbuffer.hscales)
            p.properties().set_imp(topology::IMPProperties::build().transfert_energy_by_accessible_surface(_pbuffer.hscales[i]));
            // impProperties.set_transfert_energy_by_accessible_surface(_pbuffer.hscales[i]);
            // p.properties().imp().set_transfert_energy_by_accessible_surface(_pbuffer.hscales[i]);

        _topology.add_particle(p);
    }
}

void NetCDFReader::read()
{
    try
    {
        _file = new netCDF::NcFile(_filename, netCDF::NcFile::read);
        readParticles();
        addParticlesToSpn();

        readSprings();
        addSpringsToSpn();
    }
    catch (netCDF::exceptions::NcException & e)
    {
        logging::die("%s", e.what());
        logging::die("%s not found.", _filename.c_str());
    }
}

void NetCDFReader::readParticles()
{
    readNumberOfParticles();

    netCDF::NcVar data;

    data = getNcVar("particleids");
    checkNDims(data, 1);
    checkDim(data, 0, _pbuffer.number_of_particles);
    data.getVar(_pbuffer.ids);

    data = getNcVar("resids");
    checkNDims(data, 1);
    checkDim(data, 0, _pbuffer.number_of_particles);
    data.getVar(_pbuffer.resids);

    data = getNcVar("resnames");
    checkNDims(data, 2);
    checkDim(data, 0, _pbuffer.number_of_particles);
    checkDim(data, 1, RESIDUE_NAME_LENGTH);
    data.getVar(_pbuffer.resnames);

    data = getNcVar("particlenames");
    checkNDims(data, 2);
    checkDim(data, 0, _pbuffer.number_of_particles);
    checkDim(data, 1, PARTICLE_NAME_LENGTH);
    data.getVar(_pbuffer.particlenames);

    data = getNcVar("chainnames");
    checkNDims(data, 2);
    checkDim(data, 0, _pbuffer.number_of_particles);
    checkDim(data, 1, CHAIN_NAME_LENGTH);
    data.getVar(_pbuffer.chainnames);

    data = getNcVar("dynamicstate");
    checkNDims(data, 1);
    checkDim(data, 0, _pbuffer.number_of_particles);
    data.getVar((bool *)_pbuffer.dynamic_states);

    data = getNcVar("coordinates");
    checkNDims(data, 2);
    checkDim(data, 0, _pbuffer.number_of_particles);
    checkDim(data, 1, 3);
    data.getVar(_pbuffer.coordinates);

    data = getNcVar("radii");
    checkNDims(data, 1);
    checkDim(data, 0, _pbuffer.number_of_particles);
    data.getVar(_pbuffer.radii);

    data = getNcVar("mass");
    checkNDims(data, 1);
    checkDim(data, 0, _pbuffer.number_of_particles);
    data.getVar(_pbuffer.masses);

    data = getNcVar("epsilon");
    checkNDims(data, 1);
    checkDim(data, 0, _pbuffer.number_of_particles);
    data.getVar(_pbuffer.epsilons);

    data = getNcVar("charges");
    checkNDims(data, 1);
    checkDim(data, 0, _pbuffer.number_of_particles);
    data.getVar(_pbuffer.charges);

    data = getNcVar("surfaceaccessibility", false);
    if (not data.isNull())
    {
        checkNDims(data, 1);
        checkDim(data, 0, _pbuffer.number_of_particles);
        data.getVar(_pbuffer.surface_accessibilities);
    }

    data = getNcVar("hydrophobicityscale", false);
    if (not data.isNull())
    {
        checkNDims(data, 1);
        checkDim(data, 0, _pbuffer.number_of_particles);
        data.getVar(_pbuffer.hscales);
    }
}

void NetCDFReader::readSprings()
{
    readNumberOfSprings();

    if (_sbuffer.number_of_springs > 0)
    {
        netCDF::NcVar data;

        data = getNcVar("springs");
        checkNDims(data, 2);
        checkDim(data, 0, _sbuffer.number_of_springs);
        checkDim(data, 1, 2);
        data.getVar(_sbuffer.springs);

        data = getNcVar("springsstiffness");
        checkNDims(data, 1);
        checkDim(data, 0, _sbuffer.number_of_springs);
        data.getVar(_sbuffer.springsstiffnesses);

        data = getNcVar("springsequilibrium");
        checkNDims(data, 1);
        checkDim(data, 0, _sbuffer.number_of_springs);
        data.getVar(_sbuffer.springsequilibriums);
    }
}

void NetCDFReader::readNumberOfParticles()
{
    netCDF::NcDim dim = _file->getDim("particle_number");
    size_t nparticles = 0;
    if (dim.isNull())
    {
        logging::die("Dimension \"particle_number\" not present in nc file.");
    }
    else
    {
        nparticles = dim.getSize();
        if (!nparticles)
        {
            logging::die("No particle found in topology");
        }
    }
    _pbuffer.initialize(nparticles);
}

void NetCDFReader::readNumberOfSprings()
{
    netCDF::NcDim dim = _file->getDim("spring_number");
    size_t nsprings = 0;
    if (dim.isNull())
    {
        logging::warning("Dimension \"spring_number\" not present in nc file.");
    }
    else
    {
        nsprings = dim.getSize();
        if (!nsprings)
        {
            logging::warning("No spring found in topology");
        }
    }
    _sbuffer.initialize(nsprings, _pbuffer.number_of_particles);
}

void NetCDFReader::checkNDims(const netCDF::NcVar & var, int ref)
{
    if (var.getDimCount() != ref)
    {
        logging::die("\"%s\" array has %d dimension(s) whereas it should be %d.", var.getName().c_str(),
                     var.getDimCount(), ref);
    }
}

void NetCDFReader::checkDim(const netCDF::NcVar & var, int dimid, size_t size)
{
    if (var.getDim(dimid).getSize() != size)
    {
        logging::die("Dimension %d of array \"%s\" has size %d whereas it should be %d.", dimid, var.getName().c_str(),
                     var.getDim(dimid).getSize(), size);
    }
}

netCDF::NcVar NetCDFReader::getNcVar(const char * varname, bool mandatory)
{
    netCDF::NcVar var = _file->getVar(varname);

    if (var.isNull())
    {
        if (mandatory)
        {
            logging::die("Array \"%s\" not found in topology file.", varname);
        }
        else
        {
            logging::warning("Array \"%s\" not found in topology file.", varname);
        }
    }
    return var;
}
