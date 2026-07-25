
#include "IO/NetCDFWriter.h"
#include "IO/SpnBuffer.h"
#include "Particle.h"
#include "logging.h"

#include <fstream>

using namespace netCDF;
using namespace netCDF::exceptions;

namespace logging = biospring::logging;

NcFile * NetCDFWriter::safeOpenBinary()
{
    NcFile * nc = nullptr;
    try
    {
        nc = new NcFile(_filename, NcFile::replace);
    }
    catch (NcException & e)
    {
        logging::die("cannot open output file '%s': %s", _filename.c_str(), e.what());
    }
    return nc;
}

void NetCDFWriter::writeBinary()
{
    NcFile * nc = safeOpenBinary();

    size_t particlenb = _spn->getNumberOfParticles();
    size_t springnb = _spn->getNumberOfSprings();

    // Defines the dimensions.
    NcDim spatialdim = nc->addDim("spatialdim", 3);
    NcDim pnbdim = nc->addDim("particle_number", particlenb);
    NcDim pnamelendim = nc->addDim("particlename_length", PARTICLE_NAME_LENGTH);
    NcDim chainnamelendim = nc->addDim("chainname_length", CHAIN_NAME_LENGTH);
    NcDim resnamelendim = nc->addDim("resname_length", RESIDUE_NAME_LENGTH);

    std::vector<NcDim> coordsdim(2);
    coordsdim[0] = pnbdim;
    coordsdim[1] = spatialdim;

    std::vector<NcDim> stringdim(2);

    // Defines the variables.
    NcVar coords = nc->addVar("coordinates", ncFloat, coordsdim);
    coords.putAtt("units", "angstrom");
    coords.putAtt("long_name", "Particle coordinates");

    NcVar pids = nc->addVar("particleids", ncInt, pnbdim);
    pids.putAtt("long_name", "Particle ids in source database");

    NcVar charges = nc->addVar("charges", ncFloat, pnbdim);
    charges.putAtt("units", "electron");
    charges.putAtt("long_name", "Particle charge id");

    NcVar radii = nc->addVar("radii", ncFloat, pnbdim);
    radii.putAtt("units", "A");
    radii.putAtt("long_name", "Particle radius");

    NcVar epsilon = nc->addVar("epsilon", ncFloat, pnbdim);
    epsilon.putAtt("units", "kJ.mol-1");
    epsilon.putAtt("long_name", "Particle epsilon for Lennard-Jones");

    NcVar mass = nc->addVar("mass", ncFloat, pnbdim);
    mass.putAtt("units", "Da");
    mass.putAtt("long_name", "Particle mass");

    NcVar surfacc = nc->addVar("surfaceaccessibility", ncFloat, pnbdim);
    surfacc.putAtt("units", "A2 or percent");
    surfacc.putAtt("long_name", "Particle surface accessibility");

    NcVar hscale = nc->addVar("hydrophobicityscale", ncFloat, pnbdim);
    hscale.putAtt("units", "kJ.mol-1");
    hscale.putAtt("long_name", "Particle hydrophobicity scale (transfer energy)");

    NcVar resids = nc->addVar("resids", ncInt, pnbdim);
    resids.putAtt("long_name", "Particle residue id");

    NcVar dynst = nc->addVar("dynamicstate", ncByte, pnbdim);
    dynst.putAtt("long_name", "Particle dynamic state (static 0, dynamic 1)");

    stringdim[0] = pnbdim;
    stringdim[1] = pnamelendim;
    NcVar pnames = nc->addVar("particlenames", ncChar, stringdim);
    pnames.putAtt("long_name", "Particle name");

    stringdim[1] = resnamelendim;
    NcVar resnames = nc->addVar("resnames", ncChar, stringdim);
    resnames.putAtt("long_name", "Particle residue name");

    stringdim[1] = chainnamelendim;
    NcVar chainnames = nc->addVar("chainnames", ncChar, stringdim);
    chainnames.putAtt("long_name", "Particle chain name");

    // Writes particle data to output file.
    ParticleBuffer pbuffer(particlenb);
    pbuffer.bufferize(_spn);

    coords.putVar(pbuffer.coordinates);
    pids.putVar(pbuffer.ids);
    resids.putVar(pbuffer.resids);
    dynst.putVar(pbuffer.dynamic_states);
    radii.putVar(pbuffer.radii);
    mass.putVar(pbuffer.masses);
    epsilon.putVar(pbuffer.epsilons);
    charges.putVar(pbuffer.charges);
    surfacc.putVar(pbuffer.surface_accessibilities);
    hscale.putVar(pbuffer.hscales);
    chainnames.putVar(pbuffer.chainnames);
    pnames.putVar(pbuffer.particlenames);
    resnames.putVar(pbuffer.resnames);

    if (springnb > 0)
    {
        // Defines the dimensions.
        NcDim springdim = nc->addDim("springdim", 2);
        NcDim springnbdim = nc->addDim("spring_number", springnb);
        std::vector<NcDim> springssdim(2);
        springssdim[0] = springnbdim;
        springssdim[1] = springdim;

        // Defines the variables.
        NcVar springs = nc->addVar("springs", ncInt, springssdim);
        springs.putAtt("long_name", "Spring between particles referenced by 2 particle ids");

        NcVar springsstiffness = nc->addVar("springsstiffness", ncFloat, springnbdim);
        springsstiffness.putAtt("long_name", "Spring stiffness");

        NcVar springsequilibrium = nc->addVar("springsequilibrium", ncFloat, springnbdim);
        springsequilibrium.putAtt("long_name", "Spring distance equilibrium");

        NcVar nbofsprings = nc->addVar("nbspringsperparticle", ncInt, pnbdim);
        nbofsprings.putAtt("long_name", "Number of springs per particle");

        // Writes spring data to output file.
        SpringBuffer sbuffer(springnb, particlenb);
        sbuffer.bufferize(_spn);

        springs.putVar(sbuffer.springs);
        springsequilibrium.putVar(sbuffer.springsequilibriums);
        springsstiffness.putVar(sbuffer.springsstiffnesses);
        nbofsprings.putVar(sbuffer.nbofspringsperparticle);
    }

    delete nc;
}

void NetCDFWriter::write()
{
    safeOpen();
    _writeHeaderCDL();
    _writeParticleDataCDL();
    _writeSpringDataCDL();
}

void NetCDFWriter::_writeHeaderCDL()
{
    size_t nbparticles = _spn->getNumberOfParticles();
    size_t nbsprings = _spn->getNumberOfSprings();

    _ostream << "netcdf SpringNetwork {" << std::endl;
    _ostream << std::endl;
    _ostream << "dimensions:" << std::endl;
    _ostream << "\tspatialdim = 3;" << std::endl;
    _ostream << "\tparticle_number = " << nbparticles << ";" << std::endl;
    _ostream << "\tparticlename_length = " << PARTICLE_NAME_LENGTH << ";" << std::endl;
    _ostream << "\tchainname_length = " << CHAIN_NAME_LENGTH << ";" << std::endl;
    _ostream << "\tresname_length = " << RESIDUE_NAME_LENGTH << ";" << std::endl;
    _ostream << std::endl;
    if (nbsprings > 0)
    {
        _ostream << "\tspringdim = 2; " << std::endl;
        _ostream << "\tspring_number = " << nbsprings << ";" << std::endl;
        _ostream << std::endl;
    }

    _ostream << "variables:" << std::endl;
    _ostream << "\tfloat   coordinates(particle_number, spatialdim); " << std::endl;
    _ostream << "\t        coordinates:units = \"angstrom\" ;" << std::endl;
    _ostream << "\t        coordinates:long_name = \"Particle coordinates\";" << std::endl;
    _ostream << std::endl;
    _ostream << "\tint     particleids(particle_number); " << std::endl;
    _ostream << "\t        particleids:long_name = \"Particle ids in source database\";" << std::endl;
    _ostream << std::endl;
    _ostream << "\tchar    particlenames(particle_number,particlename_length); " << std::endl;
    _ostream << "\t        particlenames:long_name = \"Particle name\";" << std::endl;
    _ostream << std::endl;
    _ostream << "\tfloat   charges(particle_number);" << std::endl;
    _ostream << "\t        charges:long_name = \"Particle charge id\";" << std::endl;
    _ostream << "\t        charges:units = \"electron\" ;" << std::endl;
    _ostream << std::endl;
    _ostream << "\tfloat   radii(particle_number);" << std::endl;
    _ostream << "\t        radii:units = \"A\" ;" << std::endl;
    _ostream << "\t        radii:long_name = \"Particle radius\";" << std::endl;
    _ostream << std::endl;
    _ostream << "\tfloat   epsilon(particle_number);" << std::endl;
    _ostream << "\t        epsilon:units = \"kJ.mol-1\" ;" << std::endl;
    _ostream << "\t        epsilon:long_name = \"Particle epsilon for Lennard-Jones\";" << std::endl;
    _ostream << std::endl;
    _ostream << "\tfloat   mass(particle_number);" << std::endl;
    _ostream << "\t        mass:units = \"Da\" ;" << std::endl;
    _ostream << "\t        mass:long_name = \"Particle mass\";" << std::endl;
    _ostream << std::endl;
    _ostream << "\tfloat   surfaceaccessibility(particle_number);" << std::endl;
    _ostream << "\t        surfaceaccessibility:units = \"A2 or percent\" ;" << std::endl;
    _ostream << "\t        surfaceaccessibility:long_name = \"Particle surface accessibility\";" << std::endl;
    _ostream << std::endl;

    _ostream << "\tfloat   hydrophobicityscale(particle_number);" << std::endl;
    _ostream << "\t        hydrophobicityscale:units = \"kJ.mol-1\" ;" << std::endl;
    _ostream << "\t        hydrophobicityscale:long_name = \"Particle hydrophobicity scale (transfer energy)\";"
             << std::endl;
    _ostream << std::endl;

    _ostream << "\tchar    resnames(particle_number,resname_length); " << std::endl;
    _ostream << "\t        resnames:long_name = \"particle residue name\";" << std::endl;
    _ostream << std::endl;
    _ostream << "\tint     resids(particle_number); " << std::endl;
    _ostream << "\t        resids:long_name = \"particle residue id\";" << std::endl;
    _ostream << std::endl;

    _ostream << "\tchar    chainnames(particle_number,chainname_length); " << std::endl;
    _ostream << "\t        chainnames:long_name = \"Chain name \";" << std::endl;
    _ostream << std::endl;

    _ostream << "\tbyte    dynamicstate(particle_number); " << std::endl;
    _ostream << "\t        dynamicstate:long_name = \"particle dynamic state (static 0 or dynamic 1)\";" << std::endl;
    _ostream << std::endl;

    _ostream << "\tint     nbspringsperparticle(particle_number); " << std::endl;
    _ostream << "\t        nbspringsperparticle:long_name = \"Number of springs per particle\";" << std::endl;
    _ostream << std::endl;

    if (nbsprings > 0)
    {
        _ostream << "\tint     springs(spring_number,springdim); " << std::endl;
        _ostream << "\t        springs:long_name = \"Spring between particle referenced by 2 particle ids\"; "
                 << std::endl;
        _ostream << std::endl;
        _ostream << "\tfloat   springsstiffness(spring_number); " << std::endl;
        _ostream << "\t        springsstiffness:long_name = \"Spring stiffness\";" << std::endl;

        _ostream << "\tfloat   springsequilibrium(spring_number); " << std::endl;
        _ostream << "\t        springsequilibrium:long_name = \"Spring distance equilibrium\";" << std::endl;
    }
}

void NetCDFWriter::_writeParticleDataCDL()
{
    size_t nbparticles = _spn->getNumberOfParticles();

    _ostream << "data:" << std::endl;

    _ostream << "\tparticleids = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << p.getExtid();
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }

    _ostream << "\tcoordinates = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        Vector3f pos = p.getPosition();
        _ostream << pos.getX() << ", " << pos.getY() << ", " << pos.getZ();
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }
    _ostream << "\tcharges = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << p.getCharge();
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }
    _ostream << "\tradii = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << p.getRadius();
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }

    _ostream << "\tepsilon = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << p.getEpsilon();
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }

    _ostream << "\tmass = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << p.getMass();
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }

    _ostream << "\tsurfaceaccessibility = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << p.getSolventAccessibilitySurface();
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }

    _ostream << "\thydrophobicityscale = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << p.getTransferEnergyByAccessibleSurface();
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }

    _ostream << "\tparticlenames = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << "\"" << p.getName() << "\"";
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }

    _ostream << "\tresnames = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << "\"" << p.getResName() << "\"";
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }

    _ostream << "\tresids = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << p.getResId();
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }
    _ostream << "\tchainnames = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << "\"" << p.getChainName() << "\"";
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }

    _ostream << "\tdynamicstate = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        if (p.isDynamic())
            _ostream << 1;
        else
            _ostream << 0;
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }
    _ostream << "\tnbspringsperparticle = " << std::endl;
    for (size_t i = 0; i < nbparticles; i++)
    {
        const biospring::spn::Particle & p = _spn->getParticle(i);
        _ostream << p.getNumberOfSprings();
        if (i == (nbparticles - 1))
            _ostream << ";" << std::endl;
        else
            _ostream << "," << std::endl;
    }
}

void NetCDFWriter::_writeSpringDataCDL()
{
    size_t nbsprings = _spn->getNumberOfSprings();

    if (nbsprings > 0)
    {
        _ostream << "\tsprings = " << std::endl;
        for (size_t i = 0; i < nbsprings; i++)
        {
            const biospring::spn::Particle & p1 = _spn->getSpring(i).getParticle1();
            const biospring::spn::Particle & p2 = _spn->getSpring(i).getParticle2();
            _ostream << p1.getId() << ", " << p2.getId();
            if (i == (nbsprings - 1))
                _ostream << ";" << std::endl;
            else
                _ostream << "," << std::endl;
        }
        _ostream << "\tspringsstiffness = " << std::endl;
        for (size_t i = 0; i < nbsprings; i++)
        {
            const biospring::spn::Spring & s = _spn->getSpring(i);
            _ostream << static_cast<float>(s.getStiffness());
            if (i == (nbsprings - 1))
                _ostream << ";" << std::endl;
            else
                _ostream << "," << std::endl;
        }
        _ostream << "\tspringsequilibrium = " << std::endl;
        for (size_t i = 0; i < nbsprings; i++)
        {
            const biospring::spn::Spring & s = _spn->getSpring(i);
            _ostream << static_cast<float>(s.getEquilibrium());
            if (i == (nbsprings - 1))
                _ostream << ";" << std::endl;
            else
                _ostream << "," << std::endl;
        }
    }
    _ostream << "}" << std::endl;
}
