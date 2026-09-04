#include "NetCDFReader.h"
#include "Vector3f.h"
#include "logging.h"
#include "topology.hpp"

using namespace biospring;

void NetCDFReader::addSpringsToSpn()
{
    for (size_t i = 0; i < _sbuffer.number_of_springs; ++i)
    {
        // springs[][] stays a plain int: it mirrors the NetCDF "springs"
        // variable's on-disk type (nc_INT), read as-is by the netCDF API.
        _topology.add_spring(static_cast<size_t>(_sbuffer.springs[i][0]), static_cast<size_t>(_sbuffer.springs[i][1]),
                             _sbuffer.springsequilibriums[i], _sbuffer.springsstiffnesses[i]);
    }
}

void NetCDFReader::addStretchSpringsToSpn()
{
    for (size_t i = 0; i < _stretchbuffer.number_of_springs; ++i)
        _topology.add_stretch_spring(_topology.get_particle(static_cast<size_t>(_stretchbuffer.springs[i][0])),
                                     _topology.get_particle(static_cast<size_t>(_stretchbuffer.springs[i][1])),
                                     _stretchbuffer.springsequilibriums[i], _stretchbuffer.springsstiffnesses[i]);
}

void NetCDFReader::addBendSpringsToSpn()
{
    for (size_t i = 0; i < _bendbuffer.number_of_springs; ++i)
        _topology.add_bend_spring(_topology.get_particle(static_cast<size_t>(_bendbuffer.springs[i][0])),
                                  _topology.get_particle(static_cast<size_t>(_bendbuffer.springs[i][1])),
                                  _bendbuffer.springsequilibriums[i], _bendbuffer.springsstiffnesses[i]);
}

void NetCDFReader::addDihedralSpringsToSpn(unsigned family)
{
    const DihedralSpringBuffer & buffer = _dihedralbuffers[family];
    for (size_t i = 0; i < buffer.number_of_springs; ++i)
        _topology
            .add_dihedral_spring(family, _topology.get_particle(static_cast<size_t>(buffer.springs[i][0])),
                                 _topology.get_particle(static_cast<size_t>(buffer.springs[i][1])),
                                 buffer.springsequilibriums[i], buffer.springsstiffnesses[i])
            .set_dc_offset(buffer.springsdcoffsets[i]);
}

void NetCDFReader::addGhostParticlesToSpn()
{
    for (size_t i = 0; i < _ghostparticlebuffer.number_of_ghostparticles; ++i)
        _topology.register_ghost_particle(
            static_cast<size_t>(_ghostparticlebuffer.ownindices[i]),
            static_cast<size_t>(_ghostparticlebuffer.anchorindices[i][0]),
            static_cast<size_t>(_ghostparticlebuffer.anchorindices[i][1]),
            static_cast<size_t>(_ghostparticlebuffer.anchorindices[i][2]), _ghostparticlebuffer.rs[i],
            _ghostparticlebuffer.thetas[i], _ghostparticlebuffer.deltas[i],
            static_cast<unsigned>(_ghostparticlebuffer.placements[i]));
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

        // Ghost (virtual-site) particles are deliberately massless and static
        // -- see spn::GhostParticle -- so this guard against an accidental
        // zero mass only applies to dynamic particles, whose mass actually
        // matters for integration.
        if (p.properties().is_dynamic() && std::abs(p.properties().mass()) < 1e-6)
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
        auto imp = topology::IMPProperties::build();
        if (_pbuffer.surface_accessibilities)
            imp.solvent_accessible_surface(_pbuffer.surface_accessibilities[i]);
        if (_pbuffer.hscales)
            imp.transfert_energy_by_accessible_surface(_pbuffer.hscales[i]);
        p.properties().set_imp(imp);

        // Hydrophobicity is a ParticleProperties member in its own right, not an
        // IMP property: it feeds the pairwise hydrophobic force, not IMPALA.
        if (_pbuffer.hydrophobicities)
            p.properties().set_hydrophobicity(_pbuffer.hydrophobicities[i]);

        _topology.add_particle(p);
    }
}

void NetCDFReader::read()
{
    try
    {
        _file = std::make_unique<netCDF::NcFile>(_filename, netCDF::NcFile::read);
        readParticles();
        addParticlesToSpn();

        readSprings();
        addSpringsToSpn();

        readDihedralSpringGroup("stretch", _stretchbuffer);
        addStretchSpringsToSpn();
        readDihedralSpringGroup("bend", _bendbuffer);
        addBendSpringsToSpn();

        for (unsigned family = 0; family < biospring::spn::SpringNetwork::DIHEDRAL_FAMILY_COUNT; ++family)
        {
            readDihedralSpringGroup(biospring::spn::SpringNetwork::DIHEDRAL_FAMILY_NAMES[family],
                                    _dihedralbuffers[family]);
            addDihedralSpringsToSpn(family);
        }

        // Ghost particles already exist as regular particles at this point
        // (added by addParticlesToSpn() above, in file order) -- this only
        // registers which ones are ghosts and their anchor/placement info,
        // so it must run after addParticlesToSpn() but has no other ordering
        // constraint relative to the spring/dihedral steps above.
        readGhostParticles();
        addGhostParticlesToSpn();

        // Axis substituent lists, which only name already-added particles --
        // so, like the ghosts, this needs addParticlesToSpn() to have run and
        // nothing else. The ordering that does matter is on the network side,
        // where the lists attach to axes the ghosts create (see
        // Topology::to_spring_network).
        readDihedralAxes();
        addDihedralAxesToSpn();
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

    // Per-particle hydrophobicity for the pairwise hydrophobic FORCE term
    // (Particle::getHydrophobicity()). Distinct from `hydrophobicityscale`
    // above, which is the IMPALA transfer energy. Optional, so a file written
    // before this variable existed still reads, with hydrophobicity left at 0.
    data = getNcVar("hydrophobicity", false);
    if (not data.isNull())
    {
        checkNDims(data, 1);
        checkDim(data, 0, _pbuffer.number_of_particles);
        data.getVar(_pbuffer.hydrophobicities);
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

void NetCDFReader::readDihedralSpringGroup(const char * prefix, DihedralSpringBuffer & buffer)
{
    // A dihedral ghost-spring family's dimension is entirely optional: an
    // older .nc file (written before dihedral support existed), or one
    // where this particular family happened to be empty, simply doesn't
    // have it -- silently leave `buffer` at zero springs rather than warn
    // (unlike the real springs' dimension, whose absence is unusual enough
    // to warrant a warning).
    const std::string number_dim = std::string(prefix) + "_number";
    netCDF::NcDim dim = _file->getDim(number_dim);
    if (dim.isNull())
    {
        buffer.initialize(0);
        return;
    }

    const size_t n = dim.getSize();
    buffer.initialize(n);
    if (n == 0)
        return;

    netCDF::NcVar data;

    data = getNcVar((std::string(prefix) + "springs").c_str());
    checkNDims(data, 2);
    checkDim(data, 0, n);
    checkDim(data, 1, 2);
    data.getVar(buffer.springs);

    data = getNcVar((std::string(prefix) + "springsstiffness").c_str());
    checkNDims(data, 1);
    checkDim(data, 0, n);
    data.getVar(buffer.springsstiffnesses);

    data = getNcVar((std::string(prefix) + "springsequilibrium").c_str());
    checkNDims(data, 1);
    checkDim(data, 0, n);
    data.getVar(buffer.springsequilibriums);

    // Optional: a .nc file written before this correction existed simply
    // doesn't have it -- buffer.springsdcoffsets stays zero-initialized
    // (initialize() already zeroes it), so those springs just report their
    // raw (uncorrected) energy, same as before this feature existed.
    data = getNcVar((std::string(prefix) + "springsdcoffset").c_str(), false);
    if (not data.isNull())
    {
        checkNDims(data, 1);
        checkDim(data, 0, n);
        data.getVar(buffer.springsdcoffsets);
    }
}

void NetCDFReader::readGhostParticles()
{
    // Entirely optional, same convention as readDihedralSpringGroup: an
    // older .nc file, or one with no ghost particles at write time, simply
    // doesn't have this dimension -- silently leave the buffer empty.
    netCDF::NcDim dim = _file->getDim("ghostparticle_number");
    if (dim.isNull())
    {
        _ghostparticlebuffer.initialize(0);
        return;
    }

    const size_t n = dim.getSize();
    _ghostparticlebuffer.initialize(n);
    if (n == 0)
        return;

    netCDF::NcVar data;

    data = getNcVar("ghostparticleownindex");
    checkNDims(data, 1);
    checkDim(data, 0, n);
    data.getVar(_ghostparticlebuffer.ownindices);

    data = getNcVar("ghostparticleanchorindices");
    checkNDims(data, 2);
    checkDim(data, 0, n);
    checkDim(data, 1, 3);
    data.getVar(_ghostparticlebuffer.anchorindices);

    data = getNcVar("ghostparticler");
    checkNDims(data, 1);
    checkDim(data, 0, n);
    data.getVar(_ghostparticlebuffer.rs);

    data = getNcVar("ghostparticletheta");
    checkNDims(data, 1);
    checkDim(data, 0, n);
    data.getVar(_ghostparticlebuffer.thetas);

    data = getNcVar("ghostparticledelta");
    checkNDims(data, 1);
    checkDim(data, 0, n);
    data.getVar(_ghostparticlebuffer.deltas);

    data = getNcVar("ghostparticleplacement");
    checkNDims(data, 1);
    checkDim(data, 0, n);
    data.getVar(_ghostparticlebuffer.placements);
}

void NetCDFReader::readDihedralAxes()
{
    // Optional in the same way and for the same reasons as the ghosts above,
    // with one extra case: a .nc written before DIHEDRALAXIS records existed
    // has ghosts but no axis lists. dihedral.distributetorque then has
    // nothing to act on and says so at setup, rather than the file failing to
    // load over a setting it may not even use.
    netCDF::NcDim dim = _file->getDim("dihedralaxis_number");
    if (dim.isNull())
    {
        _dihedralaxisbuffer.initialize(0, 0);
        return;
    }

    const size_t n = dim.getSize();
    netCDF::NcDim atomdim = _file->getDim("dihedralaxisatom_number");
    const size_t natoms = atomdim.isNull() ? 0 : atomdim.getSize();
    _dihedralaxisbuffer.initialize(n, natoms);
    if (n == 0)
        return;

    netCDF::NcVar data;

    data = getNcVar("dihedralaxisanchorindices");
    checkNDims(data, 2);
    checkDim(data, 0, n);
    checkDim(data, 1, 2);
    data.getVar(_dihedralaxisbuffer.anchorindices);

    data = getNcVar("dihedralaxissubstituentcounts");
    checkNDims(data, 2);
    checkDim(data, 0, n);
    checkDim(data, 1, 2);
    data.getVar(_dihedralaxisbuffer.counts);

    if (natoms > 0)
    {
        data = getNcVar("dihedralaxissubstituents");
        checkNDims(data, 1);
        checkDim(data, 0, natoms);
        data.getVar(_dihedralaxisbuffer.atoms);

        // Optional on its own: a .nc written before the weights existed has
        // the substituent lists but no shares, and equal shares are the right
        // reading of that. -1 is the same "none given" marker the .bi.ff's
        // unweighted form produces.
        data = getNcVar("dihedralaxissubstituentweights", /*mandatory=*/false);
        if (!data.isNull())
        {
            checkNDims(data, 1);
            checkDim(data, 0, natoms);
            data.getVar(_dihedralaxisbuffer.weights);
        }
        else
            for (size_t i = 0; i < natoms; ++i)
                _dihedralaxisbuffer.weights[i] = -1.0f;
    }
}

void NetCDFReader::addDihedralAxesToSpn()
{
    // The counts slice the flat atom list, so the offset has to run forward
    // across axes: axis i starts where axis i-1 ended. A count that overruns
    // the list would mean a truncated or hand-edited file, so it is checked
    // rather than trusted -- reading past the end would be silent corruption.
    size_t offset = 0;
    for (size_t i = 0; i < _dihedralaxisbuffer.number_of_axes; ++i)
    {
        const size_t countB = static_cast<size_t>(_dihedralaxisbuffer.counts[i][0]);
        const size_t countC = static_cast<size_t>(_dihedralaxisbuffer.counts[i][1]);
        if (offset + countB + countC > _dihedralaxisbuffer.number_of_atoms)
        {
            logging::warning("Torsion axis %zu declares %zu substituents but only %zu remain in the file: the "
                             "dihedral axis lists are inconsistent and are dropped.",
                             i, countB + countC, _dihedralaxisbuffer.number_of_atoms - offset);
            return;
        }

        std::vector<size_t> b_indices, c_indices;
        std::vector<double> b_weights, c_weights;
        for (size_t k = 0; k < countB; ++k)
        {
            b_indices.push_back(static_cast<size_t>(_dihedralaxisbuffer.atoms[offset + k]));
            b_weights.push_back(_dihedralaxisbuffer.weights[offset + k]);
        }
        for (size_t k = 0; k < countC; ++k)
        {
            c_indices.push_back(static_cast<size_t>(_dihedralaxisbuffer.atoms[offset + countB + k]));
            c_weights.push_back(_dihedralaxisbuffer.weights[offset + countB + k]);
        }
        offset += countB + countC;

        _topology.register_dihedral_axis(static_cast<size_t>(_dihedralaxisbuffer.anchorindices[i][0]),
                                         static_cast<size_t>(_dihedralaxisbuffer.anchorindices[i][1]), b_indices,
                                         c_indices, b_weights, c_weights);
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
