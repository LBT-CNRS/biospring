#ifdef MDDRIVER_SUPPORT

#include "interactor/mddriver/InteractorMDDriver.h"
#include "interactor/mddriver/CustomData.h"
#include "SpringNetwork.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "forcefield/constants.hpp"


namespace biospring
{
namespace interactor
{

InteractorMDDriver::InteractorMDDriver()
{
	_IMDdebug  = 2;
	strcpy(_IMDlogfilename,"");
	_IMDlog = NULL;
	_IMDwait           = 1;
	_IMDport           = 8888;
	_IMDmode = 1; //server
	_IMDforcescale =1.0;
	_IMDenergies.tstep  = 0.0;
	_IMDenergies.T      = 0.0;
	_IMDenergies.Etot   = 0.0;
	_IMDenergies.Epot   = 0.0;
	_IMDenergies.Evdw   = 0.0;
	_IMDenergies.Eelec  = 0.0;
	_IMDenergies.Ebond  = 0.0;
	_IMDenergies.Eangle = 0.0;
	_IMDenergies.Edihe  = 0.0;
	_IMDenergies.Eimpr  = 0.0;
	_sleepDuration = 1000;  // For example, set sleep duration to 1000 microseconds for InteractorMDDriver.

	auto initializeGrid = [](IMDGrid& grid)
	{
		grid.tstep = 0; //!< integer timestep index

		grid.Xorigin = 0.0;   //!< Grid origin X (diagonal origin)
		grid.Yorigin = 0.0;   //!< Grid origin Y (diagonal origin)
		grid.Zorigin = 0.0;   //!< Grid origin Z (diagonal origin)

		grid.Xend = 0.0;   //!< Grid end X (diagonal origin)
		grid.Yend = 0.0;   //!< Grid end Y (diagonal origin)
		grid.Zend = 0.0;   //!< Grid end Z (diagonal origin)

		grid.XdirectionX = 0.0;   //!< X component X vector direction of the grid frame (right hand frame)
		grid.YdirectionX = 0.0;   //!< Y component X vector direction of the grid frame (right hand frame)
		grid.ZdirectionX = 0.0;   //!< Z component X vector direction of the grid frame (right hand frame)

		grid.XdirectionY = 0.0;   //!< X component Y vector direction of the grid frame (right hand frame)
		grid.YdirectionY = 0.0;   //!< Y component Y vector direction of the grid frame (right hand frame)
		grid.ZdirectionY = 0.0;   //!< Z component Y vector direction of the grid frame (right hand frame)

		grid.XdirectionZ = 0.0;   //!< X component Z vector direction of the grid frame (right hand frame)
		grid.YdirectionZ = 0.0;   //!< Y component Z vector direction of the grid frame (right hand frame)
		grid.ZdirectionZ = 0.0;   //!< Z component Z vector direction of the grid frame (right hand frame)

		grid.nbcellx = 0;         //!< Number of cells along X axis
		grid.nbcelly = 0;         //!< Number of cells along Y axis
		grid.nbcellz = 0;         //!< Number of cells along Z axis

		grid.sizecellx = 0.0;     //!< Size of the cell on x axis
		grid.sizecelly = 0.0;     //!< Size of the cell on y axis
		grid.sizecellz = 0.0;     //!< Size of the cell on z axis
	};

	initializeGrid(_IMDpotentialGrid);
	initializeGrid(_IMDdensityGrid);
}

InteractorMDDriver::~InteractorMDDriver()
{

}

void InteractorMDDriver::startInteractionThread()
{
	initializeSystemState();
	Interactor::startInteractionThread();
}

/**
 * @brief Initializes the system state specific to the InteractorMDDriver.
 * 
 * This method extends the base initialization from the `Interactor` class by setting up
 * variables and resources specific to the MDDriver interaction. It prepares the system 
 * for subsequent MD driver operations and interactions.
 */
void InteractorMDDriver::initializeSystemState()
{
	Interactor::initializeSystemState(); // Set _nbpositions int defined in Interactor.

	initializeDataManager();

	CustomData::initializeDataManager(this);

	_IMDforcescale=_IMDforcescale/biospring::forcefield::AVOGADRO_NUMBER;//kcal.A-1
	_IMDforcescale=_IMDforcescale*biospring::forcefield::KCAL_TO_KJOULE*1.0E3;//J.A-1
	_IMDforcescale=_IMDforcescale/biospring::forcefield::ANGSTROM_TO_METER;//J.m-1 ou N ou kg.m.s-2
	_IMDforcescale=_IMDforcescale*biospring::forcefield::NEWTON_TO_DALTON_ANGSTROM_PER_FEMTOSECOND_2;//Da.A.fs-2
}

void InteractorMDDriver::initializeDataManager()
{
	// floatManager and intManager : data you can store here with defined size
	// refFloatManager and refIntManager : pointers to data you recieve from client

	// Position sent to client
	floatManager.add("positions", getNbPositions() * 3);

	// Store recieved forces here with fixed maximum size
	floatManager.add("forces", getNbPositions() * 3);
	// Particle ids of recieved forces from client with varying size
	refIntManager.add("particleforceids");
	// Recieved forces of particles from client with varying size
	refFloatManager.add("particleforces");
}


void InteractorMDDriver::setupIMDInteractions( InteractorMDDriver * imdl)
{
	static int fp_comm = -1;
	if ( fp_comm == -1)
	{
		biospring::spn::SpringNetwork * spn=imdl->getSpringNetwork();
		if(spn!=NULL && imdl->_IMDwait!=0)
		{
			spn->setPause(true);
		}
		imdl->_IMDlog = IIMD_init( "", &(imdl->_IMDmode),&(imdl->_IMDwait),&(imdl->_IMDport), &(imdl->_IMDdebug),imdl->_IMDlogfilename );

		IIMD_probeconnection();
		IIMD_treatprotocol();

		// Set IMDGrid at initialization
		if (spn->isElectrostaticEnabled())
		{
			biospring::grid::PotentialGrid potentialGrid = spn->getPotentialGrid();
			updateGridFromSource(imdl->_IMDpotentialGrid, potentialGrid);
		}

		if (spn->isDensityGridEnabled())
		{
			biospring::grid::PotentialGrid densityGrid = spn->getDensityGrid();
			updateGridFromSource(imdl->_IMDdensityGrid, densityGrid);
		}


		if(spn!=NULL && imdl->_IMDwait!=0)
		{
			spn->setPause(false);
		}
		fp_comm = 1;
		imdl->_isRunning.store(true, std::memory_order_release);
	}
}

void InteractorMDDriver::updateGridFromSource(IMDGrid& targetGrid, const biospring::grid::PotentialGrid& sourceGrid)
{
    targetGrid.tstep = 0; //!< integer timestep index

    targetGrid.Xorigin = sourceGrid.boundaries().origin_x();   //!< Grid origin X (diagonal origin)
    targetGrid.Yorigin = sourceGrid.boundaries().origin_y();   //!< Grid origin Y (diagonal origin)
    targetGrid.Zorigin = sourceGrid.boundaries().origin_z();   //!< Grid origin Z (diagonal origin)

    targetGrid.Xend = sourceGrid.boundaries().max_x();   //!< Grid end X (diagonal origin)
    targetGrid.Yend = sourceGrid.boundaries().max_y();   //!< Grid end Y (diagonal origin)
    targetGrid.Zend = sourceGrid.boundaries().max_z();   //!< Grid end Z (diagonal origin)

    targetGrid.XdirectionX = 1.0;   //!< X component X vector direction of the grid frame (right hand frame)
    targetGrid.YdirectionX = 1.0;   //!< Y component X vector direction of the grid frame (right hand frame)
    targetGrid.ZdirectionX = 1.0;   //!< Z component X vector direction of the grid frame (right hand frame)

    targetGrid.XdirectionY = 2.0;   //!< X component Y vector direction of the grid frame (right hand frame)
    targetGrid.YdirectionY = 2.0;   //!< Y component Y vector direction of the grid frame (right hand frame)
    targetGrid.ZdirectionY = 2.0;   //!< Z component Y vector direction of the grid frame (right hand frame)

    targetGrid.XdirectionZ = 3.0;   //!< X component Z vector direction of the grid frame (right hand frame)
    targetGrid.YdirectionZ = 3.0;   //!< Y component Z vector direction of the grid frame (right hand frame)
    targetGrid.ZdirectionZ = 3.0;   //!< Z component Z vector direction of the grid frame (right hand frame)

    targetGrid.nbcellx = sourceGrid.shape()[0];         //!< Number of cells along X axis
    targetGrid.nbcelly = sourceGrid.shape()[1];         //!< Number of cells along Y axis
    targetGrid.nbcellz = sourceGrid.shape()[2];         //!< Number of cells along Z axis

    targetGrid.sizecellx = sourceGrid.cell_size()[0];     //!< Size of the cell on x axis
    targetGrid.sizecelly = sourceGrid.cell_size()[1];     //!< Size of the cell on y axis
    targetGrid.sizecellz = sourceGrid.cell_size()[2];     //!< Size of the cell on z axis
}

void InteractorMDDriver::sendPotentialGrid(InteractorMDDriver * imdl)
{
	handleIMDWorkflow(imdl);
    IIMD_send_grid(&(imdl->_IMDpotentialGrid));
}

void InteractorMDDriver::sendDensityGrid(InteractorMDDriver * imdl)
{
	handleIMDWorkflow(imdl);
    IIMD_send_grid(&(imdl->_IMDdensityGrid));
}

int InteractorMDDriver::processIMDInteractions(InteractorMDDriver * imdl) {
    int ret = 0;
    std::lock_guard<std::mutex> lock(imdl->mutex);

    // This thread loops far faster than the simulation steps (once per
    // _sleepDuration, 1 ms by default), so without this it re-sends the same
    // buffer many times per frame. Send only what syncSystemStateData has
    // actually refreshed since last time -- which is what makes the
    // transmission rate a rate of frames and not just of memory copies.
    const unsigned long long version = imdl->_stateVersion.load(std::memory_order_acquire);
    const bool fresh = (version != imdl->_lastSentVersion);

    if (fresh)
    {
        // Send positions
        handleIMDWorkflow(imdl);
        float* positions = imdl->floatManager.get("positions").getData();
        IIMD_send_coords(&(imdl->_nbpositions), positions);

        // Send energies
        handleIMDWorkflow(imdl);
        IIMD_send_energies(&(imdl->_IMDenergies));

        imdl->_lastSentVersion = version;
    }

    // Forces and events are read every loop whatever the frame rate: a pull
    // from the client must not wait for the next frame to be produced.
    // Get forces
    handleIMDWorkflow(imdl);
    int nbforces;
    int* particleforceids = imdl->refIntManager.get("particleforceids").getPointer();
    float* particleforces = imdl->refFloatManager.get("particleforces").getPointer();
    IIMD_get_forces(&nbforces, &particleforceids, &particleforces);

    // Update local array with new forces
    updateForces(imdl, nbforces, particleforceids, particleforces);

    handleIMDEvents(imdl);

	CustomData::processCustomIntData(imdl);
	CustomData::processCustomFloatData(imdl);
	
	
    return ret;
}


void  InteractorMDDriver::handleIMDWorkflow(InteractorMDDriver * imdl)
{
    setupIMDInteractions(imdl);
    handleIMDConnection();
}

void InteractorMDDriver::handleIMDConnection()
{
	IIMD_probeconnection();
    IIMD_treatprotocol();
}

void InteractorMDDriver::handleIMDEvents(InteractorMDDriver * imdl)
{
    switch(imd_event) {
        case IMD_KILL:
			imdl->terminateInteraction();
            imd_event = -1;
            break;
        case IMD_TRATE:
            // The client is telling us how often it wants frames. Ignoring it
            // meant every client got one per step whatever it asked for, and
            // the simulation paid a full state copy for each.
            if (imd_value > 0)
                imdl->setTransmissionRate(static_cast<unsigned>(imd_value));
            imd_event = -1;
            break;
        case IMD_PAUSE:
            auto spn = imdl->getSpringNetwork();
            if(spn) {
                spn->setPause(!spn->getPause());
            }
            imd_event = -1;
    }
}

void InteractorMDDriver::updateForces(InteractorMDDriver * imdl, int nbforces, int* particleforceids, float* particleforces)
{
    float* forces = imdl->floatManager.get("forces").getData();
    imdl->floatManager.get("forces").clearData();

    float scale = imdl->_IMDforcescale;
    for(int i = 0; i < nbforces; ++i) {
        int index = particleforceids[i] * 3;
        forces[index]   = particleforces[i*3] * scale;
        forces[index+1] = particleforces[i*3+1] * scale;
        forces[index+2] = particleforces[i*3+2] * scale;
    }
}

void InteractorMDDriver::syncSystemStateData()
{
	// This runs on the simulation thread, once per step, and copies every
	// particle's position into the buffer the interaction thread sends. On a
	// 17338-particle system it measured 7.8 ms against 6.4 ms for the whole
	// rest of the step -- the transport costing more than the physics.
	//
	// The interaction thread consumes at its own pace (1 kHz by default) and
	// no viewer redraws faster than the display, so producing a fresh copy
	// every step is only useful when the client actually asked for every step.
	// That is what IMD_TRATE is for; it used to be received and dropped.
	if(_springnetwork!=NULL)
	{
		const unsigned rate = getTransmissionRate();
		if (rate > 1 && (static_cast<unsigned>(_springnetwork->getNbIterations()) % rate) != 0)
			return;
	}

	Interactor::syncSystemStateData();
	_stateVersion.fetch_add(1, std::memory_order_release);
	if(_springnetwork!=NULL)
	{
		_IMDenergies.tstep  = _springnetwork->getNbIterations(); //!< integer timestep index
		_IMDenergies.T = 0.0;          											//!< Temperature in degrees Kelvin
		_IMDenergies.Etot = _IMDenergies.Eelec+_IMDenergies.Evdw+_IMDenergies.Ebond;  //!< Total energy, in Kcal/mol
		_IMDenergies.Epot = 0.0;       //!< Potential energy, in Kcal/mol
		_IMDenergies.Evdw = _springnetwork->getStericEnergy();       //!< Van der Waals energy, in Kcal/mol
		_IMDenergies.Eelec = _springnetwork->getElectrostaticEnergy();      //!< Electrostatic energy, in Kcal/mol
		_IMDenergies.Ebond = _springnetwork->getSpringEnergy();      //!< Bond energy, Kcal/mol
		_IMDenergies.Eangle = 0.0;     //!< Angle energy, Kcal/mol
		_IMDenergies.Edihe = 0.0;      //!< Dihedral energy, Kcal/mol
		_IMDenergies.Eimpr = 0.0;      //!< Improper energy, Kcal/mol
	}
}

// Update the `_positions` array in the base Interactor class and the force 
// of the given particle based on MDDriver-related computations.
void InteractorMDDriver::syncParticleStateData(unsigned index)
{
	biospring::spn::Particle particle = _springnetwork->getParticle(index);
	float position[3];
	_springnetwork->getParticlePosition(index, position);
	// Update the `positions` array for the given particle.
	float* positions = floatManager.get("positions").getData();
	memcpy(&(positions[index * 3]), position, sizeof(float) * 3);
	// Update the force for the given particle.
	float* forces = floatManager.get("forces").getData();
	_springnetwork->setForce(index, &(forces[index * 3]));

	CustomData::syncParticleStateData(this, index);
}

} // namespace interactor
} // namespace biospring

#endif
