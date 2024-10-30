#ifdef MDDRIVER_SUPPORT

#ifndef _INTERACTORMDDRIVER_H_
#define _INTERACTORMDDRIVER_H_

#include <stdio.h>
#include <stdlib.h>

#include "interactor/Interactor.h"
#include "imd_interface.h"

#include "forcefield/constants.hpp"

#include <iostream>

#include "grid/PotentialGrid.hpp"

#define FILTERSIZE 256
#define LOGFILENAMESIZE 256

namespace biospring
{
namespace interactor
{

class InteractorMDDriver : public Interactor
{
	public : 
		InteractorMDDriver();
		virtual ~InteractorMDDriver();

		inline void setPort(unsigned port) { _IMDport = port; }
		inline int getPort() { return _IMDport; }
		inline void setWait(unsigned wait) { _IMDwait = wait; }
		inline void setDebug(unsigned debug) { _IMDdebug = debug; }
		inline void setLog(const char* logfilename) { strcpy(_IMDlogfilename, logfilename); }
		inline void setForceScale(float forcescale) { _IMDforcescale = forcescale; }
		
		virtual void startInteractionThread() override;
		virtual bool continueInteractionThread() override { return _isRunning; }
		virtual void stopInteractionThread() override { _isRunning = false; }

		virtual void syncSystemStateData() override;

		// Managers for float/int data assigned in the server side (BioSpring)
		DataArrayManager<float> floatManager;
		DataArrayManager<int> intManager;

		// Managers for float/int data assigned in the client side.
		ExternalDataRefManager<float> refFloatManager;
		ExternalDataRefManager<int> refIntManager;

		static void handleIMDWorkflow(InteractorMDDriver * imdl);
		static void handleIMDConnection();
		static void handleIMDEvents(InteractorMDDriver * imdl);

		static void updateForces(InteractorMDDriver * imdl, int nbforces, int* particleforceids, float* particleforces);

	protected : 

		char _IMDlogfilename[LOGFILENAMESIZE];
		int    _IMDdebug ;
		FILE * _IMDlog;
		int _IMDmode ; 
		int _IMDwait ;
		int _IMDport ;
		float     _IMDforcescale ;
		IMDEnergies _IMDenergies;
		int _nbforces;
		int * _particleforceids; 
		float * _particleforces;

		IMDGrid _IMDpotentialGrid;
		IMDGrid _IMDdensityGrid;
		static void updateGridFromSource(IMDGrid& targetGrid, const biospring::grid::PotentialGrid& sourceGrid);
		static void sendPotentialGrid(InteractorMDDriver * imdl);
		static void sendDensityGrid(InteractorMDDriver * imdl);
		static void sendDensityGrid();
			
		virtual void setupInteraction() override { setupIMDInteractions(this); }
    	virtual void processInteractions() override { processIMDInteractions(this); }
    	virtual void terminateInteraction() override { IIMD_terminate(); }

		static void setupIMDInteractions( InteractorMDDriver * imdl);
		static int processIMDInteractions( InteractorMDDriver * imdl);
		
		virtual void initializeSystemState() override;

		virtual void initializeDataManager() override;

		virtual void syncParticleStateData(unsigned index) override;

        bool _isRunning;
		
	};
} // namespace interactor
} // namespace biospring
#endif

#endif
