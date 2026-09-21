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
		virtual bool continueInteractionThread() override { return _isRunning.load(std::memory_order_acquire); }
		virtual void stopInteractionThread() override { _isRunning.store(false, std::memory_order_release); }

		virtual void syncSystemStateData() override;

		// The arrays the per-particle sync writes into, resolved once per step
		// rather than once per particle. DataArrayManager::get() is a lookup in
		// an unordered_map keyed by std::string, and the sync used to run six
		// of them for every particle of every step -- 223 200 string hashes per
		// step on a 37 200 bead capsid, for six pointers that cannot move: the
		// arrays are added in initializeDataManager() and never again.
		struct SyncTargets
			{
			float * positions = nullptr;
			float * forces = nullptr;
			float * sasa = nullptr;
			float * insertionforce = nullptr;   // "imsf"
			float * transferenergy = nullptr;   // "trimp"
			float * impala = nullptr;
			};
		SyncTargets synctargets;
		void resolveSyncTargets();

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

		
	};
} // namespace interactor
} // namespace biospring
#endif

#endif
