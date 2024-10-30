#include "interactor/Interactor.h"
#include "SpringNetwork.h"

Interactor::Interactor() : _springnetwork(nullptr), _nbpositions(0), _sleepDuration(1000), _isRunning(false) {}


void Interactor::initializeSystemState()
{
    if (_springnetwork != NULL)
    {
        _setNbPositions(_springnetwork->getNumberOfParticles());
    }
}

void Interactor::initializeDataManager() {}

void Interactor::syncSystemStateData()
{
    if (this->_springnetwork != NULL)
    {
        unsigned nbparticles = this->_springnetwork->getNumberOfParticles();
        for (unsigned j = 0; j < nbparticles; j++)
        {
            this->syncParticleStateData(j);
        }
    }
}


void Interactor::startInteractionThread()
{
	_isRunning = true;
	int rc;
    _thread = thread(runthread,(void *) this);
}


void* Interactor::runthread(void* userdata)
{
    Interactor* interactor = static_cast<Interactor*>(userdata);
    &interactor->_mutex.lock();
    interactor->setupInteraction();
    &interactor->_mutex.unlock();
    while (interactor->continueInteractionThread())
    {
        interactor->processInteractions();
        usleep(interactor->_sleepDuration);
    }
    interactor->terminateInteraction();
    return NULL;
}
