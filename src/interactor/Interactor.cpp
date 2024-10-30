#include "interactor/Interactor.h"
#include "SpringNetwork.h"

Interactor::Interactor() : _springnetwork(nullptr), _nbpositions(0), _sleepDuration(1000), _isRunning(false)
{
	pthread_mutex_t mutextmp1 = PTHREAD_MUTEX_INITIALIZER;
	mutex=mutextmp1;
}
	

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
	rc = pthread_create(&_thread, NULL, runthread, (void *) this);

	if (rc)
	{
		exit(-1);
	}
}


void* Interactor::runthread(void* userdata)
{
    Interactor* interactor = static_cast<Interactor*>(userdata);
    pthread_mutex_lock(&interactor->mutex);
    interactor->setupInteraction();
    pthread_mutex_unlock(&interactor->mutex);
    while (interactor->continueInteractionThread())
    {
        interactor->processInteractions();
        usleep(interactor->_sleepDuration);
    }
    interactor->terminateInteraction();
    return NULL;
}