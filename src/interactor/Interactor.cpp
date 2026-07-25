#include "interactor/Interactor.h"
#include "SpringNetwork.h"

Interactor::Interactor() : _springnetwork(nullptr), _nbpositions(0), _sleepDuration(1000), _isRunning(false) {}

Interactor::~Interactor()
{
    _isRunning.store(false, std::memory_order_release);
    if (_thread.joinable())
        _thread.join();
}

void Interactor::initializeSystemState()
{
    // _nbpositions stays a signed int: its address is passed directly to
    // MDDriver's IIMD_send_coords(const int *N, ...), a fixed external C ABI.
    if (_springnetwork != nullptr)
        _setNbPositions(static_cast<int>(_springnetwork->getNumberOfParticles()));
}

void Interactor::initializeDataManager() {}

void Interactor::syncSystemStateData()
{
    if (_springnetwork != nullptr)
    {
        const unsigned nbparticles = _springnetwork->getNumberOfParticles();
        for (unsigned j = 0; j < nbparticles; ++j)
            syncParticleStateData(j);
    }
}

void Interactor::waitForInteractionThread()
{
    if (_thread.joinable())
        _thread.join();
}

void Interactor::startInteractionThread()
{
    waitForInteractionThread();

    _isRunning.store(true, std::memory_order_release);
    _thread = std::thread(&Interactor::runthread, this);
}

void Interactor::runthread(Interactor* interactor)
{
    {
        std::lock_guard<std::mutex> lock(interactor->mutex);
        interactor->setupInteraction();
    }

    while (interactor->continueInteractionThread())
    {
        interactor->processInteractions();
        std::this_thread::sleep_for(std::chrono::microseconds(interactor->_sleepDuration));
    }
    interactor->terminateInteraction();
}
