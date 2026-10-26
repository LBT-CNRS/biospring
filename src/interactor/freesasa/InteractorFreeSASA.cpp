#ifdef FREESASA_SUPPORT
#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "interactor/freesasa/InteractorFreeSASA.h"
#include "logging.h"
#include "Vector3f.h"
#include "SpringNetwork.h"
#include "Particle.h"


namespace biospring
{
namespace interactor
{

InteractorFreeSASA::InteractorFreeSASA() : _sasaValid(false)
{
	// https://freesasa.github.io/doxygen/group__core.html
	_freeSASA_alg = "lr";
	_proberadius = 1.4;
	_lr_n = 20;
	_sr_n = 100;
	_radiiclassifier = "default_classifier";
	_nthreads = 1;
	_isDynamic = false;
}

InteractorFreeSASA::~InteractorFreeSASA() {}

void InteractorFreeSASA::startInteractionThread()
{
	initializeSystemState();

	// THE FIRST SURFACES ARE COMPUTED HERE, on the main thread, before the run
	// starts -- not left to the worker to produce whenever it gets round to it.
	//
	// The worker used to produce them, and the first steps ran before it had.
	// syncParticleStateData rightly refuses to publish a surface it does not
	// have yet -- the surface MULTIPLIES the IMPALA energy and force, so
	// publishing zero would switch the term off in silence -- so those steps
	// ran on whatever surface the particles already carried. How MANY steps
	// that was depended on how fast the worker happened to be, which made the
	// trajectory of an IMPALA run depend on timing.
	//
	// That was not hypothetical, and it had a reproducible trigger: registering
	// the MDDriver interactor makes SpringNetwork::run() sleep 100 ms before the
	// loop, which is ample for the worker to finish. Measured on 052.Trypsin,
	// three steps at 10 fs, spring energy 1.82 kJ.mol-1 with the interactor
	// registered against 0.57 without -- the same binary on the same input, a
	// transport deciding the physics. The three examples that moved when
	// MDDriver stopped being registered by default were exactly the three with
	// impala.enable = 1, which is the only term that reads a surface.
	//
	// Published to the particles here too, not merely computed: the two
	// backends call SpringNetwork::idleRun() at opposite ends of a step (see
	// SpringNetworkOpenCL::idleRun), so a surface that exists but has not been
	// copied onto the particles would still reach the CPU's first force
	// evaluation and not the device's.
	{
		std::lock_guard<std::mutex> lock(mutex);
		setupFreesasaInteractions();
		processFreesasaInteractions();
	}
	syncSystemStateData();

	// The worker exists only to REFRESH them. Without --sasa-dynamic there is
	// nothing to refresh, so no thread is started at all: one computation, on
	// the main thread, and done.
	if (isDynamic())
		Interactor::startInteractionThread();
}

void InteractorFreeSASA::initializeSystemState()
{
	Interactor::initializeSystemState(); // Set _nbpositions int defined in Interactor.

	// Sized here, on the main thread, before Interactor::startInteractionThread()
	// spawns the worker: the buffer never has to grow while the worker writes it.
	// _sasaValid is cleared on every start, not just at construction, because
	// startInteractionThread() waits for the previous thread and may be entered
	// again with a different _nbpositions.
	radii_array.resize(_nbpositions);
    coords_array.resize(_nbpositions * 3);
	_sasa.assign(static_cast<size_t>(_nbpositions), 0.0);
	_sasaValid.store(false, std::memory_order_release);
	// Cleared for the same reason as _sasaValid: a second start may come with a
	// different _nbpositions, so the radii and the parameters have to be built
	// again. See _setupDone.
	_setupDone = false;

	initializeDataManager();
}

void InteractorFreeSASA::initializeDataManager()
{
	floatManager.add("sasa", _nbpositions);
}

void InteractorFreeSASA::setupFreesasaInteractions()
{
	// ONCE per start. startInteractionThread() now runs this on the main thread,
	// and Interactor::runthread() calls setupInteraction() again when a worker
	// is started -- running it twice would re-emit every missing-radius warning
	// and leak the freesasa_parameters the first pass allocated. The flag is
	// still set, because the worker's loop condition reads it.
	if (_setupDone)
	{
		_isRunning.store(true, std::memory_order_release);
		return;
	}

	_sleepDuration = 1500;  // For example, set sleep duration to 1500 microseconds for InteractorFreeSASA.
	
    bool use_bs_radii = false;
	// Set parameters
	freesasa_parameters *params = new freesasa_parameters;

	if (getAlg() == "sr")
	{
		params->alg = FREESASA_SHRAKE_RUPLEY;
		params->shrake_rupley_n_points = get_sr_n();
	}
	else if (getAlg() == "lr")
	{
		params->alg = FREESASA_LEE_RICHARDS;
		params->lee_richards_n_slices = getFreeSASA_lr_n();
	}
	params->probe_radius = getProbeRad();
	params->n_threads = getNthreads();
	setParams(params);

	// Get atom classifier
	const freesasa_classifier* classifier = nullptr;
	std::string chosen_classifier = getRadiiClassifier();
	if (strcmp(chosen_classifier.c_str(), "protor") == 0)
	{
		// biospring::logging::info("Protor classifier");
			classifier = &freesasa_protor_classifier;
	}
	else if (strcmp(chosen_classifier.c_str(), "naccess") == 0)
	{
		// biospring::logging::info("NACCESS classifier");
			classifier = &freesasa_naccess_classifier;
	}
	else if (strcmp(chosen_classifier.c_str(), "oons") == 0)
	{
		// biospring::logging::info("OONS classifier");
			classifier = &freesasa_oons_classifier;
	}
	else if (strcmp(chosen_classifier.c_str(), "default") == 0)
	{
		// biospring::logging::info("Default classifier");
			classifier = &freesasa_default_classifier;
	}
    else if (strcmp(chosen_classifier.c_str(), "biospring") == 0)
    {
        // biospring::logging::info("BioSpring radii values");
        use_bs_radii = true;
    }
	else
	{
			biospring::logging::error("unknown FreeSASA classifier specified");
      		return;
	}
	setClassifier(classifier);

	// Get radii
	// double radii_array[getNbPositions()];

	for (int i = 0; i < _nbpositions; ++i)
	{
		biospring::spn::Particle & p = getSpringNetwork()->getParticle(i);

        if (use_bs_radii)
        {
            radii_array[i] = p.getRadius();
            continue;
        }

		/*
		* This portion of the code is responsible for retrieving radii of particles.
		* Initially, the particle's full name (p.getName()) is used to fetch the radius 
		* from freesasa_classifier_radius. If the returned radius is -1 (indicating the radius 
		* was not found), an attempt is made using p.getName().substr(1). This is because 
		* the particle name might be the result of a coarse-grained reduction procedure based on,
		* information contained in the "--grp <file>".
		* In our identification method, the first letter typically represents the one-letter 
		* residue code, helping to identify each bead by its type and the residue it belongs to.
		* Note: This reduction procedure can also apply to every atom, in the specific case where
		* the reduction results in a single bead per atom.
		*/

		double radius = freesasa_classifier_radius(classifier, p.getResName().c_str(), p.getName().c_str());

		// substr(1) on an EMPTY name throws std::out_of_range, which left the
		// run with an uncaught exception rather than a missing radius. A
		// particle with no name is not hypothetical -- anything built without a
		// topology has one -- and the answer here is the same as for a name the
		// classifier does not know: fall through to the warning below.
		if (radius == -1.0 && !p.getName().empty())
		{
			radius = freesasa_classifier_radius(classifier, p.getResName().c_str(), p.getName().substr(1).c_str());
		}

		if (radius > 0.0)
		{
			radii_array[i] = radius;
		}
		else if (radius < 0.0 && p.getName().find("H") != std::string::npos)
		{
            // Defined by default in freesasa
			radii_array[i] = 1.1;
		}
		else
		{
			radii_array[i] = radius;
			biospring::logging::warning("freesasa_classifier_radius: radius not found fo res: %s name: %s", p.getResName().c_str(), p.getName().c_str());
		}
	}
	// setRadii(radii_array);
	// Set only HERE, at the end: every error path above leaves by `return`
	// without building the radii, and a setup that failed must be retriable
	// rather than counted as done.
	_setupDone = true;
	_isRunning.store(true, std::memory_order_release);
}

void InteractorFreeSASA::processFreesasaInteractions()
{
    // double coords_array[_nbpositions * 3];

	for (int i = 0; i < _nbpositions; ++i)
	{
		biospring::spn::Particle & p = getSpringNetwork()->getParticle(i);

		Vector3f c = p.getPosition();
		coords_array[i * 3] = c.getX();       // Store the X coordinate
		coords_array[i * 3 + 1] = c.getY();   // Store the Y coordinate
		coords_array[i * 3 + 2] = c.getZ();   // Store the Z coordinate
	}

    freesasa_result * result = nullptr;
    freesasa_parameters * param = getParams();
    // param.alg = FREESASA_SHRAKE_RUPLEY;
    // param.shrake_rupley_n_points = 1000;
    // param.alg = FREESASA_SHRAKE_RUPLEY;
    // param.shrake_rupley_n_points = 10;

    result = freesasa_calc_coord(coords_array.data(), radii_array.data(), _nbpositions, param);

	// freesasa_calc_coord returns NULL on allocation failure or invalid input.
	// Reading result->sasa in that case dereferences a null pointer, so keep
	// every use of `result` inside this branch and leave _sasa untouched: the
	// previous areas, if any, stay valid and a later pass can still succeed.
	if (result == NULL)
	{
		BIOSPRING_WARN_ONCE("FreeSASA: area computation failed, keeping the previous solvent accessible surfaces");
		return;
	}

	for (int i = 0; i < _nbpositions; ++i)
	{
		_sasa[static_cast<size_t>(i)] = result->sasa[i];
	}

	// logging::info("New total sasa : %f", _total);
	_total = result->total;
    freesasa_result_free(result);

	// Publishes _sasa to the main thread, which reads it in syncParticleStateData.
	_sasaValid.store(true, std::memory_order_release);

	if (!isDynamic())
		stopInteractionThread();
}

void InteractorFreeSASA::syncSystemStateData()
{
	Interactor::syncSystemStateData();
	if(_springnetwork!=NULL)
	{
		_springnetwork->setSASATotal(_total);
		_springnetwork->isFreeSASADynamic(_isDynamic);
	}

}

// Update the `_solventaccessibilitysurface` of the given particle
// based on FreeSASA-related computations.
void InteractorFreeSASA::syncParticleStateData(unsigned index)
{
	// The areas are computed by the worker thread, but the main thread's first
	// idleRun -> syncSystemStateData -> syncParticleStateData can run before that
	// first computation completes -- _isRunning is already true by then, and
	// nothing else gates this path. Do nothing until the areas exist rather than
	// publish a surface of 0, which would silently zero the IMPALA energy and
	// force (surface is a multiplicative factor there) for the first frames. A
	// later idleRun applies the real areas once the worker has produced them.
	if (!_sasaValid.load(std::memory_order_acquire))
	{
		BIOSPRING_WARN_ONCE("FreeSASA: solvent accessible surfaces are not computed yet, "
		                    "skipping the update for this frame");
		return;
	}

	getSpringNetwork()->getParticle(index).setSolventAccessibilitySurface(_sasa[index]);
}

} // namespace interactor
} // namespace biospring

#endif