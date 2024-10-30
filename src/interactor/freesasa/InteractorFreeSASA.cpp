#ifdef FREESASA_SUPPORT
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "interactor/freesasa/InteractorFreeSASA.h"
#include "logging.h"
#include "Vector3f.h"
#include "SpringNetwork.h"
#include "Particle.h"


namespace biospring
{
namespace interactor
{

InteractorFreeSASA::InteractorFreeSASA() : _sasa(nullptr)
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

InteractorFreeSASA::~InteractorFreeSASA()
{
	delete[] _sasa;
}

void InteractorFreeSASA::startInteractionThread()
{
	initializeSystemState();
	Interactor::startInteractionThread();
}

void InteractorFreeSASA::initializeSystemState()
{
	Interactor::initializeSystemState(); // Set _nbpositions int defined in Interactor.

	radii_array.resize(_nbpositions);
    coords_array.resize(_nbpositions * 3);

	initializeDataManager();
}

void InteractorFreeSASA::initializeDataManager()
{
	floatManager.add("sasa", _nbpositions);
}

void InteractorFreeSASA::setupFreesasaInteractions()
{
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
	const freesasa_classifier* classifier;
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

	for (unsigned i=0; i<_nbpositions; ++i)
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

		if (radius == -1.0)
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
	_isRunning = true;
}

void InteractorFreeSASA::processFreesasaInteractions()
{
    // double coords_array[_nbpositions * 3];

	for (unsigned i = 0; i < _nbpositions; ++i)
	{
		biospring::spn::Particle & p = getSpringNetwork()->getParticle(i);

		Vector3f c = p.getPosition();
		coords_array[i * 3] = c.getX();       // Stockage de la coordonnée X
		coords_array[i * 3 + 1] = c.getY();   // Stockage de la coordonnée Y
		coords_array[i * 3 + 2] = c.getZ();   // Stockage de la coordonnée Z
	}

    freesasa_result * result = nullptr;
    freesasa_parameters * param = getParams();
    // param.alg = FREESASA_SHRAKE_RUPLEY;
    // param.shrake_rupley_n_points = 1000;
    // param.alg = FREESASA_SHRAKE_RUPLEY;
    // param.shrake_rupley_n_points = 10;

    result = freesasa_calc_coord(coords_array.data(), radii_array.data(), _nbpositions, param);

	if (result != NULL)
	{
		// Allocate memory for _sasa if not already allocated
    	if (_sasa == NULL)
		{
			_sasa = (double *)malloc(sizeof(double) * _nbpositions);
		}
	}
	for (unsigned i = 0; i < _nbpositions; ++i)
	{
		_sasa[i] = result->sasa[i];
	}
	
    // _sasa = result->sasa;
	// logging::info("New total sasa : %f", _total);
	_total = result->total;
    freesasa_result_free(result);

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
	getSpringNetwork()->getParticle(index).setSolventAccessibilitySurface(_sasa[index]);
}

} // namespace interactor
} // namespace biospring

#endif