#ifdef FREESASA_SUPPORT

#ifndef _INTERACTORFREESASA_H_
#define _INTERACTORFREESASA_H_

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include "freesasa.h"
#include <string>

#include "interactor/Interactor.h"

namespace biospring
{
namespace interactor
{


class InteractorFreeSASA : public Interactor
{
    public :

        InteractorFreeSASA();
        virtual ~InteractorFreeSASA();
                
        inline void setDynamic(bool isDynamic) {_isDynamic = isDynamic;};
        inline bool isDynamic() {return _isDynamic;};

        inline void setStep(unsigned step) {_step = step;};
        inline unsigned getStep() const { return _step;};

        // FreeSASA params

        inline void setAlg(std::string alg) {_freeSASA_alg = alg;};
        inline std::string getAlg() const { return _freeSASA_alg;};

        inline void set_lr_n(unsigned lr_n) {_lr_n = lr_n;};
        inline unsigned getFreeSASA_lr_n() const { return _lr_n;};

        inline void set_sr_n(unsigned sr_n) {_sr_n = sr_n;};
        inline unsigned get_sr_n() const { return _sr_n;};

        inline void setProbeRad(float probe_radius) {_proberadius = probe_radius;};
        inline float getProbeRad() const { return _proberadius;};

        inline void setRadiiClassifier(std::string radii_classifier) {_radiiclassifier = radii_classifier;};
        inline std::string getRadiiClassifier() const { return _radiiclassifier;};

        inline void setNthreads(int nthreads) {_nthreads = nthreads;};
        inline int getNthreads() const { return _nthreads;};

        inline void setParams(freesasa_parameters *params) {_params = params;};
        inline freesasa_parameters* getParams() const { return _params;};

        inline void setClassifier(const freesasa_classifier* classifier) {_classifier = classifier;};
        inline const freesasa_classifier* getClassifier() const { return _classifier;};

        // inline void setRadii(double* radii) {_radii = radii;};
        // inline double* getRadii() const { return _radii;};

        inline float getSASA_total() { return _total;};

        virtual void startInteractionThread() override;
        virtual bool continueInteractionThread() override { return _isRunning; }
        virtual void stopInteractionThread() override { _isRunning = false; }

        virtual void syncSystemStateData() override;

        virtual void initializeDataManager() override;

        DataArrayManager<float> floatManager;

    protected : 

        double *_sasa;

        bool _isDynamic;
        unsigned _step;
        std::string _freeSASA_alg;
        unsigned _lr_n;
        unsigned _sr_n;
        double _proberadius;
        std::string _radiiclassifier;
        int _nthreads;
        freesasa_parameters *_params;
        const freesasa_classifier *_classifier;

        std::vector<double> radii_array;
        std::vector<double> coords_array;

        double _total;

        virtual void setupInteraction() override { setupFreesasaInteractions(); }
        virtual void processInteractions() override { processFreesasaInteractions(); }

        void setupFreesasaInteractions();
        void processFreesasaInteractions();

        virtual void initializeSystemState() override;

        virtual void syncParticleStateData(unsigned index) override;

        bool _isRunning;

};

} // namespace interactor
} // namespace biospring

#endif

#endif
