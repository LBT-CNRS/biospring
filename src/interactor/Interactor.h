#ifndef _INTERACTOR_H_
#define _INTERACTOR_H_

#include <cstring>
#include <thread>
#include <mutex>
#include <vector>
#include <unordered_map>
#include <string>
#include "logging.h"
#include <memory>


namespace biospring
{
namespace spn
{
class SpringNetwork;
} // namespace spn
} // namespace biospring

class Interactor
{
  public:
    Interactor();

    virtual ~Interactor() {};

    inline void setSpringNetwork(biospring::spn::SpringNetwork * springnetwork) { _springnetwork = springnetwork; }
    inline biospring::spn::SpringNetwork * getSpringNetwork() const { return _springnetwork; }

    inline int getNbPositions() const { return _nbpositions; }

    virtual void startInteractionThread();
	virtual bool continueInteractionThread() = 0;
    virtual void stopInteractionThread() = 0;

    virtual void syncSystemStateData();

    template <typename T>
    static T* getInteractorInstance(const std::vector<Interactor*>& interactors)
    {
        for (Interactor* interactor : interactors) {
            T* specificInstance = dynamic_cast<T*>(interactor);
            if (specificInstance) {
                return specificInstance;
            }
        }
        return nullptr; // Return null if no instance is found.
    }


  protected:
    biospring::spn::SpringNetwork * _springnetwork;
    int _nbpositions;
    void _setNbPositions(int nbpositions) { _nbpositions = nbpositions; }

    unsigned int _sleepDuration;

    template <typename T>
    struct DataArray
    {
    private:
        std::unique_ptr<T[]> data;
        size_t size;

    public:
        DataArray(size_t initSize) : size(initSize) { reset(); }

        T* getData() const { return data.get(); }

        size_t getSize() const { return size; }

        void reset() {
            erase();
            reserveData(size);
        }

        void erase() { data.reset(); }

        void reserveData(size_t newSize) {
            data = std::make_unique<T[]>(newSize);
            memset(data.get(), 0, sizeof(T) * newSize);
        }

        void clearData() {
            if (data && size > 0)
                memset(data.get(), 0, sizeof(T) * size);
        }
    };

    template <typename T>
    class DataArrayManager {
    private:
        std::unordered_map<std::string, DataArray<T>> dataArrays;

    public:
        void add(const std::string& name, size_t size) {
            dataArrays.emplace(name, DataArray<T>(size));
        }

        DataArray<T>& get(const std::string& name) {
            auto it = dataArrays.find(name);
            if (it == dataArrays.end()) {
                biospring::logging::die("Name %s not found in dataArrays!", name.c_str());
            }
            return it->second;
        }
    };

    template <typename T>
    struct ExternalDataRef {
        private:
            T* pointer;
            char* dataName;

        public:
            ExternalDataRef() {
                dataName = (char *) malloc(sizeof(char) * 8);
            }
            ExternalDataRef(T* externalData) : pointer(externalData) {}
            ~ExternalDataRef() { if (dataName) { free(dataName); } }

            T* getPointer() const { return pointer; }
            char* getDataName() const { return dataName; }

        // Warning! This structure should not attempt to release "pointing"
        // because it doesn't own the data.
    };

    template <typename T>
    struct ExternalDataRefManager {
        private:
            std::unordered_map<std::string, ExternalDataRef<T>> references;

        public:
            // Ajoute une nouvelle référence externe (sans la créer).
            void add(const std::string& name) {
                references[name] = ExternalDataRef<T>();
            }

            // Get an external reference by name.
            ExternalDataRef<T>& get(const std::string& name) {
                return references[name];
            }

            // Deletes an external reference by name.
            void remove(const std::string& name) {
                references.erase(name);
            }
    };


    virtual void initializeSystemState();

    virtual void initializeDataManager();

    virtual void syncParticleStateData(unsigned index) = 0;

    std::thread _thread;
    std::mutex _mutex;
    bool _isRunning;

    static void* runthread(void* userdata);
    virtual void setupInteraction() = 0;        // Pure virtual method
    virtual void processInteractions() = 0;     // Pure virtual method
    virtual void terminateInteraction() {}      // Virtual method with empty default implementation

};

#endif
