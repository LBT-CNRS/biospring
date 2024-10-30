#ifdef MDDRIVER_SUPPORT
#ifndef _CUSTOMDATA_H_
#define _CUSTOMDATA_H_


#include <unordered_map>
#include <functional>

#include "interactor/mddriver/InteractorMDDriver.h"


namespace biospring
{
namespace interactor
{

class CustomData : public InteractorMDDriver
{
    public:
        static void processCustomFloatData(InteractorMDDriver* imdl);
        static void processCustomIntData(InteractorMDDriver* imdl);
        static void initializeDataManager(InteractorMDDriver* imdl);
        static void syncParticleStateData(InteractorMDDriver* imdl, unsigned index);
    
    protected:
        template <typename T>
		static void sendCustomDataToClient(const char *sendName, int *sendSize, T *sendData);
		template <typename T>
		static int getCustomDataFromClient(char *getName, int *getSize, T **getData);


};


} // namespace interactor
} // namespace biospring

#endif
#endif