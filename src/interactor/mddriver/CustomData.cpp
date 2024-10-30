#ifdef MDDRIVER_SUPPORT

#include <string>
#include <unordered_map>


#include "CustomData.h"
#include "SpringNetwork.h"

#include "rigidbody/RigidBodiesManager.h"

namespace biospring
{
namespace interactor
{

void CustomData::initializeDataManager(InteractorMDDriver * imdl)
{
    // Sasa sent to client
	imdl->floatManager.add("sasa", imdl->getNbPositions());

	// IMS particle forces sent to client
	imdl->floatManager.add("imsf", imdl->getNbPositions() * 3);

	// transfer energy sent to client
	imdl->floatManager.add("trimp", imdl->getNbPositions());

	// Recieved ids of new defined rigid particle (for rigidbody)
	imdl->floatManager.add("rigidparticlesids", imdl->getNbPositions());

	// RigidBody informations sent to client
	imdl->floatManager.add("rbinfo", 9);

    // IMPALA informations sent to client
    imdl->floatManager.add("impala", 3);

    // Insertion vector information sent to client
    imdl->floatManager.add("insvec", 2);
    imdl->intManager.add("insvec", 2);

    // Simulation informations sent to client
    imdl->floatManager.add("siminfo", 8);

    // CSV information (insertion depth & angle + nb iter) sent to client
    imdl->floatManager.add("csvinfo", 3);

    imdl->refFloatManager.add("getcustomfloat");
	imdl->refIntManager.add("getcustomint");
}

void CustomData::processCustomFloatData(InteractorMDDriver * imdl) {

    int nbfloat;
	char* customFloatDataName = imdl->refFloatManager.get("getcustomfloat").getDataName();
	float* customFloat = imdl->refFloatManager.get("getcustomfloat").getPointer();
    handleIMDWorkflow(imdl);
	int isNewCustomFloatData = getCustomDataFromClient(customFloatDataName, &nbfloat, &customFloat);

	if (isNewCustomFloatData)
    {
        // Check if customFloatDataName matches "dmou": Double Membrane Offset Upper
        if (customFloatDataName && std::strcmp(customFloatDataName, "dmou") == 0) {
            // Check if the SpringNetwork is enabled and if there is exactly one float value
            if (imdl->getSpringNetwork()->isIMPEnabled()) {
                if (nbfloat == 1) {
                    biospring::forcefield::ForceField* ff = imdl->getSpringNetwork()->getForceField();
                    logging::info("Set new upper membrane offset to: %f", customFloat[0]);
                    ff->setImpDoubleMembraneUpperMembOffset(customFloat[0]);
                }
            }
        }
        // Check if customFloatDataName matches "dmol": Double Membrane Offset Lower
        if (customFloatDataName && std::strcmp(customFloatDataName, "dmol") == 0) {
            // Check if the SpringNetwork is enabled and if there is exactly one float value
            if (imdl->getSpringNetwork()->isIMPEnabled()) {
                if (nbfloat == 1) {
                    biospring::forcefield::ForceField* ff = imdl->getSpringNetwork()->getForceField();
                    logging::info("Set new lower membrane offset to: %f", customFloat[0]);
                    ff->setImpDoubleMembraneLowerMembOffset(customFloat[0]);
                }
            }
        }
        // Check if customFloatDataName matches "dmtc" : Double Membrane Tube Curvature
        if (customFloatDataName && std::strcmp(customFloatDataName, "dmtc") == 0) {
            // Check if the SpringNetwork is enabled and if there is exactly one float value
            if (imdl->getSpringNetwork()->isIMPEnabled()) {
                if (nbfloat == 2) {
                    biospring::forcefield::ForceField* ff = imdl->getSpringNetwork()->getForceField();
                    logging::info("Recieved new double membrane tube curvature to: %f(upper memb) and %f(lower memb)", customFloat[0], customFloat[1]);
                    if (!isnan(customFloat[0]) && isnan(customFloat[1]))
                    {
                        logging::info("Set new double membrane tube curvature to: %f(upper memb) ", customFloat[0]);
                        logging::info("Keep double membrane tube curvature to: %f(lower memb) ", ff->getImpDoubleMembraneLowerMembTubeCurv());
                        ff->setImpDoubleMembraneUpperMembTubeCurv(customFloat[0]);
                    }
                      
                    if (!isnan(customFloat[1]) && isnan(customFloat[0]))
                    {
                        logging::info("Set new double membrane tube curvature to: %f(lower memb) ", customFloat[1]);
                        logging::info("Keep double membrane tube curvature to: %f(upper memb) ", ff->getImpDoubleMembraneUpperMembTubeCurv());
                        ff->setImpDoubleMembraneLowerMembTubeCurv(customFloat[1]);
                    }

                    if (!isnan(customFloat[1]) && !isnan(customFloat[0]))
                    {
                        logging::info("Set new double membrane tube curvature to: %f(upper memb) and %f(lower memb)", customFloat[0], customFloat[1]);
                        ff->setImpDoubleMembraneUpperMembTubeCurv(customFloat[0]);
                        ff->setImpDoubleMembraneLowerMembTubeCurv(customFloat[1]);
                    }
                        
                }
            }
        }
        else if (customFloatDataName && std::strcmp(customFloatDataName, "sasa") == 0)
        {
    #ifdef FREESASA_SUPPORT
            float* sasa = imdl->floatManager.get("sasa").getData();
            int nbPositions = imdl->getNbPositions();
            sendCustomDataToClient("sasa", &(nbPositions), sasa);
    #endif
        }
        // Check if client ask for ims particles forces
        else if (customFloatDataName && std::strcmp(customFloatDataName, "imsf") == 0)
        {
            float* imsf = imdl->floatManager.get("imsf").getData();
            int nbPositions = imdl->getNbPositions() * 3;
            sendCustomDataToClient("imsf", &(nbPositions), imsf);
        }
        // Check if client ask for impala transfer energies
        else if (customFloatDataName && std::strcmp(customFloatDataName, "trimp") == 0)
        {
            if (imdl->getSpringNetwork()->isIMPEnabled())
            {
                float* trimp = imdl->floatManager.get("trimp").getData();
                int nbPositions = imdl->getNbPositions();
                sendCustomDataToClient("trimp", &(nbPositions), trimp);
            }
            else
            {
                logging::info("Client asked for the transfer energies but the IMPALA option is not enabled.");
            }
            
        }
        else if (customFloatDataName && std::strcmp(customFloatDataName, "impala") == 0)
        {
            if (imdl->getSpringNetwork()->isIMPEnabled())
            {
                float* impala = imdl->floatManager.get("impala").getData();
                int sizeInfo = imdl->floatManager.get("impala").getSize();

                // Check if user send new data to update
                // Special condition here because the size is 1 so check if not
                // just ping from client.
                if (nbfloat == 1 && customFloat[0] >= 0.0)
                {
                    logging::info("Set new IMPALA scale to: %f", customFloat[0]);
                    imdl->getSpringNetwork()->getForceField()->setIMPScale(customFloat[0]);
                }

                // impala[0] <-- !!! IMPALA energy is updated in syncParticleStateData (bellow)
                // Before this energy set to 0 in the main SpringNetwork running loop
                // (avoid sending 0 value sometimes because here we're not in the main thread)
                impala[1] = imdl->getSpringNetwork()->getCentroid()[2]; // get z of centroid for insertionDepth
                impala[2] = imdl->getSpringNetwork()->getForceField()->getIMPScale();

                sendCustomDataToClient("impala", &(sizeInfo), impala);
            }
            else
            {
                logging::info("Client asked for IMPALA informations but the option is not enabled.");
            }
        }
        // Check if client ask for rigidbody informations
        else if (customFloatDataName && std::strcmp(customFloatDataName, "rbinfo") == 0)
        {
            if (imdl->getSpringNetwork()->isRigidBodyEnabled())
            {
                float* rbinfo = imdl->floatManager.get("rbinfo").getData();
                int sizeInfo = imdl->floatManager.get("rbinfo").getSize();


                rigidbody::RigidBody *rb = rigidbody::RigidBodiesManager::getCollection()[0]; // Consider first rb as the main rb every time
                Vector3f pos = rb->getPos();
                Vector3f torque = rb->getTorque();
                Vector3f omega = rb->getOmega();

                auto copyToFloatArray = [](const Vector3f& vec, float* dest) {
                    dest[0] = static_cast<float>(vec.getX());
                    dest[1] = static_cast<float>(vec.getY());
                    dest[2] = static_cast<float>(vec.getZ());
                };

                copyToFloatArray(pos, rbinfo);
                copyToFloatArray(torque, rbinfo + 3);
                copyToFloatArray(omega, rbinfo + 6);

                sendCustomDataToClient("rbinfo", &(sizeInfo), rbinfo);
            }
            else
            {
                logging::info("Client asked for RigidBody informations but the option is not enabled.");
            }
        }
        else if (customFloatDataName && std::strcmp(customFloatDataName, "insvec") == 0)
        {
            if (imdl->getSpringNetwork()->isInsertionVectorEnabled())
            {
                float* insvec = imdl->floatManager.get("insvec").getData();
                int sizeInfo = imdl->floatManager.get("insvec").getSize();

                const InsertionVector iv = imdl->getSpringNetwork()->getInsertionVector();
                insvec[0] = iv.getAngle();
                insvec[1] = iv.getRollAngle();

                sendCustomDataToClient("insvec", &(sizeInfo), insvec);
            }
            else
            {
                logging::info("Client asked for Insertion Vector informations but the option is not enabled.");
            }
        }
        // This part regroups simulation informations to send to the client:
        // - viscosity
        // - MonteCarlo values
        // - Tube curvature
        else if (customFloatDataName && std::strcmp(customFloatDataName, "siminfo") == 0)
        {
            
            float* siminfo = imdl->floatManager.get("siminfo").getData();
            int sizeInfo = imdl->floatManager.get("siminfo").getSize();

            // Send the current viscosity if the option is enabled
            if (imdl->getSpringNetwork()->isViscosityEnabled())
                siminfo[0] = imdl->getSpringNetwork()->getViscosity();
            else
                siminfo[0] = 0.0;
            
            // Senf the current MonteCarlo values if rigidbody and MonteCarlo
            // enabled
            if (imdl->getSpringNetwork()->isRigidBodyEnabled() &&
                imdl->getSpringNetwork()->isMonteCarloEnabled())
            {
                siminfo[1] = imdl->getSpringNetwork()->getMonteCarloTemperature();
                siminfo[2] = imdl->getSpringNetwork()->getMonteCarloTranslationNorm();
                siminfo[3] = imdl->getSpringNetwork()->getMonteCarloRotationNorm();
            }
            else
            {
                siminfo[1] = 0.0;
                siminfo[2] = 0.0;
                siminfo[3] = 0.0;
            }

            if (imdl->getSpringNetwork()->isIMPEnabled())
            {
                siminfo[4] = imdl->getSpringNetwork()->getForceField()->getImpDoubleMembraneUpperMembOffset();
                siminfo[5] = imdl->getSpringNetwork()->getForceField()->getImpDoubleMembraneLowerMembOffset();
                siminfo[6] = imdl->getSpringNetwork()->getForceField()->getImpDoubleMembraneUpperMembTubeCurv();
                siminfo[7] = imdl->getSpringNetwork()->getForceField()->getImpDoubleMembraneLowerMembTubeCurv();
            }
            
            sendCustomDataToClient("siminfo", &(sizeInfo), siminfo);
            
        }
        // Update viscosity
        else if (customFloatDataName && std::strcmp(customFloatDataName, "visc") == 0)
        {
            if (imdl->getSpringNetwork()->isViscosityEnabled())
            {
                // Check if user send new data to update
                // If there is only one value to update you must add a special
                // condition here because the size is 1 so check if not just 
                // ping from client to send the value. The condition to use:
                // if (nbfloat == 1 && customFloat[0] >= 0.0)
                if (nbfloat == 1 && customFloat[0] >= 0.0)
                {
                    logging::info("Set new viscosity to: %f", customFloat[0]);
                    imdl->getSpringNetwork()->setViscosity(customFloat[0]);
                }
            }
        }
        else if (customFloatDataName && std::strcmp(customFloatDataName, "monteca") == 0)
        {
            if (imdl->getSpringNetwork()->isRigidBodyEnabled() &&
                imdl->getSpringNetwork()->isMonteCarloEnabled())
            {
                if (nbfloat == 3)
                {
                    logging::info("Set new Monte Carlo values: %f°K, translation norm %f, rotation norm %f",
                        customFloat[0], customFloat[1], customFloat[2]);
                    imdl->getSpringNetwork()->setMonteCarloTemperature(customFloat[0]);
                    imdl->getSpringNetwork()->setMonteCarloTranslationNorm(customFloat[1]);
                    imdl->getSpringNetwork()->setMonteCarloRotationNorm(customFloat[2]);
                }
            }
        }


        isNewCustomFloatData = 0;
    }

    // Send data directly without waiting to recieve the ping from client
    // Send angle and insertion depth + nb inter when writing frame to csv
    // to synchronize information sent to client
    if (imdl->getSpringNetwork()->isCSVTrajectoryWriterEnabled() &&
        imdl->getSpringNetwork()->isInsertionVectorEnabled() &&
        !(imdl->getSpringNetwork()->isImpalaSamplingEnabled()) &&
        imdl->getSpringNetwork()->getNbIterations() % imdl->getSpringNetwork()->getCSVTrajectoryWriterFreq() == 0)
    {
        float* csvinfo = imdl->floatManager.get("csvinfo").getData();
        int sizeInfo = imdl->floatManager.get("csvinfo").getSize();

        csvinfo[0] = imdl->getSpringNetwork()->getNbIterations();
        const InsertionVector iv = imdl->getSpringNetwork()->getInsertionVector();
        csvinfo[1] =  iv.getAngle();
        csvinfo[2] =  iv.getInsertionDepth();

        sendCustomDataToClient("csvinfo", &(sizeInfo), csvinfo);
    }
}

void CustomData::processCustomIntData(InteractorMDDriver* imdl) {

    int nbint;
	char* customIntDataName = imdl->refIntManager.get("getcustomint").getDataName();
	int* customInt = imdl->refIntManager.get("getcustomint").getPointer();
    handleIMDWorkflow(imdl);
	int isNewCustomIntData = getCustomDataFromClient(customIntDataName, &nbint, &customInt);

	if (!isNewCustomIntData) return;

    // Rigid particles ids
    if (customIntDataName && strcmp(customIntDataName, "rpids") == 0) {
        if (imdl->getSpringNetwork()->isRigidBodyEnabled()) {
            if (nbint==1 && customInt[0]==-1)
                rigidbody::RigidBodiesManager::CleanRigidBodies();
            else {
                std::vector<unsigned int> rigidparticlesids;
                for (int i = 0; i < nbint; i++) {
                    rigidparticlesids.push_back(static_cast<unsigned int>(customInt[i]));
                }
                rigidbody::RigidBodiesManager::InitRigidBodies(imdl->getSpringNetwork(), rigidparticlesids);
            }

        } else {
            logging::warning("Recieved rigid particles ids but RigidBody dynamic is not enabled!");
        }
    }
    else if (customIntDataName && std::strcmp(customIntDataName, "paused") == 0)
    {
        int isPaused = imdl->getSpringNetwork()->getPause();
        int isPausedArray[1] = {isPaused}; 
        int nbData = 1;
        // logging::info("Send paused : %d", isPaused);
        sendCustomDataToClient("paused", &nbData, isPausedArray);
    }
    else if (customIntDataName && std::strcmp(customIntDataName, "insvec") == 0)
    {
        if (imdl->getSpringNetwork()->isInsertionVectorEnabled())
        {
            int* insvec = imdl->intManager.get("insvec").getData();
            int sizeInfo = imdl->intManager.get("insvec").getSize();

            InsertionVector iv = imdl->getSpringNetwork()->getInsertionVector();
            insvec[0] = iv.getParticle(0).getId();
            insvec[1] = iv.getParticle(1).getId();

            sendCustomDataToClient("insvec", &(sizeInfo), insvec);
        }
        else
        {
            logging::info("Client asked for Insertion Vector informations but the option is not enabled.");
        }
    }
    // If client ask for Potential grid
    else if (customIntDataName && std::strcmp(customIntDataName, "pgrid") == 0)
    {
        int isGridArray[1] = {0}; 
        int nbData = 1;
        if (imdl->getSpringNetwork()->isElectrostaticFieldEnabled())
        {
            isGridArray[0] = 1;
            sendPotentialGrid(imdl);
            sendCustomDataToClient("pgrid", &(nbData), isGridArray);

        } else {
            sendCustomDataToClient("pgrid", &(nbData), isGridArray);
        }
    }
    // If client ask for Density grid
    else if (customIntDataName && std::strcmp(customIntDataName, "dgrid") == 0)
    {
        int isGridArray[1] = {0}; 
        int nbData = 1;
        if (imdl->getSpringNetwork()->isDensityGridEnabled())
        {
            isGridArray[0] = 1;
            sendDensityGrid(imdl);
            // Send int 1 to tell the cilent that the grid has been sent
            sendCustomDataToClient("dgrid", &(nbData), isGridArray);

        } else {
            // Send int 0 otherwise
            sendCustomDataToClient("dgrid", &(nbData), isGridArray);
        }
    }
}

template <typename T>
void CustomData::sendCustomDataToClient(const char *sendName, int *sendSize, T *sendData) {
    logging::info("Send %s to client.", sendName);

    if constexpr (std::is_same_v<T, float>) {
        IIMD_send_custom_float(sendName, sendSize, sendData);
    } else if constexpr (std::is_same_v<T, int>) {
        IIMD_send_custom_int(sendName, sendSize, sendData);
    }
}

template <typename T>
int CustomData::getCustomDataFromClient(char *getName, int *getSize, T **getData) {
	
    if constexpr (std::is_same_v<T, float>) {
        return IIMD_get_custom_float(&getName, getSize, getData);
    } else if constexpr (std::is_same_v<T, int>) {
        return IIMD_get_custom_int(&getName, getSize, getData);
    }
}


void CustomData::syncParticleStateData(InteractorMDDriver* imdl, unsigned index)
{
    biospring::spn::Particle particle = imdl->getSpringNetwork()->getParticle(index);
    // Update sasa array for the given particle
	float sasa = particle.getSolventAccessibilitySurface();
	float* sasaArray = imdl->floatManager.get("sasa").getData();
	memcpy(&(sasaArray[index]), &sasa, sizeof(float));

	// Update particle force array for the given particle
	float particleForce[3];
    particle.getPreviousForce().to_array(particleForce);
	float* imsf = imdl->floatManager.get("imsf").getData();
	memcpy(&(imsf[index * 3]), particleForce, sizeof(float) * 3);

	// Update transfert energies array for the given particle
	float tre = particle.getTransferEnergyByAccessibleSurface();
	float* trimpArray = imdl->floatManager.get("trimp").getData();
	memcpy(&(trimpArray[index]), &tre, sizeof(float));

    // !!! Check processCustomFloatData to see other impala data being updated
    float* impala = imdl->floatManager.get("impala").getData();
    impala[0] = imdl->getSpringNetwork()->getIMPEnergy(); //!< IMPALA energy, KJoule/mol
}



} // namespace interactor
} // namespace biospring

#endif
