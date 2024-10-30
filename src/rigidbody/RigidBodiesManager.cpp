#include "RigidBodiesManager.h"
#include "SpringNetwork.h"
#include <vector>

namespace biospring
{
namespace rigidbody
{
    std::vector<RigidBody*> RigidBodiesManager::collection;

    void RigidBodiesManager::InitRigidBodies(spn::SpringNetwork * sp, vector<unsigned> newParticlesIds)
    {
        if (newParticlesIds.empty()) {
            // If rigidparticlesids is empty, clear the collection and return.
            CleanRigidBodies();
            return;
        }

        bool isAllNewParticle = true;
        for (auto const& rb : collection)
        {
            // Get existing ids of the rigid body found
            const auto e = rb->getParticlesIds();
            // Check if some ids in rigidparticlesids are found in existing collection of rigid bodies
            for (auto const n : newParticlesIds) {
                if (std::find(e.begin(), e.end(), n) == e.end()) {
                    // Particle found in existing rigid body ids
                    isAllNewParticle = false;
                }
            }
        }


        // If the IDs are not found in the collection, create and add a new RigidBody.
        if (isAllNewParticle) {
            unsigned newRbId = collection.size();
            RigidBody *newRb = new RigidBody(sp, newRbId, newParticlesIds);
            collection.push_back(newRb);
        }
    }

    void RigidBodiesManager::SolveRigidBodiesDynamic()
    {
#ifdef OPENMP_SUPPORT
#pragma omp parallel default(shared)
#endif
        {
#ifdef OPENMP_SUPPORT
#pragma omp for
#endif 
            for (unsigned i = 0; i < collection.size() ; i++)
            {
                //logging::info("solve rbid:%d", collection[i]->rbid);
                collection[i]->solve();
            }
        }
    }

    void RigidBodiesManager::CleanRigidBodies()
    {
        for (auto rb : collection)
        {
            rb->UpdateSpringsState(false);
            for (unsigned i = 0; i < rb->getParticlesIds().size(); i++)
            {
                spn::Particle & p = rb->getSpringNetwork()->getParticle( rb->getParticlesIds()[i]);
                p.setRigidBodyId(0);
                p.setRigid(false);
            }
            delete rb;
        }
        collection.clear();
    }

} // namespace rigidbody
} // namespace biospring
