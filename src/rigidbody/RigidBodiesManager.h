#ifndef __RIGIDBODIESMANAGER_H__
#define __RIGIDBODIESMANAGER_H__

#include "RigidBody.h"

namespace biospring
{
namespace rigidbody
{

class RigidBodiesManager
{
  public:
    static void InitRigidBodies(spn::SpringNetwork * sp, std::vector<unsigned> newParticlesIds);
    static void SolveRigidBodiesDynamic();
    static void CleanRigidBodies();

    static std::vector<RigidBody*> getCollection() { return collection; }
    
  private:
      static std::vector<RigidBody*> collection;
};

} // namespace rigidbody
} // namespace biospring

#endif