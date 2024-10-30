#ifndef __RIGIDBODY_H__
#define __RIGIDBODY_H__

#include <Vector3f.h>
#include "Quaternion.h"
#include <cmath>
#include <iostream>
#include "InsertionVector.h"
#include "matrix.h"
#include <string>
#include "random.hpp"
#include "SpringNetwork.h"

using Random = effolkronium::random_static;

#define STATE_SIZE 13
#define NUMSUBSTEPS 2

class InsertionVector;

namespace biospring
{
namespace rigidbody
{

// http://www.cs.cmu.edu/~baraff/sigcourse/notesd1.pdf

class RigidBody
{
  public:

    RigidBody()
      : pos_ini(0.0),
        cur_pos(0.0)
    {}
    
    inline biospring::spn::SpringNetwork * getSpringNetwork() const { return _spn; }

    void fetchParameters();
    RigidBody(spn::SpringNetwork * sp,
          unsigned newRbId,
          std::vector<unsigned int> particles, 
          std::vector<unsigned int> interactionParticles = std::vector<unsigned int>(0));

    void UpdateSpringsState(bool isStatic);
    unsigned rbid;

    std::vector<unsigned> getParticlesIds() const { return _particulesIds; }
    
  
    // Initialization
    void init();
    double computeTotalMass();
    Vector3f computeBarycentre() const;
    void computeAllLocalPositions();
    Matrix computeIbody(Vector3f com);
    /* ---------------------------------------------------------------------------------------------------------------*/

    // IMPALA Sampling
    double pos_ini;
    double pos_fin;
    double cur_pos;
    double _inser_angle;
    double inser_angle_ini;
    double _roll_angle;
    Vector3f iv_vec_rot;
    double min_roll_energy;

    Vector3f computeFarthestToCentroid() const;
    void initImpalaSampling();
    void updateImpalaSampling(int nbiter, double timestep);
    Vector3f getImpalaSamplingParticlePosition(spn::Particle & p, int ind);
    /* ---------------------------------------------------------------------------------------------------------------*/

    // Monte Carlo sampling
    Vector3f _montecarlo_current_position;
    Vector3f _montecarlo_current_angles;

    Vector3f _montecarlo_next_translation_vector; // Keep the current translation to apply
    Vector3f _montecarlo_next_angles;  // Track the current 2 axis (x and y) angles to apply [+/- _random_rotation_norm, +/- _random_rotation_norm]

    double _prev_impenergy;

    void initMonteCarlo(double random_translation_norm, double random_rotation_norm);
    Vector3f getRandomTranslation() { return RandomVector3f(_spn->getMonteCarloTranslationNorm()); };
    double getRandomRotation() { return Random::get({-_spn->getMonteCarloRotationNorm(), _spn->getMonteCarloRotationNorm()}); };
    Vector3f getMonteCarloParticlePosition(spn::Particle & p, int ind);
    void solveMonteCarlo();
    /* ---------------------------------------------------------------------------------------------------------------*/
    // Rigid body dynamics computation
    
    int applydepthrestraint = 0;
    double depth_restraint = 0;
    double depth_restraint_scaling = 0;
    void applyDepthRestraint();

    int applyanglerestraint = 0;
    double angle_restraint = 0;
    double angle_restraint_scaling = 0;
    void applyAngleRestraint();
    
    int withoutRotational = 0;
    Vector3f external_force;
    Vector3f external_torque;

    void applyHydrophobicRestraint();

    void solve();

    static void computeParticleForceAndTorque(spn::Particle & p);
    static void resetRigidBodiesForceAndTorque();
    static void integrateParticleVelocity(spn::Particle & p, int ind, double timestep); // Called in spn::SpringNetwork::integrateParticles

    /* ---------------------------------------------------------------------------------------------------------------*/

    Matrix createVector3Matrix();
    void updateVector3Matrix(Matrix &m, Vector3f v);

    // Get Set
    void setPos(Vector3f pos) { _pos = pos; }
    Vector3f getPos() const { return _pos; }
    Vector3f getForce() const { return _force; }
    Vector3f getTorque() const { return _torque; }
    Vector3f getOmega() const { return _omega_v; }

    double _montecarloenergy; // Energy associated to each bodie


    // Rigid Body properties
    static bool _splitchainsenabled;

    // https://www.softwaretestinghelp.com/random-number-generator-cpp/#:~:text=generate%20the%20output.-,C%2B%2B%20Random%20Number%20Between%200%20And%201,value%20either%20float%20or%20double.
    // https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
    // http://www.math.uaa.alaska.edu/~afkjm/csce211/whandouts/RandomFunctions.pdf
    Vector3f RandomVector3f(float norm)
    {
      float x = Random::get<Random::common>(-0.5, 0.5);
      float y = Random::get<Random::common>(-0.5, 0.5);
      float z = Random::get<Random::common>(-0.5, 0.5);
      Vector3f v = Vector3f(x, y, z);
      v.normalize();
      v = v * norm;
      return v;
    }


  private:
    spn::SpringNetwork * _spn;

    std::vector<unsigned> _particulesIds;
    std::vector<unsigned> _interactionParticles;
    std::vector<unsigned> _selfHydrophilicParticles;
    std::vector<unsigned> _interactionHydrophilicParticles;

    std::vector<Vector3f> _p0;
    double _dt;        
    double _dmax; // Maximum distance from barycentre to the farthest atom in the structure
    // time step size

    /* ------------------------------------------------ */

    /* Constant quantities */
    double     _mass;                       /* mass M */
    Matrix     _ibody,                      /* Ibody */
               _ibodyinv;                   /* I−1 body (inverse of Ibody) */

    /* State variables */
    Vector3f   _pos;                       /* x(t) is _com at initial state*/

    Quaternion _orientation;

    /* Derived quantities (auxiliary variables) */            
    Vector3f   _v;                         /* v(t) */
    Vector3f   _a;                         /* a(t) */
    Vector3f   _alpha;                         /* angular acceleration */
    Vector3f   _L;                         /* angular momentum */

    Vector3f   _omega_v;                    /* ω(t) */

    /* Computed quantities */
    Vector3f   _force,                     /* F(t) */
               _torque;                    /* τ(t) */
};

} // namespace rigidbody
} // namespace biospring

#endif