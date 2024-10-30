#include "RigidBody.h"
#include "logging.h"
#include "IO/CSVSampleWriter.h"
#include "SpringNetwork.h"
#include <Vector3f.h>
#include "RigidBodiesManager.h"
#include "spn/Spring.h"

namespace logging = biospring::logging;

namespace biospring
{
namespace rigidbody
{


bool RigidBody::_splitchainsenabled = false;

RigidBody::RigidBody(spn::SpringNetwork * sp, unsigned newRbId, std::vector<unsigned int> particlesids, std::vector<unsigned int> interactionParticles)
{
    this->_spn = sp;
    this->rbid = newRbId;
    this->_particulesIds = particlesids;
    this->_interactionParticles = interactionParticles;
    this->_montecarloenergy = -1.;
    this->init();

    // Add rigid body id ref to each particle
    for (unsigned i = 0; i < _particulesIds.size(); i++)
    {
        spn::Particle & p = _spn->getParticle( this->_particulesIds[i]);
        p.setRigidBodyId(this->rbid);
        p.setRigid(true);
    }

    UpdateSpringsState(true);

}

void RigidBody::UpdateSpringsState(bool isStatic)
{
    // Loop on all springs and update their status to "static" if the
    // particle pair is rigid.
    for (unsigned i = 0; i < this->_spn->getNumberOfSprings(); i++)
    {
        biospring::spn::Spring & s = _spn->getSpring(i);
        spn::Particle & p1 = s.getParticle1();
        spn::Particle & p2 = s.getParticle2();
        if (p1.isRigid() && p2.isRigid())
            _spn->updateSpringState(s.getId(), isStatic);
    }
}

// Initialization ------------------------------------------------------------------------------------------------------

void RigidBody::init()
{
    if (_spn != NULL)
    {
        for (unsigned i = 0; i < _particulesIds.size() ; i++)
        {
            spn::Particle p;
            try {
                p = _spn->getParticle(_particulesIds[i]);
            } catch (const std::out_of_range& e) {
                logging::die("Error: %s", e.what());
            }
            
            if (p.getTransferEnergyByAccessibleSurface() > 0.0)
                _selfHydrophilicParticles.push_back(p.getId());

            Vector3f intialPosition = p.getPosition();
            p.setPreviousPosition(intialPosition);
        }

        // EXPERIMENTAL
        for (unsigned i = 0; i < _interactionParticles.size() ; i++)
        {
            spn::Particle & p = _spn->getParticle(_interactionParticles[i]);
            if (p.getTransferEnergyByAccessibleSurface() > 0.0)
                _interactionHydrophilicParticles.push_back(p.getId());
        }
        // logging::info("#_selfHydrophilicParticles: %d", _selfHydrophilicParticles.size());
        // logging::info("#_interactionHydrophilicParticles: %d", _interactionHydrophilicParticles.size());

        external_force = Vector3f();
        external_torque = Vector3f();

        // Private variables and state variables
        _dt = _spn->getTimeStep();
        _mass = computeTotalMass();
        _pos = computeBarycentre();
        _orientation= Quaternion(0, 0, 0, 1);
        _force = Vector3f();
        _torque = Vector3f();
        _omega_v = Vector3f();


        // Compute all local positions
        computeAllLocalPositions();

        // For the moment I set _dmax to 1
        // _dmax = Vector3f::distance(_spn->computeBarycentre(), _spn->computeFarthestToCentroid());
        _dmax = 1;
        // Inertial matrix of a spheric object (ball) in body frame
        _ibody = (2.0/5.0) * _mass * _dmax * _dmax * Matrix::createIdentity(3);
        _ibodyinv = _ibody.inverse();

        // IMPALA automatic sampling initialization
        if (_spn->isImpalaSamplingEnabled() && !_spn->isMonteCarloEnabled())
        {
            if (RigidBodiesManager::getCollection().size() > 0)
                logging::die("Cannot initialize multiple rigid bodies for IMPALA automatic sampling.");
            initImpalaSampling();
        }
        else if (_spn->isMonteCarloEnabled())
        {
            initMonteCarlo(_spn->getMonteCarloTranslationNorm(), _spn->getMonteCarloRotationNorm());
        }

        logging::info("Rigid body index: %d | Initial bary: %s | Mass: %f", rbid, _pos.to_string().c_str(), _mass); 
    }
}


double RigidBody::computeTotalMass()
{
    double mass = 0.;
    for (unsigned i = 0; i < _particulesIds.size(); i++)
    {
        spn::Particle & p = _spn->getParticle(_particulesIds[i]);
        mass += p.getMass();
    }
    return mass;
}

Vector3f RigidBody::computeBarycentre() const
{
    Vector3f bary = Vector3f();
    double total_mass = 0.;
    for (unsigned i = 0; i < _particulesIds.size(); i++)
    {
        spn::Particle & p = _spn->getParticle(_particulesIds[i]);
        bary = bary + p.getPosition() * p.getMass();
        total_mass = total_mass + p.getMass();
    }
    if (_particulesIds.size() != 0)
        bary = bary * (1.0 / total_mass);
    return bary;
}

void RigidBody::computeAllLocalPositions()
{
    Vector3f checkCom = Vector3f();
    for (unsigned i = 0; i < _particulesIds.size() ; i++)
    {
        spn::Particle & p = _spn->getParticle(_particulesIds[i]);
        Vector3f p0 = p.getPosition() - _pos;
        _p0.push_back(p0); // = "radius"
        checkCom = checkCom + (p0) * p.getMass();
    }
    //logging::info("Check center of mass, need to be equal to (0,0,0) and has value of (%f, %f, %f)", checkCom.getX(), checkCom.getY(), checkCom.getZ());
}

// Compute inertia tensor
Matrix RigidBody::computeIbody(Vector3f com)
{
    Matrix ibody = Matrix(3,3);

    for (unsigned i = 0; i < _particulesIds.size() ; i++)
    {
        spn::Particle & p = _spn->getParticle(_particulesIds[i]);
        double x = _p0[i].getX(); double y = _p0[i].getY(); double z = _p0[i].getZ(); 
        double m = p.getMass();
        ibody(0,0) =+ m *(y*y + z*z); ibody(0,1) =+ -m * x*y;         ibody(0,2) =+ -m * x*z;
        ibody(1,0) =+ -m * x*y;       ibody(1,1) =+ -m * (x*x + z*z); ibody(1,2) =+ -m * y*z;
        ibody(2,0) =+ m * x*z;        ibody(2,1) =+ -m * y*z;         ibody(2,2) =+ -m * (x*x + y*y);

    }
    return ibody;
}

// ---------------------------------------------------------------------------------------------------------------------
// IMPALA exhaustive sampling ------------------------------------------------------------------------------------------

Vector3f RigidBody::computeFarthestToCentroid() const
{
    Vector3f farthest = _spn->getCentroid();
    Vector3f bary = _spn->getCentroid();
    double distance = 0.0;
    double distancemax = 0.0;
    for (unsigned i = 0; i < _particulesIds.size(); i++)
    {
        spn::Particle & p = _spn->getParticle(_particulesIds[i]);
        distance = Vector3f::distance(p.getPosition(), bary);
        if (distance > distancemax)
        {
            distancemax = distance;
            farthest = p.getPosition();
        }
    }
    return farthest;
}

void RigidBody::initImpalaSampling()
{
    logging::info("RigidBody: IMPALA sampling mode.");
    if (!_spn->isInsertionVectorEnabled())
        logging::die("Insertion vector must be set to compute IMPALA insertion sampling");
    

    // Set z axis insertion range.
    double dmax = Vector3f::distance(_spn->getCentroid(), computeFarthestToCentroid());
    pos_ini = 18 + dmax;
    pos_fin = -pos_ini;
    cur_pos = pos_ini;

    // Initialize insertion and roll angles
    _inser_angle = 0.;
    _roll_angle = 0.;
    min_roll_energy = 0.;

    // Initialize local insertion vector
    InsertionVector & iv = _spn->getInsertionVector();
    iv.computeVector();
    iv.computeAngle();
    iv.computeRollAngle();
    iv_vec_rot = _p0[iv.getParticle(0).getExtid()] - _p0[iv.getParticle(1).getExtid()];
    iv_vec_rot.normalize();
    inser_angle_ini = iv.getAngle();
}

void RigidBody::updateImpalaSampling(int nbiter, double timestep)
{
    if (!_spn->isInsertionVectorEnabled())
        logging::die("Insertion vector must be set to compute IMPALA insertion sampling");

    //double insertion_timestep = 1; // in milisec
    Vector3f newPos;
    if (cur_pos > pos_fin)
    {
            if (_inser_angle < 180)
            {
                //if (nbiter%(int)(insertion_timestep/timestep)==0) // Rotate every insertion_timestep
                if (true)
                {
                    if (_roll_angle == 0)
                        min_roll_energy = _spn->getIMPEnergy();


                    if (_roll_angle < 360)
                    {
                        // Find minimum roll energy
                        double new_roll_energy = _spn->getIMPEnergy();
                        if (new_roll_energy < min_roll_energy)
                            min_roll_energy = new_roll_energy;
                        _roll_angle++;
                    }
                }
                if (_roll_angle >= 360)
                {
                    // Set the minimum IMPALA energy of every roll angles at
                    // this insertion angle
                    _spn->setIMPEnergy(min_roll_energy);
                    // Write step in output csv when roll rotation is completed.
                    // IMPALA energy is then the minimum roll set just before.
                    _spn->writeNextStepNow();

                    _roll_angle = 0;
                    _inser_angle++;
                }
            }
            
            if (_inser_angle >= 180)
            {
                _inser_angle = 0.0;
                cur_pos--;
                logging::info("cur_pos: %f ; pos_fin: %f ; %f%%", cur_pos, pos_fin, (pos_ini-cur_pos)/(pos_ini-pos_fin)*100);
            }
    }
    else
        _spn->setEnd(true);
}

Vector3f RigidBody::getImpalaSamplingParticlePosition(spn::Particle & p, int ind)
{
    _pos = Vector3f(0., 0., cur_pos);

    // Rotation vector to vary insertion angle
    Vector3f iv_vec_rot_normale = iv_vec_rot^Vector3f(0,0,1);
    iv_vec_rot_normale.normalize();

    // Apply First rotation:
    // Roll rotation -> insertion vector is the rotation vector
    double roll_angle_rad = _roll_angle * (M_PI / 180);
    Vector3f newPos = Quaternion(_p0[ind], 0.).rotateVectorAboutAxisAndAngle(iv_vec_rot, roll_angle_rad).getV();
    
    // Apply second rotation:
    // Insertion angle rotation
    double inser_angle_rad = (_inser_angle + inser_angle_ini + 90) * (M_PI / 180);
    newPos = Quaternion(newPos, 0.).rotateVectorAboutAxisAndAngle(iv_vec_rot_normale, inser_angle_rad).getV();
    
    // Rescale new position and update particle position 
    // Maybe not useful now but I don't know why I needed rescale before ._.'
    // The line below seems to work good.
    // Vector3f pi = _pos + newPos * _p0[ind].norm(); /
    Vector3f pi = _pos + newPos;
    return pi;
}

// ---------------------------------------------------------------------------------------------------------------------
// Monte Carlo simulation ----------------------------------------------------------------------------------------------

void RigidBody::initMonteCarlo(double translation_norm, double rotation_norm)
{
    // Some usefull sources to implement the methods:
    // https://youtu.be/xVvUFB5Hk-g
    // https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
    // http://www.math.uaa.alaska.edu/~afkjm/csce211/handouts/RandomFunctions.pdf
    // https://github.com/effolkronium/random
    //
    logging::info("RigidBody: Monte Carlo mode.");
    if (translation_norm < 0 || rotation_norm < 0)
        logging::die("MonteCarlo sampling: Need positive or null random translation and rotation.");

    _montecarlo_current_position = _pos; // initial pos
    _montecarlo_current_angles = Vector3f(0., 0., 0.);

    _montecarlo_next_translation_vector = getRandomTranslation();
    _montecarlo_next_angles.setX(getRandomRotation());
    _montecarlo_next_angles.setY(getRandomRotation());
}

Vector3f RigidBody::getMonteCarloParticlePosition(spn::Particle & p, int ind)
{
        // Apply First rotation:
        double xangle_rad = _montecarlo_current_angles.getX() * (M_PI / 180);
        Vector3f newPos = Quaternion(_p0[ind], 0.).rotateVectorAboutAxisAndAngle(Vector3f(0,1,0), xangle_rad).getV();
        
        // Apply second rotation:
        double yangle_rad = _montecarlo_current_angles.getY() * (M_PI / 180);
        newPos = Quaternion(newPos, 0.).rotateVectorAboutAxisAndAngle(Vector3f(1,0,0), yangle_rad).getV();
        
        // Rescale new position and update particle position 
        // Maybe not useful now but I don't know why I needed rescale before ._.'
        // The line below seems to work good.
        // Vector3f pi = _montecarlo_current_position + newPos * _p0[ind].norm();
        Vector3f pi = _montecarlo_current_position + newPos;
        return pi;
}

// Monte Carlo Evaluation
void RigidBody::solveMonteCarlo()
{
    if (_prev_impenergy == -1.) // First step
        _prev_impenergy = _spn->getIMPEnergy();
    
    // Evaluation
    // Metropolis Monte Carlo ----------------------------------------------------------------------------------
    // logging::info("new ener: %f", _impenergy);
    // logging::info("pre ener: %f", rb->_prev_impenergy);
    double delta_e = _spn->getIMPEnergy() - _prev_impenergy; // < 0 if new position is more favorable
    //logging::info("delta_e: %f", delta_e);
    // p > 1 if delta_e < 0 -> accepted ; p decreases progressively when delta_e > 0 increases.
    double T = _spn->getMonteCarloTemperature();
    double p = exp(min(1., -(delta_e * 1000) / (1.380650 * T * 6.0221407))); // 
    //logging::info("delta_e: %f, p: %f", p);
    //logging::info((to_string(bf).c_str()));
    bool accepted = p >= Random::get<Random::common>(0., 1.);

    // We accept the position and we try new random movement
    if (accepted) 
    {
        _prev_impenergy = _spn->getIMPEnergy();
    }
    // We refuse the position so we go back with the previous position
    // and we try with new random movement
    else 
    {
        // logging::info("Monte Carlo: refused");
        // We go back so we reverse the next moves
        _montecarlo_current_position -= Vector3f(0., 0., -_montecarlo_next_translation_vector.getZ());
        _montecarlo_current_angles.setX(-_montecarlo_next_angles.getX());
        _montecarlo_current_angles.setY(-_montecarlo_next_angles.getY());
    }

    // and we add new random moves and positions will be  updated in next 
    // intergration
    _montecarlo_next_translation_vector = getRandomTranslation();
    _montecarlo_next_angles.setX(getRandomRotation());
    _montecarlo_next_angles.setY(getRandomRotation());

    // Update current position and angles with next random values
    // For the current position we keep the protein in the center and we only
    // modition position on z (no need to translate for a plane membrane)
    _montecarlo_current_position += Vector3f(0., 0., _montecarlo_next_translation_vector.getZ());
    _montecarlo_current_angles += _montecarlo_next_angles;
}

// ---------------------------------------------------------------------------------------------------------------------
// Rigid body dynamics computation -------------------------------------------------------------------------------------

// Numerical solver
// Update _pos, _v, _omega_v
// Reinitialize _force & _torque
void RigidBody::solve()
{
    // Angle and depth insertion restraints to apply
    if (!_spn->isImpalaSamplingEnabled())
    {
        // Apply depth and external restraints got from visualizer and MDDriver
        applyDepthRestraint();
        if (_spn->isInsertionVectorEnabled())
            applyAngleRestraint();

        //applyHydrophobicRestraint();
    }
    // Automatic sampling solving
    if (_spn->isImpalaSamplingEnabled() && !_spn->isMonteCarloEnabled())
    {
        updateImpalaSampling(_spn->getNbIterations(), _dt); // Vary insertion angle and roll angle
        return; // Finish here rigid body solve for automatic sampling
    }
    // Monte Carlo solving
    else if (!_spn->isImpalaSamplingEnabled() && _spn->isMonteCarloEnabled())
    {
        solveMonteCarlo();
        return; // Finish here rigid body solve for automatic sampling
    }
    
    // Solve dynamic of rigid body (if not sampling of monte carlo)
    _force  += external_force;
    _torque += external_torque;
    external_force  = Vector3f(0., 0., 0.);
    external_torque = Vector3f(0., 0., 0.);

    // Linear movement of a rigid body -----------------------------------------
    
    // Linear acceleration
    _a = _force / _mass;
    
    // Linear velocity
    _v += _a * _dt;

    // Rigid body position
    _pos += _v * _dt + _a * 0.5 * _dt * _dt;

    // Angular momentum
    _L += _torque * _dt;

    _omega_v = (_ibodyinv * _L).toVector3f();

    Quaternion delta_orientation = Quaternion(0, 0, 0, 1);
    delta_orientation.fromAxisAngle(_omega_v, _omega_v.norm() * _dt);
    _orientation = delta_orientation * _orientation;
    _orientation.normalize();

    _force = Vector3f();
    _torque = Vector3f();
}

void RigidBody::applyDepthRestraint()
{
    if (applydepthrestraint)
    {
        double forceScale = (_pos.getZ() - depth_restraint);
        Vector3f fp = Vector3f(0,0,-forceScale);
        fp.normalize();
        _force += (fp * depth_restraint_scaling);
    }
}

void RigidBody::applyAngleRestraint()
{
    // Angle and insertion depth restraints
    InsertionVector iv = _spn->getInsertionVector();
    spn::Particle & p1 = iv.getParticle(0);
    spn::Particle & p2 = iv.getParticle(1);

    if (_spn->isInsertionVectorEnabled() && applyanglerestraint)
    {
        Vector3f dist = p2.getPosition() - p1.getPosition();
        double forceScale = (iv.getAngle() - angle_restraint) * 0.0174533; // Degree to rad
        Vector3f fp1 = dist^(Vector3f(0,0,forceScale))^(dist);
        Vector3f fp2 = dist^(Vector3f(0,0,-forceScale))^(dist);
        fp1.normalize() ; fp2.normalize();
        p1.addForce(fp1 * angle_restraint_scaling);
        p2.addForce(fp2 * angle_restraint_scaling);
        computeParticleForceAndTorque(p1);
        computeParticleForceAndTorque(p2);
    }
}

void RigidBody::applyHydrophobicRestraint()
{
    Vector3f force = Vector3f();
    Vector3f torque = Vector3f();
#ifdef OPENMP_SUPPORT
#pragma omp parallel default(shared)
#endif
{
#ifdef OPENMP_SUPPORT
#pragma omp for schedule(static)
#endif 
    for (unsigned i = 0; i < _selfHydrophilicParticles.size() ; i++)
    {
        spn::Particle & pi = _spn->getParticle(_selfHydrophilicParticles[i]);
        // Into membrane, consider only particle with positive transfert energy
        double tri = pi.getTransferEnergyByAccessibleSurface();
        double sasai = pi.getSolventAccessibilitySurface();
        if (sasai < 10 || abs(pi.getPosition().getZ()) > 13.5) continue;

        for (unsigned j = 0; j < _interactionHydrophilicParticles.size() ; j++)
        {
            spn::Particle & pj = _spn->getParticle(_interactionHydrophilicParticles[j]);
            double trj = pj.getTransferEnergyByAccessibleSurface();
            double sasaj = pj.getSolventAccessibilitySurface();
            if (sasaj < 10 || abs(pj.getPosition().getZ()) > 13.5) continue;
            Vector3f f = pi.getPosition() - pj.getPosition();
            double distance = f.norm();
            if (distance > 5) continue;
            double hydrophobicitymodule = 0.0;
            hydrophobicitymodule = (sasai * tri * sasaj * trj) * exp(-distance);
            //hydrophobicitymodule *= ForceField::AVOGADRONUMBER * 1.0E-3;
            //hydrophobicitymodule *= 1.0E-3;
            hydrophobicitymodule *= 10;
            f = f * hydrophobicitymodule;
            //pi.addf(f);
            //Vector3f visc = pi.getVelocity() * _spn->getViscosity();
            force += f;
            Vector3f t = (pi.getPosition() - _pos) ^ f;
            torque += t;
        }
    }
    //logging::info("force:%s ; torque:%s", _force.to_string().c_str(), _torque.to_string().c_str());
    _force += force;
    _torque += torque;
}
}

// Compute resultant force and torque that acts on the rigid system.
// Update total force and torque applied to rigid body
void RigidBody::computeParticleForceAndTorque(spn::Particle & p)
{
    if (RigidBodiesManager::getCollection().size() == 0) return;
    RigidBody *rb = RigidBodiesManager::getCollection()[p.getRigidBodyId()];
    Vector3f torque = Vector3f();
    Vector3f force = p.getForce();
    if (!rb->withoutRotational)
    {
        Vector3f t = (p.getPosition() - rb->_pos) ^ force;
        torque = torque + t;
    }
    rb->_force += force;
    rb->_torque += torque;
}

void RigidBody::resetRigidBodiesForceAndTorque()
{
#ifdef OPENMP_SUPPORT
#pragma omp parallel default(shared)
#endif
{
#ifdef OPENMP_SUPPORT
#pragma omp for schedule(static)
#endif 
    for (unsigned i = 0; i < RigidBodiesManager::getCollection().size() ; i++)
    {
        //logging::info("solve rbid:%d", RigidBodiesManager::getCollection()[i]->rbid);
        RigidBodiesManager::getCollection()[i]->_force = Vector3f();
        RigidBodiesManager::getCollection()[i]->_torque = Vector3f();
    }
}
}

void RigidBody::integrateParticleVelocity(spn::Particle & p, int ind, double timestep)
{
    if (RigidBodiesManager::getCollection().size() == 0) return;
    RigidBody *rb = RigidBodiesManager::getCollection()[p.getRigidBodyId()];
    Vector3f previousPosition = p.getPosition();
    p.setPreviousPosition(previousPosition);

    // Automatic sampling solving
    if (rb->_spn->isImpalaSamplingEnabled() && !rb->_spn->isMonteCarloEnabled())
    {
        Vector3f newPos = rb->getImpalaSamplingParticlePosition(p, ind);
        p.setPosition(newPos);
        return; // Finish here rigid body solve for automatic sampling
    }
    // Monte Carlo solving
    else if (!rb->_spn->isImpalaSamplingEnabled() && rb->_spn->isMonteCarloEnabled())
    {
        Vector3f newPos = rb->getMonteCarloParticlePosition(p, ind);
        p.setPosition(newPos);
        return; // Finish here rigid body solve for automatic sampling
    }
    else
    {
        Quaternion qnewLocalPos = rb->_orientation * Quaternion(rb->_p0[ind], 0.0) * rb->_orientation.inverse();
        Vector3f newPosition = rb->_pos + qnewLocalPos.getV();
        p.setPosition(newPosition);

        // Set particle velocity (usefull when applying viscosity)
        Vector3f velocity = (newPosition - previousPosition) / timestep;
        p.setVelocity(velocity);    

        p.resetForce();
    }
}

// ---------------------------------------------------------------------------------------------------------------------

Matrix RigidBody::createVector3Matrix()
{
    Matrix m = Matrix(1, 3);
    return m;
}

void RigidBody::updateVector3Matrix(Matrix &m, Vector3f v)
{
    m(0,0) = v.getX();
    m(0,1) = v.getY();
    m(0,2) = v.getZ();
}

} // namespace rigidbody
} // namespace biospring