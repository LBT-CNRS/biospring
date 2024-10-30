

#include "Reducer.h"
#include "IO/ForceFieldReader.h"
#include "IO/ReduceRuleReader.h"
#include "ParticleGroup.h"
#include "ReduceRule.h"
#include "logging.h"

namespace biospring
{
namespace reduce
{

void Reducer::initialize_forcefield(const std::string & path)
{
    logging::status("Reading FF file %s", path.c_str());
    ForceFieldReader reader(path);
    reader.read();
    _forcefield = reader.getForceField();
    _forcefield_initialized = true;
}

void Reducer::initialize_rules(const std::string & path)
{
    logging::status("Reading GRP file %s", path.c_str());
    ReduceRuleReader reader(path);
    reader.read();
    _rules = reader.rules();
    _rules_initialized = true;
}

// =====================================================================================
//
//                              LEGACY CODE
//
// =====================================================================================

namespace legacy
{

// Returns true if a particle should be added to a grain (according to the grain's rule).
static bool particleBelongsToGrain(const spn::Particle * const p, const ReduceRule * const rule)
{
    return rule->hasAtomNamed(p->getName());
}

static bool grainHasMissingParticles(const ParticleGroup * const grain, const ReduceRule * const rule)
{
    return grain->size() < rule->getNumberOfAtoms();
}

void Reducer::reduce(void) const
{
    // Dies if one of spring network, reduce or force field is not set.
    if (!isSpringNetworkSet())
        logging::die("Reducer::reduce: spring network not initialized or empty");
    if (!isReduceSet())
        logging::die("Reducer::reduce: reduce rules not initialized or empty");
    if (!isForceFieldSet())
        logging::die("Reducer::reduce: force field not initialized or empty");

    logging::status("Reducing model using group and forcefield");

    std::vector<spn::Particle *> cgParticles;
    std::vector<const spn::Particle *> _groupByAA;
    size_t previousResid = _spn->getParticle(0).getResId();
    std::string previousResName = _spn->getParticle(0).getResName();
    std::string previousChain = _spn->getParticle(0).getChainName();

    for (size_t i = 0; i < _spn->getNumberOfParticles(); i++)
    {
        const spn::Particle & p = _spn->getParticle(i);
        if (p.getResId() == previousResid and p.getChainName() == previousChain)
        {
            _groupByAA.push_back(&p);
        }
        else
        {
            _reduceAA(_groupByAA, cgParticles, previousResName, previousResid);
            _groupByAA.clear();
            _groupByAA.push_back(&p);
            previousResid = p.getResId();
            previousResName = p.getResName();
            previousChain = p.getChainName();
        }
    }
    _reduceAA(_groupByAA, cgParticles, previousResName, previousResid);

    // Replaces spring network's particles with coarse grain particles.
    _spn->clearParticles();
    for (size_t i = 0; i < cgParticles.size(); i++)
    {
        _spn->addParticle(*cgParticles[i]);
    }
}

void Reducer::_reduceAA(std::vector<const spn::Particle *> & particles, std::vector<spn::Particle *> & cgParticles,
                        const std::string & resname, const size_t resid) const
{
    std::vector<ReduceRule *> ruleList = _reduce->getReduceRulesByResName(resname);

    // Does nothing if found no rule for this residue.
    if (ruleList.size() == 0)
    {
        logging::warning("Residue '%s' not found in reduce rules. Skipping it.", resname.c_str());
        return;
    }

    // Detects particle that we don't know about, i.e. doesn't belong to any grain
    // according to reduce rules.
    for (const spn::Particle * const p : particles)
    {
        bool belongToAGrain = false;
        for (const ReduceRule * const rule : ruleList)
        {
            if (particleBelongsToGrain(p, rule))
            {
                belongToAGrain = true;
                break;
            }
        }
        if (not belongToAGrain)
        {
            logging::warning("Residue %s:%d: Unexpected particle '%s' found in topology", resname.c_str(), resid,
                             p->getName().c_str());
        }
    }

    // For each grain defined in the rules.
    for (const ReduceRule * const rule : ruleList)
    {
        bool addGrainToModel = true;

        ParticleGroup * grain = new ParticleGroup();
        grain->setName(rule->getName());
        grain->setCGName(rule->getName());
        grain->setResName(resname);

        grain->setExtid(resid);
        grain->setResId(resid);

        // Loop over the residue's particles.
        for (const spn::Particle * const p : particles)
        {
            if (particleBelongsToGrain(p, rule))
            {
                grain->addParticle(p);
            }
        }

        // No particle in grain at all.
        if (grain->size() == 0)
        {
            logging::warning("Residue %s:%d: No particle found to create grain %s...skipping it", resname.c_str(),
                             resid, grain->getName().c_str());
            addGrainToModel = false;
        }

        // Not enough particles in grain.
        else if (grainHasMissingParticles(grain, rule))
        {
            for (const std::string & name : rule->getAtomNames())
            {
                if (not grain->contains(name))
                {
                    if (_ignoreMissingParticle)
                    {
                        logging::warning("Residue %s:%d: Particle '%s' required for grain '%s' not found in "
                                         "topology...still adding grain to model",
                                         resname.c_str(), resid, name.c_str(), grain->getName().c_str());
                    }
                    else
                    {
                        logging::warning("Residue %s:%d: Particle '%s' required for grain '%s' not found in "
                                         "topology...skipping this grain",
                                         resname.c_str(), resid, name.c_str(), grain->getName().c_str());
                        addGrainToModel = false;
                    }
                }
            }
        }

        // Too many particles in grain.
        else if (grain->size() > rule->getNumberOfAtoms())
        {
            if (_ignoreDuplicateParticles)
            {
                logging::warning("Residue %s:%d: Duplicate particles in grain %s...still adding grain to model",
                                 resname.c_str(), resid, grain->getName().c_str());
            }
            else
            {
                logging::warning("Residue %s:%d: Duplicate particles in grain %s...skipping this grain",
                                 resname.c_str(), resid, grain->getName().c_str());
                addGrainToModel = false;
            }
        }

        // Number of particles in grain is as expected.
        if (addGrainToModel)
        {
            grain->reduce();
            grain->updateFromForceField(_ff);
            cgParticles.push_back(grain);
        }
    }
}

//
// Sets force field using path to force field setup file.
//
// Reads file and initializes force field.
//
void Reducer::setForceField(const std::string & path)
{
    logging::status("Reading FF file %s", path.c_str());
    ForceFieldReader ffreader;
    ffreader.setFileName(path);
    ffreader.read();
    _ff = ffreader.getForceField();
}

//
// Sets reduce manager using reduce setup file.
//
// Reads file and initializes reduce manager.
//
void Reducer::setReduce(const std::string & path)
{
    logging::status("Reading GRP file %s", path.c_str());
    ReduceRuleReader reducereader;
    reducereader.setFileName(path);
    reducereader.read();
    _reduce = reducereader.getReduce();
}

//
// Reduce a spring network model using reduction parameters.
//
void Reducer::reduceToCoarseGrain(spn::SpringNetwork & spn, const ReductionParameters & parameters)
{
    // Dies if reduce parameter file is not set.
    if (parameters.pathForceField.empty())
        logging::die("Reducer: force field parameter file not provided");

    // Dies if force field parameter file is not set.
    if (parameters.pathGroup.empty())
        logging::die("Reducer: coarse grain parameter file not provided");

    if (not parameters.pathGroup.empty())
    {
        setSpringNetwork(&spn);
        setForceField(parameters.pathForceField);
        setReduce(parameters.pathGroup);
        setIgnoreDuplicateParticles(parameters.ignoreDuplicate);
        setIgnoreMissingParticles(parameters.ignoreMissing);
        reduce();
    }
}
} // namespace legacy
} // namespace reduce
} // namespace biospring