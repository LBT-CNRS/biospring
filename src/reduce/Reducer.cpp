

#include "Reducer.h"
#include "IO/ForceFieldReader.h"
#include "IO/ReduceRuleReader.h"
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

} // namespace reduce
} // namespace biospring