

#include <string>
#include <unordered_map>
#include <vector>

#include "argparse.hpp"
#include "topology.hpp"

namespace biospring
{

namespace editspn
{

//
// Stores a spring descriptors:
//   - p1: first particle string representation (see Particle::tostr())
//   - p2: second particle string representation (see Particle::tostr())
//   - stiffness: spring stiffness
//   - equilibrium: distance between the two particules where the spring is at equilibrium
//
struct SpringDescriptor
{
    std::string p1;
    std::string p2;
    float stiffness;
    float equilibrium;
};

//
// Handles command-line parsing as well as argument storage.
//
class CommandLineArguments : public argparse::CommandLineArgumentsBase
{
  public:
    // I/O paths.
    std::string pathTopology;
    std::vector<std::string> pathOutputList;
    std::string pathUpdateSprings;

    // User options.
    float cutoff;    // cutoff value for spring creation
    float stiffness; // stiffness for spring creation

    // Constructor (inherited from argparse::CommandLineArgumentsBase).
    CommandLineArguments(const std::string & name, const argparse::description_t & description,
                         const std::string & version = "");

    // Implements abstract methods from parent class.
    void printArgumentValues() const;
    void parseCommandLine(int argc, const char * const argv[]);

    bool useUserCutoff() const { return _parser.get_option("--cutoff").is_set(); }
    bool useUserStiffness() const { return _parser.get_option("--stiffness").is_set(); }
};

void readSpringFile(const std::string & filename, std::vector<SpringDescriptor> & sdesc);
void addSpringsToSpn(const std::vector<SpringDescriptor> & sdesc, biospring::topology::Topology & top);
void addSpringsToSpn(const std::string & path, biospring::topology::Topology & top);

void checkParticleIdFormat(const std::string & pid, const unsigned & lineid);
void checkParticleInSpn(const std::string & pid, const std::unordered_map<std::string, unsigned> & pidmap);

int main(int argc, char ** argv);

} // namespace editspn
} // namespace biospring
