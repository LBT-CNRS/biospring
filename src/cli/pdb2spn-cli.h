
#include "argparse.hpp"
#include "topology.hpp"

#include <string>
#include <vector>

namespace biospring
{

namespace pdb2spn
{

// Handles command-line arguments storage for pdb2spn.
class CommandLineArguments : argparse::CommandLineArgumentsBase
{
  public:
    // I/O paths.
    std::string pathTopology;
    std::string pathForceField;
    std::string pathGroup;
    std::string pathRigidBody;
    // Several, because the generated .bi.ff files are per polymer type
    // (protein, RNA, DNA) and a mixed system needs more than one. They must
    // NOT be concatenated by hand: each numbers its tables from 0, so a
    // `cat` silently collides and the later file's tables replace the
    // earlier's. BondedForceFieldReader shifts each file's ids instead.
    std::vector<std::string> pathsBondedInteraction;
    std::string pathStaticHydrogenBond;
    std::vector<std::string> pathOutputList;


    // User options.
    float cutoff;
    float stiffness;
    float charge;
    bool isStatic;
    bool addStaticDisulfide;
    bool ignoreDuplicates;
    bool ignoreMissing;
    bool writePdbConect;
    bool dihedral;
    bool dihedralBackbone;
    bool dihedralSidechain;
    bool dihedralPlanarity;

    // Constructor (inherited from argparse::CommandLineArgumentsBase).
    CommandLineArguments(const std::string & name, const argparse::description_t & description,
                         const std::string & version = "");

    // Implements abstract methods from parent class.
    void printArgumentValues() const;
    void parseCommandLine(int argc, const char * const argv[]);

    bool useUserCutoff() const { return _parser.get_option("--cutoff").is_set(); }
    bool useUserCharge() const { return _parser.get_option("--charge").is_set(); }
};

void reduceToCoarseGrain(topology::Topology & top, const CommandLineArguments & args);

int main(int argc, char ** argv);

} // namespace pdb2spn

} // namespace biospring
