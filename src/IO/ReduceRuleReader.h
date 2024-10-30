#ifndef _REDUCEREADER_H_
#define _REDUCEREADER_H_

#include <string>
#include <vector>

#include "IO/ReaderBase.h"
#include "reduce/Reduce.h"
#include "reduce/ReduceRule.h"
#include "reduce/ReduceRuleContainer.hpp"

namespace biospring
{
namespace reduce
{

class ReduceRuleReader : public ReaderBase
{
  protected:
    ReduceRuleContainer _rules;

  public:
    ReduceRuleReader() : ReaderBase() {}
    ReduceRuleReader(const std::string & path) : ReaderBase(path) {}
    ReduceRuleReader(const char * const path) : ReaderBase(path) {}

    const ReduceRuleContainer & rules(void) const { return _rules; }

    void read();

  protected:
    ReduceRule parse_rule(const std::string & line, size_t line_id) const;
};

namespace legacy
{

class ReduceRuleReader : public ReaderBase
{
  public:
    ReduceRuleReader() : ReaderBase(), _reduce(nullptr) {}
    ReduceRuleReader(const std::string & path) : ReaderBase(path), _reduce(nullptr) {}
    ReduceRuleReader(const char * const path) : ReaderBase(path), _reduce(nullptr) {}

    ~ReduceRuleReader() {}

    Reduce * getReduce() const { return _reduce; }
    void setReduce(Reduce * const red) { _reduce = red; }

    void read();

  private:
    Reduce * _reduce;
};

} // namespace legacy
} // namespace reduce
} // namespace biospring

#endif
