
#ifndef __COMMANDLINE_H__
#define __COMMANDLINE_H__

#include <vector>
#include <string>
#include <ostream>

// Vector of strings to char ** conversion.
struct CommandLine
{
    std::vector<std::string> cmd;
    char ** argv;

    CommandLine(std::vector<std::string> a) : cmd(a)
    {
        argv = new char *[cmd.size()];
        for (size_t i = 0; i < cmd.size(); i++)
        {
            argv[i] = new char[cmd[i].size() + 1];
            argv[i][cmd[i].size()] = '\0';
            strncpy(argv[i], cmd[i].c_str(), cmd[i].size());
        }
    }

    virtual ~CommandLine()
    {
        for (size_t i = 0; i < cmd.size(); i++)
        {
            delete argv[i];
        }
        delete[] argv;
    }

    size_t size() const { return cmd.size(); }

    friend std::ostream & operator<<(std::ostream & os, const CommandLine & cmd)
    {
        os << "Command(\"";
        for (size_t i = 0; i < cmd.size(); i++)
        {
            os << cmd.cmd[i] << " ";
        }
        os << "\")";
        return os;
    }
};

#endif