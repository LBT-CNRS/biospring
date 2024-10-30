
#ifndef __argparse_HPP__
#define __argparse_HPP__

#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../logging.h"


namespace biospring
{

namespace argparse
{

using description_t = std::vector<std::string>;

class CommandLineParser;
class Argument;

// From https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
// The code snippet above is licensed under CC0 1.0.
template<typename ... Args>
std::string string_format(const std::string & format, Args ... args)
{
    int size_s = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
    if( size_s <= 0 ){ throw std::runtime_error("Error during formatting."); }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf(new char[size]);
    std::snprintf(buf.get(), size, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}


inline bool string_to_bool(const std::string & s)
{
    if (s == "true" || s == "True" || s == "1")
        return true;
    if (s == "false" || s == "False" || s == "0")
        return false;
    throw std::runtime_error("Cannot convert string '" + s + "' to boolean.");
}


enum class NumberOfArgumentsType
{
    EXACTLY,
    AT_LEAST,
    AT_MOST,
};

struct NumberOfArguments
{
    size_t value;
    NumberOfArgumentsType type;

    NumberOfArguments(size_t value, NumberOfArgumentsType type = NumberOfArgumentsType::EXACTLY)
        : value(value), type(type)
    {
    }

    bool is_multiple() const
    {
        return value > 1 || (type == NumberOfArgumentsType::AT_LEAST);
    }
};

enum class ActionType
{
    UNSET,
    STORE,
    STORE_TRUE,
    STORE_FALSE,
    APPEND,
    APPEND_CONST,
    COUNT,
};

enum class ArgumentType
{
    STRING,
    INTEGER,
    REAL,
    BOOLEAN,
    PATH_INPUT,
    PATH_OUTPUT,
};

class Action
{
  public:
    virtual void take_action(CommandLineParser & parser, Argument & arg, size_t & position) = 0;
    virtual ~Action() {}
};


class StoreTrueAction : public Action
{
    public:
        void take_action(CommandLineParser & parser, Argument & arg, size_t & position) override;
};

class StoreFalseAction : public Action
{
    public:
        void take_action(CommandLineParser & parser, Argument & arg, size_t & position) override;
};

class StoreAction : public Action
{
    public:
        void take_action(CommandLineParser & parser, Argument & arg, size_t & position) override;
};

class AppendAction : public Action
{
    public:
        void take_action(CommandLineParser & parser, Argument & arg, size_t & position) override;
};


// Class storing information about a command line argument.
// Uses the builder pattern.
class Argument
{
  protected:
    std::string _name_short;
    std::string _name_long;
    std::string _description;
    std::string _metavar;
    std::string _default_value;
    NumberOfArguments _number_of_arguments;
    ActionType _action_type;
    ArgumentType _argument_type;
    bool _required;
    bool _is_set;
    std::vector<std::string> _values;
    std::unique_ptr<Action> _action;

  public:
    Argument():
        _name_short(""),
        _name_long(""),
        _description(""),
        _metavar(""),
        _default_value(""),
        _number_of_arguments(NumberOfArguments(0)),
        _action_type(ActionType::UNSET),
        _argument_type(ArgumentType::STRING),
        _required(false),
        _is_set(false)
    {
    }

    Argument(const Argument & source):
        _name_short(source._name_short),
        _name_long(source._name_long),
        _description(source._description),
        _metavar(source._metavar),
        _default_value(source._default_value),
        _number_of_arguments(source._number_of_arguments),
        _action_type(source._action_type),
        _argument_type(source._argument_type),
        _required(source._required),
        _is_set(source._is_set),
        _values(source._values)
    {
        this->action_type(source._action_type);
    }

    // Assignment operator.
    Argument & operator=(const Argument & source)
    {
        _name_short = source._name_short;
        _name_long = source._name_long;
        _description = source._description;
        _metavar = source._metavar;
        _default_value = source._default_value;
        _number_of_arguments = source._number_of_arguments;
        _argument_type = source._argument_type;
        _required = source._required;
        _is_set = source._is_set;
        _values = source._values;

        this->action_type(source._action_type);
        return *this;
    }

    void set_value(const std::string & value) { _values = {value}; _is_set = true; }


    // == Returns the value of the argument ============================================
    template<typename T> T get_value() const;


    // == Returns the value of the argument (array) ====================================
    template<typename T> std::vector<T> get_values() const;

    // ================================================================================



    // Appends a value to the argument.
    // When the number of argument is AT_MOST or EXACTLY, the `is_set` flag is set to true
    // when the number of values is equal to the number of arguments.
    // When the number of argument is AT_LEAST, the `Action::take_action` method
    // is responsible for setting the `is_set` flag.
    void append_value(const std::string & value)
    {
        if (is_set())
            throw std::runtime_error("Cannot append value to argument that has already been set.");
        _values.push_back(value);
        auto nargs_type = number_of_arguments().type;

        if ((nargs_type == NumberOfArgumentsType::AT_MOST || nargs_type == NumberOfArgumentsType::EXACTLY) && _values.size() == number_of_arguments().value)
            _is_set = true;
    }

    void take_action(CommandLineParser & parser, Argument & arg, size_t & position)
    {
        if (_action)
            _action->take_action(parser, arg, position);
        else
            throw std::runtime_error("No action defined for argument " + get_name());
    }

    bool is_store_type() const
    {
        return _action_type == ActionType::STORE || _action_type == ActionType::APPEND;
    }

    bool is_flag_type() const
    {
        return _action_type == ActionType::STORE_TRUE || _action_type == ActionType::STORE_FALSE;
    }

    static bool is_option(const std::string & name)
    {
        return (name != "" && name[0] == '-') || (name != "" && name.substr(0, 2) == "--");
    }

    bool is_option() const
    {
        return is_option(_name_short) || is_option(_name_long);
    }

    bool is_positional() const
    {
        return !is_option();
    }

    bool is_valid() const
    {
        if (_name_short == "" && _name_long == "")
            return false;

        if (is_option())
        {
            if (_name_short != "" && _name_short.size() < 2)
                return false;

            if (_name_long != "" && _name_long.size() < 3)
                return false;

            if (_name_long != "" && _name_long.substr(0, 2) != "--")
                return false;
        }
        return true;
    }

    std::string get_name() const
    {
        if (_name_short != "")
            return _name_short;
        return _name_long;
    }

    std::string get_full_name() const
    {
        std::string name = "";

        if (is_option())
        {
            if (_name_short != "")
                name += _name_short;
            if (_name_long != "")
            {
                if (_name_short != "")
                    name += ", ";
                name += _name_long;
            }
        }
        else
        {
            name += _name_short;
        }

        return name;
    }

    std::string get_metavar() const
    {
        if (_metavar != "")
            return _metavar;
        if (_name_long == "")
            return "";
        std::string metavar = _name_long.substr(2);
        std::replace(metavar.begin(), metavar.end(), '-', '_');
        std::transform(metavar.begin(), metavar.end(), metavar.begin(), ::toupper);
        return metavar;
    }

    std::string get_usage_string() const
    {
        std::string usage = get_name();

        if (is_option())
        {
            std::string metavar = get_metavar();
            if (is_store_type())
            {
                if (metavar != "")
                    usage += " " + metavar;

                if (_number_of_arguments.is_multiple())
                    usage += " ...";
            }

            if (!required())
                usage = "[" + usage + "]";
        }

        // positional argument
        else
        {
            if (!required())
                usage = "[" + usage + "]";
            else
                usage = "<" + usage + ">";
        }

        return usage;
    }

    std::string get_help_string() const
    {
        std::string name_buf = "";

        if (is_option())
        {
            if (_name_short != "")
                name_buf += _name_short;
            if (_name_long != "")
            {
                if (_name_short != "")
                    name_buf += ", ";
                name_buf += _name_long;
            }
        }
        else
        {
            name_buf += _name_short;
        }

        std::string type_buf;
        switch (_argument_type)
        {
            case ArgumentType::STRING:
                type_buf = "string";
                break;
            case ArgumentType::REAL:
                type_buf = "real";
                break;
            case ArgumentType::INTEGER:
                type_buf = "integer";
                break;
            case ArgumentType::BOOLEAN:
                type_buf = "flag";
                break;
            case ArgumentType::PATH_INPUT:
                type_buf = "Input";
                if (!_required) type_buf += ", Opt.";
                if (_number_of_arguments.is_multiple()) type_buf += ", Mult.";
                break;
            case ArgumentType::PATH_OUTPUT:
                type_buf = "Output";
                if (!_required) type_buf += ", Opt.";
                if (_number_of_arguments.is_multiple()) type_buf += ", Mult.";
                break;
        }

        std::string description_buf = _description;

        std::string default_buf = _default_value;

        std::string help = string_format(
            " %-20s %-20s %-15s %-s",
            name_buf.c_str(),
            type_buf.c_str(),
            default_buf.c_str(),
            description_buf.c_str());

        return help;
    }


    // Getters //////////////////////////////////////////////////////////////////

    std::string name_short() const { return _name_short; }
    std::string name_long() const { return _name_long; }
    std::string description() const { return _description; }
    std::string metavar() const { return _metavar; }
    std::string default_value() const { return _default_value; }
    NumberOfArguments number_of_arguments() const { return _number_of_arguments; }
    ActionType action_type() const { return _action_type; }
    ArgumentType argument_type() const { return _argument_type; }
    bool required() const { return _required; }
    bool is_set() const { return _is_set; }
    std::vector<std::string> values() const { return _values; }

    // Setters //////////////////////////////////////////////////////////////////

    Argument & name_short(const std::string & name_short)
    {
        _name_short = name_short;
        return *this;
    }

    Argument & name_long(const std::string & name_long)
    {
        _name_long = name_long;
        return *this;
    }

    Argument & description(const std::string & description)
    {
        _description = description;
        return *this;
    }

    Argument & metavar(const std::string & metavar)
    {
        _metavar = metavar;
        return *this;
    }

    Argument & default_value(const std::string & default_value)
    {
        _default_value = default_value;
        return *this;
    }

    // Sets the number of arguments.
    // Also adapts the action type adequatly.
    Argument & number_of_arguments(const NumberOfArguments & number_of_arguments)
    {
        _number_of_arguments = number_of_arguments;
        if (_number_of_arguments.is_multiple())
            _action_type = ActionType::APPEND;
        else if (_number_of_arguments.value == 1)
            _action_type = ActionType::STORE;
        return *this;
    }

    Argument & number_of_arguments(const size_t n)
    {
        return number_of_arguments(NumberOfArguments(n, NumberOfArgumentsType::EXACTLY));
    }

    Argument & number_of_arguments(const std::string & n)
    {
        if (n == "+")
            return number_of_arguments(NumberOfArguments(1, NumberOfArgumentsType::AT_LEAST));
        else if (n == "?")
            return number_of_arguments(NumberOfArguments(1, NumberOfArgumentsType::AT_MOST));
        else if (n == "*")
            return number_of_arguments(NumberOfArguments(0, NumberOfArgumentsType::AT_LEAST));
        return number_of_arguments(NumberOfArguments(std::stoi(n), NumberOfArgumentsType::EXACTLY));
    }

    Argument & action_type(const ActionType & action_type)
    {
        if (number_of_arguments().is_multiple() && action_type != ActionType::APPEND)
            throw std::runtime_error("Argument that accept multiple values have to be of type ActionType::APPEND.");

        _action_type = action_type;
        switch (_action_type)
        {
            case ActionType::STORE:
                _action = std::make_unique<StoreAction>();
                break;
            case ActionType::APPEND:
                _action = std::make_unique<AppendAction>();
                break;
            case ActionType::STORE_TRUE:
                _action = std::make_unique<StoreTrueAction>();
                break;
            case ActionType::STORE_FALSE:
                _action = std::make_unique<StoreFalseAction>();
                break;

            default:
                throw std::runtime_error(get_full_name() + ": Unknown action type.");
                break;
        }
        return *this;
    }

    Argument & argument_type(const ArgumentType & argument_type)
    {
        _argument_type = argument_type;
        if (_action_type == ActionType::UNSET)
        {
            switch (_argument_type)
            {
                case ArgumentType::STRING:
                    _action_type = ActionType::STORE;
                    break;
                case ArgumentType::REAL:
                    _action_type = ActionType::STORE;
                    break;
                case ArgumentType::INTEGER:
                    _action_type = ActionType::STORE;
                    break;
                case ArgumentType::BOOLEAN:
                    _action_type = ActionType::STORE_TRUE;
                    break;
                case ArgumentType::PATH_INPUT:
                    _action_type = ActionType::STORE;
                    break;
                case ArgumentType::PATH_OUTPUT:
                    _action_type = ActionType::STORE;
                    break;
            }
        }

        if ((_argument_type == ArgumentType::PATH_INPUT || _argument_type == ArgumentType::PATH_OUTPUT) && _number_of_arguments.value == 0)
            _number_of_arguments = NumberOfArguments(1, NumberOfArgumentsType::EXACTLY);
        return *this;
    }

    Argument & required(const bool & required)
    {
        _required = required;
        return *this;
    }

    Argument & is_set(const bool & is_set)
    {
        _is_set = is_set;
        return *this;
    }


  protected:

    // Returns the vector of values.
    // If the argument is not set, returns the default value.
    std::vector<std::string> _get_string_values() const
    {
        std::vector<std::string> values;
        if (!is_set())
            values.push_back(_default_value);
        else
            values = _values;
        return values;
    }
};



Argument StoreTrueArgument(const std::string & name_short, const std::string & name_long, const std::string & description);
Argument HelpArgument();
Argument VersionArgument();




class CommandLineParser
{
  protected:
    std::string _program_name;
    std::vector<std::string> _program_description;
    std::string _program_version;

    std::vector<Argument> _arguments;
    std::vector<Argument> _options;

    std::vector<std::string> _cl_arguments;

    // Special arguments.
    Argument _help_argument;
    Argument _version_argument;


    // Returns true if and argument is required from the command-line.
    // Important: does not check if the argument has been initialized.
    bool _argument_is_required(const Argument & argument) const
    {
        for (const auto & arg : _cl_arguments)
            if (arg == argument.name_short() || arg == argument.name_long())
                return true;
        return false;
    }

  public:
    CommandLineParser(const std::string & program_name, const std::vector<std::string> & program_description, const std::string & program_version = "")
        : _program_name(program_name), _program_description(program_description), _program_version(program_version)
    {
        _help_argument = HelpArgument();
        if (program_version != "")
            _version_argument = VersionArgument();
    }

    CommandLineParser(const std::string & program_name, const std::string & program_description, const std::string & program_version = "")
        : _program_name(program_name), _program_description(), _program_version(program_version)
    {
        _program_description.push_back(program_description);
    }

    // Returns an option value.
    template <typename T>
    T get_option_value(const std::string & name) const
    {
        for (const Argument & opt : _options)
            if (opt.name_short() == name || opt.name_long() == name)
                return opt.get_value<T>();
        throw std::runtime_error(string_format("Option '%s' does not exist", name.c_str()));
    }

    template<typename T> std::vector<T> get_option_values(const std::string & name) const
    {
        for (const Argument & opt : _options)
            if (opt.name_short() == name || opt.name_long() == name)
                return opt.get_values<T>();
        throw std::runtime_error(string_format("Option '%s' does not exist", name.c_str()));
    }

    // Returns true if an given option was provided (i.e. it was set).
    bool option_was_provided(const std::string & name) const
    {
        for (const Argument & opt : _options)
            if (opt.name_short() == name || opt.name_long() == name)
                return opt.is_set();
        throw std::runtime_error(string_format("Option '%s' does not exist", name.c_str()));
    }

    const std::string get_command_line_argument(const size_t & position) const
    {
        return _cl_arguments[position];
    }

    // Returns the number of command line arguments.
    const size_t get_number_of_command_line_arguments() const
    {
        return _cl_arguments.size();
    }

    // Returns true if the argument at pos + 1 exists.
    bool has_next_argument(size_t pos) const
    {
        return pos + 1 < _cl_arguments.size();
    }

    // Returns true if the help argument has been initialized.
    bool has_help_option() const
    {
        return _help_argument.is_valid();
    }

    // Returns true if the version argument has been initialized.
    bool has_version_option() const
    {
        return _version_argument.is_valid();
    }

    // Returns true if the parser has an option with the given name.
    bool has_option(const std::string & name) const
    {
        for (const Argument & opt : _options)
            if (opt.name_short() == name || opt.name_long() == name)
                return true;
        return false;
    }

    // Returns true if the parser has a positional argument with the given id.
    bool has_positional_argument(const size_t & id) const
    {
        return id < _arguments.size();
    }


    // Returns true if the argument at pos + 1 is an option.
    bool next_argument_is_option(size_t pos) const
    {
        if (!has_next_argument(pos))
            return false;
        return Argument::is_option(_cl_arguments[pos + 1]);
    }

    // Registers an argument.
    void add_argument(const Argument & argument)
    {
        if (!argument.is_valid())
            throw std::runtime_error(string_format("Invalid argument '%s'", argument.get_name().c_str()));
        if (has_been_registered(argument))
            throw std::runtime_error(string_format("Argument '%s' has already been registered", argument.get_name().c_str()));

        if (argument.is_option())
            _options.push_back(argument);
        else
            _arguments.push_back(argument);
    }

    // Returns true if the argument has been registered.
    bool has_been_registered(const Argument & arg) const
    {
        for (const Argument & argument : _options)
            if (!arg.name_short().empty() && (argument.name_short() == arg.name_short() || argument.name_long() == arg.name_long()))
                return true;
        for (const Argument & argument : _arguments)
            if (!arg.name_short().empty() && (argument.name_short() == arg.name_short() || argument.name_long() == arg.name_long()))
                return true;
        return false;
    }

    // =================================================================================

    // Creates and returns description string.
    std::string get_description_string() const
    {
        std::string desc = "";
        for (auto & line : _program_description)
            desc += line + "\n";
        return desc;
    }

    // Creates and returns usage string.
    std::string get_usage_string() const
    {
        std::string usage = "Usage: " + _program_name;

        if (has_help_option())
            usage += " " + _help_argument.get_usage_string();

        if (has_version_option())
            usage += " " + _version_argument.get_usage_string();

        // Shows options first.
        for (auto & argument : _options)
        {
            if (not argument.required())
                usage += " " + argument.get_usage_string();
        }

        // Shows required options second.
        for (auto & argument : _options)
        {
            if (argument.required())
                usage += " " + argument.get_usage_string();
        }

        // Finally shows arguments.
        for (auto & argument : _arguments)
        {
            usage += " " + argument.get_usage_string();
        }

        return usage;
    }

    // Creates and returns version string.
    std::string get_version_string() const
    {
        return _program_name + " " + _program_version;
    }

    // Creates and returns help string.
    std::string get_help_string() const
    {
        std::string help = get_description_string() + "\n\n" + get_usage_string() + "\n\n";

        help += "Option                Type                 Default         Description\n";
        help += "-----------------------------------------------------------------------------\n";

        // Shows positional arguments first.
        for (auto & argument : _arguments)
        {
            help += argument.get_help_string() + "\n";
        }

        // Then shows required options.
        for (auto & argument : _options)
        {
            if (argument.required())
                help += argument.get_help_string() + "\n";
        }

        // Then shows non-required options.
        for (auto & argument : _options)
        {
            if (not argument.required())
                help += argument.get_help_string() + "\n";
        }

        // Then show version and help arguments.
        if (has_version_option())
            help += _version_argument.get_help_string() + "\n";
        if (has_help_option())
            help += _help_argument.get_help_string() + "\n";

        return help;
    }

    // =================================================================================

    // Prints help to output string stream (default: std::cerr).
    void print_help(std::ostream & out = std::cerr) const
    {
        out << get_help_string() << std::endl;
    }

    // Prints version to output string stream (default: std::cerr).
    void print_version(std::ostream & out = std::cerr) const
    {
        out << _program_version << std::endl;
    }

    // =================================================================================

    // Parses the command line arguments.
    void parse_arguments(int argc, const char * const argv[])
    {
        // Stores command line arguments for later use.
        for (int i = 1; i < argc; i++)
            _cl_arguments.push_back(argv[i]);

        // If help is required, prints help and exits.
        if (has_help_option() && help_is_required())
        {
            print_help();
            exit(0);
        }

        // If version is required, prints version and exits.
        if (has_version_option() && version_is_required())
        {
            print_version();
            exit(0);
        }

        // Parses the command line arguments.
        size_t positional_id = 0;
        for (size_t i = 0; i < _cl_arguments.size(); i++)
        {
            const std::string arg(_cl_arguments[i]);
            if (Argument::is_option(arg))
            {
                if (!has_option(arg))
                    die(string_format("Unknown option '%s'", arg.c_str()));
                Argument & option = get_option(arg);
                option.take_action(*this, option, i);
            }
            else
            {
                if (!has_positional_argument(positional_id))
                    die(string_format("Unknown positional argument '%s'", arg.c_str()));

                Argument & positional = get_positional(positional_id++);
                positional.take_action(*this, positional, i);
            }
        }

        // Checks if all required arguments have been provided.
        for (const Argument & arg : _options)
        {
            if (arg.required() && !arg.is_set())
                die(string_format("Required option '%s' has not been provided", arg.get_name().c_str()));
        }

        for (const Argument & arg : _arguments)
        {
            if (arg.required() && !arg.is_set())
                die(string_format("Required positional argument '%s' has not been provided", arg.get_name().c_str()));
        }
    }

    // =================================================================================

    // Returns the option with the given name.
    std::vector<Argument>::reference get_option(const std::string & name)
    {
        for (Argument & opt : _options)
            if (opt.name_short() == name || opt.name_long() == name)
                return opt;
        throw std::runtime_error("Option '" + name + "' not found.");
    }

    // Returns an option with the given name (const version).
    std::vector<Argument>::const_reference get_option(const std::string & name) const
    {
        for (const Argument & opt : _options)
            if (opt.name_short() == name || opt.name_long() == name)
                return opt;
        throw std::runtime_error("Option '" + name + "' not found.");
    }

    // Return the positional argument with the given name.
    std::vector<Argument>::reference get_positional(const std::string & name)
    {
        for (Argument & arg : _arguments)
            if (arg.name_short() == name || arg.name_long() == name)
                return arg;
        throw std::runtime_error("Positional argument '" + name + "' not found.");
    }

    // Return the positional argument with the given name (const version).
    std::vector<Argument>::const_reference get_positional(const std::string & name) const
    {
        for (const Argument & arg : _arguments)
            if (arg.name_short() == name || arg.name_long() == name)
                return arg;
        throw std::runtime_error("Positional argument '" + name + "' not found.");
    }

    // Returns the positional argument with the given id.
    std::vector<Argument>::reference get_positional(const size_t & id)
    {
        if (id >= _arguments.size())
            throw std::runtime_error("Positional argument with id " + std::to_string(id) + " not found.");
        return _arguments[id];
    }

    // Returns the positional argument with the given id (const version).
    std::vector<Argument>::const_reference get_positional(const size_t & id) const
    {
        if (id >= _arguments.size())
            throw std::runtime_error("Positional argument with id " + std::to_string(id) + " not found.");
        return _arguments[id];
    }

    // =================================================================================

    // Prints error message to output string stream (default: std::cerr) and exists
    // with given status.
    void die(const std::string & msg, int status = 1, std::ostream & out = std::cerr) const
    {
        out << logging::ERROR_COLOR << logging::ERROR_PREFIX << msg << logging::RESET_COLOR << std::endl;
        exit(status);
    }

    // Returns true if help argument has been provided.
    // Important: does not check if the help argument has been initialized.
    bool help_is_required() const
    {
        return _argument_is_required(_help_argument);
    }

    // Returns true if version argument has been provided.
    // Important: does not check if the version argument has been initialized.
    bool version_is_required() const
    {
        return _argument_is_required(_version_argument);
    }
};

// Base class for command-line argument handling.
// Should be inherited by every program.
class CommandLineArgumentsBase
{
  public:
    CommandLineArgumentsBase(const std::string & name, const description_t & description, const std::string & version = "")
        : _parser(name, description, version)
    {
    }

    virtual void printArgumentValues() const = 0;
    virtual void parseCommandLine(int argc, const char * const argv[]) = 0;

  protected:
    CommandLineParser _parser;
};


} // namespace argparse
} // namespace biospring

#endif // __argparse_HPP__