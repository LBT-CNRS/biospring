#include "argparse.hpp"

namespace biospring
{

namespace argparse
{

// =====================================================================================
// Template specializations

// Returns the value of the argument.
template<> std::string Argument::get_value<std::string>() const
{
    if (!is_set())
        return _default_value;
    return _values[0];
}

template<> int Argument::get_value<int>() const
{
    if (!is_set())
        return std::stoi(_default_value);
    return std::stoi(_values[0]);
}

template<> unsigned Argument::get_value<unsigned>() const
{
    if (!is_set())
        return std::stoul(_default_value);
    return std::stoul(_values[0]);
}

template<> double Argument::get_value<double>() const
{
    if (!is_set())
        return std::stod(_default_value);
    return std::stod(_values[0]);
}

template<> float Argument::get_value<float>() const
{
    if (!is_set())
        return std::stof(_default_value);
    return std::stof(_values[0]);
}

template<> bool Argument::get_value<bool>() const
{
    if (!is_set())
        return string_to_bool(_default_value);
    return string_to_bool(_values[0]);
}

template<> std::vector<std::string> Argument::get_value<std::vector<std::string>>() const
{
    if (!is_set())
        return {_default_value};
    return _values;
}

// Returns the value of the argument (array versions).
template<> std::vector<int> Argument::get_values() const
{
    std::vector<std::string> string_values = _get_string_values();
    std::vector<int> values;
    for (auto & value : string_values)
        values.push_back(std::stoi(value));
    return values;
}

template<> std::vector<double> Argument::get_values() const
{
    std::vector<std::string> string_values = _get_string_values();
    std::vector<double> values;
    for (auto & value : string_values)
        values.push_back(std::stod(value));
    return values;
}

template<> std::vector<bool> Argument::get_values() const
{
    std::vector<std::string> string_values = _get_string_values();
    std::vector<bool> values;
    for (auto & value : string_values)
        values.push_back(string_to_bool(value));
    return values;
}

template<> std::vector<std::string> Argument::get_values() const
{
    return _get_string_values();
}

// =====================================================================================


void StoreTrueAction::take_action(CommandLineParser & parser, Argument & arg, size_t & position)
{
    arg.set_value("true");
}


void StoreFalseAction::take_action(CommandLineParser & parser, Argument & arg, size_t & position)
{
    arg.set_value("false");
}



void die_no_argument(const CommandLineParser & parser, const Argument & arg)
{
    if (arg.number_of_arguments().type == NumberOfArgumentsType::AT_LEAST)
        parser.die(string_format("Option '%s' requires at least %d argument", arg.get_name().c_str(), arg.number_of_arguments().value));
    else if (arg.number_of_arguments().type == NumberOfArgumentsType::EXACTLY)
    {
        if (arg.number_of_arguments().value == 1)
            parser.die(string_format("Option '%s' requires an argument", arg.get_name().c_str()));
        else
            parser.die(string_format("Option '%s' requires %d arguments", arg.get_name().c_str(), arg.number_of_arguments().value));
    }
    parser.die(string_format("Option '%s' requires an argument", arg.get_name().c_str()));
}


void StoreAction::take_action(CommandLineParser & parser, Argument & arg, size_t & position)
{
    std::string value;
    if (arg.is_option())
    {
        if (!parser.has_next_argument(position))
            die_no_argument(parser, arg);
        value = parser.get_command_line_argument(++position);
    }

    else
    {
        value = parser.get_command_line_argument(position);
    }

    arg.set_value(value);
}


void AppendAction::take_action(CommandLineParser & parser, Argument & arg, size_t & position)
{
    // Stores the number of values that are required for this option.
    size_t number_of_values = 0;

    if (arg.is_set())
        parser.die("duplicate option '" + arg.get_name() + "'");

    while (!arg.is_set() && position < parser.get_number_of_command_line_arguments())
    {
        std::string value;

        // If arg is an option, get the value right after the current position.
        if (arg.is_option())
        {
            if (number_of_values == 0 && !parser.has_next_argument(position))
                die_no_argument(parser, arg);
            if (parser.has_next_argument(position))
                value = parser.get_command_line_argument(position + 1);
        }
        else
        {
            value = parser.get_command_line_argument(position);
        }

        // If the next argument is an option, stops.
        if (Argument::is_option(value))
            break;

        // Else, stores the value.
        if (!value.empty())
            arg.append_value(value);

        // Increments the number of values.
        number_of_values++;

        // Increments the index.
        if (!arg.is_set())
            position++;
    }

    // If required "at most N" arguments, and the number of values is less than N, then the argument is set.
    if (arg.number_of_arguments().type == NumberOfArgumentsType::AT_MOST && number_of_values <= arg.number_of_arguments().value)
        arg.is_set(true);

    // If required "at least N" arguments, and the number of values is greater or equal to N, then the argument is set.
    if (arg.number_of_arguments().type == NumberOfArgumentsType::AT_LEAST && number_of_values >= arg.number_of_arguments().value)
        arg.is_set(true);

    // If argument is not set, the wrong number of argument has been provided.
    if (!arg.is_set())
    {
        std::string prefix(string_format("Option '%s' requires", arg.get_name().c_str()));
        std::string suffix(string_format("%lu arguments (found %lu)", arg.number_of_arguments().value, number_of_values));

        if (arg.number_of_arguments().type == NumberOfArgumentsType::AT_LEAST)
            prefix += " at least";

        parser.die(prefix + " " + suffix);
    }

    // If required "exactly N" or "at most N" arguments, and provided more than N, then the next arguments
    // will be treated as positional arguments.
}



Argument StoreTrueArgument(const std::string & name_short, const std::string & name_long, const std::string & description)
{
    return Argument()
        .name_short(name_short)
        .name_long(name_long)
        .description(description)
        .default_value("false")
        .action_type(ActionType::STORE_TRUE)
        .argument_type(ArgumentType::BOOLEAN)
        .required(false);
};

Argument HelpArgument()
{
    return StoreTrueArgument("-h", "--help", "show this help message and exit");
}

Argument VersionArgument()
{
    return StoreTrueArgument("", "--version", "show program's version number and exit");
}

} // namespace argparse
} // namespace biospring