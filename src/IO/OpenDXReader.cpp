
#include "IO/OpenDXReader.h"
#include "logging.h"
#include "utils/string.hpp"

namespace logging = biospring::logging;

std::array<size_t, 3> OpenDXReader::readSize()
{
    std::string buffer;
    std::vector<std::string> tokens;

    // Skip comment lines.
    do
    {
        std::getline(_instream, buffer);
    } while (buffer.substr(0, 1) == "#");

    // At this point, buffer should be the line containing the dx grid size.
    if (buffer.size() < 35 or buffer.substr(0, 35) != "object 1 class gridpositions counts")
    {
        logging::die("OpenDXReader: cannot read grid dimensions (found %s)", buffer.c_str());
    }

    tokens = biospring::utils::string::split(buffer.substr(35, buffer.size() - 1));
    if (tokens.size() != 3)
    {
        logging::die("OpenDXReader: cannot read grid dimensions (found '%s')", buffer.c_str());
    }

    size_t sizei = static_cast<size_t>(std::stoi(tokens[0]));
    size_t sizej = static_cast<size_t>(std::stoi(tokens[1]));
    size_t sizek = static_cast<size_t>(std::stoi(tokens[2]));

    return {sizei, sizej, sizek};
}

std::array<double, 3> OpenDXReader::readOrigin()
{
    std::string buffer;
    std::vector<std::string> tokens;

    std::getline(_instream, buffer);
    if (buffer.substr(0, 6) != "origin")
    {
        logging::die("OpenDXReader: cannot read grid origin (found '%s')", buffer.c_str());
    }

    tokens = biospring::utils::string::split(buffer.substr(7, buffer.size() - 1));
    if (tokens.size() != 3)
    {
        logging::die("OpenDXReader: cannot read grid origin (found '%s')", buffer.c_str());
    }

    float offsetx = std::stof(tokens[0]);
    float offsety = std::stof(tokens[1]);
    float offsetz = std::stof(tokens[2]);

    return {offsetx, offsety, offsetz};
}

std::array<double, 3> OpenDXReader::readScalingFactors()
{
    std::string buffer;
    std::vector<std::string> tokens;
    float delta[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    for (size_t i = 0; i < 3; i++)
    {
        std::getline(_instream, buffer);
        if (buffer.substr(0, 5) != "delta")
        {
            logging::die("OpenDXReader: cannot read scaling factors (found '%s')", buffer.c_str());
        }

        tokens = biospring::utils::string::split(buffer.substr(6, buffer.size() - 1));
        if (tokens.size() != 3)
        {
            logging::die("OpenDXReader: cannot read scaling factors (found '%s')", buffer.c_str());
        }
        delta[i][0] = std::stof(tokens[0]);
        delta[i][1] = std::stof(tokens[1]);
        delta[i][2] = std::stof(tokens[2]);
    }

    return {delta[0][0], delta[1][1], delta[2][2]};
}

void OpenDXReader::readGrid()
{
    std::string buffer;
    std::vector<std::string> tokens;

    float unityconvert = 1.0;
    size_t gx = 0, gy = 0, gz = 0;
    size_t totalsize = _grid.size();

    // OpenDX text data simply wraps the flat array at (conventionally) 3
    // values per line; the total point count has no reason to be a
    // multiple of 3 (e.g. a 97x97x97 grid is not), so the last line
    // legitimately holds fewer than 3 values. Read by total value count
    // instead of assuming every line, including the last, has exactly 3.
    size_t read_count = 0;
    while (read_count < totalsize)
    {
        if (not std::getline(_instream, buffer))
        {
            logging::die("OpenDXReader: cannot access grid data");
        }

        tokens = biospring::utils::string::split(buffer);
        if (tokens.empty() || tokens.size() > 3)
        {
            logging::die("OpenDXReader: misformatted grid data (expected 1 to 3 tokens, found '%s'", buffer.c_str());
        }

        for (const std::string & token : tokens)
        {
            if (read_count >= totalsize)
            {
                logging::die("OpenDXReader: more grid data than the declared %zu values", totalsize);
            }

            biospring::grid::discrete_coordinates cell(gx, gy, gz);

            float scalar = std::stof(token);
            _grid.at(cell).scalar = scalar * unityconvert;
            read_count++;

            gz++;
            if (gz >= _grid.shape()[2])
            {
                gz = 0;
                gy++;
                if (gy >= _grid.shape()[1])
                {
                    gy = 0;
                    gx++;
                }
            }
        }
    }
}

void OpenDXReader::read()
{
    std::string buffer;              // buffer for lines
    std::vector<std::string> tokens; // buffer for tokens

    // Opens file for reading, dies if error occurs.
    safeOpen();

    auto shape = readSize();
    auto origin = readOrigin();
    auto spacing = readScalingFactors();

    std::array<double, 6> boundaries{origin[0],
                                     origin[1],
                                     origin[2],
                                     origin[0] + shape[0] * spacing[0],
                                     origin[1] + shape[1] * spacing[1],
                                     origin[2] + shape[2] * spacing[2]};
    _grid.reshape(boundaries, spacing);

    // Skips next line ("object 2 class gridconnections counts <xn> <yn> <zn>")
    std::getline(_instream, buffer);

    // Skips next line ("object 3 class array type "double" rank 0 items <n> data follows")
    std::getline(_instream, buffer);

    readGrid();
    _grid.compute_gradient();
}
