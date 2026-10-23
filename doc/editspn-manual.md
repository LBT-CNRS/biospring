# editspn manual

`editspn` edits and creates spring networks. It reads a topology, can add
springs within a cutoff distance, update or add springs listed in a file,
remove springs listed in a file, and writes the resulting topology to one or
more files.


## Usage

```
editspn [-h] [--version] [-o OUTPUT_FILE ...] [-u SPRING_FILE] [-r SPRING_FILE]
        [-c CUTOFF] [-d STIFFNESS] -s INPUT_FILE
```

| Option                   | Argument         | Default     | Description                                                    |
|--------------------------|------------------|-------------|----------------------------------------------------------------|
| `-s`, `--topology`       | `INPUT_FILE`     | (required)  | Input topology (`.pdb`, `.pqr` or `.nc`).                      |
| `-o`, `--output`         | `OUTPUT_FILE...` | `system.nc` | One or more output files (`.pdb`, `.pqr`, `.nc` or `.cdl`).    |
| `-c`, `--cutoff`         | `CUTOFF`         | `-1.0`      | Creates springs between all particles closer than `CUTOFF` (Å). |
| `-d`, `--stiffness`      | `STIFFNESS`      | `1.0`       | Stiffness of the springs created with `-c`.                    |
| `-u`, `--update-springs` | `SPRING_FILE`    |             | Updates or adds the springs listed in the file.                |
| `-r`, `--remove-springs` | `SPRING_FILE`    |             | Removes the springs listed in the file.                        |
| `--version`              |                  |             | Prints the program version and exits.                          |
| `-h`, `--help`           |                  |             | Prints the help message and exits.                             |


## Order of operations

The operations are applied in the following order, whatever the order of the
options on the command line:

1. read the input topology (springs it already contains are kept);
2. create springs from the cutoff distance (`-c`, `-d`);
3. update or add springs from the update file (`-u`);
4. remove springs from the removal file (`-r`);
5. write the output file(s) (`-o`).

Removing springs last means springs created with `-c` or `-u` can be removed in
the same run.


## Input and output formats

The format is guessed from the file extension.

| Extension | Input | Output | Notes                                                    |
|-----------|-------|--------|----------------------------------------------------------|
| `.pdb`    | yes   | yes    | Springs are written as `CONECT` records.                 |
| `.pqr`    | yes   | yes    | Springs are not written.                                 |
| `.nc`     | yes   | yes    | NetCDF binary; stores stiffness and equilibrium distance. |
| `.cdl`    |       | yes    | NetCDF ASCII.                                            |

An output file with no extension or an unknown extension is written in NetCDF
binary format, with a warning. Use `.nc` to keep spring stiffnesses and
equilibrium distances.


## Creating springs from a cutoff (`-c`, `-d`)

With `-c CUTOFF`, a spring is created between every pair of particles closer
than `CUTOFF`. Each new spring has:

- an equilibrium distance equal to the current distance between its particles;
- a stiffness equal to `STIFFNESS` if `-d` is given, `1.0` otherwise.

Pairs of particles that are both static are skipped. Springs that already exist
between two particles (for instance read from the input file) are left
unchanged: their stiffness is not modified by `-d`.

`-d` has no effect without `-c`.


## Particle identifiers

Spring files refer to particles with identifiers of the form:

```
chain::resname::resid::name
```

For example `A::ALA::12::CA`. Fields may be empty: a particle without a chain
name is written `::ALA::12::CA`. Identifiers must match the particles of the
input topology exactly (atom names without surrounding spaces); `editspn` stops
with an error if a particle is not found.

If several particles share the same identifier, only one of them can be
referred to.


## Spring update file (`-u`)

One spring per line:

```
<particle1> <particle2> <stiffness> <equilibrium>
```

Empty lines and lines starting with `#` are ignored. All four columns are
required.

For each line:

- if a spring already exists between the two particles, its stiffness is set
  to `<stiffness>`; its equilibrium distance is set to `<equilibrium>` if it is
  positive or zero, and left unchanged otherwise;
- if no spring exists, one is created with the given stiffness. If
  `<equilibrium>` is `-1`, the equilibrium distance is the current distance
  between the two particles.

Example:

```
# chain::resname::resid::name  chain::resname::resid::name  stiffness  equilibrium
A::LYS::10::CA                 A::GLU::54::CA               5.0        -1
A::LYS::10::CA                 A::LYS::11::CA               10.0       3.8
```


## Spring removal file (`-r`)

One spring per line:

```
<particle1> <particle2>
```

Empty lines and lines starting with `#` are ignored. Extra columns are ignored,
so a spring update file can be used as a removal file. The order of the two
particles does not matter.

`editspn` stops with an error if a line has fewer than two columns or if a
particle is not found in the topology. If no spring exists between two
particles, a warning is printed and the line is skipped. The number of removed
springs is reported in the log.

Example:

```
# springs to remove
A::LYS::10::CA  A::GLU::54::CA
A::LYS::11::CA  A::LYS::10::CA
```


## Examples

Create an elastic network with an 8 Å cutoff and a stiffness of 2.5, written
in both NetCDF and PDB formats:

```
editspn -s model.pdb -c 8 -d 2.5 -o model.nc model.pdb
```

Change the stiffness of some springs of an existing network and add new ones:

```
editspn -s model.nc -u springs.txt -o model-updated.nc
```

Create a network from a cutoff, then remove a few springs from it:

```
editspn -s model.pdb -c 8 -r remove.txt -o model.nc
```
