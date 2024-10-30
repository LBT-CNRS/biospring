#include "IO/PDBTrajectoryWriter.h"

void PDBTrajectoryWriter::write()
{
    safeOpen();
	writeModel(_frameid++);
	_ostream << "END";
    _ostream.close();
}
