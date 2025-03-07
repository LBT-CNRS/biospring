# Changes in the "BioSpring" branch of UnityMolX:
# - LoadNCFromPath is a new function
# - The UMol working directory is set to the directory of this history file (PythonConsole2.cs > ExecuteScript)
LoadNCFromPath(filePath="model.nc")
clearSelections()

setStructurePositionRotation("model", Vector3(0.0000, 0.0000, 0.0000), Vector3(0.0000, 0.0000, 0.0000))
#Save parent position
setMolParentTransform( Vector3(0.1525, -0.0934, 0.1519), Vector3(0.0145, 0.0145, 0.0145), Vector3(282.9991, 180.0000, 206.8996), Vector3(0.1350, -0.1200, 0.1528) )

ShowMeshMembrane()
UpdateMeshMembraneSize(100)