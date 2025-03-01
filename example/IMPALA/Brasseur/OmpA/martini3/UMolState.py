LoadNCFromPath(filePath="model.nc")
clearSelections()

setStructurePositionRotation("model", Vector3(0.0000, 0.0000, 0.0000), Vector3(0.0000, 0.0000, 0.0000))
#Save parent position
setMolParentTransform( Vector3(0.4355, -0.2465, 0.2936), Vector3(0.0254, 0.0254, 0.0254), Vector3(284.2994, 180.0004, 165.3995), Vector3(0.4400, -0.2850, 0.2909) )


ShowMeshMembrane()
UpdateMeshMembraneSize(100)