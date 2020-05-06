// Gmsh project created on Wed May 06 13:34:33 2020
SetFactory("Built-in");
Mesh.CharacteristicLengthMin = .001;
Mesh.CharacteristicLengthMax = .001;
Geometry.LineNumbers = 0;
Geometry.SurfaceNumbers = 0;
Merge "Beam.STL";
//+
Surface Loop(1) = {1};
//+
Volume(1) = {1};
//+
