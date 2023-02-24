# Code: Offset Mesh

## Compile Instructions

###  Windows
#### Visual Studio Build Tool
	cl /EHsc /I include src\*.cpp /std:c++17 /Fe:cod.exe /O2 /openmp && del *.obj


#### MSYS2 MINGW64
	g++ -I include src/*.cpp -fopenmp -O2 -o cod


### Linux (Tested under windows WSL)
	g++ -I include src/*.cpp -fopenmp -O2 -o cod


## Usage
	cod meshfile.stl d [clean]

	* d -> distance to offset
	* [clean] - cleans invalid triangle after offset. Caution!!! Not optimised yet. Does not scale for large mesh.

## Output
	output.stl - Mesh that is offset by distance d

## Yet to Implement
	Constrained Delauny Triangulation, to fill the hoe created by deleted triangles.

### Example
	cod mesh\curve.stl 0.5 
	cod mesh\curve.stl 0.5 clean
