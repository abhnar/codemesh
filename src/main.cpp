/* Author: Abhijith
  ClassName: Mesh

Decription: Reads the command from command line. Reads mesh file name and offset distance. 
Calls the stlfile read function.
	During read, Vertex, Face information is generated, edge information is generated and winding of the triangles are corrected.
	stored in the Mesh Model object mesh

*/
#include "meshmodel.h"
#include "intersection.h"


int main(int argc, char** argv) {
	Mesh mesh;
	bool clean = false;
	double dn = 0.5;

/* This is the command selection block */
	std::vector<std::string> all_args;

	if (argc == 1) {
		std::cout << "Please provide stl filename and offset distance\n\n\tFormat: cod meshfile.stl 0.5 [clean]\n\n" << std::endl;;
		return 0;
	}
	if (argc == 2) {
		std::cout << "Please provide offset distancen\n\tFormat: cod meshfile.stl 0.5 [clean]\n\n" << std::endl;
		return 0;
	}

	if (argc == 4) {
		if (std::strcmp(argv[3], "clean")==0){
		clean = true;
		std::cout << "Warning: clean function is not optimized and may not scale for large models\n";}
		else{
			std::cout << "clean is off.\n\tFormat: cod meshfile.stl 0.5 [clean]\n\n" << std::endl;
		}
	}

	char* e;
	errno = 0;
	dn = std::strtod(argv[2], &e);

	if (*e != '\0' || errno != 0)   
	{
		std::cout << "Invalid number" << std::endl;
		return 0;
	}
	if (dn == 0) {
		std::cout << "Nothing to do." << std::endl;
		return 0;
	}

/* End ofcommand selection block */

	if(!mesh.readSTL(argv[1])){
		std::cout<<"File not found!!";
		return 0;
	}
		
	


	
	mesh.generateNormal(); 					// Generates the outward normals to triangle surface
	Eigen::MatrixXd nrm = mesh.Normals;		// Stores the normal, for further optimization of clean functionality

	
	mesh.getVetexNormals();					// Generates the outward normals to vertex. 
											//The normal vectors of triangles connected to the vertex are averaged.
	

	
	mesh.getIsoSolution();					// Generates the solution for bi-Laplace equation
	
	mesh.offsetMesh(dn);					// Move the mesh vertices by the iso-vectors obtained in above step

	std::vector<int> interSectingTriangles; 
	
	/* The following part of the program is to delte selfi=-intersecting and invalid triangles.
		This is not optimized and runs in O(n^2) time. Do not work for large mesh.*/
	if(clean){
		mesh.generateNormal();		

		
		std::cout<<"Checking intersection"<<std::endl;
		
		interSectingTriangles = selfIntersection(mesh,nrm);		// Detecting self intersection triangles

		
		std::cout<<"Checking badTriangles"<<std::endl;
		std::vector<int> badTriangles = invalidTriangles(mesh, nrm); // Detecting invalid triangles, penetrated into the model surface
		
		
		
		interSectingTriangles.insert( interSectingTriangles.end(), badTriangles.begin(), badTriangles.end() ); //Combinig both to delete later
		

		
	}
	
	mesh.writeToSTL("output.stl",interSectingTriangles);  //Write to STL Mesh, less intersecting and invalid triangles

	
	return 0;
}
