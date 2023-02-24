/* Author: Abhijith
  ClassName: Mesh

Decription: Handles all the aspect of the mesh. Reads an STL file and converts the triangle mesh information
 to vertices and vertex index format. Vertices are stored in an Nx3 matrix of double. Triangles/Faces are stored 
 in an Ntx3 matrix of integer. Further edge information and edge to triangle map is generated. Another details 
 about triangles shared by a vertex is also generated. The winding order is corrected in the Mesh class, so that 
 all the normal are pointing outward.

*/

#pragma once
#ifndef MESHMODEL_H
#define MESHMODEL_H

#include <stdio.h>
#include <stdlib.h>
#include<string>
#include<vector>
#include <cstring>
#include<iostream>
#include<map>
#include<numeric>
#include <queue>


#include "eigenhelper.h"

class Mesh {
private:
	void generateEdgeInformation();

	void correct_Winding(int ix);


	Eigen::MatrixXi edge_indices;

public:

	int ntri, nvertices, nedges;

	Eigen::MatrixXd vertices;
	Eigen::MatrixXi triangles;

	Eigen::MatrixXi tri2edge;
	Eigen::MatrixXi edge2tri;

	Eigen::MatrixXd area;

	Eigen::MatrixXd Atilde;

	std::vector<std::vector<int>> vert2tri;

	Eigen::MatrixXd Centroids;
	Eigen::MatrixXd Normals;
	Eigen::MatrixXd dvertices;

	Mesh(void);
	int readSTL(std::string name);
	int writeToSTL(std::string fname,std::vector<int>);

	void generateNormal();
	void getVetexNormals();


	void offsetMesh(double);

	
	void Laplace_smoothing();
	void getIsoSolution();


};

#endif
