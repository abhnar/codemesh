/*
	This is where the processing part of the mesh is written.
	A bi-Laplacian equation is is solved on the mesh by applying certain forcing conditions
	to ensure the mesh consistency. The field obtained as the solution of the differential equation 
	is used to move the mesh in the normal direction. Force is applied such that the vertices with
	more number of triangles are applied proportionaly less force so that they move less, while the other 
	parts of the mesh move more.

	getK() --> to generate the stiffness matrix
	getM() --> generates the mass matrix
	getD() --> This is an additional constraint matrix, which is used to generate bi-Laplacian -KDK
				Interestingly -DK is similiar the Laplace-Beltrami operator
	getf() --> Calculates the forcing term, or right hand side off the equation. 

	getIsoSolution() --> Assembles and solves the equation and modifies dvertice, which is later used for
						 modifying the vertices.

*/

#include "meshmodel.h"

/*
Generates the standard stiffness matrix.
Elemental matrixes are created of the form  Integral{(Del W),{Del W)} over each triangle
The elemental matrices are assembled to get the global stiffness matrix

*/
Eigen::SparseMatrix<double> getK(Mesh mesh){
	Eigen::SparseMatrix<double> K(mesh.nvertices, mesh.nvertices);
	double B[3][3];
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(mesh.nvertices * 10);
	for (int i = 0; i < mesh.ntri; i++) {

		double x1 = mesh.vertices(mesh.triangles(i, 0), 0);
		double x2 = mesh.vertices(mesh.triangles(i, 1), 0);
		double x3 = mesh.vertices(mesh.triangles(i, 2), 0);

		double y1 = mesh.vertices(mesh.triangles(i, 0), 1);
		double y2 = mesh.vertices(mesh.triangles(i, 1), 1);
		double y3 = mesh.vertices(mesh.triangles(i, 2), 1);

		double z1 = mesh.vertices(mesh.triangles(i, 0), 2);
		double z2 = mesh.vertices(mesh.triangles(i, 1), 2);
		double z3 = mesh.vertices(mesh.triangles(i, 2), 2);

		B[0][0] = (x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3) + (z2 - z3) * (z2 - z3); //b
		B[0][1] = (x2 - x3) * (x3 - x1) + (y2 - y3) * (y3 - y1) + (z2 - z3) * (z3 - z1);
		B[0][2] = (x2 - x3) * (x1 - x2) + (y2 - y3) * (y1 - y2) + (z2 - z3) * (z1 - z2);

		B[1][0] = (x3 - x1) * (x2 - x3) + (y3 - y1) * (y2 - y3) + (z3 - z1) * (z2 - z3);
		B[1][1] = (x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1) + (z3 - z1) * (z3 - z1); //c
		B[1][2] = (x3 - x1) * (x1 - x2) + (y3 - y1) * (y1 - y2) + (z3 - z1) * (z1 - z2);

		B[2][0] = (x1 - x2) * (x2 - x3) + (y1 - y2) * (y2 - y3) + (z1 - z2) * (z2 - z3);
		B[2][1] = (x1 - x2) * (x3 - x1) + (y1 - y2) * (y3 - y1) + (z1 - z2) * (z3 - z1);
		B[2][2] = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2); //a
		
		
		

		

		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++) {

				tripletList.push_back(T(mesh.triangles(i, j), mesh.triangles(i, k), B[j][k] / (4 * mesh.area(i))));
			}
	}
	K.setFromTriplets(tripletList.begin(), tripletList.end());
	return K;
}

/*
Mass matrix is a diagonal matrix.
Elemental matrixes are created of the form  as the one by third othe  triangle area
The elemental matrices are assembled to get the global stiffness matrix

*/
Eigen::SparseMatrix<double> getM(Mesh mesh){
	Eigen::SparseMatrix<double> M(mesh.nvertices, mesh.nvertices);
	
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;

	tripletList.reserve(mesh.nvertices * 1);
	for (int i = 0; i < mesh.ntri; i++) {
		for (int j = 0; j < 3; j++)
			

				tripletList.push_back(T(mesh.triangles(i, j), mesh.triangles(i, j), mesh.area(i)/3));
			
	
	}
	M.setFromTriplets(tripletList.begin(), tripletList.end());
	return M;
	

}

/*
This matrix is also a diagonal matrix.
These matrices are formed by iterating over each vertex, and finding the total area (Atilde) of triangles
connected to the vertex.


*/
Eigen::SparseMatrix<double>  calcD(Mesh mesh){
    Eigen::SparseMatrix<double> D(mesh.nvertices, mesh.nvertices);
    typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
    for (int i = 0; i < mesh.nvertices; i++) {
        tripletList.push_back(T(i,i, mesh.Atilde(i)/3));
    }
    D.setFromTriplets(tripletList.begin(), tripletList.end());
	
    return D;
    
}
/* Forcing term is created, by keeping in mind that the vertex with more than three triangles connected to it
may not have a unique solution for the offset position. Thus they should have less force.
By dividing by square of Atilde and square of the length of the vertex normal vector ensures that
the force is proportinally reduced when more number of traingles are connected to a vertex
*/
Eigen::VectorXd getf(Mesh mesh){
    Eigen::VectorXd f = Eigen::VectorXd::Zero(mesh.nvertices,1);
    
    for (int i= 0; i<mesh.nvertices;i++)
        f(i) = f(i) + mesh.area(i)/(mesh.Atilde(i)*mesh.Atilde(i))/(mesh.dvertices.row(i).norm()*mesh.dvertices.row(i).norm());
    return f;
}

/*
Assembles the equation (11) in https://smartech.gatech.edu/bitstream/handle/1853/61/04-12.pdf
alpha Utt − gamma1 Del^2 U' + Gamma2 Del^4 U = f
Value of alpha is obtained from trian and error, the value is small but plays a big role in
keeping the system stable. alpha provides momentum to the system so that the solution does not
overshoot. 

*/
void Mesh::getIsoSolution() {
	
	Eigen::SparseMatrix<double> K = getK(*this);
	

	Eigen::SparseMatrix<double> M = getM(*this);


	Eigen::SparseMatrix<double> D = getM(*this);

	
	Eigen::VectorXd f = getf(*this);

	Eigen::VectorXd un = Eigen::VectorXd::Zero(nvertices);
	Eigen::VectorXd utn = Eigen::VectorXd::Zero(nvertices);

	

	double alpha = 1/area.sum();
    double beta1 = 0;
    double beta2 = 0;
    double gamma1 = 1;
    double gamma2 = 1;
    double dt = 1;
    std::cout<<"alpha: "<<alpha<<std::endl;

    Eigen::SparseMatrix<double> Ktilde = gamma1*K + gamma2*K*D*K;

	Eigen::VectorXd  b = -dt*Ktilde*un+alpha*M*utn+dt*f;
    Eigen::SparseMatrix<double> A = (alpha + beta1*dt)*M + (beta2+dt)*dt*Ktilde;

	 

	std::cout << "Generated FEM Matrix of size: " << A.rows() <<"x"<<A.cols() << std::endl;
	Eigen::SparseLU< Eigen::SparseMatrix<double> >  solver;

	std::cout << "Solving Ax=b" << std::endl;
	solver.analyzePattern(A);
	solver.compute(A);
	
	un = solver.solve(b);


	// Normalizing the solution
	Eigen::MatrixXf::Index   maxIndex;

	double maxVal = un.maxCoeff();

	un = un/maxVal;

	// Modifying dvertices, to use later in moving the mesh
	for(int i=0;i<nvertices;i++)
		dvertices.row(i) = dvertices.row(i)/dvertices.row(i).norm()*un(i);
	
	std::cout << "Solved" << std::endl;
	

}