/* Author: Abhijith
  ClassName: Mesh

Decription: Handles all the aspect of the mesh. 
 Reads an STL file and converts the triangle mesh information to vertices and faces
 	1. Vertices are stored in an Nx3 matrix of double. 
 	2. Triangles/Faces are stored in an Ntx3 matrix of integer3. 
 	3. Further edge information and edge to triangle map is generated. 
 	4. Details  about triangles shared by a vertex is also generated. 
 	5. The winding order is corrected in the Mesh class, so that all the normal are pointing outward.

*/
#include "meshmodel.h"

// Class constructor
Mesh::Mesh(void) {

	edge_indices = Eigen::MatrixXi::Zero(3, 2);
	edge_indices << 0, 1, 1, 2, 2, 0;
}


//Writes the STL File. Faster in C format
int Mesh::writeToSTL(std::string fname, std::vector<int> dellist) {
	FILE* fid;
	fid = fopen(fname.c_str(), "w");
	fprintf(fid, "solid 000\n");


	if(!fid){
		std::cout<<"File not found!!"<<std::endl;
		return 01;
	}

	for (int i = 0; i < ntri; i++) {
		int  is_in = (std::count(dellist.begin(), dellist.end(), i));
		if(is_in) continue;
		fprintf(fid, "facet normal 0 0 0\n");
		fprintf(fid, "outer loop\n");
		fprintf(fid, "vertex %2.15f %2.15f %2.15f\n", vertices(triangles(i, 0), 0), vertices(triangles(i, 0), 1), vertices(triangles(i, 0), 2));
		fprintf(fid, "vertex %2.15f %2.15f %2.15f\n", vertices(triangles(i, 1), 0), vertices(triangles(i, 1), 1), vertices(triangles(i, 1), 2));
		fprintf(fid, "vertex %2.15f %2.15f %2.15f\n", vertices(triangles(i, 2), 0), vertices(triangles(i, 2), 1), vertices(triangles(i, 2), 2));
		fprintf(fid, "endloop\n");
		fprintf(fid, "endfacet\n");
	}
	fclose(fid);
	std::cout << "successully Written STL: "<<fname << std::endl;
	return 1;
}

//Reads the STL File. Faster in C format
int Mesh::readSTL(std::string fname) {
	std::vector<float> Normals = std::vector<float>();;
	std::vector<double> Triangles_ = std::vector<double>();
	std::vector<int> Triangles = std::vector<int>();
	std::vector<double> Vertex = std::vector<double>();

	

	FILE* file_pointer;
	char buffer[200];
	char str1[30];
	char str2[30];
	int roll;
	int std;
	float ax, ay, az;
	int count = 0;



	file_pointer = fopen(fname.c_str(), "r");
	if (file_pointer) {

		fscanf(file_pointer, "%[^\n]\n", buffer);


		while (1) {
			fscanf(file_pointer, "%s %s %f %f %f\n", str1, str2, &ax, &ay, &az);

			if (std::strcmp(str1, "endsolid") == 0) {
				
				break;
			}
			Normals.push_back(ax);
			Normals.push_back(ay);
			Normals.push_back(az);

			fgets(buffer, 200, file_pointer);

			fscanf(file_pointer, "%s  %f %f %f\n", str1, &ax, &ay, &az);
			Triangles_.push_back(ax);
			Triangles_.push_back(ay);
			Triangles_.push_back(az);


			fscanf(file_pointer, "%s  %f %f %f\n", str1, &ax, &ay, &az);
			Triangles_.push_back(ax);
			Triangles_.push_back(ay);
			Triangles_.push_back(az);

			fscanf(file_pointer, "%s  %f %f %f\n", str1, &ax, &ay, &az);
			Triangles_.push_back(ax);
			Triangles_.push_back(ay);
			Triangles_.push_back(az);

			fgets(buffer, 200, file_pointer);

			fgets(buffer, 200, file_pointer);

		}



		fclose(file_pointer);
		std::cout << "STL File loaded successully\n";

	}else{
		return 0;
	}

	

	/* The following part of the code converts the vertex ist to V and F format 
		Vertices V stored in Mesh::vertices.
		Faces F stored in Mesh::triangels 
		nvertices is the number of vertices
		ntri is the number of triangles.
		
		vertices[i,0] --> x of vertex i
		vertices[i,1] --> y of vertex i
		vertices[i,2] --> z of vertex i

		triangles[i,0] --> vertex 0 of triangle i
		triangles[i,1] --> vertex 1 of triangle i
		triangles[i,2] --> vertex 2 of triangle i
		*/
	int vertexCount = 0;
	float p[3];

	//Generetaes the triangle index list with nonunique vertices.
	for (int i = 0; i < Triangles_.size(); i = i + 3) {
		p[0] = Triangles_[i];
		p[1] = Triangles_[i + 1];
		p[2] = Triangles_[i + 2];

		Vertex.push_back(p[0]);						// Vertex might repeat in this Vertex list
		Vertex.push_back(p[1]);
		Vertex.push_back(p[2]);
		Triangles.push_back(vertexCount); // Each vertex is stored in Triangle vector
		vertexCount++;




	}

	Triangles_.clear();



	//Reshape vertex vector to form an Nx3 matrix of vertices (Still non-unique)
	Eigen::VectorXd vert_ = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Vertex.data(), Vertex.size());
	vertices = (Eigen::Map<Eigen::MatrixXd>(vert_.data(), 3, Vertex.size() / 3).cast<double>()).transpose();

	//Reshape triangle vector  to form an Nx3 matrix of triangles 
	Eigen::VectorXi tri_ = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(Triangles.data(), Triangles.size());
	triangles = (Eigen::Map<Eigen::MatrixXi>(tri_.data(), 3, Triangles.size() / 3).cast<int>()).transpose();


	// Get unique vertices with the indices of unique entries to map the triangle list
	Eigen::VectorXi ix = getUnique(vertices);

	//Map triangle indices to uniqur vertices
	for (int i = 0; i < triangles.rows(); i++) {
		for (int j = 0; j < triangles.cols(); j++) {
			triangles(i, j) = ix(triangles(i, j));
		}


	}
	


	ntri = triangles.rows();
	nvertices = vertices.rows();

	std::cout << "# Trangles : " << ntri << std::endl;
	std::cout << "# Vertices: " << nvertices << std::endl;

	/* Generates the list of triangles connected to a vertex.
		vert2tri[i] is the list of triangles connected to vertex i
		*/

	vert2tri.resize(nvertices);

	for (int i = 0; i < ntri; i++) {
		for (int j = 0; j < 3; j++)
			vert2tri[triangles(i, j)].push_back(i);
	}


	//Post processing mesh information

	generateEdgeInformation();
	correct_Winding(1);
	
	return 1;
}


/*Generates the edge data. Finds all edges in the triangle, maps them to each 
//		triangle and the following data are generated.
//	1.	Eigen::MatrixXi tri2edge;
//			Represents 3 edges of each triangle (ntri x 3 Matrix)
//	2.	Eigen::MatrixXi edge2tri;
//			Represents triangles shared by an edge. There could be either 1 or 
//			two triangle, depends on the existence of open surfaces.

	** Works on O(nlogn)
*/

void Mesh::generateEdgeInformation() {

	//Create edge list - non -unique. Three lsits for three edes from each triangle.
	auto ed_list1 = triangles(Eigen::all, edge_indices.row(0));
	auto ed_list2 = triangles(Eigen::all, edge_indices.row(1));
	auto ed_list3 = triangles(Eigen::all, edge_indices.row(2));

	//Merge three lists to one
	Eigen::MatrixXi edgelist(ed_list1.rows() * 3, 2);
	edgelist << ed_list1, ed_list2, ed_list3;

	/* Creates  unique id for edges. Uses Cantor Pairig function. */

	unsigned long long  temp;
	std::vector<unsigned long long> ed_id(edgelist.rows());
	VectorXl ed_id_temp(edgelist.rows());
	int counttt = 0;
	for (int i = 0; i < edgelist.rows(); i++) {
		if (edgelist(i, 0) > edgelist(i, 1)) {
			temp = edgelist(i, 0);
			edgelist(i, 0) = edgelist(i, 1);
			edgelist(i, 1) = temp;

		}
		ed_id[i] = (edgelist(i, 0) + edgelist(i, 1) + 1LL + 2LL) * (edgelist(i, 0) + edgelist(i, 1) + 2LL) / 2LL + edgelist(i, 1) + 1LL;
		ed_id_temp[i] = (edgelist(i, 0) + 2LL + edgelist(i, 1) + 1LL) * (edgelist(i, 0) + edgelist(i, 1) + 2LL) / 2LL + edgelist(i, 1) + 1LL;
		counttt = counttt + 1;
	}



	/* Get the sprted indices of Id list. */

	std::vector<size_t> iv(ed_id.size());
	std::iota(iv.begin(), iv.end(), 0);
	std::sort(iv.begin(), iv.end(),
		[&ed_id](size_t i1, size_t i2) {return ed_id[i1] < ed_id[i2]; });



	/* Map each triangle to corresponding edge.
		tri2edge(i,0) -> 1st edge of the triangle i given by triangles(i,[1,2,3])
		tri2edge(i,1) -> 2nd edge of the triangle i given by triangles(i,[1,2,3])
		tri2edge(i,2) -> 3rd edge of the triangle i given by triangles(i,[1,2,3])
		edge2tri(i,0) -> First triangle connected to edge i
		edge2tri(i,1) -> Second triangle connected to edge i, -1 if the edge has only one triangle
	*/
	std::vector<int> uix;
	std::set<int> s = {};
	uix.push_back(iv[0]);
	for (int i = 1; i < ed_id.size(); i++)
		if (ed_id[iv[i]] != ed_id[iv[i - 1]]) {
			uix.push_back(iv[i]);
		}

		VectorXl edgeids_unique(uix.size());
		for (int i = 0; i < uix.size(); i++) {
			edgeids_unique(i) = ed_id[uix[i]];
		}
		nedges = edgeids_unique.rows();
		std::map <unsigned long long, int> edge_map;
		for (int i = 0; i < edgeids_unique.rows(); i++) {
			edge_map[edgeids_unique(i)] = i;
		}

		std::vector<int> ed_id_mapped(edgelist.rows());
		for (int i = 0; i < ed_id.size(); i++) {
			ed_id_mapped[i] = edge_map[ed_id[i]];
		}

		Eigen::VectorXi  edges = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(ed_id_mapped.data(), ed_id_mapped.size());

		tri2edge = (Eigen::Map<Eigen::MatrixXi>(edges.data(), edges.size() / 3, 3).cast<int>());

		VectorXl uix_temp = Eigen::Map<VectorXl, Eigen::Unaligned>(ed_id.data(), ed_id.size());

		edge2tri = -Eigen::MatrixXi::Ones(nedges, 2);

		Eigen::VectorXi edg2tri_n = Eigen::VectorXi::Zero(nedges);
		int i, j;
		try {

			for (i = 0; i < ntri; i++)
				for (j = 0; j < 3; j++) {
					edge2tri(tri2edge(i, j), edg2tri_n(tri2edge(i, j))) = i;
					edg2tri_n(tri2edge(i, j)) = edg2tri_n(tri2edge(i, j)) + 1;


				}


			}
			catch (...) {


			}

			std::cout << "# Edges " << nedges << std::endl;

		}

/*Calculate the normal direction of the vertex, based on the direction normal 
	ofadjacent triangles. Averages the normal and finds out the vector through 
	which the vertex has to move on the offset.
 	Along with this, the sum of triangles connected to one vetex is also calculated.
		dvertices(i,[1,2,3]) --> normally outward vector of vertex i
		Atilde(i) --> Sum of areas of triangles connected to ith vertex.
*/

void Mesh::getVetexNormals() {

		dvertices = Eigen::MatrixXd::Zero(nvertices, 3);
		Eigen::MatrixXi nNormals = Eigen::MatrixXi::Zero(nvertices, 1);
		std::vector<bool> moved(nvertices);
		Atilde = Eigen::MatrixXd::Zero(triangles.rows(), 1);
		for (int i = 0; i < nvertices; i++) {

			for (int j : vert2tri[i]) {

				dvertices.row(i) = dvertices.row(i) + Normals.row(j);
				Atilde(i) = Atilde(i) + area(j);

			}


			dvertices.row(i) = dvertices.row(i) / vert2tri[i].size();



		}

		

}


/*  Corrects the winding of triangles. ENsures that all the triangles are in clockwise direction.
		Traverses through adjacent triangles based on the edge information. When one 
   edge is visited, two triangles are visited. The other edges of the two triangles
   are sted in a queue, and popped one by one so that the traversal is through 
   adjacent triangle.
   Workd in O(n)
*/

		void Mesh::correct_Winding(int ix) {

			// Start with one random edge. Technically should start with a triangle with correct order
			// If we start with triangle on the open surface, the algorithm may not converge, hence starting at 100000
			// for given mesh
			ix = 100000 - 1;
			if (ix > nedges) {
				ix = 10;
			}

			//Queue to store next edges to visit
			std::queue<int> edge_queue;

			Eigen::VectorXi visited = Eigen::VectorXi::Zero(nedges);
			Eigen::VectorXi visited1 = Eigen::VectorXi::Zero(ntri);
			Eigen::VectorXi correct = Eigen::VectorXi::Zero(ntri);
			int trcnt = 1;
			int initial_edge = ix;

			edge_queue.push(initial_edge);			//Store the edge to queue

			int t1 = edge2tri(initial_edge, 0), t2; 

			correct(t1) = 1;			//mark the first triangle of the edge correct order and visited
			visited1(t1) = 1;
			int limit = ntri;
			int ed, tri_ref, tri_edit;
			int ix1, ix2;
			Eigen::VectorXi ed1(2);
			Eigen::VectorXi ed2(2);
			bool flag = false;
			int loopcount = 0;

			//Traverse through every adjascent triangle, found by the edges
			while (!edge_queue.empty()) {


				loopcount++;
				ed = edge_queue.front();
				edge_queue.pop();
				while (visited(ed) == 1 && !edge_queue.empty()) {
					ed = edge_queue.front();
					edge_queue.pop();

				}

				visited(ed) = 1;
				t1 = edge2tri(ed, 0);  //traingle1 to the edge
				t2 = edge2tri(ed, 1);  //traingle2 to the edge

				//Already correct triangle is taken as reference for winding order
				if (correct(t1) == 1) {	
					tri_ref = t1;
					tri_edit = t2;
				}
				else if (correct(t2) == 1) {
					tri_ref = t2;
					tri_edit = t1;
				}
				

				std::vector<int> e;
				if (tri_edit != -1) {    //Ensure that there is a second triangle

					if (visited1(tri_edit) == 0) {	//If it is not visited

						auto q1 = tri2edge.row(tri_ref);

						auto q2 = tri2edge.row(tri_edit);
						ix1 = (q1(0) == ed ? 0 : (q1(1) == ed ? 1 : 2));
						ix2 = (q2(0) == ed ? 0 : (q2(1) == ed ? 1 : 2));


						ed1 = triangles(tri_ref, edge_indices.row(ix1));
						ed2 = triangles(tri_edit, edge_indices.row(ix2));


						// Reverese the vertices on the edge if the triangle winding order is reversed
						if (ed1(0) == ed2(0)) {

							triangles(tri_edit, edge_indices(ix2, 0)) = ed1(1);
							triangles(tri_edit, edge_indices(ix2, 1)) = ed1(0);
							std::vector<int> x{ 0, 1, 2 };

							x.erase(x.begin() + ix2);
							int temp = tri2edge(tri_edit, x[0]);
							tri2edge(tri_edit, x[0]) = tri2edge(tri_edit, x[1]);
							tri2edge(tri_edit, x[1]) = temp;


						}
						//Mark as visited and correct
						visited1(tri_edit) = 1;
						correct(tri_edit) = 1;
						trcnt = trcnt + 1;
					}

				}

				//Add edges of the current triangles to the queue
				if (tri_edit != -1) {
					auto e_set_1 = tri2edge.row(t1);
					auto e_set_2 = tri2edge.row(t2);
					for (int k = 0; k < 3; k++) {
						if (visited(e_set_1(k)) == 0) edge_queue.push(e_set_1(k));
						if (visited(e_set_2(k)) == 0) edge_queue.push(e_set_2(k));

					}
				}

			}




		}

/*Calculates the centroid of all triangle and the normal from the surface. 
//If the winding order is correct, the outward notmal of all triangles will
// be obtained. Area of the triangle is also calculated in the process
	Centroid(i,[0,1,2]) --> Centroid of the ith triangle
	Normal(i,[0,1,2]) --> Outward Normal from the surface of ith triangle
	area(i) --> Area of the ith triangle

*/
		void Mesh::generateNormal() {

			Eigen::MatrixXd   p1(triangles.rows(), 3);
			p1 = vertices(triangles.col(0), Eigen::all);

			Eigen::MatrixXd  p2(triangles.rows(), 3);
			p2 = vertices(triangles.col(1), Eigen::all);

			Eigen::MatrixXd  p3(triangles.rows(), 3);
			p3 = vertices(triangles.col(2), Eigen::all);

			Centroids = (p1 + p2 + p3) / 3.0;

			Eigen::MatrixX3d temp1(triangles.rows(), 3);
			temp1 = p2 - p1;
			Eigen::MatrixX3d temp2(triangles.rows(), 3);
			temp2 = p3 - p1;
			Normals = Eigen::MatrixXd::Zero(triangles.rows(), 3);
			area = Eigen::MatrixXd::Zero(triangles.rows(), 1);
			for (int i = 0; i < temp1.rows(); i++) {
				Normals.row(i) = temp1.row(i).cross(temp2.row(i));
				area(i) = Normals.row(i).norm()*0.5;
				Normals.row(i) = Normals.row(i)/(2*area(i));

			}

		}

/* Adds the vertex normal (scaled by dn) to the existing vertex to get the offset 
  vertices. This may create bad mesh. Non-smooth mesh can be smoothened by Laplace smoothing.
  */

		void Mesh::offsetMesh(double dn) {
			vertices = vertices + dn * dvertices;
	//write_csv("final_vertices.csv", vertices);
		}


/* Averages the existing vertex with, surrounding vertices.

  */

	void Mesh::Laplace_smoothing() {

		int n;
		for (int i = 0; i < nvertices; i++) {
			auto ix = vert2tri[i];
			auto v = triangles(ix, Eigen::all);
			std::set<int> s;
			for (int j = 0; j < v.rows(); j++) {
				for(int k=0;k<v.cols();k++)
					if (i != v(j,k))
						s.insert(v(j,k));
				}
				n = s.size();
				for (auto j : s) {
					vertices.row(i) = (vertices.row(j) + vertices.row(i));
				}

				vertices.row(i) = vertices.row(i) / (n+1);


			}
		}



