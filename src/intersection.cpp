#include "meshmodel.h"
#include "intersection.h"

typedef double real;

    #define ZERO_TEST(x)  (x == 0)
    //#define ZERO_TEST(x)  ((x) > -0.001 && (x) < .001)

    #include "stdio.h"

    /* function prototype for triangle-triangle intersection test */


    int tri_tri_overlap_test_3d(Eigen::MatrixXd p1, Eigen::MatrixXd p2);


    int coplanar_tri_tri3d(real  p1[3], real  q1[3], real  r1[3],
                           real  p2[3], real  q2[3], real  r2[3],
                           real  N1[3], real  N2[3]);


    int tri_tri_overlap_test_2d(real p1[2], real q1[2], real r1[2], 
                                real p2[2], real q2[2], real r2[2]);


    int tri_tri_intersection_test_3d(real p1[3], real q1[3], real r1[3], 
                                     real p2[3], real q2[3], real r2[3],
                                     int * coplanar, 
                                     real source[3],real target[3]);
     

// Perturbs the triangle, makes it smaller, compresses towards the centroid
//  to ensure the adjascent triangles are not found as self-intersecting
Eigen::MatrixXd perturb(Eigen::MatrixXd T){

       Eigen::MatrixXd C = Eigen::MatrixXd::Zero(1,T.cols());
       C.row(0) = (T.row(0) + T.row(1)+T.row(2))/3.0;
       Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(T.rows(),T.cols());
 
       Q.row(0) = T.row(0) - (T.row(0)-C)*0.1;
       Q.row(1) = T.row(1) - (T.row(1)-C)*0.1;
       Q.row(2) = T.row(2) - (T.row(2)-C)*0.1;

       return Q;
}
/*
This code tests self intersection of pair of triangles.
Refer: https://web.stanford.edu/class/cs277/resources/papers/Moller1997b.pdf
When an intersection is found, the triangles are added to the list of bad triangles.
Runs in less than O(n^2), but close to O(n^2).
Tried to make it faster by adding additional constraints, but misses many triangles.
Does not scale with mesh size. There is an algorithm to make thos faster, but
additional datatructures are needed.

*/
std::vector<int> selfIntersection(Mesh mesh, Eigen::MatrixXd nrm){
    int r;
    std::vector<int> triangles;
    //std::iota(triangles.begin(), triangles.end(), 0);
    for(int i =0;i<nrm.rows();i++){
      //  if(mesh.Normals.row(i).dot(nrm.row(i))<0.8){
            triangles.push_back(i);
      //  }
    }
    std::cout<<"Triangles under check"<<triangles.size()<<std::endl;
    std::set<int> bad_triangles;

    std::vector<int*> intersect_list;
    while(!triangles.empty()){
        int i = triangles.back();
        triangles.pop_back();
        
        std::vector<int> todelete;
        r = 0;
        bool adj = false;
        for(int j:triangles){
            if(i==j) continue;
            if(mesh.Normals.row(i).dot(mesh.Normals.row(j))>0.8) continue;
            
           
            r = tri_tri_overlap_test_3d( mesh.vertices(mesh.triangles(i,Eigen::all),Eigen::all),  mesh.vertices(mesh.triangles(j,Eigen::all),Eigen::all));
        
            if(r!=0)
            {   bad_triangles.insert(j);
                bad_triangles.insert(i);
                todelete.push_back(j);
                std::cout<<i<<","<<j<<std::endl;
            }
        }
        
        for(int j:todelete){
            std::remove(triangles.begin(), triangles.end(),j);
        }
  

    }

    
    std::vector<int> bt(bad_triangles.begin(), bad_triangles.end());

    return bt;
}

/*
    This code checks if a normal from a triangle's centroid touches any other triangle.
    If it touches and aligns with the normal of touching triangle, the triangle is inside the geometry and invalid.
    This does not get all the invalid triangles, but most of them.
    Implementation of Möller–Trumbore intersection algorithm. Vectorized implementation.
    Refer: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm


*/
std::vector<int> invalidTriangles(Mesh mesh, Eigen::MatrixXd nrm){

    double EPSILON = 0.0000001; 

    std::vector<int> invTri;
    std::vector<int> triangles;
    //std::iota(triangles.begin(), triangles.end(), 0);
    for(int i =0;i<nrm.rows();i++){
       // if(mesh.Normals.row(i).dot(nrm.row(i))<0.8){
            triangles.push_back(i);
       // }
    }

    Eigen::MatrixXd v1, v2, v0;
        v0 = mesh.vertices(mesh.triangles.col(0),Eigen::all);
        v1 = mesh.vertices(mesh.triangles.col(1),Eigen::all);
        v2 = mesh.vertices(mesh.triangles.col(2),Eigen::all);
        
        Eigen::MatrixXd edge1 = v1 - v0;
        Eigen::MatrixXd edge2 = v2 - v0;

    std::cout<<triangles.size()<<std::endl;
     #pragma omp parallel for num_threads(4)
     for(int it=0;it<triangles.size();it++){
        int j = triangles[it];
         Eigen::VectorXd rayOrigin = mesh.Centroids.row(j);
        Eigen::VectorXd rayVector = mesh.Normals.row(j);

        
        

        Eigen::MatrixXd h;
       

        Eigen::MatrixXd s;
        Eigen::MatrixXd q;
       

  
        Eigen::VectorXi res = Eigen::VectorXi::Zero(v1.rows(),1);
        
        Eigen::VectorXd a = Eigen::VectorXd::Zero(v1.rows(),1);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(v1.rows(),1);
        Eigen::VectorXd u = Eigen::VectorXd::Zero(v1.rows(),1);
        Eigen::VectorXd v = Eigen::VectorXd::Zero(v1.rows(),1);
        Eigen::VectorXd t = Eigen::VectorXd::Zero(v1.rows(),1);
       
        h = cross(rayVector,edge2);

        a = dot(h,edge1);
        
        
        
        
        res = (a.cwiseAbs().array()<EPSILON).select(1,res);
           
        
        f = a.cwiseInverse();
        s = diff(rayOrigin, v0);
        u = dot(s,h).cwiseProduct(f);
        

         res = (u.array()<0.0).select(1,res);
          res = (u.array()>1.0).select(1,res);
       
        q = crossm(s,edge1);
      
        v = dot(rayVector,q).cwiseProduct(f);

        
        res = (v.array()<0.0).select(1,res);
          res = ((u+v).array()>1.0).select(1,res);

        

            t = dot(q,edge2).cwiseProduct(f);

            res = (t.array()<EPSILON).select(1,res);
            if(res.sum()<mesh.ntri) 
                invTri.push_back(j);
            
        
        
       
    }
    
    return invTri;
}


/* The following part of the code is the implementation of Moller Algorithm for Triangle-Triangle intersection test.
    Got this from online and converted to compatible for this code, using Eigen matrix package
*/

    /* coplanar returns whether the triangles are coplanar  
     *  source and target are the endpoints of the segment of 
     *  intersection if it exists) 
     */


    /* some 3D macros */

    #define CROSS(dest,v1,v2)                       \
    dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
    dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
    dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

    #define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

    #define SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; \
    dest[1]=v1[1]-v2[1]; \
    dest[2]=v1[2]-v2[2]; 

    #define SCALAR(dest,alpha,v) dest[0] = alpha * v[0]; \
    dest[1] = alpha * v[1]; \
    dest[2] = alpha * v[2];

    #define CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) {\
    SUB(v1,p2,q1)\
    SUB(v2,p1,q1)\
    CROSS(N1,v1,v2)\
    SUB(v1,q2,q1)\
    if (DOT(v1,N1) > 0.0f) return 0;\
    SUB(v1,p2,p1)\
    SUB(v2,r1,p1)\
    CROSS(N1,v1,v2)\
    SUB(v1,r2,p1) \
    if (DOT(v1,N1) > 0.0f) return 0;\
    else return 1; }



    /* Permutation in a canonical form of T2's vertices */

    #define TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
    if (dp2 > 0.0f) { \
    if (dq2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2) \
    else if (dr2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
    else CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) }\
    else if (dp2 < 0.0f) { \
    if (dq2 < 0.0f) CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0.0f) CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
    else CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
    } else { \
    if (dq2 < 0.0f) { \
    if (dr2 >= 0.0f)  CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
    else CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2)\
    } \
    else if (dq2 > 0.0f) { \
    if (dr2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
    else  CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
    } \
    else  { \
    if (dr2 > 0.0f) CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0.0f) CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2)\
    else return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);\
    }}}



   
int tri_tri_overlap_test_3d( Eigen::MatrixXd p0, Eigen::MatrixXd q0)
    {
        Eigen::MatrixXd p = perturb(p0);
        Eigen::MatrixXd q = perturb(q0);
       
        real p1[3], q1[3], r1[3], p2[3], q2[3], r2[3];
        p1[0] = p(0,0);
        p1[1] = p(0,1);
        p1[2] = p(0,2);

        q1[0] = p(1,0);
        q1[1] = p(1,1);
        q1[2] = p(1,2);

        r1[0] = p(2,0);
        r1[1] = p(2,1);
        r1[2] = p(2,2);

        p2[0] = q(0,0);
        p2[1] = q(0,1);
        p2[2] = q(0,2);

        q2[0] = q(1,0);
        q2[1] = q(1,1);
        q2[2] = q(1,2);

        r2[0] = q(2,0);
        r2[1] = q(2,1);
        r2[2] = q(2,2);
        
        
        real dp1, dq1, dr1, dp2, dq2, dr2;
        real v1[3], v2[3];
        real N1[3], N2[3]; 

        /* Compute distance signs  of p1, q1 and r1 to the plane of
         triangle(p2,q2,r2) */


        SUB(v1,p2,r2)
        SUB(v2,q2,r2)
        CROSS(N2,v1,v2)

        SUB(v1,p1,r2)
        dp1 = DOT(v1,N2);
        SUB(v1,q1,r2)
        dq1 = DOT(v1,N2);
        SUB(v1,r1,r2)
        dr1 = DOT(v1,N2);

        if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return 0; 

        /* Compute distance signs  of p2, q2 and r2 to the plane of
         triangle(p1,q1,r1) */


        SUB(v1,q1,p1)
        SUB(v2,r1,p1)
        CROSS(N1,v1,v2)

        SUB(v1,p2,r1)
        dp2 = DOT(v1,N1);
        SUB(v1,q2,r1)
        dq2 = DOT(v1,N1);
        SUB(v1,r2,r1)
        dr2 = DOT(v1,N1);

        if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return 0;

        /* Permutation in a canonical form of T1's vertices */




        if (dp1 > 0.0f) {
            if (dq1 > 0.0f) TRI_TRI_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
                else if (dr1 > 0.0f) TRI_TRI_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)  
                    else TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
                        } else if (dp1 < 0.0f) {
                            if (dq1 < 0.0f) TRI_TRI_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
                                else if (dr1 < 0.0f) TRI_TRI_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
                                    else TRI_TRI_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
                                        } else {
                                            if (dq1 < 0.0f) {
                                                if (dr1 >= 0.0f) TRI_TRI_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
                                                    else TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
                                                        }
                                            else if (dq1 > 0.0f) {
                                                if (dr1 > 0.0f) TRI_TRI_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
                                                    else TRI_TRI_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
                                                        }
                                            else  {
                                                if (dr1 > 0.0f) TRI_TRI_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
                                                    else if (dr1 < 0.0f) TRI_TRI_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
                                                        else return 2*coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);
                                            }
                                        }
    };



    int coplanar_tri_tri3d(real p1[3], real q1[3], real r1[3],
                           real p2[3], real q2[3], real r2[3],
                           real normal_1[3], real normal_2[3]){

        real P1[2],Q1[2],R1[2];
        real P2[2],Q2[2],R2[2];

        real n_x, n_y, n_z;

        n_x = ((normal_1[0]<0)?-normal_1[0]:normal_1[0]);
        n_y = ((normal_1[1]<0)?-normal_1[1]:normal_1[1]);
        n_z = ((normal_1[2]<0)?-normal_1[2]:normal_1[2]);


        /* Projection of the triangles in 3D onto 2D such that the area of
         the projection is maximized. */


        if (( n_x > n_z ) && ( n_x >= n_y )) {
            // Project onto plane YZ

            P1[0] = q1[2]; P1[1] = q1[1];
            Q1[0] = p1[2]; Q1[1] = p1[1];
            R1[0] = r1[2]; R1[1] = r1[1]; 

            P2[0] = q2[2]; P2[1] = q2[1];
            Q2[0] = p2[2]; Q2[1] = p2[1];
            R2[0] = r2[2]; R2[1] = r2[1]; 

        } else if (( n_y > n_z ) && ( n_y >= n_x )) {
            // Project onto plane XZ

            P1[0] = q1[0]; P1[1] = q1[2];
            Q1[0] = p1[0]; Q1[1] = p1[2];
            R1[0] = r1[0]; R1[1] = r1[2]; 

            P2[0] = q2[0]; P2[1] = q2[2];
            Q2[0] = p2[0]; Q2[1] = p2[2];
            R2[0] = r2[0]; R2[1] = r2[2]; 

        } else {
            // Project onto plane XY

            P1[0] = p1[0]; P1[1] = p1[1]; 
            Q1[0] = q1[0]; Q1[1] = q1[1]; 
            R1[0] = r1[0]; R1[1] = r1[1]; 

            P2[0] = p2[0]; P2[1] = p2[1]; 
            Q2[0] = q2[0]; Q2[1] = q2[1]; 
            R2[0] = r2[0]; R2[1] = r2[1]; 
        }

        return tri_tri_overlap_test_2d(P1,Q1,R1,P2,Q2,R2);

    };



    /*
     *                                                                
     *  Three-dimensional Triangle-Triangle Intersection              
     *
     */

    /*
     This macro is called when the triangles surely intersect
     It constructs the segment of intersection of the two triangles
     if they are not coplanar.
     */

    #define CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) { \
    SUB(v1,q1,p1) \
    SUB(v2,r2,p1) \
    CROSS(N,v1,v2) \
    SUB(v,p2,p1) \
    if (DOT(v,N) > 0.0f) {\
    SUB(v1,r1,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) <= 0.0f) { \
    SUB(v2,q2,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) > 0.0f) { \
    SUB(v1,p1,p2) \
    SUB(v2,p1,r1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p1,v1) \
    SUB(v1,p2,p1) \
    SUB(v2,p2,r2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p2,v1) \
    return 1; \
    } else { \
    SUB(v1,p2,p1) \
    SUB(v2,p2,q2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p2,v1) \
    SUB(v1,p2,p1) \
    SUB(v2,p2,r2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p2,v1) \
    return 1; \
    } \
    } else { \
    return 0; \
    } \
    } else { \
    SUB(v2,q2,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) < 0.0f) { \
    return 0; \
    } else { \
    SUB(v1,r1,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) >= 0.0f) { \
    SUB(v1,p1,p2) \
    SUB(v2,p1,r1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p1,v1) \
    SUB(v1,p1,p2) \
    SUB(v2,p1,q1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p1,v1) \
    return 1; \
    } else { \
    SUB(v1,p2,p1) \
    SUB(v2,p2,q2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p2,v1) \
    SUB(v1,p1,p2) \
    SUB(v2,p1,q1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p1,v1) \
    return 1; \
    }}}} 



    #define TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
    if (dp2 > 0.0f) { \
    if (dq2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2) \
    else if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
    else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) }\
    else if (dp2 < 0.0f) { \
    if (dq2 < 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    else CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
    } else { \
    if (dq2 < 0.0f) { \
    if (dr2 >= 0.0f)  CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
    else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2)\
    } \
    else if (dq2 > 0.0f) { \
    if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
    else  CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
    } \
    else  { \
    if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
    else if (dr2 < 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2)\
    else { \
    *coplanar = 1; \
    return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);\
    } \
    }} }


    /*
     The following version computes the segment of intersection of the
     two triangles if it exists. 
     coplanar returns whether the triangles are coplanar
     source and target are the endpoints of the line segment of intersection 
     */

    int tri_tri_intersection_test_3d(real p1[3], real q1[3], real r1[3], 
                                     real p2[3], real q2[3], real r2[3],
                                     int * coplanar, 
                                     real source[3], real target[3] )

    {
        real dp1, dq1, dr1, dp2, dq2, dr2;
        real v1[3], v2[3], v[3];
        real N1[3], N2[3], N[3];
        real alpha;

        // Compute distance signs  of p1, q1 and r1 
        // to the plane of triangle(p2,q2,r2)


        SUB(v1,p2,r2)
        SUB(v2,q2,r2)
        CROSS(N2,v1,v2)

        SUB(v1,p1,r2)
        dp1 = DOT(v1,N2);
        SUB(v1,q1,r2)
        dq1 = DOT(v1,N2);
        SUB(v1,r1,r2)
        dr1 = DOT(v1,N2);

        if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return 0; 

        // Compute distance signs  of p2, q2 and r2 
        // to the plane of triangle(p1,q1,r1)


        SUB(v1,q1,p1)
        SUB(v2,r1,p1)
        CROSS(N1,v1,v2)

        SUB(v1,p2,r1)
        dp2 = DOT(v1,N1);
        SUB(v1,q2,r1)
        dq2 = DOT(v1,N1);
        SUB(v1,r2,r1)
        dr2 = DOT(v1,N1);

        if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return 0;

        // Permutation in a canonical form of T1's vertices


        //  printf("d1 = [%f %f %f], d2 = [%f %f %f]\n", dp1, dq1, dr1, dp2, dq2, dr2);
        /*
         // added by Aaron
         if (ZERO_TEST(dp1) || ZERO_TEST(dq1) ||ZERO_TEST(dr1) ||ZERO_TEST(dp2) ||ZERO_TEST(dq2) ||ZERO_TEST(dr2))
         {
         coplanar = 1;
         return 0;
         }
         */


        if (dp1 > 0.0f) {
            if (dq1 > 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
                else if (dr1 > 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)

                    else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
                        } else if (dp1 < 0.0f) {
                            if (dq1 < 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
                                else if (dr1 < 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
                                    else TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
                                        } else {
                                            if (dq1 < 0.0f) {
                                                if (dr1 >= 0.0f) TRI_TRI_INTER_3D(q1,r1,p1,p2,r2,q2,dp2,dr2,dq2)
                                                    else TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2)
                                                        }
                                            else if (dq1 > 0.0f) {
                                                if (dr1 > 0.0f) TRI_TRI_INTER_3D(p1,q1,r1,p2,r2,q2,dp2,dr2,dq2)
                                                    else TRI_TRI_INTER_3D(q1,r1,p1,p2,q2,r2,dp2,dq2,dr2)
                                                        }
                                            else  {
                                                if (dr1 > 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,q2,r2,dp2,dq2,dr2)
                                                    else if (dr1 < 0.0f) TRI_TRI_INTER_3D(r1,p1,q1,p2,r2,q2,dp2,dr2,dq2)
                                                        else {
                                                            // triangles are co-planar

                                                            *coplanar = 1;
                                                            return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);
                                                        }
                                            }
                                        }
    };





    /*
     *
     *  Two dimensional Triangle-Triangle Overlap Test    
     *
     */


    /* some 2D macros */

    #define ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))


    #define INTERSECTION_TEST_VERTEXA(P1, Q1, R1, P2, Q2, R2) {\
    if (ORIENT_2D(R2,P2,Q1) >= 0.0f)\
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0f)\
    if (ORIENT_2D(P1,P2,Q1) > 0.0f) {\
    if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) return 1; \
    else return 0;} else {\
    if (ORIENT_2D(P1,P2,R1) >= 0.0f)\
    if (ORIENT_2D(Q1,R1,P2) >= 0.0f) return 1; \
    else return 0;\
    else return 0;}\
    else \
    if (ORIENT_2D(P1,Q2,Q1) <= 0.0f)\
    if (ORIENT_2D(R2,Q2,R1) <= 0.0f)\
    if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) return 1; \
    else return 0;\
    else return 0;\
    else return 0;\
    else\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) \
    if (ORIENT_2D(Q1,R1,R2) >= 0.0f)\
    if (ORIENT_2D(P1,P2,R1) >= 0.0f) return 1;\
    else return 0;\
    else \
    if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) {\
    if (ORIENT_2D(R2,R1,Q2) >= 0.0f) return 1; \
    else return 0; }\
    else return 0; \
    else  return 0; \
    };

    #define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) {\
      if (ORIENT_2D(R2,P2,Q1) >= 0.0f)\
        if (ORIENT_2D(R2,Q2,Q1) <= 0.0f)\
          if (ORIENT_2D(P1,P2,Q1) > 0.0f) {\
            if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) return 1; \
            else return 0;} else {\
            if (ORIENT_2D(P1,P2,R1) >= 0.0f)\
              if (ORIENT_2D(Q1,R1,P2) >= 0.0f) return 1; \
              else return 0;\
            else return 0;}\
        else \
          if (ORIENT_2D(P1,Q2,Q1) <= 0.0f)\
            if (ORIENT_2D(R2,Q2,R1) <= 0.0f)\
              if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) return 1; \
              else return 0;\
            else return 0;\
          else return 0;\
      else\
        if (ORIENT_2D(R2,P2,R1) >= 0.0f) \
          if (ORIENT_2D(Q1,R1,R2) >= 0.0f)\
            if (ORIENT_2D(P1,P2,R1) >= 0.0f) return 1;\
            else return 0;\
          else \
            if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) {\
              if (ORIENT_2D(R2,R1,Q2) >= 0.0f) return 1; \
              else return 0; }\
            else return 0; \
        else  return 0; \
     };


    #define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
    if (ORIENT_2D(R2,P2,Q1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
    if (ORIENT_2D(P1,Q1,R2) >= 0.0f) return 1; \
    else return 0;} else { \
    if (ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
    if (ORIENT_2D(R1,P1,P2) >= 0.0f) return 1; else return 0;} \
    else return 0; } \
    } else {\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,R1) >= 0.0f) {\
    if (ORIENT_2D(P1,R1,R2) >= 0.0f) return 1;  \
    else {\
    if (ORIENT_2D(Q1,R1,R2) >= 0.0f) return 1; else return 0;}}\
    else  return 0; }\
    else return 0; }}



    int ccw_tri_tri_intersection_2d(real p1[2], real q1[2], real r1[2], 
                                    real p2[2], real q2[2], real r2[2]) {
        if ( ORIENT_2D(p2,q2,p1) >= 0.0f ) {
            if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
                if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) return 1;
                else INTERSECTION_TEST_EDGE(p1,q1,r1,p2,q2,r2)
                    } else {  
                        if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
                            INTERSECTION_TEST_EDGE(p1,q1,r1,r2,p2,q2)
                            else INTERSECTION_TEST_VERTEX(p1,q1,r1,p2,q2,r2)}}
        else {
            if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
                if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
                    INTERSECTION_TEST_EDGE(p1,q1,r1,q2,r2,p2)
                    else  INTERSECTION_TEST_VERTEX(p1,q1,r1,q2,r2,p2)}
            else INTERSECTION_TEST_VERTEX(p1,q1,r1,r2,p2,q2)}
    };


    int tri_tri_overlap_test_2d(real p1[2], real q1[2], real r1[2], 
                                real p2[2], real q2[2], real r2[2]) {
        if ( ORIENT_2D(p1,q1,r1) < 0.0f )
            if ( ORIENT_2D(p2,q2,r2) < 0.0f )
                return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,r2,q2);
            else
                return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,q2,r2);
            else
                if ( ORIENT_2D(p2,q2,r2) < 0.0f )
                    return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,r2,q2);
                else
                    return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,q2,r2);

    };

