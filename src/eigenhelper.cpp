

/* Author: Abhijith


Decription: Helper class for Eigen function. Implemented finding unique row in
Eigen matrix. Writing to csv file functionality is also added.

*/

#include "eigenhelper.h"

Eigen::MatrixXd crossm(Eigen::MatrixXd a,Eigen::MatrixXd b){
   // std::cout<<"In crossm\n";
    Eigen::MatrixXd res(b.rows(),3);
    //std::cout<<"In crossm - 1\n";
    //std::cout<<a.rows()<<" cross "<<a.cols()<<std::endl;
    //std::cout<<b.rows()<<" cross "<<b.cols()<<std::endl;
    //std::cout<<res.rows()<<" cross "<<res.cols()<<std::endl;
     res.col(0) = a.col(1).cwiseProduct(b.col(2)) - a.col(2).cwiseProduct(b.col(1));
    res.col(1) = a.col(2).cwiseProduct(b.col(0)) - a.col(0).cwiseProduct(b.col(2));
    res.col(2) = a.col(0).cwiseProduct(b.col(1)) - a.col(1).cwiseProduct(b.col(0));
    //std::cout<<"In crossm - 2\n";
    return res;
}

Eigen::MatrixXd cross(Eigen::VectorXd v,Eigen::MatrixXd b){
    Eigen::MatrixXd res(b.rows(),3);
    //std::cout<<v.rows()<<" x "<<v.cols()<<std::endl;
    Eigen::MatrixXd a(b.rows(),b.cols());
    a = v.replicate(1,b.rows()).transpose();

    //std::cout<<a.rows()<<" xx "<<a.cols()<<std::endl;
    res.col(0) = a.col(1).cwiseProduct(b.col(2)) - a.col(2).cwiseProduct(b.col(1));
    res.col(1) = a.col(2).cwiseProduct(b.col(0)) - a.col(0).cwiseProduct(b.col(2));
    res.col(2) = a.col(0).cwiseProduct(b.col(1)) - a.col(1).cwiseProduct(b.col(0));

    return res;
}

Eigen::VectorXd dot(Eigen::MatrixXd a,Eigen::MatrixXd b){
   
    //std::cout<<b.rows()<<" "<<b.cols()<<std::endl;
    Eigen::VectorXd res(b.rows(),1);
    res = a.col(0).cwiseProduct(b.col(0)) + a.col(1).cwiseProduct(b.col(1)) +a.col(2).cwiseProduct(b.col(2)) ;
    
    return res;
}

Eigen::VectorXd dot(Eigen::VectorXd v,Eigen::MatrixXd b){
    Eigen::VectorXd res(b.rows(),1);
    //std::cout<<v.rows()<<" x "<<v.cols()<<std::endl;
    Eigen::MatrixXd a(b.rows(),b.cols());

    a = v.replicate(1,b.rows()).transpose();

    //std::cout<<a.rows()<<" x "<<a.cols()<<std::endl;
    res = a.col(0).cwiseProduct(b.col(0)) + a.col(1).cwiseProduct(b.col(1)) +a.col(2).cwiseProduct(b.col(2)) ;
    return res;
    
}

Eigen::MatrixXd diff(Eigen::VectorXd v,Eigen::MatrixXd b){
    Eigen::MatrixXd res(b.rows(),1);
    //std::cout<<v.rows()<<" x "<<v.cols()<<std::endl;
    Eigen::MatrixXd a(b.rows(),b.cols());
    a = v.replicate(1,b.rows()).transpose();

    //std::cout<<a.rows()<<" x "<<a.cols()<<std::endl;
    res  = a - b;
    return res;
    
}

// Get unique values of eigen matrix rows, fast with indices.
Eigen::VectorXi getUnique(Eigen::MatrixXd& A_nx3)
{

    Eigen::MatrixXd B(A_nx3.rows(), A_nx3.cols() + 1);
    B.leftCols(3) = A_nx3.leftCols(3);
    B.rightCols(1) = Eigen::VectorXd::LinSpaced(A_nx3.rows(), 0, A_nx3.rows() - 1);

  
    std::vector<Eigen::VectorXd> vec;
    for (int64_t i = 0; i < B.rows(); ++i)
        vec.push_back(B.row(i));


    std::sort(vec.begin(), vec.end(), [](Eigen::VectorXd const& t1, Eigen::VectorXd const& t2) {if (t1(0) != t2(0)) return t1(0) < t2(0); else if (t1(1) != t2(1)) return t1(1) < t2(1); else return t1(2) < t2(2); });

    
    for (int64_t i = 0; i < vec.size(); ++i)
        B.row(i) = vec[i];

    A_nx3.leftCols(3) = B.leftCols(3);
    Eigen::VectorXi ix = B.col(3).cast<int>();


    std::set<Eigen::VectorXd, vecCompare> s{};
    std::vector<int> uniqueidx{};
    int k = 0;
    int t = 0;
    std::vector<int> v(A_nx3.rows());
    v[ix[0]] = 0;
    for (int64_t i = 1; i < A_nx3.rows(); ++i) {
        //if(s.insert(A_nx3.row(i)).second)
        if (A_nx3(i - 1, 2) != A_nx3(i, 2) || A_nx3(i - 1, 1) != A_nx3(i, 1) || A_nx3(i - 1, 0) != A_nx3(i, 0))
            k = k + 1;

        v[ix[i]] = k;

    }
    ix = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(v.data(), v.size());

    vec.clear();
    for (int64_t i = 0; i < A_nx3.rows(); ++i)
        vec.push_back(A_nx3.row(i));

    auto it = std::unique(vec.begin(), vec.end());
    vec.resize(std::distance(vec.begin(), it));

    A_nx3.resize(vec.size(), 3);
    for (int64_t i = 0; i < vec.size(); ++i)
        A_nx3.row(i) = vec[i];


    return ix;

}

//Writing to CSV file
void write_csv(std::string name, Eigen::MatrixXd mmm)
{
    std::ofstream file(name.c_str());
    file << mmm.format(CSVFormat);
}

void write_csv(std::string name, Eigen::VectorXi mmm)
{
    std::ofstream file(name.c_str());
    file << mmm.format(CSVFormat);
}

void write_csv(std::string name, Eigen::VectorXd mmm)
{
    std::ofstream file(name.c_str());
    file << mmm.format(CSVFormat);
}

void write_csv(std::string name, Eigen::MatrixXi mmm)
{
    std::ofstream file(name.c_str());
    file << mmm.format(CSVFormat);
}

void write_csv(std::string name, MatrixXl mmm)
{
    std::ofstream file(name.c_str());
    file << mmm.format(CSVFormat);
}
void write_csv(std::string name, VectorXl mmm)
{
    std::ofstream file(name.c_str());
    file << mmm.format(CSVFormat);
}

