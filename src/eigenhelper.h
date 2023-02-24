/* Author: Abhijith


Decription: Helper class for Eigen function. Implemented finding unique row in
Eigen matrix. Writing to csv file functionality is also added.

*/

#pragma once
#ifndef EIGENHELPER_H
#define EIGENHELPER_H

//#define eigen_assert(X) do { if(!(X)) throw std::runtime_error(#X); } while(false);

#include<fstream>
#include<string>
#include <set>
#include<iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include<Eigen/SparseLU>


#define FIND(a,elem) (a(0)==elem?0:(a(1)==elem?1:2))

typedef   Eigen::Matrix<unsigned long long, Eigen::Dynamic, Eigen::Dynamic> MatrixXl;
typedef   Eigen::Matrix<unsigned long long, Eigen::Dynamic, 1> VectorXl;


Eigen::MatrixXd crossm(Eigen::MatrixXd a,Eigen::MatrixXd b);
Eigen::MatrixXd cross(Eigen::VectorXd v,Eigen::MatrixXd b);
Eigen::VectorXd dot(Eigen::MatrixXd a,Eigen::MatrixXd b);
Eigen::VectorXd dot(Eigen::VectorXd v,Eigen::MatrixXd b);
Eigen::MatrixXd diff(Eigen::VectorXd v,Eigen::MatrixXd b);

Eigen::VectorXi getUnique(Eigen::MatrixXd& A_nx3);

const static Eigen::IOFormat CSVFormat(18, Eigen::DontAlignCols, ", ", "\n");

void write_csv(std::string name, MatrixXl mmm);
void write_csv(std::string name, Eigen::MatrixXd mmm);
void write_csv(std::string name, Eigen::MatrixXi mmm);
void write_csv(std::string name, Eigen::VectorXi mmm);
void write_csv(std::string name, Eigen::VectorXd mmm);
void write_csv(std::string name, VectorXl mmm);

struct vecCompare {
    bool operator()(const Eigen::VectorXd& v, const Eigen::VectorXd& w) const
    {
        int n = v.rows();

        for (int i = 0; i < n; ++i) {
            if (v(i) < w(i)) {
                return true;
            }
            else if (v(i) > w(i)) {
                return false;
            }
        }

        return false;
    }
};

#endif