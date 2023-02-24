#pragma once
#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "meshmodel.h"



std::vector<int> selfIntersection(Mesh mesh, Eigen::MatrixXd nrm);

std::vector<int> invalidTriangles(Mesh mesh,  Eigen::MatrixXd nrm);

Eigen::MatrixXd cross(Eigen::VectorXd v,Eigen::MatrixXd b);


#endif