/*
 * @Author: Kejie Fu
 * @Date: 2022-01-08 22:54:46
 * @LastEditTime: 2022-01-31 02:39:23
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAdaptationTensorField/src/mesh.h
 */
#pragma once

#include <vector>
#include <set>
#include <Eigen/Dense>
#include <string>
namespace MeshMetric{

#define CTE2D    2.0 / 9.0
#define CTE3D    9.0 / 32.0

class Mesh{
public:
    int normalization = 2;
    double eps;

    Eigen::MatrixXd nodes;
    Eigen::MatrixXd tets;
    double DL;
    std::vector<double> nodeValues;
    std::vector<std::set<int>> node2nodes;
    std::vector<Eigen::Vector3d> nodeValueGradients;
    // std::vector<std::array<double,6>> nodeValueHessians; // m00 m01 m02 m11 m12 m22
    std::vector<double> tetVolumes;
    std::vector<Eigen::Matrix3d> nodeValueHessians;
    std::vector<Eigen::Matrix3d> nodeValueAbsoluteHessians;
    // std::vector<Eigen::Matrix3d> nodeMetricTensorEigenvectors;
    // std::vector<Eigen::Matrix3d> nodeMetricTensorEigenvalues;
    std::vector<Eigen::Matrix3d> nodeMetricTensors;

    Mesh(){}
    void readyCalculation();
    void calculateVolumes();
    void calculateNodeValueGradients();
    void calculateNodeValueHessians();
    void calculateNodeValueHessiansDirectly();
    void calculateDL(double complexity, double norm);// N is complexity, and p is Lp norm
    void defineNodeMetricTensors(double norm);

    void generateMetricTensorField(double complexity,double norm);

    void loadMesh(const std::string &filePath);
    void loadSol(const std::string &filePath);
    void loadTensorSol(const std::string &filePath);
    
    void saveVtk(const std::string &filePath);
    void saveSol(const std::string &filePath);
    void saveMet(const std::string &filePath);


};

}