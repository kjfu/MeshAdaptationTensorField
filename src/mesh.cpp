#include "mesh.h"
#include"omp.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
namespace MeshMetric{

    void Mesh::readyCalculation(){
        node2nodes.resize(nodes.rows());
        
        const int otherIndices[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}}; 
        
        for(int i=0; i<tets.rows(); i++){
            for(int j=0; j<4; j++){
                int inode = tets(i,j);
                for (int k=0; k<3; k++){
                    node2nodes[inode].insert(tets(i,otherIndices[j][k]));
                }
            }
        }
    }

    void Mesh::calculateVolumes(){
        tetVolumes.resize(tets.rows());
        #pragma omp parallel for
        for (int i = 0; i < tets.rows(); i++){
            Eigen::Vector3d a = nodes.row(tets(i,0));
            Eigen::Vector3d b = nodes.row(tets(i,1));
            Eigen::Vector3d c = nodes.row(tets(i,2));
            Eigen::Vector3d d = nodes.row(tets(i,3));

            Eigen::Vector3d ab = b-a;
            Eigen::Vector3d ac = c-a;
            Eigen::Vector3d ad = d-a;
            tetVolumes[i] = fabs((ab.cross(ac)).dot(ad))/6.0;
        }
    }

    void Mesh::calculateNodeValueGradients(){
        nodeValueGradients.resize(nodes.rows());
        #pragma omp parallel for
        for(int i=0; i<nodes.rows(); i++){
            std::vector<double> a(6,0);
            std::vector<double> b(3,0);
            double u=nodeValues[i];
              /* Gradient: ui = u + <grad u,PPi> */

            for(auto j: node2nodes[i]){
                double u1 = nodeValues[j];
                double du = u1-u;
                /* M = At*A symmetric definite positive */
                Eigen::Vector3d da = nodes.row(j) - nodes.row(i); 
                a[0] += da.x() * da.x();
                a[1] += da.x() * da.y();
                a[2] += da.x() * da.z();
                a[3] += da.y() * da.y();
                a[4] += da.y() * da.z();
                a[5] += da.z() * da.z();

                /* b = At*du */
                b[0] += da.x() * du;
                b[1] += da.y() * du;
                b[2] += da.z() * du;
            }

            /* solution of A(3,3)*grad(1,3) = b(1,3) */ 
            double aa = a[3]*a[5] - a[4]*a[4];
            double bb = a[4]*a[2] - a[1]*a[5];
            double cc = a[1]*a[4] - a[2]*a[3];
            double du = aa*a[0] + bb*a[1] + cc*a[2];

            du = 1.0 / du;

            double dd = a[0]*a[5] - a[2]*a[2];
            double ee = a[1]*a[2] - a[0]*a[4];
            double ff = a[0]*a[3] - a[1]*a[1];

            double grd[3];
            grd[0] = (aa*b[0] + bb*b[1] + cc*b[2]) * du;
            grd[1] = (bb*b[0] + dd*b[1] + ee*b[2]) * du;
            grd[2] = (cc*b[0] + ee*b[1] + ff*b[2]) * du;

            nodeValueGradients[i] << grd[0], grd[1], grd[2];
        }
    }

    void Mesh::calculateNodeValueHessiansDirectly(){
        nodeValueHessians.resize(nodes.rows());
        nodeValueAbsoluteHessians.resize(nodes.rows());
          /* Hessian: Ui = U + <gradU,PPi> + 0.5*(tPPi.Hess.PPi) */
        #pragma omp parallel for
        for(int i = 0; i<nodes.rows(); i++){
            double u = nodeValues[i];
            Eigen::Vector3d grd = nodeValueGradients[i];
            Eigen::MatrixXd m(node2nodes[i].size(), 9);
            m.setConstant(0);
            Eigen::MatrixXd b(node2nodes[i].size(), 1);
            b.setConstant(0);
            int irow=0;
            for(auto j: node2nodes[i]){
                double u1 = nodeValues[j];
                Eigen::Vector3d da = nodes.row(j) - nodes.row(i); 
                
                double a[6];
                a[0] = 0.5 * da.x()*da.x();
                a[1] =       da.x()*da.y();
                a[2] =       da.x()*da.z();
                a[3] = 0.5 * da.y()*da.y();
                a[4] =       da.y()*da.z();
                a[5] = 0.5 * da.z()*da.z();                
                double du   = u1-u;
                m.row(irow) << da.x(), da.y(), da.z(), a[0], a[1], a[2], a[3], a[4], a[5];
                b.row(irow) << du;
                irow++;      
            }
            // Eigen::MatrixXd m(6,6);
            // m.setConstant(0);
            // Eigen::MatrixXd b(6,1);
            // b.setConstant(0);
            // for(auto j: node2nodes[i]){
            //     double u1 = nodeValues[j];
            //     Eigen::Vector3d da = nodes.row(j) - nodes.row(i); 
                
            //     double a[6];
            //     a[0] = 0.5 * da.x()*da.x();
            //     a[1] =       da.x()*da.y();
            //     a[2] =       da.x()*da.z();
            //     a[3] = 0.5 * da.y()*da.y();
            //     a[4] =       da.y()*da.z();
            //     a[5] = 0.5 * da.z()*da.z();
            //     double du   = (u1-u) - (da.x()*grd[0] + da.y()*grd[1] + da.z()*grd[2]);
            //     /* M = At*A symmetric positive definite */
            //     for (int jj=0; jj<6; jj++) {
            //         m(jj,jj) += a[jj]*a[jj];
            //         for (int l=jj+1; l<6; l++) {
            //             m(jj,l) += a[jj]*a[l];
            //             m(l,jj) += a[l]*a[jj];
            //         }
            //     }
                
            //     /* c = At*b */
            //     for (int jj=0; jj<6; jj++)
            //         b(jj,0) += a[jj]*du;        
            // }
            // Eigen::MatrixXd M = m.transpose()* m;
            // Eigen::MatrixXd B = m.transpose()* b;
            Eigen::MatrixXd grad_hessian = m.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
            Eigen::Matrix3d h;
            h  <<   grad_hessian(0+3, 0), grad_hessian(1+3, 0), grad_hessian(2+3, 0),
                    grad_hessian(1+3, 0), grad_hessian(3+3, 0), grad_hessian(4+3, 0),
                    grad_hessian(2+3, 0), grad_hessian(4+3, 0), grad_hessian(5+3, 0);
            Eigen::EigenSolver<Eigen::Matrix3d> solver(h);
            Eigen::Matrix3d D= solver.pseudoEigenvalueMatrix();
            Eigen::Matrix3d V = solver.pseudoEigenvectors();
            D(0,0) = fabs(D(0,0));
            D(1,1) = fabs(D(1,1));
            D(2,2) = fabs(D(2,2));
            // std::cout << D << std::endl;
            // std::cout <<"-----------" << std::endl;
            nodeValueHessians[i] = h;
            nodeValueAbsoluteHessians[i] = V*D*V.transpose();

        }        
    }


    void Mesh::calculateNodeValueHessians(){
        nodeValueHessians.resize(nodes.rows());
        nodeValueAbsoluteHessians.resize(nodes.rows());
          /* Hessian: Ui = U + <gradU,PPi> + 0.5*(tPPi.Hess.PPi) */
        #pragma omp parallel for
        for(int i = 0; i<nodes.rows(); i++){
            double u = nodeValues[i];
            Eigen::Vector3d grd = nodeValueGradients[i];
            Eigen::MatrixXd m(node2nodes[i].size(), 6);
            m.setConstant(0);
            Eigen::MatrixXd b(node2nodes[i].size(), 1);
            b.setConstant(0);
            int irow=0;
            for(auto j: node2nodes[i]){
                double u1 = nodeValues[j];
                Eigen::Vector3d da = nodes.row(j) - nodes.row(i); 
                
                double a[6];
                a[0] = 0.5 * da.x()*da.x();
                a[1] =       da.x()*da.y();
                a[2] =       da.x()*da.z();
                a[3] = 0.5 * da.y()*da.y();
                a[4] =       da.y()*da.z();
                a[5] = 0.5 * da.z()*da.z();                
                double du   = (u1-u) - (da.x()*grd[0] + da.y()*grd[1] + da.z()*grd[2]);
                m.row(irow) << a[0], a[1], a[2], a[3], a[4], a[5];
                b.row(irow) << du;
                irow++;      
            }
            // Eigen::MatrixXd m(6,6);
            // m.setConstant(0);
            // Eigen::MatrixXd b(6,1);
            // b.setConstant(0);
            // for(auto j: node2nodes[i]){
            //     double u1 = nodeValues[j];
            //     Eigen::Vector3d da = nodes.row(j) - nodes.row(i); 
                
            //     double a[6];
            //     a[0] = 0.5 * da.x()*da.x();
            //     a[1] =       da.x()*da.y();
            //     a[2] =       da.x()*da.z();
            //     a[3] = 0.5 * da.y()*da.y();
            //     a[4] =       da.y()*da.z();
            //     a[5] = 0.5 * da.z()*da.z();
            //     double du   = (u1-u) - (da.x()*grd[0] + da.y()*grd[1] + da.z()*grd[2]);
            //     /* M = At*A symmetric positive definite */
            //     for (int jj=0; jj<6; jj++) {
            //         m(jj,jj) += a[jj]*a[jj];
            //         for (int l=jj+1; l<6; l++) {
            //             m(jj,l) += a[jj]*a[l];
            //             m(l,jj) += a[l]*a[jj];
            //         }
            //     }
                
            //     /* c = At*b */
            //     for (int jj=0; jj<6; jj++)
            //         b(jj,0) += a[jj]*du;        
            // }
            // Eigen::MatrixXd M = m.transpose()* m;
            // Eigen::MatrixXd B = m.transpose()* b;
            Eigen::MatrixXd hessian = m.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
            // if (i==6115)std::cout << hessian << std::endl;
            Eigen::Matrix3d h;
            h  <<   hessian(0, 0), hessian(1, 0), hessian(2, 0),
                    hessian(1, 0), hessian(3, 0), hessian(4, 0),
                    hessian(2, 0), hessian(4, 0), hessian(5, 0);
            Eigen::EigenSolver<Eigen::Matrix3d> solver(h);
            Eigen::Matrix3d D= solver.pseudoEigenvalueMatrix();
            Eigen::Matrix3d V = solver.pseudoEigenvectors();
            D(0,0) = fabs(D(0,0));
            D(1,1) = fabs(D(1,1));
            D(2,2) = fabs(D(2,2));
            // std::cout << D << std::endl;
            // std::cout <<"-----------" << std::endl;
            nodeValueHessians[i] = h;
            nodeValueAbsoluteHessians[i] = V*D*V.transpose();

        }        
    }

    void Mesh::loadMesh(const std::string &filePath){
        std::cout << "* Load " << filePath << std::endl;
        std::ifstream file(filePath);
        if (file.is_open()){
            while (file)
            {
                std::string line;
                std::string keyString;
                std::getline(file, line);
                std::stringstream lineStream(line);
                lineStream >> keyString;

                if (keyString == "Vertices"){
                    std::getline(file, line);
                    std::stringstream lineStream(line);
                    int nv;
                    lineStream >> nv;
                    nodes.resize(nv,3);
                    for(int i=0; i<nv; i++){
                        std::getline(file, line);
                        std::stringstream vertexLineStream(line);
                        double x, y, z;
                        int label;
                        vertexLineStream >> x >> y >> z >> label;
                        nodes.row(i) << x,y,z;
                    }
                }
                else if (keyString == "Tetrahedra"){
                    std::getline(file, line);
                    std::stringstream lineStream(line);
                    int nt;
                    lineStream >> nt;
                    tets.resize(nt,4);
                    for(int i=0; i<nt;  i++){
                        std::getline(file, line);
                        std::stringstream tetLineStream(line);
                        int n0, n1, n2, n3, label;
                        tetLineStream >> n0 >> n1 >> n2 >> n3 >> label;
                        tets.row(i) << n0-1, n1-1, n2-1, n3-1;
                    }
                }
                else if (keyString == "End"){
                    break;
                }
            }
            file.close();
        }
    }

    void Mesh::loadTensorSol(const std::string &filePath){
        std::cout << "* Load " << filePath << std::endl;
        std::ifstream file(filePath);
        if (file.is_open()){
            while(file){
                std::string line;
                std::getline(file, line);
                std::stringstream lineStream(line);
                std::string keyString;
                lineStream >> keyString;
                if (keyString=="SolAtVertices"){
                    std::getline(file, line);
                    std::stringstream subStream(line);
                    int nv;
                    subStream >> nv;
                    std::getline(file, line);
                    nodeMetricTensors.resize(nv);
                    for (int i=0; i<nv; i++){
                        std::string solLine;
                        std::getline(file, solLine);
                        std::stringstream solLineStream(solLine);
                        double m11, m12, m22, m13, m23, m33;
                        solLineStream >> m11 >> m12 >> m22 >> m13 >> m23 >> m33;
                        nodeMetricTensors[i] << m11, m12, m13, 
                                                m12, m22, m23,
                                                m13, m23, m33;
                    }
                    break;
                }
            }
            file.close();
        }        
    }


    void Mesh::loadSol(const std::string &filePath){
        std::cout << "* Load " << filePath << std::endl;
        std::ifstream file(filePath);
        if (file.is_open()){
            while(file){
                std::string line;
                std::getline(file, line);
                std::stringstream lineStream(line);
                std::string keyString;
                lineStream >> keyString;
                if (keyString=="SolAtVertices"){
                    std::getline(file, line);
                    std::stringstream subStream(line);
                    int nv;
                    subStream >> nv;
                    std::getline(file, line);
                    for (int i=0; i<nv; i++){
                        std::string solLine;
                        std::getline(file, solLine);
                        std::stringstream solLineStream(solLine);
                        double s;
                        solLineStream >> s;
                        nodeValues.push_back(s);
                    }
                    break;
                }
            }
            file.close();
        }
    }

    void Mesh::saveSol(const std::string &filePath){
        std::cout << "* Save " << filePath << std::endl;
        std::ofstream file(filePath);
        file << "MeshVersionFormatted 2\n" << std::endl;
        file << "Dimension 3\n" << std::endl;
        file << "SolAtVertices" << std::endl;
        file << nodes.rows() << std::endl;
        file << "1 3" << std::endl;
        for(int i = 0; i < nodes.rows(); i++){
            Eigen::Matrix3d &m = nodeMetricTensors[i];
            file << m(0, 0) << " " 
            << m(0, 1) << " " 
            << m(1, 1) << " " 
            << m(0, 2) << " "
            << m(1, 2) << " "
            << m(2, 2) << std::endl;
        }
        file << "\nEnd" <<std::endl;
        file.close();
    }

    void Mesh::saveMet(const std::string &filePath){
        std::cout << "* Save " << filePath << std::endl;
        std::ofstream file(filePath);
        for (size_t i = 0; i < nodes.rows(); i++){
            Eigen::Matrix3d &m = nodeMetricTensors[i];
            file << nodes(i, 0) << " "
            << nodes(i,1) << " "
            << nodes(i,2) << " "
            << m(0, 0) << " " 
            << m(0, 1) << " "        
            << m(0, 2) << " " 
            << m(1, 1) << " " 
            << m(1, 2) << " "
            << m(2, 2) << std::endl;
        }
        file.close();
        
    }

    void Mesh::calculateDL(double N, double p){
        if (N<0){
            N = nodes.rows();
        }


        std:: vector<double> hessianNorms(nodes.rows());
        #pragma omp parallel for
        for (int i=0; i<nodes.rows(); i++){
            hessianNorms[i] = nodeValueAbsoluteHessians[i].norm();
        }

        double sum;
        for(int i=0; i<tets.rows();i++){
            double averageNorm =0.25*(hessianNorms[tets(i,0)]+hessianNorms[tets(i,1)]+hessianNorms[tets(i,2)]+hessianNorms[tets(i,3)]);
            sum+=pow(averageNorm, p/(2*p+3)) * tetVolumes[i];
        }
        DL = pow(N, 2.0/3.0) * pow(sum, -2.0/3.0);
        std::cout << "DL = " << DL << std::endl;
    }

    void Mesh::defineNodeMetricTensors(double p){
        nodeMetricTensors.resize(nodes.rows());
        #pragma omp parallel for
        for(int i = 0; i <nodes.rows(); i++){
            nodeMetricTensors[i] = DL * pow(nodeValueAbsoluteHessians[i].norm(), -1.0/(2*p+3)) * nodeValueAbsoluteHessians[i];
        }
    }

    void Mesh::generateMetricTensorField(double N, double p){
        readyCalculation();
        calculateNodeValueGradients();
        calculateNodeValueHessians();
        // calculateNodeValueHessiansDirectly();
        calculateVolumes();
        calculateDL(N, p);
        defineNodeMetricTensors(p);
    }

    void Mesh::saveVtk(const std::string &filePath){
        std::cout << "* Save " << filePath << std::endl;
        std::ofstream file(filePath);
        file << "# vtk DataFile Version 3.0" << std::endl;
        file << "vtk output"                 << std::endl;
        file << "ASCII"                      << std::endl;
        file << "DATASET UNSTRUCTURED_GRID"  << std::endl;
        file << "POINTS " << nodes.rows() << " double" << std::endl;
        for(int i=0; i<nodes.rows(); i++){
            file <<nodes(i, 0) <<  "  " << nodes(i, 1) << "  " << nodes(i, 2) << std::endl;
        }

        file << "CELLS " << tets.rows() << "  " << 5* tets.rows() << std::endl;
        for(int i=0; i<tets.rows(); i++){
            file << "4 " << tets(i, 0)
            << "  " << tets(i, 1)
            << "  " << tets(i, 2)
            << "  " << tets(i, 3)
            <<std::endl;
        }
        file << "CELL_TYPES " << tets.rows() << std::endl;
        for(int i=0; i<tets.rows(); i++){
            file << 10 << std::endl;
        }

        file << "POINT_DATA " << nodes.rows() << std::endl;            
        file << "SCALARS  det_metric double 1" << std::endl;
        file << "LOOKUP_TABLE default" << std::endl;
        for(int i=0; i<nodes.rows(); i++){
            file << nodeMetricTensors[i].norm() << std::endl;
        }


        file.close();
    }



}