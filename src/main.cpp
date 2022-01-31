/*
 * @Author: Kejie Fu
 * @Date: 2022-01-08 22:48:27
 * @LastEditTime: 2022-01-31 03:52:50
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAdaptationTensorField/src/main.cpp
 */
#include "mesh.h"
#include <string>
#include "parser.h"
#include <iostream>

int main(int argc, char* argv[]){
    Parser parser;
    parser.parseCMDLine(argc, argv);
    MeshMetric::Mesh mesh;
    
    if (parser.translateFile){
        mesh.loadMesh(parser.meshFile);
        mesh.loadTensorSol(parser.solFile);
        mesh.saveSol(parser.outSolFile);
        mesh.saveMet(parser.outMetFile);
        mesh.saveVtk(parser.outVtkFile);        
        return 0;
    }

    mesh.loadMesh(parser.meshFile);
    mesh.loadSol(parser.solFile);



    mesh.generateMetricTensorField(parser.complexity, parser.normNumber, parser.hmin, parser.hmax);
    
    mesh.saveSol(parser.outSolFile);
    mesh.saveMet(parser.outMetFile);
    mesh.saveVtk(parser.outVtkFile);
    std::cout << "Finished!" << std::endl;
    return 1;
}