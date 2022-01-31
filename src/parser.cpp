/*
 * @Author: Kejie Fu
 * @Date: 2022-01-29 22:29:53
 * @LastEditTime: 2022-01-31 03:50:25
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAdaptationTensorField/src/parser.cpp
 */
#include "parser.h"
#include <iostream>
bool Parser::parseFilePathWithDot(const std::string &filePath, std::string &prefix, std::string &postfix){
    size_t pos = filePath.find_last_of('.');
    if (pos !=std::string::npos){
        prefix = filePath.substr(0, pos);
        postfix = filePath.substr(pos+1, filePath.length()-pos);
    }
    else{
        return false;
    }

    return true;
}

bool Parser::parseCMDLine(int argc, char **argv){
    for(int i=1; i<argc; i++){
        std::string key(argv[i]);
        if (key=="-in"){
            i++;
            filePath = argv[i];
        }
        else if (key=="-outSol"){
            i++;
            outSolFile = argv[i];
        }
        else if (key=="-outMet"){
            i++;
            outMetFile = argv[i];
        }
        else if (key=="-sol"){
            i++;
            solFile = argv[i];
        }
        else if (key=="-N"){
            i++;
            complexity = std::stod(argv[i]);
        }
        else if (key=="-p"){
            i++;
            normNumber = std::stod(argv[i]);
        }
        else if (key=="-t"){
            translateFile = true;
        }
        else if (key=="-hmax"){
            i++;
            hmax = std::stod(argv[i]);
        }
        else if (key=="-hmin"){
            i++;
            hmin = std::stod(argv[i]);            
        }
        else if(key=="-h"||key=="--help"){
            std::cout << "$MeshAdaptationTensorField -in input.mesh -sol input.sol -N complexity -p norm" << std::endl;
            exit(1);
        }
        // else if (key=="-hmax"){
        //     i++;
        //     hmax = std::stod(argv[i]);
        // }
        // else if (key=="-hmin"){
        //     i++;
        //     hmin = std::stod(argv[i]);
        // }

    }
        
    if (filePath.empty()){
        std::cout <<"ERROR: No FilePath specified" << std::endl;
        return false;
    }
    std::string prefix, postfix;
    if (parseFilePathWithDot(filePath, prefix, postfix)){
        meshFile = filePath;
    }
    else{
        prefix = filePath;
        meshFile = filePath + ".mesh";
    }
    if (solFile.empty()){
        solFile = prefix + ".sol";            
    }
    
    outMetFile = prefix + ".o.met";
    if (outSolFile.empty()){
        outSolFile = prefix + ".o.sol";
    }

    if (outVtkFile.empty()){
        outVtkFile = prefix + ".o.vtk";
    }

    if (outMetFile.empty()){
        outMetFile = prefix + ".o.met";
    }

    return true;    
}