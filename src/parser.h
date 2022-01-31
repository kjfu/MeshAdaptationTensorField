/*
 * @Author: Kejie Fu
 * @Date: 2022-01-29 22:27:20
 * @LastEditTime: 2022-01-31 03:50:54
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshAdaptationTensorField/src/parser.h
 */
#pragma once
#include <string>
class Parser{
public:
    Parser(){}
    bool translateFile = false;
    std::string filePath;
    std::string meshFile;
    std::string solFile;
    std::string outVtkFile;
    std::string outSolFile;
    std::string outMetFile;
    double hmax;
    double hmin;
    double normNumber = 2;
    double complexity = -1;
    bool parseFilePathWithDot(const std::string &filePath, std::string &prefix, std::string &postfix);

    bool parseCMDLine(int argc, char **argv);


};