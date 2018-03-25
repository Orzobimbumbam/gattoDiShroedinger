//
//  IOUtils.hpp
//  Codice
//
//  Created by Alberto Campi on 24/03/2018.
//  Copyright © 2018 Alberto Campi. All rights reserved.
//

#ifndef IOUtils_hpp
#define IOUtils_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

//declare IO global functions
//can have templetized types
template <class T, class D> std::ostream& writeMap(const std::map<T, D>& inputMap, std::ostream& wStream, bool header)
{
    if (header)
        wStream << "Key" << "\t" << "Value" << std::endl;
    
    for (const auto& it : inputMap)
        wStream << it.first << "\t" << it.second << std::endl;
    
    return wStream;
}

template <class T, class D> std::ifstream& readMap(std::map<T, D>& outputMap, std::ifstream& rStream, bool header)
{
    if (header)
    {
        std::string headerLine;
        std::getline(rStream, headerLine);
    }
    
    while(true)
    {
        T tempKey;
        D tempValue;
        rStream >> tempKey >> tempValue;
        outputMap.insert(std::make_pair(tempKey, tempValue));
        
        if (rStream.eof())
            break;
        
    }
    
    return rStream;
}

template <class T> std::ostream& writeRowVector(const std::vector<T>& inputVector, std::ostream& wStream)
{
    for (const auto& it : inputVector)
        wStream << it << "\t" ; //tab delimiter
    
    return wStream;
}

template <class T> std::ostream& writeMatrix(const std::vector<std::vector<T>>& inputMatrix, std::ostream& wStream, bool header, const std::vector<std::string> headerLabels = std::vector<std::string>())
{
    if (header && !headerLabels.empty())
        writeRowVector(headerLabels, wStream);
    
    for (const auto& rows : inputMatrix)
        writeRowVector(rows, wStream) << std::endl;
    
    return wStream;
}

template <class T> std::istream& readMatrix(std::vector<std::vector<T>>& outputMatrix, std::istream& rStream, bool header)
{
    if (header)
    {
        std::string headerLine;
        std::getline(rStream, headerLine);
    }
    
    std::string line;
    unsigned int rowIndex = 0;
    while (std::getline(rStream, line))
    {
        
        std::istringstream sLine(line);
        std::string value;
        while (sLine >> value)
            outputMatrix[rowIndex].push_back(std::stod(value));
        
        ++rowIndex;
     }
    return rStream;
}





#endif /* IOUtils_hpp */
