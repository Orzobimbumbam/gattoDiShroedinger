#ifndef IOUtils_hpp
#define IOUtils_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

//declare IO global functions
//can have templetized types and be friended with relevant classes for I/O operator overloading
//NB the following functions are tested to work on happy-path cases. Corner case behaviour must be specified or handled somewhere else.
template <class T, class D> std::ostream& writeMap(const std::map<T, D>& inputMap, std::ostream& wStream, bool header)
{
    if (header) // header is the first row of the file used for titles columns. If present set boolean true, else false.
        wStream << "Key" << "\t" << "Value" << std::endl;
    
    for (const auto& it : inputMap)
        wStream << it.first << "\t" << it.second << std::endl;
    
    return wStream;
}

template <class T, class D> std::istream& readMap(std::map<T, D>& outputMap, std::istream& rStream, bool header)
{
    if (header) // header is the first row of the file used for titles columns. If present set boolean true, else false.
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
        wStream << it << "\t" ; //tab delimiter, all values on same row
    
    return wStream;
}

template <class T> std::ostream& writeMatrix(const std::vector<std::vector<T>>& inputMatrix, std::ostream& wStream, bool header, const std::vector<std::string> headerLabels = std::vector<std::string>())
{
    if (header && !headerLabels.empty()) // header is the first row of the file used for titles columns. If present set boolean true, else false.
        writeRowVector(headerLabels, wStream);
    
    for (const auto& rows : inputMatrix)
        writeRowVector(rows, wStream) << std::endl;
    
    return wStream;
}

template <class T> std::istream& readMatrix(std::vector<std::vector<T>>& outputMatrix, std::istream& rStream, bool header)
{
    if (header) // header is the first row of the file used for titles columns. If present set boolean true, else false.
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
            outputMatrix[rowIndex].push_back(std::stod(value)); //convert string to double
        
        ++rowIndex;
     }
    return rStream;
}




#endif /* IOUtils_hpp */
