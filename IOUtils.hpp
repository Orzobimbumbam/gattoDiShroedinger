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
#include <map>
#include <string>

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

template <class T, class D> std::istream& readMap(const std::map<T, D>& outputMap, std::istream& rStream, bool header)
{
    if (header)
    {
        std::string headerLine;
        std::getline(rStream, headerLine);
    }
    
    for (const auto& it : outputMap)
        rStream >> it.first >> "\t" >> it.second >> std::endl;
    
    return rStream;
}

template <class T> std::ostream& writeMatrix(const std::vector<std::vector<T>>& inputMatrix, std::ostream& wStream, bool header)
{
    if (header)
    {
        //need to overload << for vectors to do this
    }
    
    for (const auto& rows : inputMatrix)
        for (const auto& elements : rows)
            //need to overload << for vectors to do this
            
    return wStream;
}


#endif /* IOUtils_hpp */
