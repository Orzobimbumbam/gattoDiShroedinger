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

//declare IO global functions
//can have templetized types
template <class T, class D> std::ostream& writeMap(const std::map<T, D>& inputMap, std::ostream& wStream)
{
    wStream << "Key" << "\t" << "Value" << std::endl;
    for (const auto& it : inputMap)
        wStream << it.first << "\t" << it.second << std::endl;
    
    return wStream;
}

#endif /* IOUtils_hpp */
