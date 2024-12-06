//
// MIT License
// Copyright (c) 2020 Jonathan R. Madsen
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED
// "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
// LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//

#pragma once

#include <cctype>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>  // IWYU pragma: keep
#include <string>
#include <tuple>
#include <utility>

namespace PTL
{
//--------------------------------------------------------------------------------------//
//  use this function to get an environment variable setting +
//  a default if not defined, e.g.
//      int num_threads =
//          GetEnv<int>("FORCENUMBEROFTHREADS",
//                          std::thread::hardware_concurrency());
//
template <typename Tp>
Tp
GetEnv(const std::string& env_id, Tp _default = Tp())
{
    char* env_var = std::getenv(env_id.c_str());
    if(env_var)
    {
        std::string        str_var = std::string(env_var);
        std::istringstream iss(str_var);
        Tp                 var = Tp();
        iss >> var;
        return var;
    }

    // return default if not specified in environment
    return _default;
}

//--------------------------------------------------------------------------------------//
//  overload for boolean
//
template <>
inline bool
GetEnv(const std::string& env_id, bool _default)
{
    char* env_var = std::getenv(env_id.c_str());
    if(env_var)
    {
        std::string var = std::string(env_var);
        bool        val = true;
        if(var.find_first_not_of("0123456789") == std::string::npos)
            val = (bool) atoi(var.c_str());
        else
        {
            for(auto& itr : var)
                itr = tolower(itr);
            if(var == "off" || var == "false")
                val = false;
        }
        return val;
    }

    // return default if not specified in environment
    return _default;
}

//--------------------------------------------------------------------------------------//

}  // namespace PTL
