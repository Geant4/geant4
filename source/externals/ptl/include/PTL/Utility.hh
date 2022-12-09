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
// Global utility functions
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
// use this function to get rid of "unused parameter" warnings
//
template <typename... Args>
void
ConsumeParameters(Args&&...)
{}

//--------------------------------------------------------------------------------------//
// a non-string environment option with a string identifier
template <typename Tp>
using EnvChoice = std::tuple<Tp, std::string, std::string>;

//--------------------------------------------------------------------------------------//
// list of environment choices with non-string and string identifiers
template <typename Tp>
using EnvChoiceList = std::set<EnvChoice<Tp>>;

//--------------------------------------------------------------------------------------//

class EnvSettings
{
public:
    using mutex_t    = std::mutex;
    using string_t   = std::string;
    using env_map_t  = std::multimap<string_t, string_t>;
    using env_pair_t = std::pair<string_t, string_t>;

public:
    static EnvSettings* GetInstance()
    {
        static EnvSettings* _instance = new EnvSettings();
        return _instance;
    }

public:
    template <typename Tp>
    void insert(const std::string& env_id, Tp val)
    {
        std::stringstream ss;
        ss << std::boolalpha << val;
        m_mutex.lock();
        if(m_env.find(env_id) != m_env.end())
        {
            for(const auto& itr : m_env)
                if(itr.first == env_id && itr.second == ss.str())
                {
                    m_mutex.unlock();
                    return;
                }
        }
        m_env.insert(env_pair_t(env_id, ss.str()));
        m_mutex.unlock();
    }

    template <typename Tp>
    void insert(const std::string& env_id, EnvChoice<Tp> choice)
    {
        Tp&          val      = std::get<0>(choice);
        std::string& str_val  = std::get<1>(choice);
        std::string& descript = std::get<2>(choice);

        std::stringstream ss, ss_long;
        ss << std::boolalpha << val;
        ss_long << std::boolalpha << std::setw(8) << std::left << val << " # (\""
                << str_val << "\") " << descript;
        m_mutex.lock();
        if(m_env.find(env_id) != m_env.end())
        {
            for(const auto& itr : m_env)
                if(itr.first == env_id && itr.second == ss.str())
                {
                    m_mutex.unlock();
                    return;
                }
        }
        m_env.insert(env_pair_t(env_id, ss_long.str()));
        m_mutex.unlock();
    }

    const env_map_t& get() const { return m_env; }
    mutex_t&         mutex() const { return m_mutex; }

    friend std::ostream& operator<<(std::ostream& os, const EnvSettings& env)
    {
        std::stringstream filler;
        filler.fill('#');
        filler << std::setw(90) << "";
        std::stringstream ss;
        ss << filler.str() << "\n# Environment settings:\n";
        env.mutex().lock();
        for(const auto& itr : env.get())
        {
            ss << "# " << std::setw(35) << std::right << itr.first << "\t = \t"
               << std::left << itr.second << "\n";
        }
        env.mutex().unlock();
        ss << filler.str();
        os << ss.str() << std::endl;
        return os;
    }

private:
    env_map_t       m_env;
    mutable mutex_t m_mutex;
};

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
        // record value defined by environment
        EnvSettings::GetInstance()->insert<Tp>(env_id, var);
        return var;
    }
    // record default value
    EnvSettings::GetInstance()->insert<Tp>(env_id, _default);

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
                itr = (char)std::tolower(itr);
            if(var == "off" || var == "false")
                val = false;
        }
        // record value defined by environment
        EnvSettings::GetInstance()->insert<bool>(env_id, val);
        return val;
    }
    // record default value
    EnvSettings::GetInstance()->insert<bool>(env_id, false);

    // return default if not specified in environment
    return _default;
}

//--------------------------------------------------------------------------------------//
//  overload for GetEnv + message when set
//
template <typename Tp>
Tp
GetEnv(const std::string& env_id, Tp _default, const std::string& msg)
{
    char* env_var = std::getenv(env_id.c_str());
    if(env_var)
    {
        std::string        str_var = std::string(env_var);
        std::istringstream iss(str_var);
        Tp                 var = Tp();
        iss >> var;
        std::cout << "Environment variable \"" << env_id << "\" enabled with "
                  << "value == " << var << ". " << msg << std::endl;
        // record value defined by environment
        EnvSettings::GetInstance()->insert<Tp>(env_id, var);
        return var;
    }
    // record default value
    EnvSettings::GetInstance()->insert<Tp>(env_id, _default);

    // return default if not specified in environment
    return _default;
}

//--------------------------------------------------------------------------------------//
//  use this function to get an environment variable setting from set of choices
//
//      EnvChoiceList<int> choices =
//              { EnvChoice<int>(NN,     "NN",     "nearest neighbor interpolation"),
//                EnvChoice<int>(LINEAR, "LINEAR", "bilinear interpolation"),
//                EnvChoice<int>(CUBIC,  "CUBIC",  "bicubic interpolation") };
//
//      int eInterp = GetEnv<int>("INTERPOLATION", choices, CUBIC);
//
template <typename Tp>
Tp
GetEnv(const std::string& env_id, const EnvChoiceList<Tp>& _choices, Tp _default)
{
    auto asupper = [](std::string var) {
        for(auto& itr : var)
            itr = (char)std::toupper(itr);
        return var;
    };

    char* env_var = std::getenv(env_id.c_str());
    if(env_var)
    {
        std::string str_var = std::string(env_var);
        std::string upp_var = asupper(str_var);
        Tp          var     = Tp();
        // check to see if string matches a choice
        for(const auto& itr : _choices)
        {
            if(asupper(std::get<1>(itr)) == upp_var)
            {
                // record value defined by environment
                EnvSettings::GetInstance()->insert(env_id, itr);
                return std::get<0>(itr);
            }
        }
        std::istringstream iss(str_var);
        iss >> var;
        // check to see if string matches a choice
        for(const auto& itr : _choices)
        {
            if(var == std::get<0>(itr))
            {
                // record value defined by environment
                EnvSettings::GetInstance()->insert(env_id, itr);
                return var;
            }
        }
        // the value set in env did not match any choices
        std::stringstream ss;
        ss << "\n### Environment setting error @ " << __FUNCTION__ << " (line "
           << __LINE__ << ")! Invalid selection for \"" << env_id
           << "\". Valid choices are:\n";
        for(const auto& itr : _choices)
            ss << "\t\"" << std::get<0>(itr) << "\" or \"" << std::get<1>(itr) << "\" ("
               << std::get<2>(itr) << ")\n";
        std::cerr << ss.str() << std::endl;
        abort();
    }

    std::string _name = "???";
    std::string _desc = "description not provided";
    for(const auto& itr : _choices)
        if(std::get<0>(itr) == _default)
        {
            _name = std::get<1>(itr);
            _desc = std::get<2>(itr);
            break;
        }

    // record default value
    EnvSettings::GetInstance()->insert(env_id, EnvChoice<Tp>(_default, _name, _desc));

    // return default if not specified in environment
    return _default;
}

//--------------------------------------------------------------------------------------//

template <typename Tp>
Tp
GetChoice(const EnvChoiceList<Tp>& _choices, const std::string& str_var)
{
    auto asupper = [](std::string var) {
        for(auto& itr : var)
            itr = (char)std::toupper(itr);
        return var;
    };

    std::string upp_var = asupper(str_var);
    Tp          var     = Tp();
    // check to see if string matches a choice
    for(const auto& itr : _choices)
    {
        if(asupper(std::get<1>(itr)) == upp_var)
        {
            // record value defined by environment
            return std::get<0>(itr);
        }
    }
    std::istringstream iss(str_var);
    iss >> var;
    // check to see if string matches a choice
    for(const auto& itr : _choices)
    {
        if(var == std::get<0>(itr))
        {
            // record value defined by environment
            return var;
        }
    }
    // the value set in env did not match any choices
    std::stringstream ss;
    ss << "\n### Environment setting error @ " << __FUNCTION__ << " (line " << __LINE__
       << ")! Invalid selection \"" << str_var << "\". Valid choices are:\n";
    for(const auto& itr : _choices)
        ss << "\t\"" << std::get<0>(itr) << "\" or \"" << std::get<1>(itr) << "\" ("
           << std::get<2>(itr) << ")\n";
    std::cerr << ss.str() << std::endl;
    abort();
}

//--------------------------------------------------------------------------------------//

inline void
PrintEnv(std::ostream& os = std::cout)
{
    os << (*EnvSettings::GetInstance());
}

//--------------------------------------------------------------------------------------//

struct ScopeDestructor
{
    template <typename FuncT>
    ScopeDestructor(FuncT&& _func)
    : m_functor(std::forward<FuncT>(_func))
    {}

    // delete copy operations
    ScopeDestructor(const ScopeDestructor&) = delete;
    ScopeDestructor& operator=(const ScopeDestructor&) = delete;

    // allow move operations
    ScopeDestructor(ScopeDestructor&& rhs) noexcept
    : m_functor(std::move(rhs.m_functor))
    {
        rhs.m_functor = []() {};
    }

    ScopeDestructor& operator=(ScopeDestructor&& rhs) noexcept
    {
        if(this != &rhs)
        {
            m_functor     = std::move(rhs.m_functor);
            rhs.m_functor = []() {};
        }
        return *this;
    }

    ~ScopeDestructor() { m_functor(); }

private:
    std::function<void()> m_functor = []() {};
};

//--------------------------------------------------------------------------------------//

}  // namespace PTL
