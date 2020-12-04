//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4Profiler
//
// Class description:
//
// Class providing the internal profiling interface for Geant4.

// Author: Jonathan Madsen, LBNL - November 2020
// --------------------------------------------------------------------
#ifndef G4Profiler_hh
#define G4Profiler_hh 1

// Fundamental definitions
#ifndef G4GMAKE
#  include "G4GlobalConfig.hh"
#endif

// for meta-programming stuff
#include "G4Profiler.icc"

#if defined(GEANT4_USE_TIMEMORY)
#  include <timemory/utility/argparse.hpp>
#endif

#include "globals.hh"

#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <type_traits>
#include <tuple>
#include <vector>
#include <array>

//----------------------------------------------------------------------------//

class G4Run;
class G4Event;
class G4Track;
class G4Step;

struct G4ProfileType
{
  enum : size_t
  {
    Run = 0,
    Event,
    Track,
    Step,
    User,
    TypeEnd
  };
};

struct G4ProfileOp
{
  enum : int
  {
    Query = 0,
    Label,
    Tool
  };
};

//----------------------------------------------------------------------------//

class G4Profiler
{
 public:
  using array_type = std::array<bool, G4ProfileType::TypeEnd>;

#if defined(GEANT4_USE_TIMEMORY)
  using ArgumentParser = tim::argparse::argument_parser;
#else
  struct ArgumentParser
  {
    explicit ArgumentParser(std::string) {}
  };
#endif

  static void Configure(const std::vector<std::string>& args);
  static void Configure(int argc, char** argv);
  static void Configure(ArgumentParser&, const std::vector<std::string>& args);
  static void Configure(ArgumentParser&, int argc, char** argv);
  static void Finalize();

  static bool GetEnabled(size_t v) { return GetEnabled().at(v); }
  static void SetEnabled(size_t v, bool val) { GetEnabled().at(v) = val; }

  static bool GetPerEvent() { return GetPerEventImpl(); }
  static void SetPerEvent(bool val) { GetPerEventImpl() = val; }

 private:
  static array_type& GetEnabled();
  static bool& GetPerEventImpl()
  {
    static bool _value = false;
    return _value;
  }
};

//----------------------------------------------------------------------------//
//  maps enumerations to types
//
template <size_t Category>
struct G4ProfilerObject
{
  using type = void;
};

template <size_t Category>
using G4ProfilerObject_t = typename G4ProfilerObject<Category>::type;

template <>
struct G4ProfilerObject<G4ProfileType::Run>
{
  using type = const G4Run*;
};

template <>
struct G4ProfilerObject<G4ProfileType::Event>
{
  using type = const G4Event*;
};

template <>
struct G4ProfilerObject<G4ProfileType::Track>
{
  using type = const G4Track*;
};

template <>
struct G4ProfilerObject<G4ProfileType::Step>
{
  using type = const G4Step*;
};

template <>
struct G4ProfilerObject<G4ProfileType::User>
{
  using type = const std::string&;
};

//----------------------------------------------------------------------------//
// default set of profiler args
template <size_t Category>
struct G4ProfilerArgs
{
  // this resolves to the G4ProfilerObject type-trait above, e.g.
  // "const G4Step*" when category is G4ProfileType::Step
  using value_type = G4ProfilerObject_t<Category>;

  // two-dimensional type-list where each inner type-list is a set
  // of arguments to support creating a profiler type from
  using type = G4TypeList<G4PROFILER_ARG_SET(value_type)>;
  // so above means there are functors in use which apply:
  //
  //  G4StepProfiler _profiler(const G4Step*);
  //
};

template <size_t Category>
using G4ProfilerArgs_t = typename G4ProfilerArgs<Category>::type;

//----------------------------------------------------------------------------//

template <size_t Category, typename RetT, typename CommonT = G4CommonTypeList<>>
struct G4ProfilerFunctors
{
  using type = G4Impl::Functors_t<RetT, CommonT, G4ProfilerArgs_t<Category>>;
};

template <size_t Category, typename RetT, typename... CommonT>
using G4ProfilerFunctors_t =
  typename G4ProfilerFunctors<Category, RetT,
                              G4CommonTypeList<CommonT...>>::type;

//----------------------------------------------------------------------------//

#ifdef GEANT4_USE_TIMEMORY

// Pre-declare the timemory component that will be used
namespace tim
{
  namespace component
  {
    template <size_t, typename Tag>
    struct user_bundle;
  }  // namespace component

  template <typename... Types>
  class auto_tuple;

  template <typename Tag, typename... Types>
  class auto_bundle;
}  // namespace tim

namespace g4tim
{
  using namespace tim;
  using tim::component::user_bundle;

  struct G4api : public tim::concepts::api
  {};

  using ProfilerArgparser = argparse::argument_parser;

}  // namespace g4tim

using G4RunProfiler   = g4tim::user_bundle<G4ProfileType::Run, G4ProfileType>;
using G4EventProfiler = g4tim::user_bundle<G4ProfileType::Event, G4ProfileType>;
using G4TrackProfiler = g4tim::user_bundle<G4ProfileType::Track, G4ProfileType>;
using G4StepProfiler  = g4tim::user_bundle<G4ProfileType::Step, G4ProfileType>;
using G4UserProfiler  = g4tim::user_bundle<G4ProfileType::User, G4ProfileType>;

template <typename... Types>
using G4ProfilerBundle = g4tim::auto_bundle<g4tim::G4api, Types...>;

#else

namespace g4tim
{
  struct ProfilerArgparser
  {};

  /// this provides a dummy wrapper for the profiling
  template <typename... Types>
  struct handler
  {
    template <typename... Args>
    handler(Args&&...)
    {}
    ~handler() = default;
    handler(const handler&) = default;
    handler(handler&&) = default;
    handler& operator=(const handler&) = default;
    handler& operator=(handler&&) = default;

    void record() {}
    template <typename... Args>
    void start(Args&&...)
    {}
    template <typename... Args>
    void stop(Args&&...)
    {}
    void push() {}
    void pop() {}
    void reset() {}
    void report_at_exit(bool) {}
    template <typename... Args>
    void mark_begin(Args&&...)
    {}
    template <typename... Args>
    void mark_end(Args&&...)
    {}
    friend std::ostream& operator<<(std::ostream& os, const handler&)
    {
      return os;
    }
  };

  template <size_t Idx, typename Tp>
  struct user_bundle
  {
    template <typename... Args>
    user_bundle(Args&&...)
    {}

    template <typename... Types, typename... Args>
    static void configure(Args&&...)
    {}

    static void reset() {}
  };

}  // namespace g4tim

using G4RunProfiler = g4tim::handler<>;
using G4EventProfiler = g4tim::handler<>;
using G4TrackProfiler = g4tim::handler<>;
using G4StepProfiler = g4tim::handler<>;
using G4UserProfiler = g4tim::handler<>;

template <typename... Types>
using G4ProfilerBundle = g4tim::handler<Types...>;

#endif

//----------------------------------------------------------------------------//

/// @brief G4ProfilerConfig
/// This class is used to determine whether to activate profiling in the code
///
/// @example extended/parallel/ThreadsafeScorers/ts_scorers.cc
/// The main for this file contains an example for configuring the G4Profiler
/// for G4Track and G4Step
///
template <size_t Category>
class G4ProfilerConfig
{
 public:
  using type = G4ProfilerBundle<g4tim::user_bundle<Category, G4ProfileType>>;
  using this_type = G4ProfilerConfig<Category>;

  static constexpr size_t arg_sets =
    G4TypeListSize<G4ProfilerArgs_t<Category>>::value;

  using QueryFunc_t = G4ProfilerFunctors_t<Category, bool>;
  using LabelFunc_t = G4ProfilerFunctors_t<Category, std::string>;
  using ToolFunc_t  = std::tuple<std::function<type*(const std::string&)>>;

 public:
  // when constructed with no args, should call operator()
  G4ProfilerConfig() = default;

  // constructor calls Query(...), Label(...), Tool(...)
  template <typename Arg, typename... Args>
  G4ProfilerConfig(Arg, Args...);

  // will delete m_bundle is allocated
  ~G4ProfilerConfig();

  // default the move-constructor and move-assignment
  G4ProfilerConfig(G4ProfilerConfig&&) = default;
  G4ProfilerConfig& operator=(G4ProfilerConfig&&) = default;

  // do not allow copy-construct and copy-assign bc of raw pointer
  G4ProfilerConfig(const G4ProfilerConfig&) = delete;
  G4ProfilerConfig& operator=(const G4ProfilerConfig&) = delete;

  // if constructed without args, this function should be called
  template <typename... Args>
  bool operator()(Args...);

  // provide nullptr check
  operator bool() const { return (m_bundle != nullptr); }

 private:
  type* m_bundle = nullptr;

  static QueryFunc_t& GetQueries()
  {
    static QueryFunc_t _value;
    return _value;
  }

  static LabelFunc_t& GetLables()
  {
    static LabelFunc_t _value;
    return _value;
  }

  static ToolFunc_t& GetTools()
  {
    static ToolFunc_t _value;
    return _value;
  }

 public:
  // invokes the functor that determines whether to enable profiling
  template <typename... Args>
  static bool Query(Args... _args);

  // invokes the functor for generating a Label when profiling is enabled
  template <typename... Args>
  static std::string Label(Args... _args);

  // invokes the functor for configuring a Tool instance
  template <typename... Args>
  static type* Tool(const std::string&);

  using QueryHandler_t = FuncHandler<this_type, QueryFunc_t, bool>;
  using LabelHandler_t = FuncHandler<this_type, LabelFunc_t, std::string>;
  using ToolHandler_t  = FuncHandler<this_type, ToolFunc_t, type*>;

  static QueryHandler_t GetQueryFunctor();
  static QueryHandler_t GetFallbackQueryFunctor();

  static LabelHandler_t GetLabelFunctor();
  static LabelHandler_t GetFallbackLabelFunctor();

  static ToolHandler_t GetToolFunctor();
  static ToolHandler_t GetFallbackToolFunctor();

 private:
  template <bool B, typename Lhs, typename Rhs>
  using conditional_t = typename std::conditional<B, Lhs, Rhs>::type;

  // this provides the global statics for the functors
  template <int Idx>
  struct PersistentSettings
  {
    // determine the functor type
    using functor_type = conditional_t<
      Idx == G4ProfileOp::Query, QueryFunc_t,
      conditional_t<Idx == G4ProfileOp::Label, LabelFunc_t, ToolFunc_t>>;

    PersistentSettings()  = default;
    ~PersistentSettings() = default;
    // default member initialization
    functor_type m_functor;
  };

  template <int Idx>
  static PersistentSettings<Idx>& GetPersistentFallback();

  template <int Idx>
  static PersistentSettings<Idx>& GetPersistent();
};

//----------------------------------------------------------------------------//

// alias for getting type
template <size_t Category>
using G4ProfilerConfig_t = typename G4ProfilerConfig<Category>::type;

// ----------------------------------------------------------------------

template <size_t Cat>
template <typename Arg, typename... Args>
G4ProfilerConfig<Cat>::G4ProfilerConfig(Arg _arg, Args... _args)
{
  this->operator()(_arg, _args...);
}

// ----------------------------------------------------------------------

template <size_t Cat>
template <typename... Args>
bool G4ProfilerConfig<Cat>::operator()(Args... _args)
{
  if(Query(_args...))
  {
    m_bundle = Tool(Label(_args...));
    if(m_bundle)
      m_bundle->start(_args...);
    return (m_bundle != nullptr);
  }
  return false;
}

// ----------------------------------------------------------------------

// invokes the functor that determines whether to enable profiling
template <size_t Cat>
template <typename... Args>
bool G4ProfilerConfig<Cat>::Query(Args... _args)
{
  return QueryHandler_t{ GetPersistent<G4ProfileOp::Query>().m_functor }(
    _args...);
}

//----------------------------------------------------------------------------//

// invokes the functor for generating a label when profiling is enabled
template <size_t Cat>
template <typename... Args>
std::string G4ProfilerConfig<Cat>::Label(Args... _args)
{
  return LabelHandler_t{ GetPersistent<G4ProfileOp::Label>().m_functor }(
    _args...);
}

//----------------------------------------------------------------------------//

// invokes the functor for configuring a tool instance
template <size_t Cat>
template <typename... Args>
typename G4ProfilerConfig<Cat>::type* G4ProfilerConfig<Cat>::Tool(
  const std::string& _args)
{
  return ToolHandler_t{ GetPersistent<G4ProfileOp::Tool>().m_functor }(_args);
}

//----------------------------------------------------------------------------//

// ensure that any implicit instantiations of the G4ProfilerConfig
// are not done in another translation units because we will
// explicitly instantiate G4ProfilerConfig in the .cc file

// line breaks make this much harder to read
// clang-format off
extern template class G4ProfilerConfig<G4ProfileType::Run>;
extern template class G4ProfilerConfig<G4ProfileType::Event>;
extern template class G4ProfilerConfig<G4ProfileType::Track>;
extern template class G4ProfilerConfig<G4ProfileType::Step>;
extern template class G4ProfilerConfig<G4ProfileType::User>;
extern template G4ProfilerConfig<G4ProfileType::Run>::G4ProfilerConfig(const G4Run*);
extern template G4ProfilerConfig<G4ProfileType::Event>::G4ProfilerConfig(const G4Event*);
extern template G4ProfilerConfig<G4ProfileType::Track>::G4ProfilerConfig(const G4Track*);
extern template G4ProfilerConfig<G4ProfileType::Step>::G4ProfilerConfig(const G4Step*);
extern template G4ProfilerConfig<G4ProfileType::User>::G4ProfilerConfig(const std::string&);
// clang-format on

//----------------------------------------------------------------------------//

#ifndef GEANT4_USE_TIMEMORY

#  include <ostream>
#  include <string>

#endif

#if defined(GEANT4_USE_TIMEMORY)
// two macros below create a unique variable name based on the line number
#  define G4USER_PROFILER_VAR_JOIN(X, Y) X##Y
#  define G4USER_PROFILER_VAR(Y) G4USER_PROFILER_VAR_JOIN(g4user_profiler_, Y)

// inserts just the string
#  define G4USER_SCOPED_PROFILE(...)                                           \
    G4ProfilerConfig<G4ProfileType::User> G4USER_PROFILER_VAR(__LINE__)(       \
      TIMEMORY_JOIN("", __VA_ARGS__))

// inserts the function
#  define G4USER_SCOPED_PROFILE_FUNC(...)                                      \
    G4ProfilerConfig<G4ProfileType::User> G4USER_PROFILER_VAR(__LINE__)(       \
      TIMEMORY_JOIN("", __FUNCTION__, "/", __VA_ARGS__))

// inserts the function and file
#  define G4USER_SCOPED_PROFILE_FUNC_FILE(...)                                 \
    G4ProfilerConfig<G4ProfileType::User> G4USER_PROFILER_VAR(__LINE__)(       \
      TIMEMORY_JOIN("", __FUNCTION__, '@', __FILE__, '/', __VA_ARGS__))

// inserts the function, file, and line number
#  define G4USER_SCOPED_PROFILE_FUNC_FILE_LINE(...)                            \
    G4ProfilerConfig<G4ProfileType::User> G4USER_PROFILER_VAR(__LINE__)(       \
      TIMEMORY_JOIN("", __FUNCTION__, '@', __FILE__, ':', __LINE__, '/',       \
                    __VA_ARGS__))
#else
#  define G4USER_SCOPED_PROFILE(...)
#  define G4USER_SCOPED_PROFILE_FUNC(...)
#  define G4USER_SCOPED_PROFILE_FUNC_FILE(...)
#  define G4USER_SCOPED_PROFILE_FUNC_FILE_LINE(...)
#endif

#endif  // G4Profiler_hh
