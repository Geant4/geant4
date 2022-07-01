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
// G4Profiler class implementation
//
// Author: J.Madsen (LBL), Aug 2020
// --------------------------------------------------------------------

#include "G4Profiler.hh"
#include "G4TiMemory.hh"

#if defined(GEANT4_USE_TIMEMORY)
#  include <timemory/runtime/configure.hpp>
#endif

#include <regex>
#include <thread>

//---------------------------------------------------------------------------//

typename G4Profiler::array_type& G4Profiler::GetEnabled()
{
  static array_type _instance = []() {
    array_type _tmp{};
    _tmp.fill(false);
    // environment settings introduce the defaults
    // but will be overridden by messenger and/or command line
    _tmp.at(G4ProfileType::Run)   = G4GetEnv<bool>("G4PROFILE_RUN", true);
    _tmp.at(G4ProfileType::Event) = G4GetEnv<bool>("G4PROFILE_EVENT", false);
    _tmp.at(G4ProfileType::Track) = G4GetEnv<bool>("G4PROFILE_TRACK", false);
    _tmp.at(G4ProfileType::Step)  = G4GetEnv<bool>("G4PROFILE_STEP", false);
    _tmp.at(G4ProfileType::User)  = G4GetEnv<bool>("G4PROFILE_USER", true);
    return _tmp;
  }();
  return _instance;
}

//---------------------------------------------------------------------------//

void G4Profiler::Configure(int argc, char** argv)
{
#if defined(GEANT4_USE_TIMEMORY)
  tim::timemory_init(argc, argv);
  std::vector<std::string> _args;
  for(int i = 0; i < argc; ++i)
    _args.push_back(argv[i]);
  Configure(_args);
#else
  G4Impl::consume_parameters(argc, argv);
#endif
}
//---------------------------------------------------------------------------//

void G4Profiler::Configure(ArgumentParser& parser, int argc, char** argv)
{
#if defined(GEANT4_USE_TIMEMORY)
  tim::timemory_init(argc, argv);
  std::vector<std::string> _args;
  for(int i = 0; i < argc; ++i)
    _args.push_back(argv[i]);
  Configure(parser, _args);
#else
  G4Impl::consume_parameters(parser, argc, argv);
#endif
}

//---------------------------------------------------------------------------//

void G4Profiler::Configure(const std::vector<std::string>& args)
{
#if defined(GEANT4_USE_TIMEMORY)
  if(args.empty())
    return;
  ArgumentParser p{ args.at(0) };
  Configure(p, args);
#else
  G4Impl::consume_parameters(args);
#endif
}

//---------------------------------------------------------------------------//

void G4Profiler::Configure(ArgumentParser& parser,
                           const std::vector<std::string>& args)
{
#if defined(GEANT4_USE_TIMEMORY)
  using parser_t     = ArgumentParser;
  using parser_err_t = typename parser_t::result_type;

  if(args.empty())
    return;

  static std::mutex mtx;
  std::unique_lock<std::mutex> lk(mtx);

  static auto tid = std::this_thread::get_id();
  if(std::this_thread::get_id() != tid)
    return;

  // std::cout << "Arguments:";
  // for(auto itr : args)
  //  std::cout << " " << itr;
  // std::cout << std::endl;

  auto help_action = [](parser_t& p) {
    p.print_help();
    exit(EXIT_FAILURE);
  };

  parser.enable_help();
  parser.on_error([=](parser_t& p, parser_err_t _err) {
    std::cerr << _err << std::endl;
    help_action(p);
  });

  auto get_bool = [](const std::string& _str, bool _default) {
    if(_str.empty())
      return _default;
    using namespace std::regex_constants;
    const auto regex_config = egrep | icase;
    if(std::regex_match(_str, std::regex("on|true|yes|1", regex_config)))
      return true;
    else if(std::regex_match(_str, std::regex("off|false|no|0", regex_config)))
      return false;
    return _default;
  };
  //
  //  Collection control
  //
  parser.add_argument()
    .names({ "-p", "--profile" })
    .description("Profiler modes")
    .choices({ "run", "event", "track", "step", "user" })
    .action([&](parser_t& p) {
      using namespace std::regex_constants;
      const auto regex_config = egrep | icase;
      for(auto&& itr : p.get<std::vector<std::string>>("profile"))
      {
        if(std::regex_match(itr, std::regex("run", regex_config)))
          G4Profiler::SetEnabled(G4ProfileType::Run, true);
        else if(std::regex_match(itr, std::regex("event", regex_config)))
          G4Profiler::SetEnabled(G4ProfileType::Event, true);
        else if(std::regex_match(itr, std::regex("track", regex_config)))
          G4Profiler::SetEnabled(G4ProfileType::Track, true);
        else if(std::regex_match(itr, std::regex("step", regex_config)))
          G4Profiler::SetEnabled(G4ProfileType::Step, true);
        else if(std::regex_match(itr, std::regex("user", regex_config)))
          G4Profiler::SetEnabled(G4ProfileType::User, true);
      }
    });
  //
  //  Component controls
  //
  parser.add_argument()
    .names({ "-r", "--run-components" })
    .description("Components for run profiling (see 'timemory-avail -s')")
    .action([&](parser_t& p) {
      G4RunProfiler::reset();
      tim::configure<G4RunProfiler>(
        p.get<std::vector<std::string>>("run-components"));
    });
  parser.add_argument()
    .names({ "-e", "--event-components" })
    .description("Components for event profiling (see 'timemory-avail -s')")
    .action([&](parser_t& p) {
      G4EventProfiler::reset();
      tim::configure<G4EventProfiler>(
        p.get<std::vector<std::string>>("event-components"));
    });
  parser.add_argument()
    .names({ "-t", "--track-components" })
    .description("Components for track profiling (see 'timemory-avail -s')")
    .action([&](parser_t& p) {
      G4TrackProfiler::reset();
      tim::configure<G4TrackProfiler>(
        p.get<std::vector<std::string>>("track-components"));
    });
  parser.add_argument()
    .names({ "-s", "--step-components" })
    .description("Components for step profiling (see 'timemory-avail -s')")
    .action([&](parser_t& p) {
      G4StepProfiler::reset();
      tim::configure<G4StepProfiler>(
        p.get<std::vector<std::string>>("step-components"));
    });
  parser.add_argument()
    .names({ "-u", "--user-components" })
    .description("Components for user profiling (see 'timemory-avail -s')")
    .action([&](parser_t& p) {
      G4StepProfiler::reset();
      tim::configure<G4UserProfiler>(
        p.get<std::vector<std::string>>("user-components"));
    });
  //
  //  Display controls
  //
  parser.add_argument()
    .names({ "-H", "--hierarchy", "--tree" })
    .description("Display the results as a call-stack hierarchy.")
    .max_count(1)
    .action([&](parser_t& p) {
      tim::settings::flat_profile() =
        !get_bool(p.get<std::string>("tree"), true);
    });
  parser.add_argument()
    .names({ "-F", "--flat" })
    .description("Display the results as a flat call-stack")
    .max_count(1)
    .action([&](parser_t& p) {
      tim::settings::flat_profile() =
        get_bool(p.get<std::string>("flat"), true);
    });
  parser.add_argument()
    .names({ "-T", "--timeline" })
    .description(
      "Do not merge duplicate entries at the same call-stack position")
    .max_count(1)
    .action([&](parser_t& p) {
      tim::settings::timeline_profile() =
        get_bool(p.get<std::string>("timeline"), true);
    });
  parser.add_argument()
    .names({ "--per-thread" })
    .description(
      "Display the results for each individual thread (default: aggregation)")
    .max_count(1)
    .action([&](parser_t& p) {
      tim::settings::flat_profile() =
        get_bool(p.get<std::string>("per-thread"), true);
    });
  parser.add_argument()
    .names({ "--per-event" })
    .description("Each G4Event is a unique entry")
    .max_count(1)
    .action([&](parser_t& p) {
      tim::settings::timeline_profile() =
        get_bool(p.get<std::string>("per-event"), true);
    });
  //
  //  Output controls
  //
  parser.add_argument()
    .names({ "-D", "--dart" })
    .description("Enable Dart output (CTest/CDash data tracking)")
    .max_count(1)
    .action([&](parser_t& p) {
      tim::settings::dart_output() = get_bool(p.get<std::string>("dart"), true);
    });
  parser.add_argument()
    .names({ "-J", "--json" })
    .description("Enable JSON output")
    .max_count(1)
    .action([&](parser_t& p) {
      tim::settings::json_output() = get_bool(p.get<std::string>("json"), true);
    });
  parser.add_argument()
    .names({ "-P", "--plot" })
    .description("Plot the JSON output")
    .max_count(1)
    .action([&](parser_t& p) {
      tim::settings::plot_output() = get_bool(p.get<std::string>("plot"), true);
    });
  parser.add_argument()
    .names({ "-X", "--text" })
    .description("Enable TEXT output")
    .max_count(1)
    .action([&](parser_t& p) {
      tim::settings::text_output() = get_bool(p.get<std::string>("text"), true);
    });
  parser.add_argument()
    .names({ "-C", "--cout" })
    .description("Enable output to terminal")
    .max_count(1)
    .action([&](parser_t& p) {
      tim::settings::cout_output() = get_bool(p.get<std::string>("cout"), true);
    });
  parser.add_argument()
    .names({ "-O", "--output-path" })
    .description("Set the output directory")
    .count(1)
    .action([&](parser_t& p) {
      tim::settings::output_path() = p.get<std::string>("output-path");
    });
  parser.add_argument()
    .names({ "-W", "--hw-counters" })
    .description(
      "Set the hardware counters to collect (see 'timemory-avail -H')")
    .action([&](parser_t& p) {
      tim::settings::papi_events() = p.get<std::string>("hw-counters");
    });

  tim::settings::time_output() = true;
  tim::settings::time_format() = "%F_%I.%M_%p";

  auto err = parser.parse(args);
  if(err)
  {
    std::cout << "Error! " << err << std::endl;
    help_action(parser);
  }

#else
  G4Impl::consume_parameters(parser, args);
#endif
}

//---------------------------------------------------------------------------//

void G4Profiler::Finalize()
{
#if defined(GEANT4_USE_TIMEMORY)
  // generally safe to call multiple times but just in case
  // it is not necessary to call from worker threads, calling
  // once on master will merge any threads that accumulated
  // results. If a thread is destroyed, before finalize is
  // called on master thread, then it will merge itself into
  // the master thread before terminating. The biggest issue
  // is when finalize is never called on any threads (including
  // the master) so finalization happens after main exits. When
  // this happens, the order that data gets deleted is very
  // hard to predict so it commonly results in a segfault. Thus,
  // if the user does not call this function and does not
  // delete G4RunManager, a segfault is quite likely.
  static thread_local bool _once = false;
  if(!_once)
    tim::timemory_finalize();
  _once = true;
#endif
}

//---------------------------------------------------------------------------//

template <size_t Ct>
template <int Idx>
typename G4ProfilerConfig<Ct>::template PersistentSettings<Idx>&
G4ProfilerConfig<Ct>::GetPersistentFallback()
{
  // G4RunManager::ConfigureProfilers() assigns defaults to these lambdas
  // based on environment variables. The users should (generally) not modify
  // these.
  static PersistentSettings<Idx> _instance = PersistentSettings<Idx>{};
  return _instance;
}

//---------------------------------------------------------------------------//

template <size_t Ct>
template <int Idx>
typename G4ProfilerConfig<Ct>::template PersistentSettings<Idx>&
G4ProfilerConfig<Ct>::GetPersistent()
{
  //
  //  The pointers here automatically initialize on the first
  //  invocation on a thread and a reference is returned so they
  //  can be cleanly "leaked" at the application termination. Assignment
  //  is thread-safe and this scheme avoids having to deal with getters
  //  and setters. It is designed such that we can set defaults but
  //  the user can override them at any time.
  //
  // the first global instance copies from the fallback rountines, which are
  // assigned defaults in G4RunManager::ConfigureProfilers() but if the user
  // assigns new settings before creating G4RunManager, those defaults will
  // be ignored
  static auto* _instance =
    new PersistentSettings<Idx>(GetPersistentFallback<Idx>());
  static thread_local PersistentSettings<Idx>* _tlinstance = [=]() {
    static std::mutex mtx;
    std::unique_lock<std::mutex> lk(mtx);
    // If this is the primary thread, just return the global instance from
    // above. Modifying that instance will result in all threads created later
    // to inherit those settings
    static bool _first = true;
    if(_first)
    {  // below uses comma operator to assign to false before return
      return ((_first = false), _instance);
    }
    // if not first, make a copy from the primary thread
    return new PersistentSettings<Idx>(*_instance);
  }();
  return *_tlinstance;
}

// ----------------------------------------------------------------------

template <size_t Cat>
G4ProfilerConfig<Cat>::~G4ProfilerConfig()
{
  delete m_bundle;
}

//----------------------------------------------------------------------------//
// primary (utilized) versions
//
template <size_t Cat>
typename G4ProfilerConfig<Cat>::QueryHandler_t
G4ProfilerConfig<Cat>::GetQueryFunctor()
{
  return QueryHandler_t{ GetPersistent<G4ProfileOp::Query>().m_functor };
}

//----------------------------------------------------------------------------//

template <size_t Cat>
typename G4ProfilerConfig<Cat>::LabelHandler_t
G4ProfilerConfig<Cat>::GetLabelFunctor()
{
  return LabelHandler_t{ GetPersistent<G4ProfileOp::Label>().m_functor };
}

//----------------------------------------------------------------------------//

template <size_t Cat>
typename G4ProfilerConfig<Cat>::ToolHandler_t
G4ProfilerConfig<Cat>::GetToolFunctor()
{
  return ToolHandler_t{ GetPersistent<G4ProfileOp::Tool>().m_functor };
}

//----------------------------------------------------------------------------//
// fallback versions
//
template <size_t Cat>
typename G4ProfilerConfig<Cat>::QueryHandler_t
G4ProfilerConfig<Cat>::GetFallbackQueryFunctor()
{
  return QueryHandler_t{
    GetPersistentFallback<G4ProfileOp::Query>().m_functor
  };
}

//----------------------------------------------------------------------------//

template <size_t Cat>
typename G4ProfilerConfig<Cat>::LabelHandler_t
G4ProfilerConfig<Cat>::GetFallbackLabelFunctor()
{
  return LabelHandler_t{
    GetPersistentFallback<G4ProfileOp::Label>().m_functor
  };
}

//----------------------------------------------------------------------------//

template <size_t Cat>
typename G4ProfilerConfig<Cat>::ToolHandler_t
G4ProfilerConfig<Cat>::GetFallbackToolFunctor()
{
  return ToolHandler_t{ GetPersistentFallback<G4ProfileOp::Tool>().m_functor };
}

//---------------------------------------------------------------------------//

class G4Run;
class G4Event;
class G4Track;
class G4Step;

using G4ProfType      = G4ProfileType;
using RunProfConfig   = G4ProfilerConfig<G4ProfType::Run>;
using EventProfConfig = G4ProfilerConfig<G4ProfType::Event>;
using TrackProfConfig = G4ProfilerConfig<G4ProfType::Track>;
using StepProfConfig  = G4ProfilerConfig<G4ProfType::Step>;
using UserProfConfig  = G4ProfilerConfig<G4ProfType::User>;

//---------------------------------------------------------------------------//

// force instantiations
template class G4ProfilerConfig<G4ProfType::Run>;
template class G4ProfilerConfig<G4ProfType::Event>;
template class G4ProfilerConfig<G4ProfType::Track>;
template class G4ProfilerConfig<G4ProfType::Step>;
template class G4ProfilerConfig<G4ProfType::User>;
template G4ProfilerConfig<G4ProfType::Run>::G4ProfilerConfig(const G4Run*);
template G4ProfilerConfig<G4ProfType::Event>::G4ProfilerConfig(const G4Event*);
template G4ProfilerConfig<G4ProfType::Track>::G4ProfilerConfig(const G4Track*);
template G4ProfilerConfig<G4ProfType::Step>::G4ProfilerConfig(const G4Step*);
template G4ProfilerConfig<G4ProfileType::User>::G4ProfilerConfig(
  const std::string&);

#define G4PROFILE_INSTANTIATION(TYPE)                                          \
  template TYPE::PersistentSettings<G4ProfileOp::Query>&                       \
  TYPE::GetPersistentFallback<G4ProfileOp::Query>();                           \
  template TYPE::PersistentSettings<G4ProfileOp::Label>&                       \
  TYPE::GetPersistentFallback<G4ProfileOp::Label>();                           \
  template TYPE::PersistentSettings<G4ProfileOp::Tool>&                        \
  TYPE::GetPersistentFallback<G4ProfileOp::Tool>();                            \
                                                                               \
  template TYPE::PersistentSettings<G4ProfileOp::Query>&                       \
  TYPE::GetPersistent<G4ProfileOp::Query>();                                   \
  template TYPE::PersistentSettings<G4ProfileOp::Label>&                       \
  TYPE::GetPersistent<G4ProfileOp::Label>();                                   \
  template TYPE::PersistentSettings<G4ProfileOp::Tool>&                        \
  TYPE::GetPersistent<G4ProfileOp::Tool>();

G4PROFILE_INSTANTIATION(RunProfConfig)
G4PROFILE_INSTANTIATION(EventProfConfig)
G4PROFILE_INSTANTIATION(TrackProfConfig)
G4PROFILE_INSTANTIATION(StepProfConfig)
G4PROFILE_INSTANTIATION(UserProfConfig)

//---------------------------------------------------------------------------//

#if !defined(GEANT4_USE_TIMEMORY)
#  define TIMEMORY_WEAK_PREFIX
#  define TIMEMORY_WEAK_POSTFIX
#endif

//---------------------------------------------------------------------------//

extern "C"
{
  // this allows the default setup to be overridden by linking
  // in an custom extern C function into the application
  TIMEMORY_WEAK_PREFIX
  void G4ProfilerInit(void) TIMEMORY_WEAK_POSTFIX;

  // this gets executed when the library gets loaded
  void G4ProfilerInit(void)
  {
#ifdef GEANT4_USE_TIMEMORY

    // guard against re-initialization
    static bool _once = false;
    if(_once)
      return;
    _once = true;

    puts(">>> G4ProfilerInit <<<");

    //
    // the default settings
    //
    // large profiles can take a very long time to plot
    tim::settings::plot_output() = false;
    // large profiles can take quite a bit of console space
    tim::settings::cout_output() = false;
    // this creates a subdirectory with the timestamp of the run
    tim::settings::time_output() = false;
    // see `man 3 strftime` for formatting keys
    tim::settings::time_format() = "%F_%I.%M_%p";
    // set the default precision for timing
    tim::settings::timing_precision() = 6;
    // set the minimum width for outputs
    tim::settings::width() = 12;
    // set dart reports (when enabled) to only print the first entry
    tim::settings::dart_count() = 1;
    // set dart reports (when enabled) to use the component label
    // instead of the string identifer of the entry, e.g.
    // >>> G4Run/0 ... peak_rss ... 50 MB would report
    // 'peak_rss 50 MB' not 'G4Run/0 50 MB'
    tim::settings::dart_label() = true;

    // allow environment overrides of the defaults
    tim::settings::parse();
#endif
  }
}  // extern "C"

//---------------------------------------------------------------------------//

#ifdef GEANT4_USE_TIMEMORY
namespace
{
  static bool profiler_is_initialized = (G4ProfilerInit(), true);
}
#endif

//---------------------------------------------------------------------------//
