#ifndef STATE_PROCESS_MAP_H
#define STATE_PROCESS_MAP_H 1

#include <string>
#include <map>

#include "RCObjectHandle.hh"
#include "SAXStateProcess.hh"

class StateProcessMap
{
public:
  typedef RCObjectHandle<SAXStateProcess>                 Process;
  typedef std::map<std::string,Process>                   ProcessMap;

public:
  StateProcessMap();
  ~StateProcessMap();
  
  void Initialize();
  void Reset();
  
  void AddProcess( const std::string& tag, Process process );
  void AddProcess( char*        tag, Process process );
  bool Check( const std::string& tag );
  bool Check( char*        tag );
  Process GetProcess( const std::string& tag );
  Process GetProcess( char*        tag );

private:
  ProcessMap fMap;
};

#endif // STATE_PROCESS_MAP_H

