#ifndef G4WATCHER_GUN_HH
#define G4WATCHER_GUN_HH

#include "G4NuclWatcher.hh"
#include "g4std/vector"

class G4WatcherGun {

public:

  G4WatcherGun();
  void setWatchers();

  G4std::vector<G4NuclWatcher> getWatchers() const { 
    return watchers; 
  };

private: 

  G4int verboseLevel;
  G4std::vector<G4NuclWatcher> watchers;

};        

#endif // G4WATCHER_GUN_HH 
