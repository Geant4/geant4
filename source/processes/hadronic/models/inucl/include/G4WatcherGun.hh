#ifndef G4WATCHER_GUN_HH
#define G4WATCHER_GUN_HH

#include "G4NuclWatcher.hh"

#include "vector"

class G4WatcherGun {

public:

G4WatcherGun() {};

void setWatchers();

vector<G4NuclWatcher> getWatchers() const { return watchers; };

private: 

vector<G4NuclWatcher> watchers;

};        

#endif // G4WATCHER_GUN_HH 
