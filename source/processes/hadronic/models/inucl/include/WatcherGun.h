#ifndef WATCHER_GUN_H
#define WATCHER_GUN_H

#include "NuclWatcher.h"

#include "vector"

class WatcherGun {

public:

WatcherGun() {};

void setWatchers();

vector<NuclWatcher> getWatchers() const { return watchers; };

private: 
vector<NuclWatcher> watchers;
};        

#endif // WATCHER_GUN_H 
