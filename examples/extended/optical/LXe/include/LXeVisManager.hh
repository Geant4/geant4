
#ifndef LXeVisManager_h
#define LXeVisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

class LXeVisManager: public G4VisManager 
{
public:
  
  LXeVisManager ();

private:

  void RegisterGraphicsSystems ();
};

#endif

#endif
