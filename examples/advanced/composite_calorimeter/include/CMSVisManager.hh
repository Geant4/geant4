// $Id: CMSVisManager.hh,v 1.1 2002-10-01 14:38:56 arce Exp $
// John Allison 24th January 1998.

// Example Visualization Manager implementing virtual function
//   RegisterGraphicsSystems.  Exploits C-pre-processor variables
//   G4VIS_USE_DAWN, etc., which are set by the GNUmakefiles if
//   environment variables of the same name are set.

// So all you have to do is set environment variables and compile and
//   instantiate this in your main().

// Alternatively, you can implement an empty function here and just
//   register the systems you want in your main(), e.g.:
//   G4VisManager* myVisManager = new CMSVisManager;
//   myVisManager -> RegisterGraphicsSystem (new MyGraphicsSystem);

#ifndef CMSEXAMPLEVISMANAGER_HH
#define CMSEXAMPLEVISMANAGER_HH

#ifdef G4VIS_USE
#include "G4VisManager.hh"

class CMSVisManager: public G4VisManager {

public:

  CMSVisManager (G4int verboseLevel = 0);
  // Controls initial verbose level of VisManager and VisMessenger.
  // Can be changed by /vis/set/verbose.

private:

  virtual void RegisterGraphicsSystems ();

};

#endif
#endif
