// $Id: A01VisManager.hh,v 1.1 2002-11-13 07:19:07 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef A01VisManager_h
#define A01VisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

class A01VisManager: public G4VisManager
{
  public:
    A01VisManager ();

  private:
    void RegisterGraphicsSystems ();

};

#endif

#endif
