// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVAVisManager.hh,v 1.3 2001-02-07 17:31:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#ifndef TstVAVisManager_h
#define TstVAVisManager_h 1

#include "G4VisManager.hh"

class TstVAVisManager: public G4VisManager
{

public:

  TstVAVisManager ();

private:

  void RegisterGraphicsSystems ();

};

#endif
