// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIxvtMessenger.hh,v 1.2 1999-04-13 01:25:01 yhajime Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// The original code (written by Makoto Asai) has been adapted for use      //
// with the XVT Geant GUI                                                   //
//                                                                          //
// By: Simon Prior 11/08/97                                                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////
#ifndef G4UIxvtMessenger_h
#define G4UIxvtMessenger_h 1

#if defined(G4UI_BUILD_XVT_SESSION) || defined(G4UI_USE_XVT)

#include "G4UImessenger.hh"
#include "G4UImanager.hh"

class G4UIxvt;

class G4UIxvtMessenger : public G4UImessenger 
{
  public:
      G4UIxvtMessenger(G4UIxvt*);
      void SetNewValue(G4UIcommand * command,G4String newValue);
  private:
      G4UIxvt* xvtptr;
};

#endif

#endif

