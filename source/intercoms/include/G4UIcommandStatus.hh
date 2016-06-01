// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcommandStatus.hh,v 2.0 1998/07/02 17:08:01 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4UIcommandStatus_h
#define G4UIcommandStatus_h 1

enum G4UIcommandStatus
{
  fCommandSucceeded         =   0,
  fCommandNotFound          = 100,
  fIllegalApplicationState  = 200,
  fParameterOutOfRange      = 300,
  fParameterUnreadable      = 400,
  fParameterOutOfCandidates = 500
};

#endif
