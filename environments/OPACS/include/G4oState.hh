//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4oState.hh,v 1.4.4.1 2001/06/28 19:06:33 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
#ifndef G4oState_h
#define G4oState_h 

#include <G4VStateDependent.hh>
#include <G4UIsession.hh>

class G4oState : public G4VStateDependent, public G4UIsession {
public:
  G4oState(G4String);
  ~G4oState();
  G4bool Notify(G4ApplicationState);
  void PauseSessionStart(G4String);
private:
  G4String name;
};

#endif

