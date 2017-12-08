//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4NeutronTrackingCut.hh 107562 2017-11-22 15:39:24Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronTrackingCut
//
// Author: Nov 2006 G.Folger
//
//
//----------------------------------------------------------------------------
//

#ifndef G4NeutronTrackingCut_h
#define G4NeutronTrackingCut_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

class G4NeutronKiller;

class G4NeutronTrackingCut : public G4VPhysicsConstructor
{
public: 
  G4NeutronTrackingCut(G4int ver=0);
  G4NeutronTrackingCut(const G4String& name,G4int ver=0);
  virtual ~G4NeutronTrackingCut();

    // This method will be invoked in the Construct() method. 
    // each particle type will be instantiated
  virtual void ConstructParticle();
 
    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type 
  virtual void ConstructProcess();

  inline void SetTimeLimit(G4double);
  inline void SetKineticEnergyLimit(G4double);

private:

  G4double timeLimit;
  G4double kineticEnergyLimit;
  
  G4int    verbose;
};

inline void G4NeutronTrackingCut::SetTimeLimit(G4double val)
{
  timeLimit = val;
}

inline void G4NeutronTrackingCut::SetKineticEnergyLimit(G4double val)
{
  kineticEnergyLimit = val;
}

#endif








