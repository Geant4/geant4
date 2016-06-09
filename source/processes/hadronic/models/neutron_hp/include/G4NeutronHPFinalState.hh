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
//
// $Id: G4NeutronHPFinalState.hh,v 1.13 2007/06/06 12:45:13 ahoward Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
#ifndef G4NeutronHPFinalState_h
#define G4NeutronHPFinalState_h

#include "G4Material.hh"
#include "G4FastVector.hh"
#include "G4HadFinalState.hh"
#include "G4NeutronHPNames.hh"
#include "G4NeutronHPVector.hh"
#include "G4HadProjectile.hh"

class G4NeutronHPFinalState
{
public:

  G4NeutronHPFinalState()
  { 
    hasFSData = true; 
    hasXsec = true;
    hasAnyData = true;
    theBaseZ = 0;
    theBaseA = 0;
  };
  
  virtual ~G4NeutronHPFinalState(){};

  virtual void Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType) = 0;
  virtual G4HadFinalState * ApplyYourself(const G4HadProjectile & ) 
  {
    throw G4HadronicException(__FILE__, __LINE__, "G4HadFinalState * ApplyYourself(const G4HadProjectile & theTrack) needs implementation");
    return 0;
  }
  
  // of course this would better be Done templating G4NeutronHPChannel..., 
  // but due to peculiarities of the DEC compiler, this way it
  // is easier to maintain.
  virtual G4NeutronHPFinalState * New() = 0;
  
  G4bool HasXsec() {return hasXsec;}
  G4bool HasFSData() {return hasFSData;}
  G4bool HasAnyData() {return hasAnyData;}
  
  virtual G4double GetXsec(G4double ) { return 0; }
  virtual G4NeutronHPVector * GetXsec() { return 0; }
  
  void     SetA_Z(G4double anA, G4double aZ) {theBaseA = anA; theBaseZ = aZ; }
  G4double GetZ() { return theBaseZ; }
  G4double GetN() { return theBaseA; }
  
  protected:
  
  G4bool hasXsec;
  G4bool hasFSData;
  G4bool hasAnyData;
  G4NeutronHPNames theNames;
  
  G4HadFinalState theResult;
  
  G4double theBaseA;
  G4double theBaseZ;

  private:
  
};
#endif
