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
// $Id: G4EmStandardPhysics72.hh,v 1.2 2006/11/23 15:30:19 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysics72
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 05.12.2005 V.Ivanchenko add controlled verbosity
// 19.12.2005 V.Ivanchenko rename 71 -> 72
//
//----------------------------------------------------------------------------
//

#ifndef G4EmStandardPhysics72_h
#define G4EmStandardPhysics72_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

// for lib list detection....
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4EmProcessOptions.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmStandardPhysics72 : public G4VPhysicsConstructor
{
public:
  G4EmStandardPhysics72(const G4String& name = "EMstandard72", G4int ver = 1);
  virtual ~G4EmStandardPhysics72();

public:
  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4int  verbose;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif






