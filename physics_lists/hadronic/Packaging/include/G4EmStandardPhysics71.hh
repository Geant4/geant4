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
// $Id: G4EmStandardPhysics71.hh,v 1.1 2005-11-30 18:00:38 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysics71
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
//
//----------------------------------------------------------------------------
//

#ifndef G4EmStandardPhysics71_h
#define G4EmStandardPhysics71_h 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

// for lib list detection....
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmStandardPhysics71 : public G4VPhysicsConstructor
{
public:
  G4EmStandardPhysics71(const G4String& name = "EMstandard71", G4bool msc=true);
  virtual ~G4EmStandardPhysics71();

public:
  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4bool mscStepLimit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif






