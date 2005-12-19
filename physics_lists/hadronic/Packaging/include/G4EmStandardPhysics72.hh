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
// $Id: G4EmStandardPhysics72.hh,v 1.1 2005-12-19 12:18:40 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  G4EmStandardPhysics72(const G4String& name = "EMstandard72", G4int ver = 1,
			G4bool msc=true);
  virtual ~G4EmStandardPhysics72();

public:
  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4int  verbose;
  G4bool mscStepLimit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif






