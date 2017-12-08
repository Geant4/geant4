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
// $Id: G4EmStandardPhysics_option2.hh 105735 2017-08-16 12:59:43Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmStandardPhysics_option2
//
// Author:      V.Ivanchenko 09.11.2005
//
// Modified:
// 05.12.2005 V.Ivanchenko add controlled verbosity
// 19.12.2005 V.Ivanchenko rename 71 -> 72
// 15.05.2007 V.Ivanchenko rename to _option2
//
//----------------------------------------------------------------------------
//
// This class provides construction of EM standard physics using set of options
// allowing to utilize sub-cutoff option for ionisation processes and higher
// production threshold than in default EM physics.
//

#ifndef G4EmStandardPhysics_option2_h
#define G4EmStandardPhysics_option2_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4EmParticleList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmStandardPhysics_option2 : public G4VPhysicsConstructor
{
public:

  explicit G4EmStandardPhysics_option2(G4int ver=1, const G4String& name="");

  virtual ~G4EmStandardPhysics_option2();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  G4int  verbose;
  G4EmParticleList partList;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif






