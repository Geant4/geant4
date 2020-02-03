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
/// \file exoticphysics/monopole/include/G4Monopole.hh
/// \brief Definition of the G4Monopole class
//
//

#ifndef G4Monopole_h
#define G4Monopole_h 1

//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
//
// ClassName:   G4Monopole
//
// Authors:   21.03.05  V.Ivanchenko 
//
// Modified:
//
//  12.07.10  S.Burdin (changed the magnetic and electric charge variables from integer to double)
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include "CLHEP/Units/SystemOfUnits.h"

// ######################################################################
// ###                       Monopole                                 ###
// ######################################################################

class G4Monopole : public G4ParticleDefinition
{
private:

  static G4Monopole*  theMonopole;

  G4Monopole(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable );

  virtual ~G4Monopole();

public: 
  
  static G4Monopole* MonopoleDefinition(G4double mass = 100.*CLHEP::GeV, 
                                        G4double magCharge = 1.0, 
                                        G4double elCharge  = 0.0);

  static G4Monopole* Monopole();

  G4double MagneticCharge() const;

private:

  static G4double magCharge;
};

#endif
