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
/// \file electromagnetic/TestEm17/include/MuCrossSections.hh
/// \brief Definition of the MuCrossSections class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MuCrossSections_h
#define MuCrossSections_h 1

#include "globals.hh"

class G4Material;
class G4Element;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MuCrossSections
{
  public:
    MuCrossSections();
   ~MuCrossSections();

  public:
    G4double CR_Macroscopic (const G4String&, const G4Material*, 
                             G4double, G4double);   
    G4double CR_PerAtom     (const G4String&, const G4Element*, 
                             G4double, G4double);
                       
  private:
    G4double CRB_Mephi (G4double, G4double, G4double, G4double);
    G4double CRK_Mephi (G4double, G4double, G4double, G4double);
    G4double CRN_Mephi (G4double, G4double, G4double, G4double);
    G4double CRP_Mephi (G4double, G4double, G4double, G4double);    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
