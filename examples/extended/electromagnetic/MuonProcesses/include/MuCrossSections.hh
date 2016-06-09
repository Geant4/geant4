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
// $Id: MuCrossSections.hh,v 1.2 2004/08/27 11:06:23 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
    G4double CR_Macroscopic (const G4String&, G4Material*, G4double, G4double);   
    G4double CR_PerAtom     (const G4String&, G4Element* , G4double, G4double);
                       
  private:
    G4double CRB_Mephi (G4double, G4double, G4double, G4double);
    G4double CRK_Mephi (G4double, G4double, G4double, G4double);
    G4double CRN_Mephi (G4double, G4double, G4double, G4double);
    G4double CRP_Mephi (G4double, G4double, G4double, G4double);    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
