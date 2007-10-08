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
// -------------------------------------------------------------------
// $Id: G4DNACrossSectionPolicyChargeDecreaseDingfelder.hh,v 1.1 2007-10-08 09:18:42 sincerti Exp $
// -------------------------------------------------------------------
//

#ifndef G4DNACrossSectionPolicyChargeDecreaseDingfelder_HH
#define G4DNACrossSectionPolicyChargeDecreaseDingfelder_HH 1

#include "G4DNACrossSectionDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DNACrossSectionPolicyChargeDecreaseDingfelder 
{
  protected:
   G4DNACrossSectionPolicyChargeDecreaseDingfelder() {}
   ~G4DNACrossSectionPolicyChargeDecreaseDingfelder() {}

   G4double CrossSection(const G4Track& track);
   G4int RandomizePartialCrossSection(G4double k);
   void BuildTotalCrossSection(void) const {}

  private:

   G4double                             PartialCrossSection(G4double k, G4int index);

   G4double lowEnergyLimit;
   G4double highEnergyLimit;
   G4double zeroBelowLowEnergyLimit;
   G4double zeroAboveHighEnergyLimit;
  
   G4int NumberOfPartialCrossSections;

   G4double f0[2];
   G4double a0[2];
   G4double a1[2];
   G4double b0[2];
   G4double b1[2];
   G4double c0[2];
   G4double d0[2];
   G4double x0[2];
   G4double x1[2];

   // Hides default constructor and assignment operator as private
   G4DNACrossSectionPolicyChargeDecreaseDingfelder(const G4DNACrossSectionPolicyChargeDecreaseDingfelder & copy);
   G4DNACrossSectionPolicyChargeDecreaseDingfelder & operator=(const G4DNACrossSectionPolicyChargeDecreaseDingfelder & right);
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4DNACrossSectionPolicyChargeDecreaseDingfelder.icc"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
