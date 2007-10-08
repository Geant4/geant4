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
// $Id: G4DNACrossSectionPolicyScreenedRutherford.hh,v 1.1 2007-10-08 09:18:43 sincerti Exp $
// -------------------------------------------------------------------
//

#ifndef G4DNACrossSectionPolicyScreenedRutherford_HH
#define G4DNACrossSectionPolicyScreenedRutherford_HH 1

#include "G4DNACrossSectionDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DNACrossSectionPolicyScreenedRutherford 
{
  protected:
   G4DNACrossSectionPolicyScreenedRutherford() {}
   ~G4DNACrossSectionPolicyScreenedRutherford() {}

   G4double CrossSection(const G4Track&);

  private:
   
   G4double RutherfordCrossSection(G4double energy, G4double z);

   G4double ScreeningFactor(G4double energy, G4double z);

   G4double lowEnergyLimit;
   G4double highEnergyLimit;
   G4double zeroBelowLowEnergyLimit;
   G4double zeroAboveHighEnergyLimit;

   // Hides default constructor and assignment operator as private
   G4DNACrossSectionPolicyScreenedRutherford(const G4DNACrossSectionPolicyScreenedRutherford & copy);
   G4DNACrossSectionPolicyScreenedRutherford & operator=(const G4DNACrossSectionPolicyScreenedRutherford & right);
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4DNACrossSectionPolicyScreenedRutherford.icc"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
