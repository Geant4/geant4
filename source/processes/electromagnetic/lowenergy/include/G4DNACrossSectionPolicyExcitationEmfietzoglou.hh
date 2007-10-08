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
// $Id: G4DNACrossSectionPolicyExcitationEmfietzoglou.hh,v 1.1 2007-10-08 09:18:42 sincerti Exp $
// -------------------------------------------------------------------
//

#ifndef G4DNACrossSectionPolicyExcitationEmfietzoglou_HH
#define G4DNACrossSectionPolicyExcitationEmfietzoglou_HH 1

#include "G4DNACrossSectionDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DNACrossSectionPolicyExcitationEmfietzoglou 
{
  public:
   G4DNACrossSectionPolicyExcitationEmfietzoglou() {}
   ~G4DNACrossSectionPolicyExcitationEmfietzoglou() {}

   G4double CrossSection(const G4Track&);

   G4double EnergyConstant(G4int excitationLevelIndex);
   G4double PartialCrossSection(G4double T, G4int excitationLevelIndex);
   G4int RandomizePartialCrossSection(G4double k, G4int z, const G4ParticleDefinition* particle);
 
  private:
 
   G4String name;  
   G4double lowEnergyLimit;
   G4double highEnergyLimit;
   G4double zeroBelowLowEnergyLimit;
   G4double zeroAboveHighEnergyLimit;

   // Hides default constructor and assignment operator as private
   G4DNACrossSectionPolicyExcitationEmfietzoglou(const G4DNACrossSectionPolicyExcitationEmfietzoglou & copy);
   G4DNACrossSectionPolicyExcitationEmfietzoglou & operator=(const G4DNACrossSectionPolicyExcitationEmfietzoglou & right);
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4DNACrossSectionPolicyExcitationEmfietzoglou.icc"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
