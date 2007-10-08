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
// $Id: G4DNACrossSectionPolicyIonizationBorn.hh,v 1.1 2007-10-08 09:18:42 sincerti Exp $
// -------------------------------------------------------------------
//

#ifndef G4DNACrossSectionPolicyIonizationBorn_HH
#define G4DNACrossSectionPolicyIonizationBorn_HH 1

#include "G4DNACrossSectionDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DNACrossSectionPolicyIonizationBorn 
{
  protected:
   G4DNACrossSectionPolicyIonizationBorn() {}
   ~G4DNACrossSectionPolicyIonizationBorn() {}

   G4double CrossSection(const G4Track& track);
   G4int RandomizePartialCrossSection(const G4Track& track, G4double k);
   G4int NumberOfPartialCrossSections(void);
   void BuildTotalCrossSection(void);

  private:

   void                                 Free(void);
   G4DNACrossSectionDataSet *           dataset;
   G4double *                           valuesBuffer;

   G4double lowEnergyLimit;
   G4double highEnergyLimit;
   G4double zeroBelowLowEnergyLimit;
   G4double zeroAboveHighEnergyLimit;
   G4double dataFileEnergyUnit;
   G4double dataFileCrossSectionUnit;
   G4String dataFileName;
 
   // Hides default constructor and assignment operator as private
   G4DNACrossSectionPolicyIonizationBorn(const G4DNACrossSectionPolicyIonizationBorn & copy);
   G4DNACrossSectionPolicyIonizationBorn & operator=(const G4DNACrossSectionPolicyIonizationBorn & right);
 };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4DNACrossSectionPolicyIonizationBorn.icc"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
