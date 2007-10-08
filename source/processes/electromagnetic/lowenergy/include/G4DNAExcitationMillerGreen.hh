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
// $Id: G4DNAExcitationMillerGreen.hh,v 1.1 2007-10-08 09:18:43 sincerti Exp $
// -------------------------------------------------------------------
//

#ifndef G4DNAExcitationMillerGreen_HH
#define G4DNAExcitationMillerGreen_HH 1

#include "G4DNAExcitationInWater.hh"
#include "G4DNACrossSectionPolicyExcitationMillerGreen.hh"
#include "G4DNAFinalStatesPolicyExcitationMillerGreen.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DNAExcitationMillerGreen :
public G4DNAExcitationInWater
<G4DNACrossSectionPolicyExcitationMillerGreen, 
G4DNAFinalStatesPolicyExcitationMillerGreen>
{
 public:
 
   G4DNAExcitationMillerGreen(const G4String & name = "G4DNAExcitationMillerGreen")
   :G4DNAExcitationInWater
   <G4DNACrossSectionPolicyExcitationMillerGreen, 
   G4DNAFinalStatesPolicyExcitationMillerGreen > (name) {}
   
   virtual ~G4DNAExcitationMillerGreen() {}
   
   G4bool IsApplicable(const G4ParticleDefinition& aParticleDefinition) 
   {
     G4DNAGenericIonsManager *instance;
     instance = G4DNAGenericIonsManager::Instance();

     return 
     ( 
       &aParticleDefinition == G4Proton::Proton()
       ||
       &aParticleDefinition == instance->GetIon("helium") 
       ||
       &aParticleDefinition == instance->GetIon("alpha+") 
       ||
       &aParticleDefinition == instance->GetIon("alpha++") 
     ); 
   } 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
