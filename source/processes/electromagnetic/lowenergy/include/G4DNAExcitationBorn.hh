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
// $Id: G4DNAExcitationBorn.hh,v 1.1 2007-10-08 09:18:43 sincerti Exp $
// -------------------------------------------------------------------
//

#ifndef G4DNAExcitationBorn_HH
#define G4DNAExcitationBorn_HH 1

#include "G4DNAExcitationInWater.hh"
#include "G4DNACrossSectionPolicyExcitationBorn.hh"
#include "G4DNAFinalStatesPolicyExcitationBorn.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DNAExcitationBorn :
public G4DNAExcitationInWater
<G4DNACrossSectionPolicyExcitationBorn, 
G4DNAFinalStatesPolicyExcitationBorn>
{
 public:
 
   G4DNAExcitationBorn(const G4String & name = "G4DNAExcitationBorn")
   :G4DNAExcitationInWater
   <G4DNACrossSectionPolicyExcitationBorn, 
   G4DNAFinalStatesPolicyExcitationBorn > (name) {}
   
   virtual ~G4DNAExcitationBorn() {}
   
   G4bool IsApplicable(const G4ParticleDefinition& aParticleDefinition) 
   {
    return 
    ( 
      &aParticleDefinition == G4Proton::Proton()
     ); 
   } 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
