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
// -------------------------------------------------------------
//
//      ---------- TestBertiniStopping -------
//                 Julia Yarba, FNAL
//                 Oct.18, 2011 
// 
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "TestBertiniStopping.hh"

#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"

#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleChange.hh"

// #include "G4HadronicInteraction.hh"
#include "G4CascadeInterface.hh"
#include "G4Nucleus.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TestBertiniStopping::TestBertiniStopping( const G4String& pName, G4ProcessType pType  )
 : G4VRestProcess(pName,pType),
   fTargetNucleus(0), fNSec(0)
{

   fModel = new G4CascadeInterface();
   fModel->SetMinEnergy(0.*GeV);
   fModel->SetMaxEnergy(15.*GeV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TestBertiniStopping::~TestBertiniStopping()
{

  if ( fModel ) delete fModel;
  if ( fTargetNucleus ) delete fTargetNucleus;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TestBertiniStopping::InitTarget( G4Material* mat )
{

/*
  if( mat->GetName() == "D") {
    Z = 1;
    N = 2;
  }
*/
/*
   const G4ElementVector* ev = mat->GetElementVector();
   G4int Z = (G4int)(((*ev)[0])->GetZ() + 0.5);
   G4int N = (G4int)(((*ev)[0])->GetN() + 0.5);
   fTargetNucleus = new G4Nucleus( ((*ev)[0]->GetN()+0.5), ((*ev)[0]->GetZ()+0.5) );
*/

   fTargetNucleus = new G4Nucleus( mat );

   return;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TestBertiniStopping::UsePreCompound()
{
   fModel->usePreCompoundDeexcitation();
   return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* TestBertiniStopping::AtRestDoIt( const G4Track& track,
                                                    const G4Step& )
{

  if ( !fTargetNucleus ) return 0;
  
  fNSec = 0;
  
  G4HadProjectile proj(track);
  G4HadFinalState* result = fModel->ApplyYourself( proj, *fTargetNucleus );
  result->SetTrafoToLab( proj.GetTrafoToLab() );
  
  ClearNumberOfInteractionLengthLeft();

  fPartChange.Initialize(track);

  fNSec = result->GetNumberOfSecondaries();
  G4int nb = fNSec;
  if( result->GetStatusChange() == isAlive ) nb++;
  
  fPartChange.ProposeTrackStatus(fStopAndKill);
  fPartChange.SetNumberOfSecondaries(nb);

  for(G4int i=0; i<fNSec; i++) {
    G4Track* tr = new G4Track(result->GetSecondary(i)->GetParticle(),
                              track.GetGlobalTime(),
	                      track.GetPosition());
    fPartChange.AddSecondary(tr);
  }

  if(result->GetStatusChange() == isAlive) {
    G4DynamicParticle* dp = new G4DynamicParticle(*(track.GetDynamicParticle()));
    G4Track* tr = new G4Track(dp,track.GetGlobalTime(),track.GetPosition());
    tr->SetKineticEnergy(result->GetEnergyChange());
    tr->SetMomentumDirection(result->GetMomentumChange());
    fPartChange.AddSecondary(tr);
  }

  result->Clear();

  return &fPartChange;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
