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
#include "ProcessWrapper.hh"

ProcessWrapper::~ProcessWrapper()
{
  
   // need to see how desctructors of other data members work !!!...
   // ... to avoid potential memory leaks...
   //

}


G4double ProcessWrapper::PostStepGetPhysicalInteractionLength(const G4Track&,
                        	     		                    G4double,
			                                            G4ForceCondition* condition)
{

  *condition = NotForced;
  G4double z = DBL_MAX;

  return z;

}

G4VParticleChange* ProcessWrapper::PostStepDoIt( const G4Track& track, const G4Step& )
{

  G4Material* mat = track.GetMaterial();
  G4Element*  elm = (G4Element*)mat->GetElement(0); // terrible trick - cast away const...
   
  G4Nucleus* tNucleus = GetTargetNucleusPointer();
  tNucleus->ChooseParameters( mat );
    
  G4HadProjectile proj(track);

  G4HadronicInteraction* hint = ChooseHadronicInteraction( proj, *tNucleus, mat, elm );
  
  if ( !hint ) return 0;
  
  G4HadFinalState* result = hint->ApplyYourself( proj, *tNucleus );

  result->SetTrafoToLab( proj.GetTrafoToLab() );
  
  ClearNumberOfInteractionLengthLeft();

  fPartChange.Initialize(track);

  G4int NSec = result->GetNumberOfSecondaries();
  G4int nb = NSec;
  if( result->GetStatusChange() == isAlive ) nb++;
  
  fPartChange.ProposeTrackStatus(fStopAndKill);
  fPartChange.SetNumberOfSecondaries(nb);

  for(G4int i=0; i<NSec; i++) {
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

