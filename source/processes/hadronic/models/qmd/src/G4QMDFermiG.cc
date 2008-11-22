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
#include "G4QMDFermiG.hh"

#include "G4QMDFermiParticipants.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QGSMFragmentation.hh"
#include "G4BinaryCascade.hh"
#include "G4PreCompoundModel.hh"
#include "G4ParticleTable.hh"

G4QMDFermiG::G4QMDFermiG()
{

   theHEG = new G4QGSModel< G4QMDFermiParticipants >;
   G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay(new G4QGSMFragmentation);
   theHEG->SetFragmentationModel(theStringDecay);

   theTransport = new G4BinaryCascade;
   G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(new G4ExcitationHandler);

   theTransport->SetDeExcitation(thePreEquilib);

}



G4QMDFermiG::~G4QMDFermiG()
{
   delete theHEG;
   delete theTransport;
}



G4HadFinalState* G4QMDFermiG::ApplyYourself( const G4HadProjectile &aTrack , G4Nucleus & targetNucleus )
{

   // init particle change
   theParticleChange.Clear();
   theParticleChange.SetStatusChange(stopAndKill);

   G4int z_tot = int ( aTrack.GetDefinition()->GetPDGCharge()/eplus  + targetNucleus.GetZ() ); 
   G4int a_tot = aTrack.GetDefinition()->GetBaryonNumber() + int ( targetNucleus.GetN() ); 

   G4LorentzVector p4_tot = aTrack.Get4Momentum() + G4LorentzVector( G4ThreeVector(0) , targetNucleus.AtomicMass( targetNucleus.GetN() , targetNucleus.GetZ() ) ); 

   //G4cout << "High Energy " << G4endl;

   G4DynamicParticle aTemp(const_cast<G4ParticleDefinition *>(aTrack.GetDefinition()),
                          aTrack.Get4Momentum().vect());
   const G4DynamicParticle * aPart = &aTemp;

   G4KineticTrackVector* theInitialResult = theHEG->Scatter( targetNucleus , *aPart );
  
/*
   for ( G4KineticTrackVector::iterator it = theInitialResult->begin() ; it != theInitialResult->end() ; it++ )
   {

      G4cout << "G4KineticTrack" 
             << " " << (*it)->GetDefinition()->GetParticleName()
             << " " << (*it)->GetDefinition()->GetPDGCharge()/eplus
             << " " << (*it)->GetDefinition()->GetBaryonNumber()
             << " " << (*it)->Get4Momentum()
             << " " << (*it)->GetPosition()/fermi
             << G4endl;
   }
*/

   G4DecayStrongResonances theDecay;

   G4ReactionProductVector * theTransportResult = NULL;
   G4int hitCount = 0;
   const std::vector<G4Nucleon *> & they = theHEG->GetWoundedNucleus()->GetNucleons();
   for(size_t them=0; them<they.size(); them++)
   {
     if(they[them]->AreYouHit()) hitCount ++;
   }
   if(hitCount != theHEG->GetWoundedNucleus()->GetMassNumber() )
   {
     theTransportResult =
                theTransport->Propagate(theInitialResult, theHEG->GetWoundedNucleus());
     if ( !theTransportResult ) {
        G4cout << "G4QMDFermiG: null ptr from transport propagate " << G4endl;
        throw G4HadronicException(__FILE__, __LINE__, "Null ptr from transport propagate");
     }
   }
   else
   {
     theTransportResult = theDecay.Propagate(theInitialResult, theHEG->GetWoundedNucleus());
     if ( !theTransportResult ) {
        G4cout << "G4QMDFermiG: null ptr from decay propagate " << G4endl;
        throw G4HadronicException(__FILE__, __LINE__, "Null ptr from decay propagate");
     }
   }

   G4int z_sum = 0; 
   G4int a_sum = 0; 
   G4LorentzVector p4_sum(0); 
   for ( G4ReactionProductVector::iterator it = theTransportResult->begin() ; it != theTransportResult->end() ; it++ )
   {

      z_sum += int ( (*it)->GetDefinition()->GetPDGCharge()/eplus ); 
      a_sum += (*it)->GetDefinition()->GetBaryonNumber();
      p4_sum += G4LorentzVector( (*it)->GetMomentum() , (*it)->GetDefinition()->GetPDGMass() + (*it)->GetKineticEnergy() );

/*
      G4cout << "G4ReactionProduct" 
             << " " << (*it)->GetDefinition()->GetParticleName()
             << " " << (*it)->GetDefinition()->GetPDGCharge()/eplus
             << " " << (*it)->GetDefinition()->GetBaryonNumber()
             << " " << (*it)->GetKineticEnergy()/GeV
             << " " << (*it)->GetMomentum().x()
             << " " << (*it)->GetMomentum().y()
             << " " << (*it)->GetMomentum().z()
             << G4endl;
*/

     theParticleChange.AddSecondary( new G4DynamicParticle ( (*it)->GetDefinition() , (*it)->GetTotalEnergy() , (*it)->GetMomentum()  ) );
     delete *it;
   } 

/*
   G4cout << "TK tot " << z_tot << " " << a_tot << " " << p4_tot << G4endl;
   G4cout << "TK sub " << z_sum << " " << a_sum << " " << p4_sum << G4endl;
   G4cout << "TK Dif " << z_tot - z_sum << " " << a_tot - a_sum << " " << p4_tot - p4_sum << G4endl;
   G4cout << "TK ExE " << (p4_tot - p4_sum).m() - targetNucleus.AtomicMass( a_tot - a_sum , z_tot - z_sum ) << G4endl;
*/
   G4double exe = (p4_tot - p4_sum).m() - targetNucleus.AtomicMass( a_tot - a_sum , z_tot - z_sum );
   if ( exe < 0 ) exe = 0.0;
   // neglect excitation energy TK 
   exe = 0.0;

   if ( z_tot - z_sum > 0 && a_tot - a_sum ) 
      theParticleChange.AddSecondary( new G4DynamicParticle ( G4ParticleTable::GetParticleTable()->GetIon( z_tot - z_sum , a_tot - a_sum , exe ) , ( p4_tot - p4_sum ).v() ) );

   return &theParticleChange;
}

