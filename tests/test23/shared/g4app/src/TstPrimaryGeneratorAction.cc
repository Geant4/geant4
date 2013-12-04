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
#include "TstPrimaryGeneratorAction.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Proton.hh"

#include "TstReader.hh"
#include "TstTarget.hh"

#include "G4RunManager.hh"
#include "G4HadronCrossSections.hh"

TstPrimaryGeneratorAction::~TstPrimaryGeneratorAction()
{

   if ( fPartGun ) delete fPartGun;

}

void TstPrimaryGeneratorAction::InitBeam( TstReader* pset )
{

   if ( fIsInit ) return;
   
   fIsInit=true;
   
   fConfigPtr=pset;

   G4int NPart = 1; 
   fPartGun = new G4ParticleGun( NPart );
      
   G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition* partDef = partTable->FindParticle(pset->GetBeamParticle()); 
   // G4ParticleDefinition* partDef = (G4ParticleTable::GetParticleTable())->FindParticle(pset->GetBeamParticle());
   G4double partMass = partDef->GetPDGMass();
   G4double partMom = pset->GetBeamMomentum();
   G4double partEnergy = std::sqrt( partMom*partMom + partMass*partMass );
   
   fPartGun->SetParticleDefinition( partDef );
   fPartGun->SetParticleEnergy( partEnergy );
   fPartGun->SetParticlePosition( pset->GetPosition() ); // it's already in mm (from TstReader

   // in principle, this should be calculated on the event-by-event basis
   // in case beam direction varies
   
   const TstTarget* target = dynamic_cast<const TstTarget*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
   
   const G4Element* elm = target->GetCurrentMaterial()->GetElement(0);
   G4int A = (G4int)(elm->GetN()+0.5);
   G4int Z = (G4int)(elm->GetZ()+0.5);
   G4double amass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, A);
   
   G4DynamicParticle dParticle( partDef, pset->GetDirection(), partEnergy);  
   fXSecOnTarget = (G4HadronCrossSections::Instance())->GetInelasticCrossSection( &dParticle, Z, A );

       // if under an agnle, then like this:
       //labv = G4LorentzVector(mom.x()/CLHEP::GeV, mom.y()/CLHEP::GeV, 
       //			      mom.z()/CLHEP::GeV, (e0+mass+amass)/CLHEP::GeV);
   
   fLabV.setX(0.);
   fLabV.setY(0.);
   fLabV.setZ( std::sqrt( partEnergy*(partEnergy+2.0*partMass) )/GeV );
   fLabV.setT( (partEnergy+partMass+amass)/GeV );
   fLabP.setX(0.);
   fLabP.setY(0.);
   fLabP.setZ( std::sqrt(partEnergy*(partEnergy+2.0*partMass))/GeV );
   fLabP.setT( (partEnergy+partMass+G4Proton::Proton()->GetPDGMass())/GeV );
   
   return;

}

void TstPrimaryGeneratorAction::GeneratePrimaries( G4Event* evt ) 
{

   // assert( fConfigPtr );

   // G4ThreeVector dir( fConfigPtr->GetDirection() );
   fPartGun->SetParticleMomentumDirection( fConfigPtr->GetDirection()  );
  
   fPartGun->GeneratePrimaryVertex( evt );
     
   return;

}


