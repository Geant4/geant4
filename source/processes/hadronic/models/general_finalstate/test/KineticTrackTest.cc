//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: KineticTrackTest.cc,v 1.2 2003-10-08 12:19:48 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"

#include "g4templates.hh"

#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4NeutrinoMu.hh"

#include "G4DeltaMinus.hh"
#include "G4AntiDeltaMinus.hh"
#include "G4DeltaZero.hh"
#include "G4AntiDeltaZero.hh"
#include "G4DeltaPlus.hh"
#include "G4AntiDeltaPlus.hh"
#include "G4DeltaPlusPlus.hh"
#include "G4AntiDeltaPlusPlus.hh"

#include "G4Lambda.hh"

#include "G4ParticleDefinition.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"



int main()
{
    G4ParticleDefinition* theGamma = G4Gamma::GammaDefinition();
    G4ParticleDefinition* theMuonPlus = G4MuonPlus::MuonPlusDefinition();
    G4ParticleDefinition* theMuonMinus = G4MuonMinus::MuonMinusDefinition();
    G4ParticleDefinition* theNeutrinoMu = G4NeutrinoMu::NeutrinoMuDefinition();
    G4ParticleDefinition* theProton = G4Proton::ProtonDefinition();
    G4ParticleDefinition* theNeutron = G4Neutron::NeutronDefinition();   
    G4ParticleDefinition* thePionPlus = G4PionPlus::PionPlusDefinition();
    G4ParticleDefinition* thePionMinus = G4PionMinus::PionMinusDefinition();
    G4ParticleDefinition* thePionZero = G4PionZero::PionZeroDefinition();
   
    G4ParticleDefinition* theLambda = G4Lambda::LambdaDefinition();

//
//     New particle (resonance) definition
//

    G4ParticleDefinition* theDeltaMinus =
                          G4DeltaMinus::DeltaMinusDefinition();
    G4ParticleDefinition* theAntiDeltaMinus =
                          G4AntiDeltaMinus::AntiDeltaMinusDefinition();
    G4ParticleDefinition* theDeltaZero =
                          G4DeltaZero::DeltaZeroDefinition();
    G4ParticleDefinition* theAntiDeltaZero =
                          G4AntiDeltaZero::AntiDeltaZeroDefinition();
    G4ParticleDefinition* theDeltaPlus =
                          G4DeltaPlus::DeltaPlusDefinition();
    G4ParticleDefinition* theAntiDeltaPlus =
                          G4AntiDeltaPlus::AntiDeltaPlusDefinition();
    G4ParticleDefinition* theDeltaPlusPlus =
                          G4DeltaPlusPlus::DeltaPlusPlusDefinition();
    G4ParticleDefinition* theAntiDeltaPlusPlus =
                          G4AntiDeltaPlusPlus::AntiDeltaPlusPlusDefinition();

//===========================================================================//
    
    G4double         thePionZeroFormationTime = 33.0*ns;
    G4ThreeVector    thePionZeroInitialPosition 
                     (12.0*cm, 62.0*cm, 65.0*cm);
    G4LorentzVector  thePionZero4Momentum 
                     (0.0*MeV, 0.0*MeV, 0.0*MeV, 134.98*MeV);
    
    G4double         thePionPlusFormationTime = 33.0*ns;
    G4ThreeVector    thePionPlusInitialPosition 
                     (12.0*cm, 62.0*cm, 65.0*cm);
    G4LorentzVector  thePionPlus4Momentum 
                     (0.0*MeV, 0.0*MeV, 0.0*MeV, 139.57*MeV);

    G4double         theProtonFormationTime = 3.0*ns;
    G4ThreeVector    theProtonInitialPosition 
                     (1.0*cm, 6.0*cm, 65.0*cm);
    G4LorentzVector  theProton4Momentum 
                     (0.0*MeV, 0.0*MeV, 0.0*MeV, 938.27*MeV);

    G4double         theDeltaPlusPlusFormationTime = 12.0*ns;
    G4ThreeVector    theDeltaPlusPlusInitialPosition 
                     (12.0*cm, 8.0*cm, 62.0*cm);
    G4LorentzVector  theDeltaPlusPlus4Momentum 
                     (0.0*MeV, 0.0*MeV, 0.0*MeV, 1232.0*MeV);

    G4double         theLambdaFormationTime = 69.0*ns;
    G4ThreeVector    theLambdaInitialPosition 
                     (8.0*cm, 4.0*cm, 25.0*cm);
    G4LorentzVector  theLambda4Momentum 
                     (0.0*MeV, 0.0*MeV, 0.0*MeV, 1115.68*MeV);

    G4KineticTrack   aKT1;
    G4KineticTrack   theProtonKT (theProton, 
                                  theProtonFormationTime,
                                  theProtonInitialPosition,
                                  theProton4Momentum);
                               
    G4KineticTrack   theDeltaPlusPlusKT (theDeltaPlusPlus,
                                         theDeltaPlusPlusFormationTime,
                                         theDeltaPlusPlusInitialPosition,
                                         theDeltaPlusPlus4Momentum);

    G4KineticTrack   theLambdaKT (theLambda,
                                  theLambdaFormationTime,
                                  theLambdaInitialPosition,
                                  theLambda4Momentum);

    G4KineticTrack   aKT2 = theDeltaPlusPlusKT;
    
    G4KineticTrack   thePionZeroKT (thePionZero,
                                    thePionZeroFormationTime,
                                    thePionZeroInitialPosition,
                                    thePionZero4Momentum);
    
    G4KineticTrack   thePionPlusKT (thePionPlus,
                                    thePionPlusFormationTime,
                                    thePionPlusInitialPosition,
                                    thePionPlus4Momentum);
                                    
///////////
// DEBUG //
///////////


   G4ParticleDefinition* theDefinition;
   G4String              theParticleName;
   G4double              theFormationTime;
   G4ThreeVector         thePosition;
   G4ThreeVector         theMomentum;
   G4LorentzVector       theQuadriMomentum, theNewQuadriMomentum,
     theLambdaNew4Momentum1(0.0*MeV, 0.0*MeV, 0.0*MeV, 1500.0*MeV),
     theLambdaNew4Momentum2(10.0*MeV, 20.0*MeV, 30.0*MeV, 1500.47*MeV);
   G4double              theTotalEnergy, theNewTotalEnergy;
   G4double              theMass;
   G4int                 theDecayChannels;
   G4double*             theActualWidth;
   G4double              theTotalActualWidth;
   G4double              theResidualLifetime;
   G4int                 index, i;
   const G4KineticTrackVector* theDecayProductList;
   G4KineticTrack*       theKineticTrack;
   
   theDefinition = theProtonKT.GetDefinition();
   G4cout << G4endl << "========================================================";   
   theDefinition->DumpTable();
   for (i = 1; i <= 10; i++)
      {   
       theDecayProductList = theProtonKT.Decay();
       G4cout << "N. of decay products " << theDecayProductList->entries() <<
               G4endl;
      }      

   theDefinition = theLambdaKT.GetDefinition();
   G4cout << G4endl << "========================================================";   
   theDefinition->DumpTable();
   for (i = 1; i <= 10; i++)
      {
       theDecayProductList = theLambdaKT.Decay();
       for (G4int j = 0; j < theDecayProductList->length(); j++)
          {
           theKineticTrack = theDecayProductList->at(j);
           theDefinition = theKineticTrack->GetDefinition();
           theParticleName = theDefinition->GetParticleName(); 
           theQuadriMomentum = theKineticTrack->Get4Momentum();
           theMomentum = theQuadriMomentum.vect();
           theMass = sqrt((theQuadriMomentum.e() * theQuadriMomentum.e()) -
                          theMomentum.mag2());
           G4cout << G4endl << theParticleName << "   " << theQuadriMomentum
                << "   |p| = " << theMomentum.mag()
                << "   M = " << theMass;
          }
       G4cout << G4endl;
      }
   G4cout << G4endl;

   G4cout << G4endl << "========================================================";   

   theLambdaKT.Set4Momentum(theLambdaNew4Momentum1);   
   for (i = 1; i <= 10; i++)
      {
       theDecayProductList = theLambdaKT.Decay();
       for (G4int j = 0; j < theDecayProductList->length(); j++)
          {
           theKineticTrack = theDecayProductList->at(j);
           theDefinition = theKineticTrack->GetDefinition();
           theParticleName = theDefinition->GetParticleName(); 
           theQuadriMomentum = theKineticTrack->Get4Momentum();
           theMomentum = theQuadriMomentum.vect();
           theMass = sqrt((theQuadriMomentum.e() * theQuadriMomentum.e()) -
                          theMomentum.mag2());
           G4cout << G4endl << theParticleName << "   " << theQuadriMomentum
                << "   |p| = " << theMomentum.mag()
                << "   M = " << theMass;
          }
       G4cout << G4endl;
      }
   G4cout << G4endl;

   G4cout << G4endl << "========================================================";   

   theLambdaKT.Set4Momentum(theLambdaNew4Momentum2);   
   for (i = 1; i <= 10; i++)
      {
       theDecayProductList = theLambdaKT.Decay();
       for (G4int j = 0; j < theDecayProductList->length(); j++)
          {
           theKineticTrack = theDecayProductList->at(j);
           theDefinition = theKineticTrack->GetDefinition();
           theParticleName = theDefinition->GetParticleName(); 
           theQuadriMomentum = theKineticTrack->Get4Momentum();
           theMomentum = theQuadriMomentum.vect();
           theMass = sqrt((theQuadriMomentum.e() * theQuadriMomentum.e()) -
                          theMomentum.mag2());
           G4cout << G4endl << theParticleName << "   " << theQuadriMomentum
                << "   |p| = " << theMomentum.mag()
                << "   M = " << theMass;
          }
       G4cout << G4endl;
      }
   G4cout << G4endl;



}
