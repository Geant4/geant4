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
// $Id: StringDecay.cc,v 1.2 2003-10-08 12:19:48 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include <stdio.h> 
#include <stdlib.h>
#include "G4LundStringFragmentation.hh"
#include "G4QGSMFragmentation.hh"
#include <fstream> 
#include "G4ios.hh" 
  
#include "G4ExcitedString.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleChange.hh"
#include "G4VShortLivedParticle.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"

#include "G4DiQuarks.hh"
#include "G4Quarks.hh"
#include "G4Gluons.hh"

//******************************************************************************************

void ConstructParticle()
  {
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  //Bosons
  G4Gamma::GammaDefinition();
  G4OpticalPhoton::OpticalPhotonDefinition();

  // leptons
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();

 //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();

  G4DMesonPlus::DMesonPlusDefinition();
  G4DMesonMinus::DMesonMinusDefinition();
  G4DMesonZero::DMesonZeroDefinition();
  G4AntiDMesonZero::AntiDMesonZeroDefinition();
  G4DsMesonPlus::DsMesonPlusDefinition();
  G4DsMesonMinus::DsMesonMinusDefinition();
  G4JPsi::JPsiDefinition();

  G4BMesonPlus::BMesonPlusDefinition();
  G4BMesonMinus::BMesonMinusDefinition();
  G4BMesonZero::BMesonZeroDefinition();
  G4AntiBMesonZero::AntiBMesonZeroDefinition();
  G4BsMesonZero::BsMesonZeroDefinition();
  G4AntiBsMesonZero::AntiBsMesonZeroDefinition();
  G4RhoZero::RhoZeroDefinition();

//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
  G4Lambda::LambdaDefinition();
  G4AntiLambda::AntiLambdaDefinition();
  G4SigmaZero::SigmaZeroDefinition();
  G4AntiSigmaZero::AntiSigmaZeroDefinition();
  G4SigmaPlus::SigmaPlusDefinition();
  G4AntiSigmaPlus::AntiSigmaPlusDefinition();
  G4SigmaMinus::SigmaMinusDefinition();
  G4AntiSigmaMinus::AntiSigmaMinusDefinition();
  G4XiZero::XiZeroDefinition();
  G4AntiXiZero::AntiXiZeroDefinition();
  G4XiMinus::XiMinusDefinition();
  G4AntiXiMinus::AntiXiMinusDefinition();
  G4OmegaMinus::OmegaMinusDefinition();
  G4AntiOmegaMinus::AntiOmegaMinusDefinition();

  G4LambdacPlus::LambdacPlusDefinition();
  G4SigmacPlusPlus::SigmacPlusPlusDefinition();
  G4SigmacPlus::SigmacPlusDefinition();
  G4SigmacZero::SigmacZeroDefinition();
  G4XicPlus::XicPlusDefinition();
  G4XicZero::XicZeroDefinition();
  G4OmegacZero::OmegacZeroDefinition();
  G4AntiLambdacPlus::AntiLambdacPlusDefinition();
  G4AntiSigmacPlusPlus::AntiSigmacPlusPlusDefinition();
  G4AntiSigmacPlus::AntiSigmacPlusDefinition();
  G4AntiSigmacZero::AntiSigmacZeroDefinition();
  G4AntiXicPlus::AntiXicPlusDefinition();
  G4AntiXicZero::AntiXicZeroDefinition();
  G4AntiOmegacZero::AntiOmegacZeroDefinition();
  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();
  }
//***************************************************************************
//*****************************************************************************************************
/*
int main(int ArgC, char* ArgS[])
   {
   if (ArgC != 2)
      ArgS[1] = "String.dat";

   //at first Construct particles
   ConstructParticle();
   std::ifstream In;
   In.open(ArgS[1]);

   G4int Encoding;
   In>>Encoding;
   G4Parton Left  = G4Parton(Encoding);
   
   G4double Px,Py,Pz;   
   In >> Px >> Py >> Pz;
   G4ThreeVector Pxyz1(Px*GeV, Py*GeV, Pz*GeV);
   G4double Mass = G4ParticleTable::GetParticleTable()->FindParticle(Encoding)->GetPDGMass();
   G4double E    = sqrt(Pxyz1.mag2() + Mass*Mass);
   G4LorentzVector Mom1(Pxyz1, E);
   Left.Set4Momentum(Mom1);
   
   G4double X, Y, Z;
   In >> X >> Y >> Z;
   G4ThreeVector Pos1(X*fermi, Y*fermi, Z*fermi);
   Left.SetPosition(Pos1);
    
   In>>Encoding;
   G4Parton Right  = G4Parton(Encoding);
   
   In >> Px >> Py >> Pz;
   G4ThreeVector Pxyz2(Px*GeV, Py*GeV, Pz*GeV);
   Mass = G4ParticleTable::GetParticleTable()->FindParticle(Encoding)->GetPDGMass();
   E    = sqrt(Pxyz2.mag2() + Mass*Mass);
   G4LorentzVector Mom2(Pxyz2, E);
   Right.Set4Momentum(Mom2);
   
   In >> X >> Y >> Z;
   G4ThreeVector Pos2(X*fermi, Y*fermi, Z*fermi);
   Right.SetPosition(Pos2);

   G4ExcitedString theString(&Left, &Right);
   G4cout<<"String 3-position [fermi]"<<G4endl;
   G4cout<<theString.GetPosition()<<G4endl;
   G4cout<<"String 4-momentum [MeV]"<<G4endl;
   G4cout<<theString.Get4Momentum().vect()<<" "<<theString.Get4Momentum().e()<<G4endl;
   G4cout<<"Mass = "<<theString.Get4Momentum().mag()<<G4endl;
   for(G4int c2 = 0; c2 < 100; c2++)
      {
      G4ExcitedString* PtrString = &theString;
      G4LundStringFragmentation StringDecay;
      G4QGSMFragmentation QGSMFragmentationTest;
      G4KineticTrackVector* Output = StringDecay.FragmentString(PtrString);
      printf("Number of hadron: %d\n", (int)(Output->length()));
      G4double px, py, pz, e;
      px = py =pz = e  = 0;
      for(int c1 = 0; c1 < Output->length(); c1++) 
	 {  
	 G4KineticTrack* Trk = Output->at(c1);
	 printf("Enc = %7d Px = %10g Py = %10g Pz = %10g E = %10f\n", 
	     Trk->GetDefinition()->GetPDGEncoding(),
	     Trk->Get4Momentum().px(), 
	     Trk->Get4Momentum().py(),
	     Trk->Get4Momentum().pz(),
	     Trk->Get4Momentum().e());
	 printf("x = %10g y = %10g z = %10g\n", 
	     Trk->GetPosition().x(), 
	     Trk->GetPosition().y(),
	     Trk->GetPosition().z()); 
	 px += Trk->Get4Momentum().px(); 
	 py += Trk->Get4Momentum().py();
	 pz += Trk->Get4Momentum().pz();
	 e  += Trk->Get4Momentum().e ();
	 }   
	 
      printf("Px = %10f Py = %10f Pz = %10f E = %10f\n", px, py, pz, e);   	   
      Output->clearAndDestroy();
      }
  }
*/
//*****************************************************************************************    


G4int LeftEncoding[] = { -1, -2, -3, 1103, 3103, 3303 }; 
G4int RightEncoding[]  = {1, 2,  3 };

int main(int ArgC, char* ArgS[])
   {
   G4double PzMax = 10000;
   if (ArgC == 2)
      PzMax = atof(ArgS[1]);
   ConstructParticle();
   for(G4int c1 = 0; c1 < sizeof(LeftEncoding)/sizeof(G4int); c1++)
      for(G4int c2 = 0; c2 < sizeof(RightEncoding)/sizeof(G4int); c2++)
          {
          G4Parton Left  = G4Parton(LeftEncoding[c1]);
          G4Parton Right = G4Parton(RightEncoding[c2]);
          G4double PzMin;
          if (c1 < 3)
              {
              PzMin = 375;
              if (c1 < 2 && c2 < 2)
                 PzMin = 150;
              if (c1 == 2 && c2 == 2)
                 PzMin = 625;
             }   
          else
              {
              PzMin = 850; 
              if ((c1-3) < 2 && c2 < 2)
                 PzMin = 770; 
              if (c2 == 2)   
                 if (c1 == 4)
                     PzMin = 925;
                 else
                     PzMin = 1000;
             }   
          if (PzMax < PzMin)
              PzMax = PzMin + 201;
          for(G4double Pz = PzMin; Pz < PzMax; Pz += 200)
              {
	      if (c1 == 2 && c2 ==1 && Pz == PzMin)
	          G4cout << "Kyky";
	      G4ThreeVector Pxyz1(0, 0, Pz);
	      G4ParticleDefinition* pLeftDef = G4ParticleTable::GetParticleTable()->FindParticle(LeftEncoding[c1]);
	      G4double Mass = pLeftDef->GetPDGMass();
	      G4double E    = sqrt(Pxyz1.mag2() + Mass*Mass);
	      G4LorentzVector Mom1(Pxyz1, E);
	      Left.Set4Momentum(Mom1);
	      G4ThreeVector Pos(0*fermi, 0*fermi, 0*fermi);
	      Left.SetPosition(Pos);

	      G4ThreeVector Pxyz2(0, 0, -Pz);
	      G4ParticleDefinition* pRightDef = G4ParticleTable::GetParticleTable()->FindParticle(RightEncoding[c2]);
	      Mass = pRightDef->GetPDGMass();
	      E    = sqrt(Pxyz2.mag2() + Mass*Mass);
	      G4LorentzVector Mom2(Pxyz2, E);
	      Right.Set4Momentum(Mom2);
	      Right.SetPosition(Pos);
	      G4ExcitedString theString(&Left, &Right);

              G4cout << "Pz min = "<<PzMin<<"  Pz = "<<Pz<<G4endl;
              G4cout<<G4endl<<"String: " <<LeftEncoding[c1]<<" "<<RightEncoding[c2]<<G4endl;
              G4cout<<theString.Get4Momentum().vect()<<" "<<theString.Get4Momentum().e()<<G4endl;
              G4cout<<"Mass = "<<theString.Get4Momentum().mag()<<G4endl;
              G4cout<<"Charge = "<<pLeftDef->GetPDGCharge() + pRightDef->GetPDGCharge()<<G4endl;
              G4cout<<"Barion charge = "<<pLeftDef->GetBaryonNumber() + pRightDef->GetBaryonNumber()<<G4endl;
              G4LundStringFragmentation StringDecay;
              G4KineticTrackVector* Output = StringDecay.FragmentString(&theString);
              G4cout<<"Number of hadron: " << Output->length()<<G4endl;
	      G4double px, py, pz, e;
	      G4double Charge, BaryonNumber;
	      px = py =pz = e = 0;
	      Charge = BaryonNumber = 0;
	      for(G4int c3 = 0; c3 < Output->length(); c3++) 
		 {  
		 G4KineticTrack* Trk = Output->at(c3);
		 px += Trk->Get4Momentum().px(); 
		 py += Trk->Get4Momentum().py();
		 pz += Trk->Get4Momentum().pz();
		 e  += Trk->Get4Momentum().e ();
		 Charge += Trk->GetDefinition()->GetPDGCharge();
		 BaryonNumber += Trk->GetDefinition()->GetBaryonNumber();
		 }   
	      
	      G4LorentzVector Mom(px, py, pz, e);
	      if (abs(theString.Get4Momentum().mag()- Mom.mag()) > 10E-3 ||
	         abs(pLeftDef->GetPDGCharge() + pRightDef->GetPDGCharge() - Charge)> 10E-3 ||	 
                 (c1/3 != BaryonNumber))
                   {
                   G4cout<<" Error:"<<G4endl;
                   G4cout<<"  Momentum: "<<Mom.vect()<< " E = "<<Mom.e()<<G4endl;
                   G4cout<<"  Charge = "<<Charge << "  BaryonNumber = "<<BaryonNumber<<G4endl;
		   for(G4int c4 = 0; c4 < Output->length(); c4++) 
		      {  
		      G4KineticTrack* Trk = Output->at(c4);
		      G4ParticleDefinition* pD = Trk->GetDefinition();
		      G4cout<<"     Encoding = "<<pD->GetPDGEncoding();
		      G4cout<<"     Mass = "<<pD->GetPDGMass();
		      G4cout<<"     Charge = "<<pD->GetPDGCharge();
		      G4cout<<"     BaryonNumber = "<<pD->GetBaryonNumber()<<G4endl;
		      }   
                   Output->clearAndDestroy();
                   delete Output; 
                   return 0;
                   }
              Output->clearAndDestroy();
              delete Output; 
              }
          }
  }

//*****************************************************************************************    
 
