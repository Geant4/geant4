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
#include <stdio.h>

#include "G4ParticleTable.hh"
#include "G4LeptonConstructor.hh" 
#include "G4BarionConstructor.hh" 
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4QuarkAnnihilator.hh"

#include "G4LundStringFragmentation.hh"
#include "G4QGSMFragmentation.hh"


//**********************************************************************************************************************

G4double MaxUniformRand();

void ConstructParticle()
  {
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4LeptonConstructor Leptons;
  Leptons.ConstructParticle();

  G4BosonConstructor Bosons;
  Bosons.ConstructParticle();

  G4MesonConstructor Mesons;
  Mesons.ConstructParticle();

  G4BarionConstructor Barions;
  Barions.ConstructParticle();

  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();
  }

//**********************************************************************************************************************


struct
   {
   double x;
   double y;
   double z;
   double t;
   double Px;
   double Py;
   double Pz;
   double e;
   double mass;
   int    Encoding;
   int    nHadron;
   int    nEvent;
   int    cHadron;
   } Buff;

int main(int ArgC, char* ArgS[])
   {
   ConstructParticle();
   G4int LeftEncoding  =  2;
   G4int RightEncoding = -2;
   G4int nEvent = 1000;
   G4double Pz = 10*GeV;
   
   if (ArgC == 2)
      Pz = atof(ArgS[1])*GeV;
   if (ArgC == 3)
      {
      LeftEncoding  = atoi(ArgS[1]);
      RightEncoding = atoi(ArgS[2]);
      }
   if (ArgC == 4)
      {
      LeftEncoding  = atoi(ArgS[1]);
      RightEncoding = atoi(ArgS[2]);
      Pz = atof(ArgS[3])*GeV;
      }
   G4ParticleDefinition* ptr = G4ParticleTable::GetParticleTable()->FindParticle(LeftEncoding);
   if (ptr == NULL)
       G4Exception("Check your particle table");
   G4double Mass = ptr->GetPDGMass();
   G4ThreeVector   Pos(0,0,0);
   G4LorentzVector Mom1(0, 0,  Pz, sqrt(Mass*Mass + Pz*Pz));
   G4Parton* Left = new G4Parton(LeftEncoding);
   Left->SetPosition (Pos);
   Left->Set4Momentum(Mom1);
    
   ptr = G4ParticleTable::GetParticleTable()->FindParticle(RightEncoding);
   if (ptr == NULL)
       G4Exception("Check your particle table");
   Mass = ptr->GetPDGMass();
   G4LorentzVector Mom2(0, 0, -Pz, sqrt(Mass*Mass + Pz*Pz));
   G4Parton* Right = new G4Parton(RightEncoding);
   Right->SetPosition(Pos);
   Right->Set4Momentum(Mom2);
   

   G4ExcitedString theString(Left, Right);
//   G4LundStringFragmentation StringDecay;
   G4QGSMFragmentation StringDecay;
   FILE* out = fopen("kyky.dat", "wb");
   if (!out) 
       {
       G4Exception("I cannot open output file");
       return 0;
       }
   for(G4int cEvent = 0; cEvent < nEvent; cEvent++)
       {
       G4KineticTrackVector* OutputDecay   = StringDecay.FragmentString(&theString);
  //     G4KineticTrackVector* OutputDecay = StringDecay.DecayResonans(Output);
       G4int nHadron = OutputDecay->length();
       for(G4int c1 = 0; c1 < nHadron; c1++)
          {
          G4KineticTrack*  KT = OutputDecay->at(c1);
          G4ThreeVector    Pos = KT->GetPosition();
          G4LorentzVector  Mom = KT->Get4Momentum();
	  Buff.x  = Pos.x();
	  Buff.y  = Pos.y();
	  Buff.z  = Pos.z();
	  Buff.t  = KT->GetFormationTime();
	  Buff.Px = Mom.px();
	  Buff.Py = Mom.py();
	  Buff.Pz = Mom.pz();
	  Buff.e  = Mom.e ();
	  Buff.mass     = KT->GetDefinition()->GetPDGMass();
	  Buff.Encoding = KT->GetDefinition()->GetPDGEncoding();
	  Buff.nEvent   = cEvent;
          Buff.nHadron  = nHadron;
          Buff.cHadron  = c1;
          
          fwrite(&Buff, sizeof(Buff), 1, out);
          }       
//       Output->clearAndDestroy();
//       delete Output;
       OutputDecay->clearAndDestroy();
       delete OutputDecay;
       }
  fclose(out);     
  }

//***************************************************************************************************************************    
 































