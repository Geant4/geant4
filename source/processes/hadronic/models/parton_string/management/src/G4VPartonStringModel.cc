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
// $Id: G4VPartonStringModel.cc 102023 2016-12-16 14:39:26Z gcosmo $
//
//// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4VPartonStringModel ----------------
//             by Gunter Folger, May 1998.
//      abstract class for all Parton String Models
// ------------------------------------------------------------
// debug switch
//#define debug_PartonStringModel
// ------------------------------------------------------------

#include "G4VPartonStringModel.hh"
#include "G4ios.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

G4VPartonStringModel::G4VPartonStringModel(const G4String& modelName)
    : G4VHighEnergyGenerator(modelName),
      stringFragmentationModel(0),
      theThis(0)
{
//  Make shure Shotrylived particles are constructed.
	G4ShortLivedConstructor ShortLived;
	ShortLived.ConstructParticle();
}

G4VPartonStringModel::~G4VPartonStringModel()
{
}

G4KineticTrackVector * G4VPartonStringModel::Scatter(const G4Nucleus &theNucleus, 
                                                const G4DynamicParticle &aPrimary)
{  
  G4ExcitedStringVector * strings = NULL;

  G4DynamicParticle thePrimary=aPrimary;

  G4LorentzVector SumStringMom(0.,0.,0.,0.);  // Uzhi 21 Sept. 2016

  G4KineticTrackVector * theResult = 0;       // Uzhi 21 Sept. 2016

  G4Nucleon * theNuclNucleon(0);              // Uzhi 11 Oct. 2016

#ifdef debug_PartonStringModel
  G4cout<<G4endl;
  G4cout<<"-----------------------Parton-String model is runnung ------------"<<G4endl;
  G4cout<<"Projectile Name Mass "<<thePrimary.GetDefinition()->GetParticleName()<<" "
                                 <<thePrimary.GetMass()<<G4endl;
  G4cout<<"           Momentum  "<<thePrimary.Get4Momentum()<<G4endl;
  G4cout<<"Target nucleus   A Z "<<theNucleus.GetA_asInt()<<" "
                                 <<theNucleus.GetZ_asInt()<<G4endl<<G4endl;	
  G4int Bsum=thePrimary.GetDefinition()->GetBaryonNumber() + theNucleus.GetA_asInt();
  G4int Qsum=thePrimary.GetDefinition()->GetPDGCharge() + theNucleus.GetZ_asInt();
  G4cout<<"Initial baryon number "<<Bsum<<G4endl;
  G4cout<<"Initial charge        "<<Qsum<<G4endl;
  G4cout<<"-------------- Parton-String model:  Generation of strings -------"<<G4endl<<G4endl;

  Bsum -= theNucleus.GetA_asInt();  Qsum -= theNucleus.GetZ_asInt();
  if(theThis->GetProjectileNucleus()) {
    Bsum -= thePrimary.GetDefinition()->GetBaryonNumber();
    Qsum -= thePrimary.GetDefinition()->GetPDGCharge();
  }

  G4int QsumSec(0), BsumSec(0);                 // Uzhi 21 Sept. 2016
  G4LorentzVector SumPsecondr(0.,0.,0.,0.);     // Uzhi 21 Sept. 2016
#endif

  G4LorentzRotation toZ;
  G4LorentzVector Ptmp=thePrimary.Get4Momentum();
  toZ.rotateZ(-1*Ptmp.phi());
  toZ.rotateY(-1*Ptmp.theta());
  thePrimary.Set4Momentum(toZ*Ptmp);
  G4LorentzRotation toLab(toZ.inverse());

  G4bool Success=true;                                                         // Uzhi 21 Sept. 2016
  G4int attempts = 0, maxAttempts=1000;                                        // Uzhi 12 Oct.  2016
  do                                                                           // Uzhi 21 Sept. 2016
  {
  	if (attempts++ > maxAttempts ) 
  	{
            theThis->Init(theNucleus,thePrimary);  // To put a nucleus into ground state
                                                   // But marks of hitted nucleons are left. They must be erased.

    	    G4V3DNucleus * ResNucleus=theThis->GetWoundedNucleus(); 
    	    theNuclNucleon = ResNucleus->StartLoop() ? ResNucleus->GetNextNucleon() : NULL;
            while( theNuclNucleon )
            {
             if(theNuclNucleon->AreYouHit()) theNuclNucleon->Hit(nullptr);
             theNuclNucleon = ResNucleus->GetNextNucleon();
            }

            G4V3DNucleus * ProjResNucleus=theThis->GetProjectileNucleus();
            if(ProjResNucleus != 0)
            {
             theNuclNucleon = ProjResNucleus->StartLoop() ? ProjResNucleus->GetNextNucleon() : NULL;
             while( theNuclNucleon )
             {
              if(theNuclNucleon->AreYouHit()) theNuclNucleon->Hit(nullptr);
              theNuclNucleon = ProjResNucleus->GetNextNucleon();
             }
            }

            G4ExceptionDescription ed;
            ed << "Projectile Name Mass " <<thePrimary.GetDefinition()->GetParticleName()
               << " " << thePrimary.GetMass()<< G4endl;
            ed << "           Momentum  " << thePrimary.Get4Momentum() <<G4endl;
            ed << "Target nucleus   A Z " << theNucleus.GetA_asInt() << " "
                                          << theNucleus.GetZ_asInt() <<G4endl;
            ed << "Initial states of projectile and target nucleus will be returned!"<<G4endl;

            G4Exception( "G4VPartonStringModel::Scatter(): fails to generate or fragment strings ",
                         "HAD_PARTON_STRING_001", JustWarning, ed );

            G4ThreeVector Position(0.,0.,2*ResNucleus->GetOuterRadius());
	    G4KineticTrack* Hadron = new G4KineticTrack(aPrimary.GetParticleDefinition(),
                                                            0., Position,
                                                        aPrimary.Get4Momentum());

            if(theResult == 0) theResult = new G4KineticTrackVector();
            theResult->push_back(Hadron);

            return theResult;
  	}

   Success=true;

   theThis->Init(theNucleus,thePrimary);

   strings = GetStrings();

   if(strings == 0) {Success=false; continue;}     // Uzhi 28 Sept. 2016

   G4double stringEnergy(0);
   SumStringMom=G4LorentzVector(0.,0.,0.,0.);      // Uzhi 21 Sept. 2016

#ifdef debug_PartonStringModel
  G4cout<<"------------ Parton-String model: Number of produced strings ---- "<<strings->size()<<G4endl;
#endif

   for ( unsigned int astring=0; astring < strings->size(); astring++)
   {
    //    rotate string to lab frame, models have it aligned to z
    if((*strings)[astring]->IsExcited())
    {
     stringEnergy += (*strings)[astring]->GetLeftParton()->Get4Momentum().t();
     stringEnergy += (*strings)[astring]->GetRightParton()->Get4Momentum().t();
     (*strings)[astring]->LorentzRotate(toLab);
     SumStringMom+=(*strings)[astring]->Get4Momentum();
#ifdef debug_PartonStringModel
G4cout<<"String No "<<astring+1<<" "<<(*strings)[astring]->Get4Momentum()<<" "
                                    <<(*strings)[astring]->Get4Momentum().mag()
     <<" Partons   "<<(*strings)[astring]->GetLeftParton()->GetDefinition()->GetPDGEncoding()
     <<"          "<<(*strings)[astring]->GetRightParton()->GetDefinition()->GetPDGEncoding()<<G4endl;
#endif
    }
    else
    {
     stringEnergy += (*strings)[astring]->GetKineticTrack()->Get4Momentum().t();
     (*strings)[astring]->LorentzRotate(toLab);
     SumStringMom+=(*strings)[astring]->GetKineticTrack()->Get4Momentum();
#ifdef debug_PartonStringModel
     G4cout<<"A track No "<<astring+1<<" "
           <<(*strings)[astring]->GetKineticTrack()->Get4Momentum()<<" "
           <<(*strings)[astring]->GetKineticTrack()->Get4Momentum().mag()<<" "
           <<(*strings)[astring]->GetKineticTrack()->GetDefinition()->GetParticleName()<<G4endl;
#endif
    }
   }

#ifdef debug_PartonStringModel
   G4cout<<G4endl<<"SumString4Mom "<<SumStringMom<<G4endl;

   G4LorentzVector TargetResidual4Momentum(0.,0.,0.,0.);
   G4LorentzVector ProjectileResidual4Momentum(0.,0.,0.,0.);
   G4int hitsT(0), charged_hitsT(0);
   G4int hitsP(0), charged_hitsP(0);
   G4double ExcitationEt(0.), ExcitationEp(0.);
#endif

   G4V3DNucleus * ProjResNucleus=theThis->GetProjectileNucleus();

   if(ProjResNucleus != 0)
   {
    theNuclNucleon = ProjResNucleus->StartLoop() ? ProjResNucleus->GetNextNucleon() : NULL;
    while( theNuclNucleon )
    {
     if(theNuclNucleon->AreYouHit())
     {
      G4LorentzVector tmp=toLab*theNuclNucleon->Get4Momentum();
      #ifdef debug_PartonStringModel
         ProjectileResidual4Momentum += tmp;
         hitsP++;
         if ( theNuclNucleon->GetDefinition() == G4Proton::Proton() )  ++charged_hitsP;
         ExcitationEp +=theNuclNucleon->GetBindingEnergy();
      #endif
      theNuclNucleon->SetMomentum(tmp);
     }
     theNuclNucleon = ProjResNucleus->GetNextNucleon();
    }
#ifdef debug_PartonStringModel
    G4cout<<"Projectile residual A, Z and E* "
          <<thePrimary.GetDefinition()->GetBaryonNumber() - hitsP<<" "
          <<thePrimary.GetDefinition()->GetPDGCharge()    - charged_hitsP<<" "
          <<ExcitationEp<<G4endl;
    G4cout<<"Projectile residual 4 momentum  "<<ProjectileResidual4Momentum<<G4endl;
#endif
   }

   G4V3DNucleus * ResNucleus=theThis->GetWoundedNucleus(); 

   // loop over wounded nucleus
   theNuclNucleon = ResNucleus->StartLoop() ? ResNucleus->GetNextNucleon() : NULL;
   while( theNuclNucleon )
   {
     if(theNuclNucleon->AreYouHit())
     {
      G4LorentzVector tmp=toLab*theNuclNucleon->Get4Momentum();
      #ifdef debug_PartonStringModel
         TargetResidual4Momentum += tmp;
         hitsT++;
         if ( theNuclNucleon->GetDefinition() == G4Proton::Proton() )  ++charged_hitsT;
         ExcitationEt +=theNuclNucleon->GetBindingEnergy();
      #endif
      theNuclNucleon->SetMomentum(tmp);
     }
     theNuclNucleon = ResNucleus->GetNextNucleon();
   }
#ifdef debug_PartonStringModel
   G4cout<<"Target residual A, Z and E* "
         <<theNucleus.GetA_asInt() - hitsT<<" "
         <<theNucleus.GetZ_asInt() - charged_hitsT<<" "
         <<ExcitationEt<<G4endl;
   G4cout<<"Target residual 4 momentum "<<TargetResidual4Momentum<<G4endl;

   Bsum+=(        hitsT +         hitsP);
   Qsum+=(charged_hitsT + charged_hitsP);

   G4cout<<"Hitted # of nucleons of projectile and target "<<hitsP<<" "<<hitsT<<G4endl;
   G4cout<<"Hitted # of protons of projectile and target  "
         <<charged_hitsP<<" "<<charged_hitsT<<G4endl<<G4endl;
   G4cout<<"Bsum Qsum "<<Bsum<<" "<<Qsum<<G4endl<<G4endl;
#endif

//=========================================================================================
//                              Fragment strings
#ifdef debug_PartonStringModel
  G4cout<<"---------------- Attempt to fragment strings ------------- "<<attempts<<G4endl;
#endif

   G4double InvMass=SumStringMom.mag();
   G4double SumMass(0.); 

#ifdef debug_PartonStringModel
  QsumSec=0; BsumSec=0;                              // Uzhi 21 Sept. 2016 
  SumPsecondr=G4LorentzVector(0.,0.,0.,0.);          // Uzhi 21 Sept. 2016
#endif

   if(theResult != 0)
   {
    std::for_each(theResult->begin(), theResult->end(), DeleteKineticTrack());
    delete theResult;
   }

   theResult = stringFragmentationModel->FragmentStrings(strings);

#ifdef debug_PartonStringModel
  G4cout<<"String fragmentation success (OK - #0, Not OK - 0) : "<<theResult<<G4endl;
#endif

   if(theResult == 0) {Success=false; continue;}           //  Uzhi 21 Sept. 2016

#ifdef debug_PartonStringModel
  G4cout<<"Parton-String model: Number of produced particles "<<theResult->size()<<G4endl;
  SumPsecondr = G4LorentzVector(0.,0.,0.,0.);
  QsumSec = 0; BsumSec = 0;
#endif

   SumMass=0.;

   for ( unsigned int i=0; i < theResult->size(); i++)
   {
     SumMass+=(*theResult)[i]->Get4Momentum().mag();
     #ifdef debug_PartonStringModel
       G4cout<<i<<" : "<<(*theResult)[i]->GetDefinition()->GetParticleName()<<" "
                       <<(*theResult)[i]->Get4Momentum()<<" "
                       <<(*theResult)[i]->Get4Momentum().mag()<<" "
                       <<(*theResult)[i]->GetDefinition()->GetPDGMass()<<G4endl;
       SumPsecondr+=(*theResult)[i]->Get4Momentum();
       BsumSec += (*theResult)[i]->GetDefinition()->GetBaryonNumber();
       QsumSec += (*theResult)[i]->GetDefinition()->GetPDGCharge();
     #endif
   }

#ifdef debug_PartonStringModel
   G4cout<<G4endl<<"-----------------------Parton-String model: balances -------------"<<G4endl;
   if(Qsum != QsumSec) {
     G4cout<<"Charge is not conserved!!! ----"<<G4endl; 
     G4cout<<" Qsum != QsumSec "<<Qsum<<" "<<QsumSec<<G4endl;
   }
   if(Bsum != BsumSec) {
     G4cout<<"Baryon number is not conserved!!!"<<G4endl;
     G4cout<<" Bsum != BsumSec "<<Bsum<<" "<<BsumSec<<G4endl;
   }
#endif

   if((SumMass > InvMass)||(SumMass == 0.)) {Success=false;}      // Uzhi 21 Sept. 2016
   std::for_each(strings->begin(), strings->end(), DeleteString() );
   delete strings;
  } while(!Success); // End of do                              // Uzhi 21 Sept. 2016
#ifdef debug_PartonStringModel
  G4cout        <<"Baryon number balance "<<Bsum-BsumSec<<G4endl;
  G4cout        <<"Charge balance        "<<Qsum-QsumSec<<G4endl;
  G4cout        <<"4 momentum balance    "<<SumStringMom-SumPsecondr<<G4endl; 
  G4cout<<"---------------------End of Parton-String model work -------------"<<G4endl<<G4endl;
#endif

  return theResult;
}

void G4VPartonStringModel::ModelDescription(std::ostream& outFile) const
{
	outFile << GetModelName() << " has no description yet.\n";
}

G4V3DNucleus * G4VPartonStringModel::GetProjectileNucleus() const 
{ return 0;}

