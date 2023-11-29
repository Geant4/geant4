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
//#define debug_heavyHadrons
// ------------------------------------------------------------

#include "G4VPartonStringModel.hh"
#include "G4ios.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

G4VPartonStringModel::G4VPartonStringModel(const G4String& modelName)
    : G4VHighEnergyGenerator(modelName),
      stringFragmentationModel(nullptr)
{
  //  Make shure Shotrylived particles are constructed.
  //  VI: should not instantiate particles by any model
  /*
  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();
  */
}

G4VPartonStringModel::~G4VPartonStringModel()
{
}

G4KineticTrackVector * G4VPartonStringModel::Scatter(const G4Nucleus &theNucleus, 
                                                const G4DynamicParticle &aPrimary)
{  
  G4ExcitedStringVector * strings = nullptr;
  G4DynamicParticle thePrimary=aPrimary;
  G4LorentzVector SumStringMom(0.,0.,0.,0.);
  G4KineticTrackVector * theResult = 0;
  G4Nucleon * theNuclNucleon(nullptr);

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
  if(GetProjectileNucleus()) {
    Bsum -= thePrimary.GetDefinition()->GetBaryonNumber();
    Qsum -= thePrimary.GetDefinition()->GetPDGCharge();
  }
  G4int QsumSec(0), BsumSec(0);
  G4LorentzVector SumPsecondr(0.,0.,0.,0.);
  #endif

  G4LorentzRotation toZ;
  G4LorentzVector Ptmp=thePrimary.Get4Momentum();
  toZ.rotateZ(-1*Ptmp.phi());
  toZ.rotateY(-1*Ptmp.theta());
  thePrimary.Set4Momentum(toZ*Ptmp);
  G4LorentzRotation toLab(toZ.inverse());

  G4bool Success=true;
  G4int attempts = 0, maxAttempts=1000;
  do
  {
    if (attempts++ > maxAttempts ) 
    {
      Init(theNucleus,thePrimary);  // To put a nucleus into ground state
                                             // But marks of hitted nucleons are left. They must be erased.
      G4V3DNucleus * ResNucleus = GetWoundedNucleus(); 
      theNuclNucleon = ResNucleus->StartLoop() ? ResNucleus->GetNextNucleon() : nullptr;
      while( theNuclNucleon )
      {
        if(theNuclNucleon->AreYouHit()) theNuclNucleon->Hit(nullptr);
        theNuclNucleon = ResNucleus->GetNextNucleon();
      }

      G4V3DNucleus * ProjResNucleus = GetProjectileNucleus();
      if(ProjResNucleus != 0)
      {
        theNuclNucleon = ProjResNucleus->StartLoop() ? ProjResNucleus->GetNextNucleon() : nullptr;
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
      G4KineticTrack* Hadron = new G4KineticTrack(aPrimary.GetParticleDefinition(), 0.,
                                                  Position, aPrimary.Get4Momentum());
      if(theResult == nullptr) theResult = new G4KineticTrackVector();
      theResult->push_back(Hadron);
      return theResult;
    }

    Success=true;

    Init(theNucleus,thePrimary);

    strings = GetStrings();

    if (strings->empty()) { Success=false; continue; }

    // G4double stringEnergy(0);
    SumStringMom=G4LorentzVector(0.,0.,0.,0.);

    #ifdef debug_PartonStringModel
    G4cout<<"------------ Parton-String model: Number of produced strings ---- "<<strings->size()<<G4endl;
    #endif

    #ifdef debug_heavyHadrons
    // Check charm and bottom numbers of the projectile:
    G4int count_charm_projectile  = thePrimary.GetDefinition()->GetQuarkContent( 4 ) -
                                    thePrimary.GetDefinition()->GetAntiQuarkContent( 4 );
    G4int count_bottom_projectile = thePrimary.GetDefinition()->GetQuarkContent( 5 ) -
                                    thePrimary.GetDefinition()->GetAntiQuarkContent( 5 );
    G4int count_charm_strings = 0, count_bottom_strings = 0;
    G4int count_charm_hadrons = 0, count_bottom_hadrons = 0;
    #endif
    
    for ( unsigned int astring=0; astring < strings->size(); astring++)
    {
      //    rotate string to lab frame, models have it aligned to z
      if((*strings)[astring]->IsExcited())
      {
        // stringEnergy += (*strings)[astring]->GetLeftParton()->Get4Momentum().t();
        // stringEnergy += (*strings)[astring]->GetRightParton()->Get4Momentum().t();
        (*strings)[astring]->LorentzRotate(toLab);
        SumStringMom+=(*strings)[astring]->Get4Momentum();
        #ifdef debug_PartonStringModel
        G4cout<<"String No "<<astring+1<<" "<<(*strings)[astring]->Get4Momentum()<<" "
                            <<(*strings)[astring]->Get4Momentum().mag()
              <<" Partons   "<<(*strings)[astring]->GetLeftParton()->GetDefinition()->GetPDGEncoding()
              <<"          "<<(*strings)[astring]->GetRightParton()->GetDefinition()->GetPDGEncoding()<<G4endl;
        #endif

        #ifdef debug_heavyHadrons
	G4int left_charm =      (*strings)[astring]->GetLeftParton()->GetDefinition()->GetQuarkContent( 4 );
	G4int left_anticharm =  (*strings)[astring]->GetLeftParton()->GetDefinition()->GetAntiQuarkContent( 4 );
	G4int right_charm =     (*strings)[astring]->GetRightParton()->GetDefinition()->GetQuarkContent( 4 );
	G4int right_anticharm = (*strings)[astring]->GetRightParton()->GetDefinition()->GetAntiQuarkContent( 4 );
	G4int left_bottom =      (*strings)[astring]->GetLeftParton()->GetDefinition()->GetQuarkContent( 5 );
	G4int left_antibottom =  (*strings)[astring]->GetLeftParton()->GetDefinition()->GetAntiQuarkContent( 5 );
	G4int right_bottom =     (*strings)[astring]->GetRightParton()->GetDefinition()->GetQuarkContent( 5 );
	G4int right_antibottom = (*strings)[astring]->GetRightParton()->GetDefinition()->GetAntiQuarkContent( 5 );
	if ( left_charm  != 0  ||  left_anticharm  != 0  ||  right_charm  != 0  ||  right_anticharm  != 0  ||
	     left_bottom != 0  ||  left_antibottom != 0  ||  right_bottom != 0  ||  right_antibottom != 0 ) {
	  count_charm_strings  += left_charm  - left_anticharm  + right_charm  - right_anticharm;
	  count_bottom_strings += left_bottom - left_antibottom + right_bottom - right_antibottom;
	  G4cout << "G4VPartonStringModel::Scatter : string #" << astring << " ("
		 << (*strings)[astring]->GetLeftParton()->GetDefinition()->GetParticleName() << " , "
		 << (*strings)[astring]->GetRightParton()->GetDefinition()->GetParticleName() << ")" << G4endl;
	}
        #endif	
      }
      else
      {
        // stringEnergy += (*strings)[astring]->GetKineticTrack()->Get4Momentum().t();
        (*strings)[astring]->LorentzRotate(toLab);
        SumStringMom+=(*strings)[astring]->GetKineticTrack()->Get4Momentum();
        #ifdef debug_PartonStringModel
        G4cout<<"A track No "<<astring+1<<" "
              <<(*strings)[astring]->GetKineticTrack()->Get4Momentum()<<" "
              <<(*strings)[astring]->GetKineticTrack()->Get4Momentum().mag()<<" "
              <<(*strings)[astring]->GetKineticTrack()->GetDefinition()->GetParticleName()<<G4endl;
        #endif

        #ifdef debug_heavyHadrons
	G4int charm =      (*strings)[astring]->GetKineticTrack()->GetDefinition()->GetQuarkContent( 4 );
	G4int anticharm =  (*strings)[astring]->GetKineticTrack()->GetDefinition()->GetAntiQuarkContent( 4 );
	G4int bottom =     (*strings)[astring]->GetKineticTrack()->GetDefinition()->GetQuarkContent( 5 );
	G4int antibottom = (*strings)[astring]->GetKineticTrack()->GetDefinition()->GetAntiQuarkContent( 5 );
        if ( charm != 0  ||  anticharm != 0  ||  bottom != 0  || antibottom != 0 ) {
	  count_charm_strings +=  charm  - anticharm;
	  count_bottom_strings += bottom - antibottom;
	  G4cout << "G4VPartonStringModel::Scatter : track #" << astring << "\t"
                 << (*strings)[astring]->GetKineticTrack()->GetDefinition()->GetParticleName() << G4endl;
	}
        #endif
      }
    }

    #ifdef debug_heavyHadrons
    if ( count_charm_projectile != count_charm_strings ) {
      G4cout << "G4VPartonStringModel::Scatter : CHARM VIOLATION in String formation ! #projectile="
	     << count_charm_projectile << " ; #strings=" << count_charm_strings << G4endl;
    }
    if ( count_bottom_projectile != count_bottom_strings ) {
      G4cout << "G4VPartonStringModel::Scatter : BOTTOM VIOLATION in String formation ! #projectile="
	     << count_bottom_projectile << " ; #strings=" << count_bottom_strings << G4endl;
    }
    #endif
    
    #ifdef debug_PartonStringModel
    G4cout<<G4endl<<"SumString4Mom "<<SumStringMom<<G4endl;
    G4LorentzVector TargetResidual4Momentum(0.,0.,0.,0.);
    G4LorentzVector ProjectileResidual4Momentum(0.,0.,0.,0.);
    G4int hitsT(0), charged_hitsT(0);
    G4int hitsP(0), charged_hitsP(0);
    G4double ExcitationEt(0.), ExcitationEp(0.);
    #endif

    // We assume that the target nucleus is never a hypernucleus, whereas
    // the projectile nucleus can be a light hypernucleus or anti-hypernucleus.
    
    G4V3DNucleus * ProjResNucleus = GetProjectileNucleus();

    G4int numberProtonProjectileResidual( 0 ), numberNeutronProjectileResidual( 0 );
    G4int numberLambdaProjectileResidual( 0 );
    if(ProjResNucleus != 0)
    {
      theNuclNucleon = ProjResNucleus->StartLoop() ? ProjResNucleus->GetNextNucleon() : nullptr;
      G4int numberProtonProjectileHits( 0 ), numberNeutronProjectileHits( 0 );
      G4int numberLambdaProjectileHits( 0 );
      while( theNuclNucleon )
      {
        if(theNuclNucleon->AreYouHit())
        {
          G4LorentzVector tmp=toLab*theNuclNucleon->Get4Momentum();
	  const G4ParticleDefinition* def = theNuclNucleon->GetDefinition();
          #ifdef debug_PartonStringModel
          ProjectileResidual4Momentum += tmp;
          hitsP++;
          if ( def == G4Proton::Definition() || def == G4AntiProton::Definition() ) ++charged_hitsP;
          ExcitationEp +=theNuclNucleon->GetBindingEnergy();
          #endif
          theNuclNucleon->SetMomentum(tmp);
          if ( def == G4Proton::Definition()  || def == G4AntiProton::Definition() )  ++numberProtonProjectileHits;
          if ( def == G4Neutron::Definition() || def == G4AntiNeutron::Definition() ) ++numberNeutronProjectileHits;
          if ( def == G4Lambda::Definition()  || def == G4AntiLambda::Definition() )  ++numberLambdaProjectileHits;
        }
        theNuclNucleon = ProjResNucleus->GetNextNucleon();
      }
      G4int numberLambdaProjectile = 0;
      if ( thePrimary.GetDefinition()->IsHypernucleus() ) {
	numberLambdaProjectile = thePrimary.GetDefinition()->GetNumberOfLambdasInHypernucleus();
      } else if ( thePrimary.GetDefinition()->IsAntiHypernucleus() ) {
	numberLambdaProjectile = thePrimary.GetDefinition()->GetNumberOfAntiLambdasInAntiHypernucleus();
      }
      #ifdef debug_PartonStringModel
      G4cout<<"Projectile residual A, Z (numberOfLambdasOrAntiLambdas) and E* "
            <<thePrimary.GetDefinition()->GetBaryonNumber() - hitsP<<" "
            <<thePrimary.GetDefinition()->GetPDGCharge()    - charged_hitsP<<" ("
	    << numberLambdaProjectile - numberLambdaProjectileHits << ") "
            <<ExcitationEp<<G4endl;
      G4cout<<"Projectile residual 4 momentum  "<<ProjectileResidual4Momentum<<G4endl;
      #endif
      numberProtonProjectileResidual = std::max( std::abs( G4int( thePrimary.GetDefinition()->GetPDGCharge() ) ) -
						 numberProtonProjectileHits, 0 );
      numberLambdaProjectileResidual = std::max( numberLambdaProjectile - numberLambdaProjectileHits, 0 );
      numberNeutronProjectileResidual = std::max( std::abs( thePrimary.GetDefinition()->GetBaryonNumber() ) -
                                                  std::abs( G4int( thePrimary.GetDefinition()->GetPDGCharge() ) ) -
						  numberLambdaProjectile - numberNeutronProjectileHits, 0 );
    }

    G4V3DNucleus * ResNucleus = GetWoundedNucleus(); 

    // loop over wounded nucleus
    theNuclNucleon = ResNucleus->StartLoop() ? ResNucleus->GetNextNucleon() : nullptr;
    G4int numberProtonTargetHits( 0 ), numberNeutronTargetHits( 0 );
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
        if ( theNuclNucleon->GetDefinition() == G4Proton::Proton() )   ++numberProtonTargetHits;
        if ( theNuclNucleon->GetDefinition() == G4Neutron::Neutron() ) ++numberNeutronTargetHits;
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

    // Re-sample in the case of unphysical nuclear residual:
    // 1 (H), 2 (2He), and 3 (3Li) protons alone without neutrons can exist, but not more;
    // no bound states of 2 or more neutrons without protons can exist.
    G4int numberProtonTargetResidual = theNucleus.GetZ_asInt() - numberProtonTargetHits;
    G4int numberNeutronTargetResidual = theNucleus.GetA_asInt() - theNucleus.GetZ_asInt() - numberNeutronTargetHits;
    G4bool unphysicalResidual = false;
    if ( ( numberProtonTargetResidual > 3   &&  numberNeutronTargetResidual == 0 ) ||
         ( numberProtonTargetResidual == 0  &&  numberNeutronTargetResidual > 1  ) ) {
      unphysicalResidual = true;
      //G4cout << "***UNPHYSICAL TARGET RESIDUAL*** Z=" << numberProtonTargetResidual 
      //       << " ; N=" << numberNeutronTargetResidual;
    }
    // The projectile residual can be a hypernucleus or anti-hypernucleus:
    // only the following combinations are currently allowed in Geant4:
    // p-n-lambda (hypertriton), p-n-n-lambda (hyperH4), p-p-n-lambda (hyperAlpha),
    // p-p-n-n-lambda (hyperHe5), n-n-lambda-lambda (doublehyperdoubleneutron),
    // p-n-lambda-lambda (doubleHyperH4)
    if ( ( numberProtonProjectileResidual > 3   &&  numberNeutronProjectileResidual == 0 ) ||
         ( numberProtonProjectileResidual == 0  &&  numberNeutronProjectileResidual > 1  &&
	   numberLambdaProjectileResidual == 0 ) ||
	 ( numberProtonProjectileResidual == 0  &&  numberNeutronProjectileResidual <= 1  &&
	   numberLambdaProjectileResidual > 0 ) ||
	 ( numberProtonProjectileResidual == 0  &&  numberNeutronProjectileResidual > 2  &&
	   numberLambdaProjectileResidual > 0 ) ||
	 ( numberLambdaProjectileResidual > 2 ) ||
	 ( numberProtonProjectileResidual > 0  &&  numberNeutronProjectileResidual == 0  &&
	   numberLambdaProjectileResidual > 0 ) ||
	 ( numberProtonProjectileResidual > 1  &&  numberNeutronProjectileResidual > 1  &&
	   numberLambdaProjectileResidual > 1 )
	 ) {
      unphysicalResidual = true;
      //G4cout << "***UNPHYSICAL PROJECTILE RESIDUAL*** Z=" << numberProtonProjectileResidual
      //       << " ; N=" << numberNeutronProjectileResidual
      //       << " ; L=" << numberLambdaProjectileResidual;
    }
    if ( unphysicalResidual ) {
      //G4cout << " -> REJECTING COLLISION because of unphysical residual !" << G4endl;
      Success = false;
      continue;
    }

    //=========================================================================================
    //                              Fragment strings
    #ifdef debug_PartonStringModel
    G4cout<<"---------------- Attempt to fragment strings ------------- "<<attempts<<G4endl;
    #endif

    G4double InvMass=SumStringMom.mag();
    G4double SumMass(0.); 

    #ifdef debug_PartonStringModel
    QsumSec=0; BsumSec=0;
    SumPsecondr=G4LorentzVector(0.,0.,0.,0.);
    #endif

    if(theResult != nullptr)
    {
      std::for_each(theResult->begin(), theResult->end(), DeleteKineticTrack());
      delete theResult;
    }

    theResult = stringFragmentationModel->FragmentStrings(strings);

    #ifdef debug_PartonStringModel
    G4cout<<"String fragmentation success (OK - #0, Not OK - 0) : "<<theResult<<G4endl;
    #endif

    if(theResult == 0) {Success=false; continue;}

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

      #ifdef debug_heavyHadrons
      G4int charm =      (*theResult)[i]->GetDefinition()->GetQuarkContent( 4 );
      G4int anticharm =  (*theResult)[i]->GetDefinition()->GetAntiQuarkContent( 4 );
      G4int bottom =     (*theResult)[i]->GetDefinition()->GetQuarkContent( 5 );
      G4int antibottom = (*theResult)[i]->GetDefinition()->GetAntiQuarkContent( 5 );
      if ( charm != 0  ||  anticharm != 0  ||  bottom != 0  || antibottom != 0 ) {
        count_charm_hadrons +=  charm - anticharm;          
        count_bottom_hadrons += bottom - antibottom;
	G4cout << "G4VPartonStringModel::Scatter : hadron #" << i << "\t"
               << (*theResult)[i]->GetDefinition()->GetParticleName() << G4endl;
      }
      #endif  
    }

    #ifdef debug_heavyHadrons
    if ( count_charm_projectile != count_charm_hadrons ) {
      G4cout << "G4VPartonStringModel::Scatter : CHARM VIOLATION in String hadronization ! #projectile="
	     << count_charm_projectile << " ; #hadrons=" << count_charm_hadrons << G4endl;
    }
    if ( count_bottom_projectile != count_bottom_hadrons ) {
      G4cout << "G4VPartonStringModel::Scatter : BOTTOM VIOLATION in String hadronization ! #projectile="
	     << count_bottom_projectile << " ; #hadrons=" << count_bottom_hadrons << G4endl;
    }
    #endif
    
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

    if((SumMass > InvMass)||(SumMass == 0.)) {Success=false;}
    std::for_each(strings->begin(), strings->end(), DeleteString() );
    delete strings;

  } while(!Success);

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
{ return nullptr; }

