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
// $Id: G4DiffractiveExcitation.cc,v 1.12 2009-08-03 13:14:19 vuzhinsk Exp $
// ------------------------------------------------------------
//      GEANT 4 class implemetation file
//
//      ---------------- G4DiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//        diffractive Excitation used by strings models
//             Take a projectile and a target
//             excite the projectile and target
//  Essential changed by V. Uzhinsky in November - December 2006
//  in order to put it in a correspondence with original FRITIOF
//  model. Variant of FRITIOF with nucleon de-excitation is implemented.
//  Other changes by V.Uzhinsky in May 2007 were introduced to fit
//  meson-nucleon interactions. Additional changes by V. Uzhinsky
//  were introduced in December 2006. They treat diffraction dissociation
//  processes more exactly.
// ---------------------------------------------------------------------


#include "globals.hh"
#include "Randomize.hh"

#include "G4DiffractiveExcitation.hh"
#include "G4LorentzRotation.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh" 
#include "G4VSplitableHadron.hh"
#include "G4ExcitedString.hh"
#include "G4FTFParameters.hh"
#include "G4ParticleTable.hh"
//#include "G4ios.hh"
//#include "UZHI_diffraction.hh"

G4DiffractiveExcitation::G4DiffractiveExcitation()
{
}

// ---------------------------------------------------------------------
G4bool G4DiffractiveExcitation::
  ExciteParticipants(G4VSplitableHadron *projectile, 
                     G4VSplitableHadron *target,
                     G4FTFParameters    *theParameters) const
{
// -------------------- Projectile parameters -----------------------
     G4LorentzVector Pprojectile=projectile->Get4Momentum();

     if(Pprojectile.z() < 0.)
     {
       target->SetStatus(2);
       return false;
     } 

     G4double ProjectileRapidity = Pprojectile.rapidity();

     G4int    ProjectilePDGcode=projectile->GetDefinition()->GetPDGEncoding();
     G4ParticleDefinition * ProjectileDefinition = projectile->GetDefinition();
     G4int    absProjectilePDGcode=std::abs(ProjectilePDGcode);

     G4bool PutOnMassShell(false);
//   G4double M0projectile=projectile->GetDefinition()->GetPDGMass(); // With de-excitation
     G4double M0projectile = Pprojectile.mag();                       // Without de-excitation

     if(M0projectile < projectile->GetDefinition()->GetPDGMass())
     {
      PutOnMassShell=true;
      M0projectile=projectile->GetDefinition()->GetPDGMass();
     }

     G4double M0projectile2 = M0projectile * M0projectile;

     G4double ProjectileDiffStateMinMass=theParameters->GetProjMinDiffMass();
     G4double ProjectileNonDiffStateMinMass=theParameters->GetProjMinNonDiffMass();
     G4double ProbProjectileDiffraction=theParameters->GetProbabilityOfProjDiff();

// -------------------- Target parameters -------------------------
     G4int    TargetPDGcode=target->GetDefinition()->GetPDGEncoding();

     G4LorentzVector Ptarget=target->Get4Momentum();
     G4double M0target = Ptarget.mag();
     G4double TargetRapidity = Ptarget.rapidity();

     if(M0target < target->GetDefinition()->GetPDGMass())
     {
      PutOnMassShell=true;
      M0target=target->GetDefinition()->GetPDGMass();
     }

     G4double M0target2 = M0target * M0target; 
 
     G4double TargetDiffStateMinMass=theParameters->GetTarMinDiffMass();    
     G4double TargetNonDiffStateMinMass=theParameters->GetTarMinNonDiffMass();    
     G4double ProbTargetDiffraction=theParameters->GetProbabilityOfTarDiff();

     G4double AveragePt2=theParameters->GetAveragePt2();

// Kinematical properties of the interactions --------------
     G4LorentzVector Psum;      // 4-momentum in CMS
     Psum=Pprojectile+Ptarget;
     G4double S=Psum.mag2(); 

// Transform momenta to cms and then rotate parallel to z axis;
     G4LorentzRotation toCms(-1*Psum.boostVector());

     G4LorentzVector Ptmp=toCms*Pprojectile;
     if ( Ptmp.pz() <= 0. )
     {
       target->SetStatus(2); 
   // "String" moving backwards in  CMS, abort collision !!
       return false;
     }

     toCms.rotateZ(-1*Ptmp.phi());
     toCms.rotateY(-1*Ptmp.theta());

     G4LorentzRotation toLab(toCms.inverse());

     Pprojectile.transform(toCms);
     Ptarget.transform(toCms);

     G4double PZcms2, PZcms;

     G4double SqrtS=std::sqrt(S);

     if(absProjectilePDGcode > 1000 && SqrtS < 2200*MeV)
     {target->SetStatus(2);  return false;}  // The model cannot work for
                                             // p+p-interactions
                                             // at Plab < 1.3 GeV/c.

     if(( absProjectilePDGcode == 211 || ProjectilePDGcode ==  111) && SqrtS < 1600*MeV)
     {target->SetStatus(2);  return false;}  // The model cannot work for
                                             // Pi+p-interactions
                                             // at Plab < 1. GeV/c.

     if(( (absProjectilePDGcode == 321) || (ProjectilePDGcode == -311)   ||
          (absProjectilePDGcode == 311) || (absProjectilePDGcode == 130) ||
          (absProjectilePDGcode == 310)) && SqrtS < 1600*MeV)
     {target->SetStatus(2);  return false;}  // The model cannot work for
                                             // K+p-interactions
                                             // at Plab < ??? GeV/c.  ???

     PZcms2=(S*S+M0projectile2*M0projectile2+M0target2*M0target2-
             2*S*M0projectile2 - 2*S*M0target2 - 2*M0projectile2*M0target2)
             /4./S;

     if(PZcms2 < 0)
     {target->SetStatus(2);  return false;}   // It can be in an interaction with 
                                              // off-shell nuclear nucleon

     PZcms = std::sqrt(PZcms2);

     if(PutOnMassShell)
     {
      if(Pprojectile.z() > 0.)
      {
       Pprojectile.setPz( PZcms);
       Ptarget.setPz(    -PZcms);
      }
      else
      {
       Pprojectile.setPz(-PZcms);
       Ptarget.setPz(     PZcms);
      };

      Pprojectile.setE(std::sqrt(M0projectile2                  +
                                 Pprojectile.x()*Pprojectile.x()+
                                 Pprojectile.y()*Pprojectile.y()+
                                 PZcms2));
      Ptarget.setE(std::sqrt(M0target2              +
                             Ptarget.x()*Ptarget.x()+
                             Ptarget.y()*Ptarget.y()+
                             PZcms2));
     }

     G4double maxPtSquare; // = PZcms2;

// Charge exchange can be possible for baryons -----------------
//        (G4UniformRand() < 0.5))
//        (G4UniformRand() < 0.5*std::exp(-0.5*(ProjectileRapidity - TargetRapidity))))
//
     if((ProjectilePDGcode != TargetPDGcode) &&
        ((ProjectilePDGcode > 1000) && (TargetPDGcode > 1000)) &&  
        (G4UniformRand() < 1.0*std::exp(-0.5*(ProjectileRapidity - TargetRapidity))))
       {
        projectile->SetDefinition(target->GetDefinition());
        target->SetDefinition(ProjectileDefinition);      
       }
//
// ------------------------------------------------------------------
//  In the case of the projectile pion an absorption is possible ----
     if((absProjectilePDGcode < 1000) &&       // Absorption  Uzhi 7.07.09
        (G4UniformRand() < 1.0*std::exp(-0.5*(ProjectileRapidity - TargetRapidity)))) 
     {

      G4int ProjQ1=  absProjectilePDGcode/ 100;
      G4int ProjQ2= (absProjectilePDGcode %100)/10;
	   
      G4int anti= 1 -2 * ( std::max( ProjQ1, ProjQ2 ) % 2 );
      if (ProjectilePDGcode < 0 ) anti *=-1;
	    
      ProjQ1 *= anti;
      ProjQ2 *= -1 * anti;

      G4int TargQ1 = TargetPDGcode           / 1000;
      G4int TargQ2 = (TargetPDGcode % 1000)  / 100;
      G4int TargQ3 = (TargetPDGcode % 100)   / 10;	    

      G4int CombinedSystemQ1(0), CombinedSystemQ2(0), CombinedSystemQ3(0);

      if(ProjQ1 < 0)                     // Q1 of the projectile is anti-quark
      {
       CombinedSystemQ1 = ProjQ2; 
       if(ProjQ1 + TargQ1 == 0)          // Pr Q1 annihilates with Tr Q1
       {
        CombinedSystemQ2 = TargQ2;
        CombinedSystemQ3 = TargQ3;
       } else 
       {
         if(ProjQ1 + TargQ2 == 0)        // Pr Q1 annihilates with Tr Q2
         {
          CombinedSystemQ2 = TargQ1;
          CombinedSystemQ3 = TargQ3;
         } else
         {
          if(ProjQ1 + TargQ3 == 0)       // Pr Q1 annihilates with Tr Q3
          {
           CombinedSystemQ2 = TargQ1;
           CombinedSystemQ3 = TargQ2;
          }
         }
       }
      } else                             // Q2 of the projectile is anti-quark
      {
       CombinedSystemQ1 = ProjQ1; 
       if(ProjQ2 + TargQ1 == 0)          // Pr Q2 annihilates with Tr Q1
       {
        CombinedSystemQ2 = TargQ2;
        CombinedSystemQ3 = TargQ3;
       } else 
       {
         if(ProjQ2 + TargQ2 == 0)        // Pr Q2 annihilates with Tr Q2
         {
          CombinedSystemQ2 = TargQ1;
          CombinedSystemQ3 = TargQ3;
         } else
         {
          if(ProjQ2 + TargQ3 == 0)       // Pr Q2 annihilates with Tr Q3
          {
           CombinedSystemQ2 = TargQ1;
           CombinedSystemQ3 = TargQ2;
          }
         }
       }     
      };
// Odering of the quarks

      G4int TmpQ;
      if( CombinedSystemQ3 > CombinedSystemQ2 ) 
      {
       TmpQ = CombinedSystemQ2;
       CombinedSystemQ2 = CombinedSystemQ3;
       CombinedSystemQ3 = TmpQ;
      } else
      {
       if( CombinedSystemQ3 > CombinedSystemQ1 )
       {
        TmpQ = CombinedSystemQ1;
        CombinedSystemQ1 = CombinedSystemQ3;
        CombinedSystemQ3 = TmpQ;
       }
      };

      if( CombinedSystemQ2 > CombinedSystemQ1 ) 
      {
       TmpQ = CombinedSystemQ1;
       CombinedSystemQ1 = CombinedSystemQ2;
       CombinedSystemQ2 = TmpQ;
      }; 
//===========================

      G4int CombinedCode = CombinedSystemQ1*1000 + 
                           CombinedSystemQ2* 100 +
                           CombinedSystemQ3*  10 + 2; 
      if( (CombinedSystemQ1 == CombinedSystemQ2) && 
          (CombinedSystemQ1 == CombinedSystemQ3) ) CombinedCode +=2; // For Delta isobar

      projectile->SetDefinition(
       G4ParticleTable::GetParticleTable()->FindParticle(CombinedCode)); 

// Calculation of the creation time ---------------------
      projectile->SetTimeOfCreation(target->GetTimeOfCreation());
      projectile->SetPosition(target->GetPosition());
// Creation time and position of target nucleon were determined at
// ReggeonCascade() of G4FTFModel
// ------------------------------------------------------

      projectile->Set4Momentum(Psum);

      theParameters->SetProjMinDiffMass(theParameters->GetTarMinDiffMass()/GeV);
      theParameters->SetProjMinNonDiffMass(theParameters->GetTarMinNonDiffMass()/GeV);
      theParameters->SetProbabilityOfProjDiff(theParameters->GetProbabilityOfTarDiff());      

      projectile->SetStatus(3);
      target->SetStatus(3);

      projectile->IncrementCollisionCount(1);
      target->IncrementCollisionCount(1);

      return true;
     };        // End Absorption =============================

// ------------------- Ordinary job of the Fritiof model -----------

     G4double ProbOfDiffraction=ProbProjectileDiffraction +
                                ProbTargetDiffraction;

     if(ProbOfDiffraction!=0.)
     {
      ProbProjectileDiffraction/=ProbOfDiffraction;
     }
     else
     {
      ProbProjectileDiffraction=0.;
     }

     G4double ProjectileDiffStateMinMass2    = ProjectileDiffStateMinMass    *
                                               ProjectileDiffStateMinMass;
     G4double ProjectileNonDiffStateMinMass2 = ProjectileNonDiffStateMinMass *
                                               ProjectileNonDiffStateMinMass;

     G4double TargetDiffStateMinMass2        = TargetDiffStateMinMass        *
                                               TargetDiffStateMinMass;
     G4double TargetNonDiffStateMinMass2     = TargetNonDiffStateMinMass     *
                                               TargetNonDiffStateMinMass;


     G4double Pt2;
     G4double ProjMassT2, ProjMassT;
     G4double TargMassT2, TargMassT;
     G4double PMinusMin, PMinusMax;
//   G4double PPlusMin , PPlusMax;
     G4double TPlusMin , TPlusMax;
     G4double PMinusNew, PPlusNew, TPlusNew, TMinusNew;

     G4LorentzVector Qmomentum;
     G4double Qminus, Qplus;

     G4int whilecount=0;
//   Choose a process ---------------------------

     if(G4UniformRand() < ProbOfDiffraction)
       {
        if(G4UniformRand() < ProbProjectileDiffraction)
        { //-------- projectile diffraction ---------------
         do {
//             Generate pt
//             if (whilecount++ >= 500 && (whilecount%100)==0)
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants possibly looping"
//	   	 << ", loop count/ maxPtSquare : "
//           	 << whilecount << " / " << maxPtSquare << G4endl;
             if (whilecount > 1000 )
             {
              Qmomentum=G4LorentzVector(0.,0.,0.,0.);
              target->SetStatus(2);  return false;    //  Ignore this interaction
             };
// --------------- Check that the interaction is possible -----------
             ProjMassT2=ProjectileDiffStateMinMass2;
             ProjMassT =ProjectileDiffStateMinMass;

             TargMassT2=M0target2;
             TargMassT =M0target;

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                     2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                    /4./S;

             if(PZcms2 < 0 ) 
             {
               target->SetStatus(2);  
               return false;
             }
             maxPtSquare=PZcms2;

             Qmomentum=G4LorentzVector(GaussianPt(AveragePt2,maxPtSquare),0);
             Pt2=G4ThreeVector(Qmomentum.vect()).mag2();

             ProjMassT2=ProjectileDiffStateMinMass2+Pt2;
             ProjMassT =std::sqrt(ProjMassT2);

             TargMassT2=M0target2+Pt2;
             TargMassT =std::sqrt(TargMassT2);

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                     2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                    /4./S;

             if(PZcms2 < 0 ) continue;
             PZcms =std::sqrt(PZcms2);

             PMinusMin=std::sqrt(ProjMassT2+PZcms2)-PZcms;
             PMinusMax=SqrtS-TargMassT;

             PMinusNew=ChooseP(PMinusMin, PMinusMax);
// PMinusNew=1./sqrt(1./PMinusMin-G4UniformRand()*(1./PMinusMin-1./PMinusMax));

             TMinusNew=SqrtS-PMinusNew;
             Qminus=Ptarget.minus()-TMinusNew;
             TPlusNew=TargMassT2/TMinusNew;
             Qplus=Ptarget.plus()-TPlusNew;

             Qmomentum.setPz( (Qplus-Qminus)/2 );
             Qmomentum.setE(  (Qplus+Qminus)/2 );
          } while (
((Pprojectile+Qmomentum).mag2() <  ProjectileDiffStateMinMass2) ||  //No without excitation
((Ptarget    -Qmomentum).mag2() <  M0target2                  ));
        }
        else
        { // -------------- Target diffraction ----------------
         do {
//             Generate pt
//             if (whilecount++ >= 500 && (whilecount%100)==0)
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants possibly looping"
//	   	 << ", loop count/ maxPtSquare : "
//           	 << whilecount << " / " << maxPtSquare << G4endl;
             if (whilecount > 1000 )
             {
              Qmomentum=G4LorentzVector(0.,0.,0.,0.);
              target->SetStatus(2);  return false;    //  Ignore this interaction
             };
// --------------- Check that the interaction is possible -----------
             ProjMassT2=M0projectile2;
             ProjMassT =M0projectile;

             TargMassT2=TargetDiffStateMinMass2;
             TargMassT =TargetDiffStateMinMass;

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                     2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                    /4./S;

             if(PZcms2 < 0 ) 
             {
               target->SetStatus(2);  
               return false;
             }
             maxPtSquare=PZcms2;

             Qmomentum=G4LorentzVector(GaussianPt(AveragePt2,maxPtSquare),0);
             Pt2=G4ThreeVector(Qmomentum.vect()).mag2();

             ProjMassT2=M0projectile2+Pt2;
             ProjMassT =std::sqrt(ProjMassT2);

             TargMassT2=TargetDiffStateMinMass2+Pt2;
             TargMassT =std::sqrt(TargMassT2);

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                     2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                    /4./S;

             if(PZcms2 < 0 ) continue;
             PZcms =std::sqrt(PZcms2);

             TPlusMin=std::sqrt(TargMassT2+PZcms2)-PZcms;
             TPlusMax=SqrtS-ProjMassT;

             TPlusNew=ChooseP(TPlusMin, TPlusMax);

//TPlusNew=TPlusMax;

             PPlusNew=SqrtS-TPlusNew;
             Qplus=PPlusNew-Pprojectile.plus();
             PMinusNew=ProjMassT2/PPlusNew;
             Qminus=PMinusNew-Pprojectile.minus();

             Qmomentum.setPz( (Qplus-Qminus)/2 );
             Qmomentum.setE(  (Qplus+Qminus)/2 );

          } while (
 ((Pprojectile+Qmomentum).mag2() <  M0projectile2          ) ||  //No without excitation
 ((Ptarget    -Qmomentum).mag2() <  TargetDiffStateMinMass2));
         }
        }
        else  //----------- Non-diffraction process ------------
        {
         do {
//             Generate pt
//             if (whilecount++ >= 500 && (whilecount%100)==0)
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants possibly looping"
//	   	 << ", loop count/ maxPtSquare : "
//           	 << whilecount << " / " << maxPtSquare << G4endl;
             if (whilecount > 1000 )
             {
              Qmomentum=G4LorentzVector(0.,0.,0.,0.);
              target->SetStatus(2);  return false;    //  Ignore this interaction
             };
// --------------- Check that the interaction is possible -----------
             ProjMassT2=ProjectileNonDiffStateMinMass2;
             ProjMassT =ProjectileNonDiffStateMinMass;

             TargMassT2=TargetNonDiffStateMinMass2;
             TargMassT =TargetNonDiffStateMinMass;

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                    2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                   /4./S;

             if(PZcms2 < 0 ) 
             {
               target->SetStatus(2);  
               return false;
             }
             maxPtSquare=PZcms2;
             Qmomentum=G4LorentzVector(GaussianPt(AveragePt2,maxPtSquare),0);
             Pt2=G4ThreeVector(Qmomentum.vect()).mag2();

             ProjMassT2=ProjectileNonDiffStateMinMass2+Pt2;
             ProjMassT =std::sqrt(ProjMassT2);

             TargMassT2=TargetNonDiffStateMinMass2+Pt2;
             TargMassT =std::sqrt(TargMassT2);

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                    2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                   /4./S;

             if(PZcms2 < 0 ) continue;
             PZcms =std::sqrt(PZcms2);

             PMinusMin=std::sqrt(ProjMassT2+PZcms2)-PZcms;
             PMinusMax=SqrtS-TargMassT;

             PMinusNew=ChooseP(PMinusMin, PMinusMax);

             Qminus=PMinusNew-Pprojectile.minus();

             TPlusMin=std::sqrt(TargMassT2+PZcms2)-PZcms;
//           TPlusMax=SqrtS-PMinusNew;                      
             TPlusMax=SqrtS-ProjMassT;      

             TPlusNew=ChooseP(TPlusMin, TPlusMax);

             Qplus=-(TPlusNew-Ptarget.plus());

             Qmomentum.setPz( (Qplus-Qminus)/2 );
             Qmomentum.setE(  (Qplus+Qminus)/2 );

       } while (
 ((Pprojectile+Qmomentum).mag2() <  ProjectileNonDiffStateMinMass2) || //No double Diffraction
 ((Ptarget    -Qmomentum).mag2() <  TargetNonDiffStateMinMass2    ));
         }

	   Pprojectile += Qmomentum;
	   Ptarget     -= Qmomentum;

// Transform back and update SplitableHadron Participant.
	   Pprojectile.transform(toLab);
	   Ptarget.transform(toLab);

// Calculation of the creation time ---------------------
      projectile->SetTimeOfCreation(target->GetTimeOfCreation());
      projectile->SetPosition(target->GetPosition());
// Creation time and position of target nucleon were determined at
// ReggeonCascade() of G4FTFModel
// ------------------------------------------------------

	   projectile->Set4Momentum(Pprojectile);
	   target->Set4Momentum(Ptarget);

           projectile->IncrementCollisionCount(1);
           target->IncrementCollisionCount(1);

	   return true;
}

// ---------------------------------------------------------------------
G4ExcitedString * G4DiffractiveExcitation::
           String(G4VSplitableHadron * hadron, G4bool isProjectile) const
{
	hadron->SplitUp();
	G4Parton *start= hadron->GetNextParton();
	if ( start==NULL)
	{ G4cout << " G4FTFModel::String() Error:No start parton found"<< G4endl;
	  return NULL;
	}
	G4Parton *end  = hadron->GetNextParton();
	if ( end==NULL)
	{ G4cout << " G4FTFModel::String() Error:No end parton found"<< G4endl;
	  return NULL;
	}

	G4ExcitedString * string;
	if ( isProjectile )
	{
		string= new G4ExcitedString(end,start, +1);
	} else {
		string= new G4ExcitedString(start,end, -1);
	}

	string->SetTimeOfCreation(hadron->GetTimeOfCreation());
	string->SetPosition(hadron->GetPosition());

// momenta of string ends
//
        G4double Momentum=hadron->Get4Momentum().vect().mag();
        G4double  Plus=hadron->Get4Momentum().e() + Momentum;
        G4double Minus=hadron->Get4Momentum().e() - Momentum;

        G4ThreeVector tmp;
        if(Momentum > 0.) 
        {
         tmp.set(hadron->Get4Momentum().px(),
                 hadron->Get4Momentum().py(),
                 hadron->Get4Momentum().pz());
         tmp/=Momentum;
        }
        else
        {
          tmp.set(0.,0.,1.);
        };

        G4LorentzVector Pstart(tmp,0.);
        G4LorentzVector   Pend(tmp,0.);

        if(isProjectile)
        {
         Pstart*=(-1.)*Minus/2.;
         Pend  *=(+1.)*Plus /2.;
        } 
        else
        {
          Pstart*=(+1.)*Plus/2.;
          Pend  *=(-1.)*Minus/2.;
        };

        Momentum=-Pstart.mag();
        Pstart.setT(Momentum);  // It is assumed that quark has m=0.

        Momentum=-Pend.mag();
        Pend.setT(Momentum);    // It is assumed that di-quark has m=0.

	start->Set4Momentum(Pstart);
	end->Set4Momentum(Pend);

#ifdef G4_FTFDEBUG
	G4cout << " generated string flavors          " 
               << start->GetPDGcode() << " / " 
               << end->GetPDGcode()                    << G4endl;
	G4cout << " generated string momenta:   quark " 
               << start->Get4Momentum() << "mass : " 
               <<start->Get4Momentum().mag()           << G4endl;
	G4cout << " generated string momenta: Diquark " 
               << end ->Get4Momentum() 
               << "mass : " <<end->Get4Momentum().mag()<< G4endl;
	G4cout << " sum of ends                       " << Pstart+Pend << G4endl;
	G4cout << " Original                          " << hadron->Get4Momentum() << G4endl;
#endif

	return string;
}


// --------- private methods ----------------------

// ---------------------------------------------------------------------
G4double G4DiffractiveExcitation::ChooseP(G4double Pmin, G4double Pmax) const
{
// choose an x between Xmin and Xmax with P(x) ~ 1/x
//  to be improved...

	G4double range=Pmax-Pmin;                    

	if ( Pmin <= 0. || range <=0. )
	{
		G4cout << " Pmin, range : " << Pmin << " , " << range << G4endl;
		throw G4HadronicException(__FILE__, __LINE__, "G4DiffractiveExcitation::ChooseP : Invalid arguments ");
	}

	G4double P=Pmin * std::pow(Pmax/Pmin,G4UniformRand()); 
	return P;
}

// ---------------------------------------------------------------------
G4ThreeVector G4DiffractiveExcitation::GaussianPt(G4double AveragePt2, 
                                                  G4double maxPtSquare) const
{            //  @@ this method is used in FTFModel as well. Should go somewhere common!

	G4double Pt2;
        Pt2 = -AveragePt2 * std::log(1. + G4UniformRand() * 
                           (std::exp(-maxPtSquare/AveragePt2)-1.));

	G4double Pt=std::sqrt(Pt2);
	G4double phi=G4UniformRand() * twopi;
	return G4ThreeVector (Pt*std::cos(phi), Pt*std::sin(phi), 0.);
}

// ---------------------------------------------------------------------
G4DiffractiveExcitation::G4DiffractiveExcitation(const G4DiffractiveExcitation &)
{
	throw G4HadronicException(__FILE__, __LINE__, "G4DiffractiveExcitation copy contructor not meant to be called");
}


G4DiffractiveExcitation::~G4DiffractiveExcitation()
{
}


const G4DiffractiveExcitation & G4DiffractiveExcitation::operator=(const G4DiffractiveExcitation &)
{
	throw G4HadronicException(__FILE__, __LINE__, "G4DiffractiveExcitation = operator meant to be called");
	return *this;
}


int G4DiffractiveExcitation::operator==(const G4DiffractiveExcitation &) const
{
	throw G4HadronicException(__FILE__, __LINE__, "G4DiffractiveExcitation == operator meant to be called");
	return false;
}

int G4DiffractiveExcitation::operator!=(const G4DiffractiveExcitation &) const
{
	throw G4HadronicException(__FILE__, __LINE__, "G4DiffractiveExcitation != operator meant to be called");
	return true;
}
