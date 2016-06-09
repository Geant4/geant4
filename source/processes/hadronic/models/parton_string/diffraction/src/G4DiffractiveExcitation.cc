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
// $Id: G4DiffractiveExcitation.cc,v 1.7 2008/12/18 13:01:58 gunter Exp $
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
#include "G4FTFParameters.hh"                            // Uzhi 19.04.08
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
     G4bool PutOnMassShell=0;

// -------------------- Projectile parameters -----------------------

     G4LorentzVector Pprojectile=projectile->Get4Momentum();
//   G4double M0projectile=projectile->GetDefinition()->GetPDGMass(); // With de-excitation
     G4double M0projectile = Pprojectile.mag();                       // Without de-excitation
/*
G4cout<<"ExciteParticipants-------------------"<<G4endl;
G4cout<<"Mom "<<Pprojectile<<" mass "<<M0projectile<<G4endl;
*/
     if(M0projectile < projectile->GetDefinition()->GetPDGMass())
     {
      PutOnMassShell=1;
      M0projectile=projectile->GetDefinition()->GetPDGMass();
     }

     G4double M0projectile2 = M0projectile * M0projectile;

     G4int    PDGcode=projectile->GetDefinition()->GetPDGEncoding();
     G4int    absPDGcode=std::abs(PDGcode);

     G4double ProjectileDiffStateMinMass=theParameters->GetProjMinDiffMass();
     G4double ProjectileNonDiffStateMinMass=theParameters->GetProjMinNonDiffMass();
     G4double ProbProjectileDiffraction=theParameters->GetProbabilityOfProjDiff();
/*
G4cout<<ProjectileDiffStateMinMass<<" "<<ProjectileNonDiffStateMinMass<<" "<<ProbProjectileDiffraction<<G4endl;
*/
// -------------------- Target paraExciteParticipantsmeters -------------------------
     G4LorentzVector Ptarget=target->Get4Momentum();
     G4double M0target = Ptarget.mag();

//G4cout<<"Mom "<<Ptarget<<" mass "<<M0target<<G4endl;

     if(M0target < target->GetDefinition()->GetPDGMass())
     {
      PutOnMassShell=1;
      M0target=target->GetDefinition()->GetPDGMass();
     }

     G4double M0target2 = M0target * M0target;             //Ptarget.mag2();
                                                           // for AA-inter.
     G4double TargetDiffStateMinMass=theParameters->GetTarMinDiffMass();    
     G4double TargetNonDiffStateMinMass=theParameters->GetTarMinNonDiffMass();    
     G4double ProbTargetDiffraction=theParameters->GetProbabilityOfTarDiff();
/*
G4cout<<TargetDiffStateMinMass<<" "<<TargetNonDiffStateMinMass<<" "<<ProbTargetDiffraction<<G4endl;
*/
     G4double AveragePt2=theParameters->GetAveragePt2();

// Kinematical properties of the interactions --------------
     G4LorentzVector Psum;      // 4-momentum in CMS
     Psum=Pprojectile+Ptarget;
     G4double S=Psum.mag2(); 

//G4cout<<" sqrt(s) "<<std::sqrt(S)<<G4endl;

// ------------------------------------------------------------------

//ProbProjectileDiffraction=1.;
//ProbTargetDiffraction    =1.;
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
//   ProbTargetDiffraction    /=ProbOfDiffraction;

//G4cout<<"ProbOfDiffraction "<<ProbOfDiffraction<<"ProbProjectileDiffraction "<<ProbProjectileDiffraction<<G4endl;   // Vova

     G4double ProjectileDiffStateMinMass2    = ProjectileDiffStateMinMass    *
                                               ProjectileDiffStateMinMass;
     G4double ProjectileNonDiffStateMinMass2 = ProjectileNonDiffStateMinMass *
                                               ProjectileNonDiffStateMinMass;

     G4double TargetDiffStateMinMass2        = TargetDiffStateMinMass        *
                                               TargetDiffStateMinMass;
     G4double TargetNonDiffStateMinMass2     = TargetNonDiffStateMinMass     *
                                               TargetNonDiffStateMinMass;

// Transform momenta to cms and then rotate parallel to z axis;

//	   G4LorentzVector Psum;
//	   Psum=Pprojectile+Ptarget;

     G4LorentzRotation toCms(-1*Psum.boostVector());

     G4LorentzVector Ptmp=toCms*Pprojectile;
     if ( Ptmp.pz() <= 0. )
        {
	   // "String" moving backwards in  CMS, abort collision !!
           //G4cout << " abort Collision!! " << G4endl;
          return false;
	 }

     toCms.rotateZ(-1*Ptmp.phi());
     toCms.rotateY(-1*Ptmp.theta());

     G4LorentzRotation toLab(toCms.inverse());

     Pprojectile.transform(toCms);
     Ptarget.transform(toCms);

     G4double Pt2;
     G4double ProjMassT2, ProjMassT;
     G4double TargMassT2, TargMassT;
     G4double PZcms2, PZcms;
     G4double PMinusMin, PMinusMax;
//     G4double PPlusMin , PPlusMax;
     G4double TPlusMin , TPlusMax;
     G4double PMinusNew, PPlusNew, TPlusNew, TMinusNew;

//   G4double S=Psum.mag2();
     G4double SqrtS=std::sqrt(S);

     if(absPDGcode > 1000 && SqrtS < 2200*MeV)
     {return false;}                         // The model cannot work for
                                             // p+p-interactions
                                             // at Plab < 1.3 GeV/c.

     if(( absPDGcode == 211 || PDGcode ==  111) && SqrtS < 1600*MeV)
     {return false;}                         // The model cannot work for
                                             // Pi+p-interactions
                                             // at Plab < 1. GeV/c.

     if(( absPDGcode == 321 || PDGcode == -311) && SqrtS < 1600*MeV)
     {return false;}                         // The model cannot work for
                                             // K+p-interactions
                                             // at Plab < ??? GeV/c.  ???

     PZcms2=(S*S+M0projectile2*M0projectile2+M0target2*M0target2-
             2*S*M0projectile2 - 2*S*M0target2 - 2*M0projectile2*M0target2)
             /4./S;

     if(PZcms2 < 0)
     {return false;}   // It can be in an interaction with off-shell nuclear nucleon

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
/*
G4cout << "Pprojectile aft boost : " << Pprojectile <<" "<<Pprojectile.mag()<< G4endl;
G4cout << "Ptarget aft boost : " << Ptarget <<" "<<Ptarget.mag()<< G4endl;
G4cout << "cms aft boost : " << (Pprojectile+ Ptarget) << G4endl;
G4cout << " Projectile Xplus / Xminus : " <<
            Pprojectile.plus() << " / " << Pprojectile.minus() << G4endl;
G4cout << " Target Xplus / Xminus : " <<           Ptarget.plus() << " / " << Ptarget.minus() << G4endl;
G4cout<<"maxPtSquare "<<maxPtSquare<<G4endl;
*/
     G4LorentzVector Qmomentum;
     G4double Qminus, Qplus;

     G4int whilecount=0;
//  Choose a process

     if(G4UniformRand() < ProbOfDiffraction)
       {
        if(G4UniformRand() < ProbProjectileDiffraction)
        { //-------- projectile diffraction ---------------
//G4cout<<"   Projectile diffraction"<<G4endl;
//Uzhi_projectilediffraction++;
         do {
//             Generate pt
//             if (whilecount++ >= 500 && (whilecount%100)==0)
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants possibly looping"
//	   	 << ", loop count/ maxPtSquare : "
//           	 << whilecount << " / " << maxPtSquare << G4endl;
             if (whilecount > 1000 )
             {
              Qmomentum=G4LorentzVector(0.,0.,0.,0.);
              return false;    //  Ignore this interaction
             };
// --------------- Check that the interaction is possible -----------
             ProjMassT2=ProjectileDiffStateMinMass2;
             ProjMassT =ProjectileDiffStateMinMass;

             TargMassT2=M0target2;
             TargMassT =M0target;

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                     2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                    /4./S;
//G4cout<<" Pt2 Mpt Mtt Pz2 "<<Pt2<<" "<<ProjMassT<<" "<<TargMassT<<" "<<PZcms2<<G4endl;

if(PZcms2 < 0 ) 
{
/*
G4cout<<"whilecount "<<whilecount<<" "<<Pt2<<" "<<ProjMassT<<" "<<TargMassT<<" "<<PZcms2<<G4endl;
G4int Uzhi; G4cin>>Uzhi;
*/
return false;
};
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
//G4cout<<" Pt2 Mpt Mtt Pz2 "<<Pt2<<" "<<ProjMassT<<" "<<TargMassT<<" "<<PZcms2<<G4endl;

//             if(PZcms2 < 0 ) {PZcms2=0;};
             if(PZcms2 < 0 ) continue;
             PZcms =std::sqrt(PZcms2);

             PMinusMin=std::sqrt(ProjMassT2+PZcms2)-PZcms;
             PMinusMax=SqrtS-TargMassT;
//G4cout<<" SqrtS P+mim max "<<SqrtS<<" "<<PMinusMin<<" "<<PMinusMax<<G4endl;

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
//G4cout<<"   Target difraction"<<G4endl;
//Uzhi_targetdiffraction++;
         do {
//             Generate pt
//             if (whilecount++ >= 500 && (whilecount%100)==0)
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants possibly looping"
//	   	 << ", loop count/ maxPtSquare : "
//           	 << whilecount << " / " << maxPtSquare << G4endl;
             if (whilecount > 1000 )
             {
              Qmomentum=G4LorentzVector(0.,0.,0.,0.);
              return false;    //  Ignore this interaction
             };
// --------------- Check that the interaction is possible -----------
             ProjMassT2=M0projectile2;
             ProjMassT =M0projectile;

             TargMassT2=TargetDiffStateMinMass2;
             TargMassT =TargetDiffStateMinMass;

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                     2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                    /4./S;
//G4cout<<" Pt2 Mpt Mtt Pz2 "<<Pt2<<" "<<ProjMassT<<" "<<TargMassT<<" "<<PZcms2<<G4endl;

if(PZcms2 < 0 ) 
{
/*
G4cout<<"whilecount "<<whilecount<<" "<<Pt2<<" "<<ProjMassT<<" "<<TargMassT<<" "<<PZcms2<<G4endl;
G4int Uzhi; G4cin>>Uzhi;
*/
return false;
};
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
/*
if(PZcms2 < 0 ) 
{
G4cout<<"whilecount "<<whilecount<<" "<<Pt2<<" "<<ProjMassT<<" "<<TargMassT<<" "<<PZcms2<<G4endl;
G4int Uzhi; G4cin>>Uzhi;
return false;
};
*/
             if(PZcms2 < 0 ) continue;
             PZcms =std::sqrt(PZcms2);

             TPlusMin=std::sqrt(TargMassT2+PZcms2)-PZcms;
             TPlusMax=SqrtS-ProjMassT;

//G4cout<<" Tmin max "<<TPlusMin<<" "<<TPlusMax<<G4endl;

             TPlusNew=ChooseP(TPlusMin, TPlusMax);

//TPlusNew=TPlusMax;
//G4cout<<"T+new "<<TPlusNew<<G4endl;

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
//G4cout<<"   Non-difraction"<<G4endl;
         do {
//             Generate pt
//             if (whilecount++ >= 500 && (whilecount%100)==0)
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants possibly looping"
//	   	 << ", loop count/ maxPtSquare : "
//           	 << whilecount << " / " << maxPtSquare << G4endl;
             if (whilecount > 1000 )
             {
              Qmomentum=G4LorentzVector(0.,0.,0.,0.);
              return false;    //  Ignore this interaction
             };
// --------------- Check that the interaction is possible -----------
             ProjMassT2=ProjectileNonDiffStateMinMass2;
             ProjMassT =ProjectileNonDiffStateMinMass;

             TargMassT2=TargetNonDiffStateMinMass2;
             TargMassT =TargetNonDiffStateMinMass;

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                    2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                   /4./S;
//G4cout<<" Pt2 Mpt Mtt Pz2 "<<Pt2<<" "<<ProjMassT<<" "<<TargMassT<<" "<<PZcms2<<G4endl;

if(PZcms2 < 0 ) 
{
/*
G4cout<<"whilecount "<<whilecount<<" "<<Pt2<<" "<<ProjMassT<<" "<<TargMassT<<" "<<PZcms2<<G4endl;
G4int Uzhi; G4cin>>Uzhi;
*/
return false;
};
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
/*
G4cout<<"ProjectileNonDiffStateMinMass2 "<<ProjectileNonDiffStateMinMass2<<G4endl;
G4cout<<"TargetNonDiffStateMinMass2     "<<TargetNonDiffStateMinMass2<<G4endl;
G4cout<<"Mt "<<ProjMassT<<" "<<TargMassT<<" "<<Pt2<<" "<<PZcms2<<G4endl<<G4endl;
*/
//             if(PZcms2 < 0 ) {PZcms2=0;};
             if(PZcms2 < 0 ) continue;
             PZcms =std::sqrt(PZcms2);

             PMinusMin=std::sqrt(ProjMassT2+PZcms2)-PZcms;
             PMinusMax=SqrtS-TargMassT;

             PMinusNew=ChooseP(PMinusMin, PMinusMax);
// PMinusNew=1./sqrt(1./PMinusMin-G4UniformRand()*(1./PMinusMin-1./PMinusMax));

//G4cout<<"Proj "<<PMinusMin<<" "<<PMinusMax<<" "<<PMinusNew<<G4endl;

//PMinusNew=PMinusMax; //+++++++++++++++++++++++++++++++++++ Vova

             Qminus=PMinusNew-Pprojectile.minus();

             TPlusMin=std::sqrt(TargMassT2+PZcms2)-PZcms;
//             TPlusMax=SqrtS-PMinusNew;                      // Vova
             TPlusMax=SqrtS-ProjMassT;                        // Vova

             TPlusNew=ChooseP(TPlusMin, TPlusMax);

//G4cout<<"Targ "<<TPlusMin<<" "<<TPlusMax<<" "<<TPlusNew<<G4endl;
//G4cout<<PMinusNew<<" "<<TPlusNew<<G4endl;

             Qplus=-(TPlusNew-Ptarget.plus());

             Qmomentum.setPz( (Qplus-Qminus)/2 );
             Qmomentum.setE(  (Qplus+Qminus)/2 );
/*
G4cout << "Qplus / Qminus " << Qplus << " / " << Qminus<<G4endl;
G4cout << "pt2" << pt2 << G4endl;
G4cout << "Qmomentum " << Qmomentum << G4endl;
G4cout << " Masses (P/T) : " << (Pprojectile+Qmomentum).mag() <<
           " / " << (Ptarget-Qmomentum).mag() << G4endl;   // mag()
G4cout<<"Mprojectile "<<std::sqrt(M0projectile2)<<G4endl;
G4cout<<"Mtarget     "<<std::sqrt(M0target2    )<<G4endl;
G4cout<<"ProjectileDiffStateMinMass "<<std::sqrt(ProjectileDiffStateMinMass2)<<G4endl;
G4cout<<"TargetDiffStateMinMass "<<std::sqrt(TargetDiffStateMinMass2)<<G4endl;
*/
       } while (
 ((Pprojectile+Qmomentum).mag2() <  ProjectileNonDiffStateMinMass2) || //No double Diffraction
 ((Ptarget    -Qmomentum).mag2() <  TargetNonDiffStateMinMass2    ));
         }

//G4int Uzhiinp; G4cin>>Uzhiinp;    // Vova

	   Pprojectile += Qmomentum;
	   Ptarget     -= Qmomentum;
/*
G4cout << "Pprojectile with Q : " << Pprojectile << G4endl;
G4cout << "Ptarget with Q : " << Ptarget << G4endl;
G4cout << "Target	 mass  " <<  Ptarget.mag() << G4endl;
G4cout << "Projectile mass  " <<  Pprojectile.mag() << G4endl;
//
//G4cout << "Projectile back: " << toLab * Pprojectile << G4endl;
//G4cout << "Target back: " << toLab * Ptarget << G4endl;
*/
//-------------- Flip if projectale moves in backward direction ------------
//G4bool Flip=Pprojectile.pz()< 0.;


// Transform back and update SplitableHadron Participant.
	   Pprojectile.transform(toLab);
	   Ptarget.transform(toLab);

//G4cout << "Pprojectile with Q M: " << Pprojectile<<" "<<  Pprojectile.mag() << G4endl;
//G4cout << "Ptarget     with Q M: " << Ptarget    <<" "<<  Ptarget.mag()     << G4endl;
//G4cout << "Target	 mass  " <<  Ptarget.mag() << G4endl;
//G4cout << "Projectile mass  " <<  Pprojectile.mag() << G4endl;

/*
if(!Flip){
	   projectile->Set4Momentum(Pprojectile);
	   target->Set4Momentum(Ptarget);
         }
else     {
           G4ParticleDefinition * t_Definition=projectile->GetDefinition();
           projectile->SetDefinition(target->GetDefinition());
           projectile->Set4Momentum(Ptarget);
           target->SetDefinition(t_Definition);
           target->Set4Momentum(Pprojectile);
          }
*/
//
/*
if(G4UniformRand() < 1.) {
           G4ParticleDefinition * t_Definition=projectile->GetDefinition();
           projectile->SetDefinition(target->GetDefinition());
           target->SetDefinition(t_Definition);
}
*/ // For flip, for HARP

           G4double ZcoordinateOfCurrentInteraction = target->GetPosition().z();
// It is assumed that nucleon z-coordinates are ordered on increasing -----------

           G4double betta_z=projectile->Get4Momentum().pz()/projectile->Get4Momentum().e();

           G4double ZcoordinateOfPreviousCollision=projectile->GetPosition().z();
           if(projectile->GetSoftCollisionCount()==0) {
              projectile->SetTimeOfCreation(0.);
              target->SetTimeOfCreation(0.);
              ZcoordinateOfPreviousCollision=ZcoordinateOfCurrentInteraction;
           }
           
           G4ThreeVector thePosition(projectile->GetPosition().x(),
                                     projectile->GetPosition().y(),
                                     ZcoordinateOfCurrentInteraction);
           projectile->SetPosition(thePosition);

           G4double TimeOfPreviousCollision=projectile->GetTimeOfCreation();
           G4double TimeOfCurrentCollision=TimeOfPreviousCollision+ 
                    (ZcoordinateOfCurrentInteraction-ZcoordinateOfPreviousCollision)/betta_z; 

           projectile->SetTimeOfCreation(TimeOfCurrentCollision);
           target->SetTimeOfCreation(TimeOfCurrentCollision);

	   projectile->Set4Momentum(Pprojectile);
	   target->Set4Momentum(Ptarget);

           projectile->IncrementCollisionCount(1);
           target->IncrementCollisionCount(1);

//
//G4cout<<"Out of Excitation --------------------"<<G4endl;
//G4int Uzhiinp; G4cin>>Uzhiinp;    // Vova

	   return true;
}

// ---------------------------------------------------------------------
G4ExcitedString * G4DiffractiveExcitation::
           String(G4VSplitableHadron * hadron, G4bool isProjectile) const
{

//G4cout<<"G4DiffractiveExcitation::String isProj"<<isProjectile<<G4endl;

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
// Uzhi
//G4cout<<"G4ExcitedString * G4DiffractiveExcitation::String"<<G4endl;
//G4cout<<hadron->GetTimeOfCreation()<<" "<<hadron->GetPosition()/fermi<<G4endl;

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
//
/* Uzhi
	G4double ptSquared= hadron->Get4Momentum().perp2();
	G4double transverseMassSquared= hadron->Get4Momentum().plus()
				    *	hadron->Get4Momentum().minus();


	G4double maxAvailMomentumSquared=
		 sqr( std::sqrt(transverseMassSquared) - std::sqrt(ptSquared) );

        G4double widthOfPtSquare = 0.25*GeV*GeV;       // Uzhi 11.07 <Pt^2>=0.25 ??????????????????
	G4ThreeVector pt=GaussianPt(widthOfPtSquare,maxAvailMomentumSquared);

	G4LorentzVector Pstart(G4LorentzVector(pt,0.));
	G4LorentzVector Pend;
	Pend.setPx(hadron->Get4Momentum().px() - pt.x());
	Pend.setPy(hadron->Get4Momentum().py() - pt.y());

	G4double tm1=hadron->Get4Momentum().minus() +
	  ( Pend.perp2()-Pstart.perp2() ) / hadron->Get4Momentum().plus();

	G4double tm2= std::sqrt( std::max(0., sqr(tm1) -
	     4. * Pend.perp2() * hadron->Get4Momentum().minus()
	      /  hadron->Get4Momentum().plus() ));

	G4int Sign= isProjectile ? -1 : 1;

	G4double endMinus  = 0.5 * (tm1 + Sign*tm2);
	G4double startMinus= hadron->Get4Momentum().minus() - endMinus;

	G4double startPlus= Pstart.perp2() /  startMinus;
	G4double endPlus  = hadron->Get4Momentum().plus() - startPlus;

	Pstart.setPz(0.5*(startPlus - startMinus));
	Pstart.setE(0.5*(startPlus + startMinus));

	Pend.setPz(0.5*(endPlus - endMinus));
	Pend.setE(0.5*(endPlus + endMinus));
*/ // Uzhi
	start->Set4Momentum(Pstart);
	end->Set4Momentum(Pend);
/*
G4cout<<"G4DiffractiveExcitation::String hadro"<<hadron->Get4Momentum()<<" "<<hadron->Get4Momentum().mag2()<<G4endl;

G4cout<<"G4DiffractiveExcitation::String start"<<start->Get4Momentum()<<" "<<start->GetPDGcode()<<G4endl;

G4cout<<"G4DiffractiveExcitation::String end  "<<  end->Get4Momentum()<<" "<<  end->GetPDGcode()<<G4endl;
G4int Uzhi; G4cin>>Uzhi;
*/
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
G4double G4DiffractiveExcitation::ChooseP(G4double Pmin, G4double Pmax) const // Uzhi
{
// choose an x between Xmin and Xmax with P(x) ~ 1/x
//  to be improved...

	G4double range=Pmax-Pmin;                                         // Uzhi

	if ( Pmin <= 0. || range <=0. )
	{
		G4cout << " Pmin, range : " << Pmin << " , " << range << G4endl;
		throw G4HadronicException(__FILE__, __LINE__, "G4DiffractiveExcitation::ChooseP : Invalid arguments ");
	}

	G4double P;
/*                                                                          // Uzhi
	do {
	    x=Xmin + G4UniformRand() * range;
	}  while ( Xmin/x < G4UniformRand() );
*/                                                                          // Uzhi

        P=Pmin * std::pow(Pmax/Pmin,G4UniformRand());                       // Uzhi

//debug-hpw	cout << "DiffractiveX "<<x<<G4endl;
	return P;
}

// ---------------------------------------------------------------------
G4ThreeVector G4DiffractiveExcitation::GaussianPt(G4double AveragePt2, 
                                                  G4double maxPtSquare) const // Uzhi
{            //  @@ this method is used in FTFModel as well. Should go somewhere common!

	G4double Pt2;
/*                                                                          // Uzhi
	do {
	    pt2=widthSquare * std::log( G4UniformRand() );
	} while ( pt2 > maxPtSquare);
*/                                                                          // Uzhi

        Pt2 = -AveragePt2 * std::log(1. + G4UniformRand() * 
                           (std::exp(-maxPtSquare/AveragePt2)-1.));// Uzhi

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
