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
// $Id: G4DiffractiveExcitation.cc,v 1.14 2009-09-17 18:24:30 vuzhinsk Exp $
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
#include "G4FTFParameters.hh"
#include "G4ElasticHNScattering.hh"

#include "G4LorentzRotation.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh" 
#include "G4VSplitableHadron.hh"
#include "G4ExcitedString.hh"
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
                     G4FTFParameters    *theParameters,
                     G4ElasticHNScattering *theElastic) const  // Uzhi 03.09.09
{
// -------------------- Projectile parameters -----------------------
     G4LorentzVector Pprojectile=projectile->Get4Momentum();
//G4cout<<Pprojectile<<G4endl;

     if(Pprojectile.z() < 0.)
     {
       target->SetStatus(2);
       return false;
     } 

     G4double ProjectileRapidity = Pprojectile.rapidity();

     G4int    ProjectilePDGcode=projectile->GetDefinition()->GetPDGEncoding();
//     G4ParticleDefinition * ProjectileDefinition = projectile->GetDefinition();
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

//G4cout<<"PDGs "<<ProjectilePDGcode<<" "<<TargetPDGcode<<G4endl;

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

     G4double ProbOfDiffraction=ProbProjectileDiffraction +
                                ProbTargetDiffraction;

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
//G4cout<<" SqrtS "<<SqrtS<<" 2200 "<<G4endl;
     if(absProjectilePDGcode > 1000 && SqrtS < 2300*MeV)
     {target->SetStatus(2);  return false;}  // The model cannot work for
                                             // p+p-interactions
                                             // at Plab < 1.62 GeV/c.

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

//G4cout<<1./std::exp(-2.0*(ProjectileRapidity - TargetRapidity))<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;

     G4double MagQuarkExchange        =theParameters->GetMagQuarkExchange();
     G4double SlopeQuarkExchange      =theParameters->GetSlopeQuarkExchange();
     G4double DeltaProbAtQuarkExchange=theParameters->GetDeltaProbAtQuarkExchange();

     if(G4UniformRand() < MagQuarkExchange*std::exp(-SlopeQuarkExchange*
              (ProjectileRapidity - TargetRapidity)))
     {    
//G4cout<<"Exchange "<<G4endl;
G4bool StandardExcitation(false); // ================================= 
      G4int ProjQ1(0), ProjQ2(0), ProjQ3(0);

      if(absProjectilePDGcode < 1000 )
      {    // projectile is meson ----------------- 
       UnpackMeson(ProjectilePDGcode, ProjQ1, ProjQ2);  // Uzhi 7.09.09
      } else 
      {    // projectile is baryon ----------------
       UnpackBaryon(ProjectilePDGcode, ProjQ1, ProjQ2, ProjQ3);  // Uzhi 7.09.09
      } // End of the hadron's unpacking ----------

//  Target unpacking ------------------------------
      G4int TargQ1(0), TargQ2(0), TargQ3(0);
      UnpackBaryon(TargetPDGcode, TargQ1, TargQ2, TargQ3);  // Uzhi 7.09.09

// Sampling of exchanged quarks -------------------
      G4int ProjExchangeQ(0);
      G4int TargExchangeQ(0);

//G4cout<<"Targ "<<TargQ1<<" "<<TargQ2<<" "<<TargQ3<<G4endl;
      if(absProjectilePDGcode < 1000 )
      {    // projectile is meson ----------------- 
//G4cout<<"Proj Q1 Q2 "<<ProjQ1<<" "<<ProjQ2<<G4endl;

       if(ProjQ1 > 0 ) // ProjQ1 is quark
       {  
        ProjExchangeQ = ProjQ1;
        if(ProjExchangeQ != TargQ1)
        {
         TargExchangeQ = TargQ1; TargQ1=ProjExchangeQ; ProjQ1=TargExchangeQ;
        } else
        if(ProjExchangeQ != TargQ2)
        {
         TargExchangeQ = TargQ2; TargQ2=ProjExchangeQ; ProjQ1=TargExchangeQ;
        } else          
        {
         TargExchangeQ = TargQ3;  TargQ3=ProjExchangeQ; ProjQ1=TargExchangeQ;
        }
//G4cout<<" Pr Tr "<<ProjQ1<<" "<<TargQ1<<" "<<TargQ2<<" "<<TargQ3<<G4endl;
       } else          // ProjQ2 is quark
       {  
        ProjExchangeQ = ProjQ2;
        if(ProjExchangeQ != TargQ1)
        {
         TargExchangeQ = TargQ1; TargQ1=ProjExchangeQ; ProjQ2=TargExchangeQ;
        } else
        if(ProjExchangeQ != TargQ2)
        {
         TargExchangeQ = TargQ2; TargQ2=ProjExchangeQ; ProjQ2=TargExchangeQ;
        } else 
        {
         TargExchangeQ = TargQ3;  TargQ3=ProjExchangeQ; ProjQ2=TargExchangeQ;
        }

       } // End of if(ProjQ1 > 0 ) // ProjQ1 is quark
        G4int NewProjCode(0);
        G4int aProjQ1=std::abs(ProjQ1);
        G4int aProjQ2=std::abs(ProjQ2);
        if(aProjQ1 == aProjQ2)          {NewProjCode = 111;} // Pi0-meson
        else  // |ProjQ1| # |ProjQ2|
        {
         if(aProjQ1 > aProjQ2)          {NewProjCode = aProjQ1*100+aProjQ2*10+1;}
         else                           {NewProjCode = aProjQ2*100+aProjQ1*10+1;}
        }
//       } // End of if(ProjQ1 > 0 ) // ProjQ1 is quark
        projectile->SetDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle(NewProjCode)); 

        G4int NewTargCode = NewNucleonId(TargQ1, TargQ2, TargQ3);
 
        if( (TargQ1 == TargQ2) && (TargQ1 == TargQ3) )  {NewTargCode +=2;} // For Delta isobar
        else if((G4UniformRand() < DeltaProbAtQuarkExchange) &&
                (SqrtS > 1100.))                        {NewTargCode +=2;} // For Delta isobar

        target->SetDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle(NewTargCode)); 
//G4cout<<"New PDGs "<<NewProjCode<<" "<<NewTargCode<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
       //}
      } else 
      {    // projectile is baryon ----------------

       G4int NewProjCode(0);
       G4int NewTargCode(0);

G4double Same=0.; //0.3; //0.5;
G4bool ProjDeltaHasCreated(false);
G4bool TargDeltaHasCreated(false);
       G4double Ksi=G4UniformRand();
       if(G4UniformRand() < 0.5) 
       {   // Sampling exchanged quark from the projectile ---
//G4cout<<"Proj Exc"<<G4endl;
        if( Ksi < 0.333333 ) 
        {ProjExchangeQ = ProjQ1;}
        else if( (0.333333 <= Ksi) && (Ksi < 0.666667))
        {ProjExchangeQ = ProjQ2;}
        else
        {ProjExchangeQ = ProjQ3;}
//G4cout<<"ProjExchangeQ "<<ProjExchangeQ<<G4endl;

        if((ProjExchangeQ != TargQ1)||(G4UniformRand()<Same)) // Vova
        {
         TargExchangeQ = TargQ1; TargQ1=ProjExchangeQ; ProjExchangeQ=TargExchangeQ;
        } else
        if((ProjExchangeQ != TargQ2)||(G4UniformRand()<Same)) // Vova
        {
         TargExchangeQ = TargQ2; TargQ2=ProjExchangeQ; ProjExchangeQ=TargExchangeQ;
        } else 
        {
         TargExchangeQ = TargQ3;  TargQ3=ProjExchangeQ; ProjExchangeQ=TargExchangeQ;
        }
//G4cout<<"ProjExchangeQ "<<ProjExchangeQ<<G4endl;
        if( Ksi < 0.333333 ) 
        {ProjQ1=ProjExchangeQ;}
        else if( (0.333333 <= Ksi) && (Ksi < 0.666667))
        {ProjQ2=ProjExchangeQ;}
        else
        {ProjQ3=ProjExchangeQ;}

        NewProjCode = NewNucleonId(ProjQ1, ProjQ2, ProjQ3);
        if( (ProjQ1 == ProjQ2) && (ProjQ1 == ProjQ3) )  
        {NewProjCode +=2; ProjDeltaHasCreated = true;}

        NewTargCode = NewNucleonId(TargQ1, TargQ2, TargQ3);
        if( (TargQ1 == TargQ2) && (TargQ1 == TargQ3) )  
        {NewTargCode +=2; TargDeltaHasCreated = true;}

        if(!ProjDeltaHasCreated && !TargDeltaHasCreated)
        {NewProjCode +=2; ProjDeltaHasCreated = true;}

        if(!(ProjDeltaHasCreated && TargDeltaHasCreated))
        {
         if((G4UniformRand() < DeltaProbAtQuarkExchange) && (SqrtS >2600.)) 
         { 
//G4cout<<"2 Deltas "<<G4endl;
          if(!ProjDeltaHasCreated)
               {NewProjCode +=2; ProjDeltaHasCreated = true;}
          else {NewTargCode +=2; TargDeltaHasCreated = true;} // For Delta isobar
         }
        }

       } else
       {   // Sampling exchanged quark from the target -------
//G4cout<<"Targ Exc"<<G4endl;
        if( Ksi < 0.333333 ) 
        {TargExchangeQ = TargQ1;}
        else if( (0.333333 <= Ksi) && (Ksi < 0.666667)) 
        {TargExchangeQ = TargQ2;}
        else
        {TargExchangeQ = TargQ3;}
//G4cout<<"TargExchangeQ "<<TargExchangeQ<<G4endl;
        if((TargExchangeQ != ProjQ1)||(G4UniformRand()<Same)) // Vova
        {
         ProjExchangeQ = ProjQ1; ProjQ1=TargExchangeQ; TargExchangeQ=ProjExchangeQ;
        } else
        if((TargExchangeQ != ProjQ2)||(G4UniformRand()<Same)) // Vova
        {
         ProjExchangeQ = ProjQ2; ProjQ2=TargExchangeQ; TargExchangeQ=ProjExchangeQ;
        } else 
        {
         ProjExchangeQ = ProjQ3;  ProjQ3=TargExchangeQ; TargExchangeQ=ProjExchangeQ;
        }
//G4cout<<"TargExchangeQ "<<TargExchangeQ<<G4endl;
        if( Ksi < 0.333333 ) 
        {TargQ1=TargExchangeQ;}
        else if( (0.333333 <= Ksi) && (Ksi < 0.666667)) 
        {TargQ2=TargExchangeQ;}
        else
        {TargQ3=TargExchangeQ;}

        NewProjCode = NewNucleonId(ProjQ1, ProjQ2, ProjQ3);
        if( (ProjQ1 == ProjQ2) && (ProjQ1 == ProjQ3) )  
        {NewProjCode +=2; ProjDeltaHasCreated = true;}

        NewTargCode = NewNucleonId(TargQ1, TargQ2, TargQ3);
        if( (TargQ1 == TargQ2) && (TargQ1 == TargQ3) )  
        {NewTargCode +=2; TargDeltaHasCreated = true;}

        if(!ProjDeltaHasCreated && !TargDeltaHasCreated)
        {NewTargCode +=2; TargDeltaHasCreated = true;}

        if(!(ProjDeltaHasCreated && TargDeltaHasCreated))
        {
         if((G4UniformRand() < DeltaProbAtQuarkExchange) && (SqrtS >2600.)) 
         { 
//G4cout<<"2 Deltas "<<G4endl;
          if(!ProjDeltaHasCreated)
               {NewProjCode +=2; ProjDeltaHasCreated = true;}
          else {NewTargCode +=2; TargDeltaHasCreated = true;}
         }
        }
       } 
// Forming baryons --------------------------------------------------

       projectile->SetDefinition(
       G4ParticleTable::GetParticleTable()->FindParticle(NewProjCode)); 

       target->SetDefinition(
       G4ParticleTable::GetParticleTable()->FindParticle(NewTargCode)); 

       if( ProjDeltaHasCreated && TargDeltaHasCreated ) {StandardExcitation = true;}

//G4cout<<"Ne PDG "<<NewProjCode<<" "<<NewTargCode<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
      } // End of if projectile is baryon ---------------------------


// If we assume that final state hadrons after the charge exchange will be
// in the ground states, we have to put ----------------------------------
      M0projectile=projectile->GetDefinition()->GetPDGMass();
      M0projectile2 = M0projectile * M0projectile;

      M0target=target->GetDefinition()->GetPDGMass();
      M0target2 = M0target * M0target;

      PZcms2=(S*S+M0projectile2*M0projectile2+M0target2*M0target2-
             2*S*M0projectile2 - 2*S*M0target2 - 2*M0projectile2*M0target2)
             /4./S;
//G4cout<<M0projectile<<" "<<M0target<<G4endl;
//      if(PZcms2 < 0)
//      {target->SetStatus(2);  return false;}   // It can be in an interaction with 
                                              // off-shell nuclear nucleon
      PZcms = std::sqrt(PZcms2);

      Pprojectile.setPz( PZcms);
      Pprojectile.setE(std::sqrt(M0projectile2+PZcms2));

      Ptarget.setPz(    -PZcms);
      Ptarget.setE(std::sqrt(M0target2+PZcms2));

//G4cout<<"Proj "<<Pprojectile<<" "<<Pprojectile.mag()<<G4endl;
//G4cout<<"Targ "<<Ptarget<<" "<<Ptarget.mag()<<G4endl;

//      if(!StandardExcitation)
      {
       Pprojectile.transform(toLab);
       Ptarget.transform(toLab);

       projectile->SetTimeOfCreation(target->GetTimeOfCreation());
       projectile->SetPosition(target->GetPosition());

       projectile->Set4Momentum(Pprojectile);
       target->Set4Momentum(Ptarget);

       G4bool Result= theElastic->ElasticScattering (projectile,target,theParameters);
//G4cout<<"1 Delta result "<<Result<<"**********"<<G4endl;
       return Result;
      } 
/*
else
      {
       if(M0projectile > ProjectileDiffStateMinMass) 
       { ProjectileDiffStateMinMass +=200.*MeV; ProjectileNonDiffStateMinMass +=200.*MeV;}

       if(M0target > TargetDiffStateMinMass) 
       { TargetDiffStateMinMass +=200.*MeV; TargetNonDiffStateMinMass +=200.*MeV;}

       ProbOfDiffraction         = 1.;
       ProbProjectileDiffraction =0.5;
       ProbTargetDiffraction     =0.5;
G4cout<<"Exc DD "<<M0projectile<<" "<<M0target<<" --------------------"<<G4endl;
//       After that standard FTF excitation
      }
*/
     }  // End of charge exchange part ------------------------------

// ------------------------------------------------------------------
//     G4double ProbOfDiffraction=ProbProjectileDiffraction +
//                                ProbTargetDiffraction;
//G4cout<<ProbOfDiffraction<<" "<<ProbProjectileDiffraction<<" "<<ProbTargetDiffraction<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
     if(ProbOfDiffraction!=0.)
     {
      ProbProjectileDiffraction/=ProbOfDiffraction;
     }
     else
     {
      ProbProjectileDiffraction=0.;
     }

//ProbOfDiffraction=1.; //0.5; // Vova 

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
//G4cout<<"Proj difr"<<G4endl;
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
//G4cout<<"PZcms2 < 0 false "<<PZcms2<<G4endl;
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
//G4cout<<"Targ difr"<<G4endl;
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
//G4cout<<"Non difr"<<G4endl;
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

//G4cout<<"Masses "<<Pprojectile.mag()<<" "<<Ptarget.mag()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;

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
//G4ExcitedString * G4DiffractiveExcitation::
//           String(G4VSplitableHadron * hadron, G4bool isProjectile) const
void G4DiffractiveExcitation::CreateStrings(G4VSplitableHadron * hadron, 
                                            G4bool isProjectile,
                                            G4ExcitedString * &FirstString, 
                                            G4ExcitedString * &SecondString,
                                            G4FTFParameters *theParameters) const
{
	hadron->SplitUp();
	G4Parton *start= hadron->GetNextParton();
	if ( start==NULL)
	{ G4cout << " G4FTFModel::String() Error:No start parton found"<< G4endl;
          FirstString=0; SecondString=0;
	  return;
	}
	G4Parton *end  = hadron->GetNextParton();
	if ( end==NULL)
	{ G4cout << " G4FTFModel::String() Error:No end parton found"<< G4endl;
          FirstString=0; SecondString=0;
	  return;
	}
        G4LorentzVector Phadron=hadron->Get4Momentum();
/*
G4cout<<"Create strings had "<<hadron->GetDefinition()->GetParticleName()<<" "<<Phadron<<" "<<Phadron.mag()<<G4endl;
G4cout<<"isProjectile "<<isProjectile<<G4endl;
G4cout<<"start Q "<<start->GetDefinition()->GetPDGEncoding()<<G4endl;
G4cout<<"end   Q "<<end->GetDefinition()->GetPDGEncoding()<<G4endl;
*/
        G4LorentzVector Pstart(0.,0.,0.,0.);
        G4LorentzVector Pend(0.,0.,0.,0.);
        G4LorentzVector Pkink(0.,0.,0.,0.);
        G4LorentzVector PkinkQ1(0.,0.,0.,0.);
        G4LorentzVector PkinkQ2(0.,0.,0.,0.);

        G4int PDGcode_startQ = std::abs(start->GetDefinition()->GetPDGEncoding());
        G4int PDGcode_endQ   = std::abs(  end->GetDefinition()->GetPDGEncoding());

//--------------------------------------------------------------------------------        
        G4double Wmin(0.);
        if(isProjectile)
        {
          Wmin=theParameters->GetProjMinDiffMass();
        } else
        {
          Wmin=theParameters->GetTarMinDiffMass();
        } // end of if(isProjectile)

        G4double W = hadron->Get4Momentum().mag();
        G4double W2=W*W;
//G4cout<<"W Wmin "<<W<<" "<<Wmin<<G4endl;
        G4double Pt(0.), x1(0.), x2(0.), x3(0.);
        G4bool Kink=false;

        if(W > Wmin)
        {                                        // Kink is possible
          G4double Pt2kink=theParameters->GetPt2Kink();
          Pt = std::sqrt(Pt2kink*(std::pow(W2/16./Pt2kink+1.,G4UniformRand()) - 1.));
// Pt=0.;
//G4cout<<"Pt2kink Pt Pt2 "<<Pt2kink<<" "<<Pt<<" "<<Pt*Pt<<G4endl;
          if(Pt > 500.*MeV)
          {
             G4double Ymax = std::log(W/2./Pt + std::sqrt(W2/4./Pt/Pt - 1.));
             G4double Y=Ymax*(1.- 2.*G4UniformRand());
//G4cout<<"Ymax Y "<<Ymax<<" "<<Y<<G4endl;
             x1=1.-Pt/W*std::exp( Y);
             x3=1.-Pt/W*std::exp(-Y);
             x2=2.-x1-x3;
//G4cout<<"X1 X2 X3 "<<x1<<" "<<x2<<" "<<x3<<G4endl;
             G4double Mass_startQ = 650.*MeV;
             if(PDGcode_startQ <  3) Mass_startQ =  325.*MeV;
             if(PDGcode_startQ == 3) Mass_startQ =  500.*MeV;
             if(PDGcode_startQ == 4) Mass_startQ = 1600.*MeV;

             G4double Mass_endQ = 650.*MeV;
             if(PDGcode_endQ <  3) Mass_endQ =  325.*MeV;
             if(PDGcode_endQ == 3) Mass_endQ =  500.*MeV;
             if(PDGcode_endQ == 4) Mass_endQ = 1600.*MeV;

             G4double P2_1=W2*x1*x1/4.-Mass_endQ  *Mass_endQ;
             G4double P2_3=W2*x3*x3/4.-Mass_startQ*Mass_startQ;
     
             G4double P2_2=sqr((2.-x1-x3)*W/2.);
//G4cout<<"P2_1 P2_2 P2_3 "<<P2_1<<" "<<P2_2<<" "<<P2_3<<G4endl;
             if((P2_1 < 0.) || (P2_3 <0.))
             { Kink=false;}
             else
             {
               G4double P_1=std::sqrt(P2_1);
               G4double P_2=std::sqrt(P2_2);
               G4double P_3=std::sqrt(P2_3);
//G4cout<<"P_1 P_2 P_3 "<<P_1<<" "<<P_2<<" "<<P_3<<G4endl;
               G4double CosT12=(P2_3-P2_1-P2_2)/(2.*P_1*P_2);
               G4double CosT13=(P2_2-P2_1-P2_3)/(2.*P_1*P_3);
//G4cout<<"CosT12 CosT13 "<<CosT12<<" "<<CosT13<<" "<<std::acos(CosT12)*180./3.14159<<" "<<std::acos(CosT13)*180./3.14159<<G4endl;
               if((std::abs(CosT12) >1.) || (std::abs(CosT13) > 1.)) 
               { Kink=false;}
               else
               { 
                 Kink=true; 
                 Pstart.setPx(-Pt); Pstart.setPy(0.); Pstart.setPz(P_3*CosT13); 
                 Pend.setPx(   0.); Pend.setPy(  0.); Pend.setPz(          P_1); 
                 Pkink.setPx(  Pt); Pkink.setPy( 0.); Pkink.setPz(  P_2*CosT12);
                 Pstart.setE(x3*W/2.);                
                 Pkink.setE(Pkink.vect().mag());
                 Pend.setE(x1*W/2.);

                 G4double XkQ=GetQuarkFractionOfKink(0.,1.);
                 if(Pkink.getZ() > 0.) 
                 {
                   if(XkQ > 0.5) {PkinkQ1=XkQ*Pkink;} else {PkinkQ1=(1.-XkQ)*Pkink;}
                 } else {
                   if(XkQ > 0.5) {PkinkQ1=(1.-XkQ)*Pkink;} else {PkinkQ1=XkQ*Pkink;}
                 }

                 PkinkQ2=Pkink - PkinkQ1;
//------------------------- Minimizing Pt1^2+Pt3^2 ------------------------------
/*
G4cout<<"Pstart "<<Pstart<<G4endl;
G4cout<<"Pkink  "<<Pkink <<G4endl;
G4cout<<"Pkink1 "<<PkinkQ1<<G4endl;
G4cout<<"Pkink2 "<<PkinkQ2<<G4endl;
G4cout<<"Pend   "<<Pend  <<G4endl;
*/
                 G4double Cos2Psi=(sqr(x1) -sqr(x3)+2.*sqr(x3*CosT13))/
                          std::sqrt(sqr(sqr(x1)-sqr(x3)) + sqr(2.*x1*x3*CosT13));
                 G4double Psi=std::acos(Cos2Psi);
//G4cout<<"Psi "<<Psi<<" "<<Psi*180./twopi<<G4endl;

G4LorentzRotation Rotate;
if(isProjectile) {Rotate.rotateY(Psi);}
else             {Rotate.rotateY(pi-Psi);}                   
Rotate.rotateZ(twopi*G4UniformRand());

//G4cout<<"Rotate"<<G4endl;

Pstart*=Rotate;
Pkink*=Rotate;
PkinkQ1*=Rotate;
PkinkQ2*=Rotate;
Pend*=Rotate;
/*
G4cout<<"Pstart "<<Pstart<<G4endl;
G4cout<<"Pkink1 "<<PkinkQ1 <<G4endl;
G4cout<<"Pkink2 "<<PkinkQ2 <<G4endl;
G4cout<<"Pend   "<<Pend  <<G4endl;
*/
               }
             }      // end of if((P2_1 < 0.) || (P2_3 <0.))
          }         // end of if(Pt > 500.*MeV)
        }           // end of if(W > Wmin) Check for a kink

//--------------------------------------------------------------------------------

//G4cout<<"Kink "<<Kink<<G4endl;

        if(Kink)
        {                                        // Kink is possible
          std::vector<G4double> QuarkProbabilitiesAtGluonSplitUp =
              theParameters->GetQuarkProbabilitiesAtGluonSplitUp();

          G4int QuarkInGluon(1); G4double Ksi=G4UniformRand();
          for(unsigned int Iq=0; Iq <3; Iq++)
          {
//G4cout<<Iq<<" "<<QuarkProbabilitiesAtGluonSplitUp[Iq]<<G4endl;
if(Ksi > QuarkProbabilitiesAtGluonSplitUp[Iq]) QuarkInGluon++;}

//G4cout<<"Gquark "<<QuarkInGluon<<G4endl;
          G4Parton * Gquark = new G4Parton(QuarkInGluon);
          G4Parton * Ganti_quark = new G4Parton(-QuarkInGluon);
/*
G4cout<<Gquark->GetDefinition()->GetParticleName()<<" "<<Gquark->GetDefinition()->GetPDGEncoding()<<" "<<Gquark->GetDefinition()->GetPDGMass()<<G4endl;
G4cout<<Ganti_quark->GetDefinition()->GetParticleName()<<" "<<Ganti_quark->GetDefinition()->GetPDGEncoding()<<" "<<Ganti_quark->GetDefinition()->GetPDGMass()<<G4endl;
*/

//G4int Uzhi; G4cin>>Uzhi;

//-------------------------------------------------------------------------------
/*
DefineMomentumInZ(G4double aLightConeMomentum, G4bool aDirection); 
      void Set4Momentum(const G4LorentzVector & aMomentum);
      G4int PDGencoding;
      G4ParticleDefinition * theDefinition;
      G4LorentzVector theMomentum;
      G4ThreeVector   thePosition;
      
      G4int theColour;
      G4double theIsoSpinZ;
      G4double theSpinZ;
      
      G4double theX;
*/
//-------------------------------------------------------------------------------

//G4cout<<"Phadron "<<Phadron<<" mass "<<Phadron.mag()<<G4endl;
          G4LorentzRotation toCMS(-1*Phadron.boostVector());
//G4cout<<"To lab"<<G4endl;
          G4LorentzRotation toLab(toCMS.inverse());

          Pstart.transform(toLab);  start->Set4Momentum(Pstart);
          PkinkQ1.transform(toLab);
          PkinkQ2.transform(toLab);
          Pend.transform(toLab);    end->Set4Momentum(Pend);
/*
G4cout<<"Pstart "<<start->GetDefinition()->GetPDGEncoding()<<Pstart<<G4endl;
G4cout<<"Pkink1 "<<PkinkQ1 <<G4endl;
G4cout<<"Pkink2 "<<PkinkQ2 <<G4endl;
G4cout<<"Pend   "<<end->GetDefinition()->GetPDGEncoding()<<Pend  <<G4endl;
*/
// !!!	  G4ExcitedString * FirstString(0); G4ExcitedString * SecondString(0);
          G4int absPDGcode=std::abs(hadron->GetDefinition()->GetPDGEncoding());

//G4cout<<"absPDGcode "<<absPDGcode<<G4endl;

          if(absPDGcode < 1000)
          {                                // meson
	    if ( isProjectile )
	    {                              // Projectile
              if(end->GetDefinition()->GetPDGEncoding() > 0 )  // A quark on the end
              {                            // Quark on the end
                FirstString = new G4ExcitedString(end   ,Ganti_quark, +1);
                SecondString= new G4ExcitedString(Gquark,start      ,+1);
                Ganti_quark->Set4Momentum(PkinkQ1);
                Gquark->Set4Momentum(PkinkQ2);
/*
G4cout<<" Proj Meson end Q"<<G4endl;
G4cout<<"First string  ============ "<<G4endl;
G4cout<<"end  Q "<<end->GetDefinition()->GetPDGEncoding()<<" "<<end->Get4Momentum()<<G4endl;
G4cout<<"G antiQ"<<Ganti_quark->GetDefinition()->GetPDGEncoding()<<" "<<Ganti_quark->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Ganti_quark->Get4Momentum()+end->Get4Momentum())<<G4endl;
G4cout<<"Secondstring   ============ "<<G4endl;
G4cout<<"G    Q "<<Gquark->GetDefinition()->GetPDGEncoding()<<" "<<Gquark->Get4Momentum()<<G4endl;
G4cout<<"startQ "<<start->GetDefinition()->GetPDGEncoding()<<" "<<start->Get4Momentum()<<G4endl;

G4cout<<"Sum P  "<<(Gquark->Get4Momentum()+start->Get4Momentum())<<G4endl;
*/
              } else
              {                            // Anti_Quark on the end
                FirstString = new G4ExcitedString(end        ,Gquark, +1);
                SecondString= new G4ExcitedString(Ganti_quark,start ,+1);
                Gquark->Set4Momentum(PkinkQ1);
                Ganti_quark->Set4Momentum(PkinkQ2);
/*
G4cout<<" Proj Meson end Qbar"<<G4endl;
G4cout<<"First string  ============ "<<G4endl;
G4cout<<"end  Q "<<end->GetDefinition()->GetPDGEncoding()<<" "<<end->Get4Momentum()<<G4endl;
G4cout<<"G     Q"<<Gquark->GetDefinition()->GetPDGEncoding()<<" "<<Gquark->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Gquark->Get4Momentum()+end->Get4Momentum())<<G4endl;
G4cout<<"Secondstring   ============ "<<G4endl;
G4cout<<"G antQ "<<Ganti_quark->GetDefinition()->GetPDGEncoding()<<" "<<Ganti_quark->Get4Momentum()<<G4endl;
G4cout<<"startQ "<<start->GetDefinition()->GetPDGEncoding()<<" "<<start->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Ganti_quark->Get4Momentum()+start->Get4Momentum())<<G4endl;
*/
              }   // end of if(end->GetPDGcode() > 0)
            } else {                      // Target
              if(end->GetDefinition()->GetPDGEncoding() > 0 )  // A quark on the end
              {                           // Quark on the end
                FirstString = new G4ExcitedString(Ganti_quark,end   ,-1);
                SecondString= new G4ExcitedString(start      ,Gquark,-1);
                Ganti_quark->Set4Momentum(PkinkQ2);
                Gquark->Set4Momentum(PkinkQ1);
/*
G4cout<<" Targ Meson end Q"<<G4endl;
G4cout<<"First string   ============ "<<G4endl;
G4cout<<"G antiQ"<<Ganti_quark->GetDefinition()->GetPDGEncoding()<<" "<<Ganti_quark->Get4Momentum()<<G4endl;
G4cout<<"end  Q "<<end->GetDefinition()->GetPDGEncoding()<<" "<<end->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Ganti_quark->Get4Momentum()+end->Get4Momentum())<<G4endl;
G4cout<<"Secondstring   ============ "<<G4endl;
G4cout<<"startQ "<<start->GetDefinition()->GetPDGEncoding()<<" "<<start->Get4Momentum()<<G4endl;
G4cout<<"G    Q "<<Gquark->GetDefinition()->GetPDGEncoding()<<" "<<Gquark->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Gquark->Get4Momentum()+start->Get4Momentum())<<G4endl;
*/
              } else
              {                            // Anti_Quark on the end
                FirstString = new G4ExcitedString(Gquark,end        ,-1);
                SecondString= new G4ExcitedString(start ,Ganti_quark,-1);
                Gquark->Set4Momentum(PkinkQ2);
                Ganti_quark->Set4Momentum(PkinkQ1);
/*
G4cout<<" Targ Meson end Qbar"<<G4endl;
G4cout<<"First string   ============ "<<G4endl;
G4cout<<"G     Q"<<Gquark->GetDefinition()->GetPDGEncoding()<<" "<<Gquark->Get4Momentum()<<G4endl;
G4cout<<"end  Q "<<end->GetDefinition()->GetPDGEncoding()<<" "<<end->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Gquark->Get4Momentum()+end->Get4Momentum())<<G4endl;
G4cout<<"Secondstring   ============ "<<G4endl;
G4cout<<"startQ "<<start->GetDefinition()->GetPDGEncoding()<<" "<<start->Get4Momentum()<<G4endl;
G4cout<<"G antQ "<<Ganti_quark->GetDefinition()->GetPDGEncoding()<<" "<<Ganti_quark->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Ganti_quark->Get4Momentum()+start->Get4Momentum())<<G4endl;
*/
              }   // end of if(end->GetPDGcode() > 0)
	    }     // end of if ( isProjectile )
          } else  // if(absPDGCode < 1000)
          {                             // Baryon/AntiBaryon
	    if ( isProjectile )
	    {                              // Projectile
              if((end->GetDefinition()->GetParticleType() == "diquarks") &&
                 (end->GetDefinition()->GetPDGEncoding() > 0           )   ) 
              {                            // DiQuark on the end
                FirstString = new G4ExcitedString(end        ,Gquark, +1);
                SecondString= new G4ExcitedString(Ganti_quark,start ,+1);
                Gquark->Set4Momentum(PkinkQ1);
                Ganti_quark->Set4Momentum(PkinkQ2);
/*
G4cout<<" Proj baryon end QQ"<<G4endl;
G4cout<<"First string   ============ "<<G4endl;
G4cout<<"end OQ "<<end->GetDefinition()->GetPDGEncoding()<<" "<<end->Get4Momentum()<<G4endl;
G4cout<<"G    Q"<<Gquark->GetDefinition()->GetPDGEncoding()<<" "<<Gquark->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Gquark->Get4Momentum()+end->Get4Momentum())<<G4endl;
G4cout<<"Secondstring   ============ "<<G4endl;
G4cout<<"G  Qbar"<<Ganti_quark->GetDefinition()->GetPDGEncoding()<<" "<<Ganti_quark->Get4Momentum()<<G4endl;
G4cout<<"startQ "<<start->GetDefinition()->GetPDGEncoding()<<" "<<start->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Ganti_quark->Get4Momentum()+start->Get4Momentum())<<G4endl;
*/
              } else
              {                            // Anti_DiQuark on the end or quark
                FirstString = new G4ExcitedString(end   ,Ganti_quark, +1);
                SecondString= new G4ExcitedString(Gquark,start      ,+1);
                Ganti_quark->Set4Momentum(PkinkQ1);
                Gquark->Set4Momentum(PkinkQ2);
/*
G4cout<<" Proj baryon end Q"<<G4endl;
G4cout<<"First string   ============ "<<G4endl;
G4cout<<"end OQ "<<end->GetDefinition()->GetPDGEncoding()<<" "<<end->Get4Momentum()<<G4endl;
G4cout<<"G antQ"<<Ganti_quark->GetDefinition()->GetPDGEncoding()<<" "<<Ganti_quark->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Ganti_quark->Get4Momentum()+end->Get4Momentum())<<G4endl;
G4cout<<"Secondstring   ============ "<<G4endl;
G4cout<<"G  Q   "<<Gquark->GetDefinition()->GetPDGEncoding()<<" "<<Gquark->Get4Momentum()<<G4endl;
G4cout<<"startQ "<<start->GetDefinition()->GetPDGEncoding()<<" "<<start->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Gquark->Get4Momentum()+start->Get4Momentum())<<G4endl;
*/
              }   // end of if(end->GetPDGcode() > 0)
            } else {                      // Target

              if((end->GetDefinition()->GetParticleType() == "diquarks") &&
                 (end->GetDefinition()->GetPDGEncoding() > 0           )   ) 
              {                            // DiQuark on the end
                FirstString = new G4ExcitedString(Gquark,end        ,-1);

                SecondString= new G4ExcitedString(start ,Ganti_quark,-1);
                Gquark->Set4Momentum(PkinkQ1);
                Ganti_quark->Set4Momentum(PkinkQ2);
/*
G4cout<<" Targ baryon end QQ"<<G4endl;
G4cout<<"First string   ============ "<<G4endl;
G4cout<<"G  Q   "<<Gquark->GetDefinition()->GetPDGEncoding()<<" "<<Gquark->Get4Momentum()<<G4endl;
G4cout<<"end OQ "<<end->GetDefinition()->GetPDGEncoding()<<" "<<end->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Gquark->Get4Momentum()+end->Get4Momentum())<<G4endl;
G4cout<<"Secondstring   ============ "<<G4endl;
G4cout<<"startQ "<<start->GetDefinition()->GetPDGEncoding()<<" "<<start->Get4Momentum()<<G4endl;
G4cout<<"G  Qbar"<<Ganti_quark->GetDefinition()->GetPDGEncoding()<<" "<<Ganti_quark->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Ganti_quark->Get4Momentum()+start->Get4Momentum())<<G4endl;
*/
              } else
              {                            // Anti_DiQuark on the end or Q
                FirstString = new G4ExcitedString(Ganti_quark,end   ,-1);
                SecondString= new G4ExcitedString(start      ,Gquark,-1);
                Gquark->Set4Momentum(PkinkQ2);
                Ganti_quark->Set4Momentum(PkinkQ1);
/*
G4cout<<" Targ baryon end Q"<<G4endl;
G4cout<<"First string   ============ "<<G4endl;
G4cout<<"G  Qbar"<<Ganti_quark->GetDefinition()->GetPDGEncoding()<<" "<<Ganti_quark->Get4Momentum()<<G4endl;
G4cout<<"end O "<<end->GetDefinition()->GetPDGEncoding()<<" "<<end->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Ganti_quark->Get4Momentum()+end->Get4Momentum())<<G4endl;
G4cout<<"Secondstring   ============ "<<G4endl;
G4cout<<"startQ "<<start->GetDefinition()->GetPDGEncoding()<<" "<<start->Get4Momentum()<<G4endl;
G4cout<<"G  Q   "<<Gquark->GetDefinition()->GetPDGEncoding()<<" "<<Gquark->Get4Momentum()<<G4endl;
G4cout<<"Sum P  "<<(Gquark->Get4Momentum()+start->Get4Momentum())<<G4endl;
*/
              }   // end of if(end->GetPDGcode() > 0)
	    }     // end of if ( isProjectile )
          }  // end of if(absPDGcode < 1000)

	  FirstString->SetTimeOfCreation(hadron->GetTimeOfCreation());
	  FirstString->SetPosition(hadron->GetPosition());

	  SecondString->SetTimeOfCreation(hadron->GetTimeOfCreation());
	  SecondString->SetPosition(hadron->GetPosition());

// -------------------------------------------------------------------------  
        } else                                   // End of kink is possible
        {                                        // Kink is impossible
	  if ( isProjectile )
	  {
		FirstString= new G4ExcitedString(end,start, +1);
	  } else {
		FirstString= new G4ExcitedString(start,end, -1);
	  }

          SecondString=0;

	  FirstString->SetTimeOfCreation(hadron->GetTimeOfCreation());
	  FirstString->SetPosition(hadron->GetPosition());

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
          }

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
          }

          Momentum=-Pstart.mag();
          Pstart.setT(Momentum);  // It is assumed that quark has m=0.

          Momentum=-Pend.mag();
          Pend.setT(Momentum);    // It is assumed that di-quark has m=0.

	  start->Set4Momentum(Pstart);
	  end->Set4Momentum(Pend);
          SecondString=0;
        }            // End of kink is impossible 

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

        return;
   
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
G4double G4DiffractiveExcitation::GetQuarkFractionOfKink(G4double zmin, G4double zmax) const
    {
       G4double z, yf;
       do {
           z = zmin + G4UniformRand()*(zmax-zmin);
           yf = z*z +sqr(1 - z);	
           } 
       while (G4UniformRand() > yf); 
       return z;
    }
// ---------------------------------------------------------------------
void G4DiffractiveExcitation::UnpackMeson(const G4int IdPDG, G4int &Q1, G4int &Q2) const // Uzhi 7.09.09
    {
       G4int absIdPDG = std::abs(IdPDG);
       Q1=  absIdPDG/ 100;
       Q2= (absIdPDG %100)/10;
	   
       G4int anti= 1 -2 * ( std::max( Q1, Q2 ) % 2 );

       if (IdPDG < 0 ) anti *=-1;
       Q1 *= anti;
       Q2 *= -1 * anti;
       return;
    }
// ---------------------------------------------------------------------
void G4DiffractiveExcitation::UnpackBaryon(G4int IdPDG, 
                              G4int &Q1, G4int &Q2, G4int &Q3) const // Uzhi 7.09.09
    {
       Q1 =  IdPDG           / 1000;
       Q2 = (IdPDG % 1000)  / 100;
       Q3 = (IdPDG % 100)   / 10;
       return;
    }
// ---------------------------------------------------------------------
G4int G4DiffractiveExcitation::NewNucleonId(G4int Q1, G4int Q2, G4int Q3) const // Uzhi 7.09.09
    {
       G4int TmpQ(0);
       if( Q3 > Q2 ) 
       {
        TmpQ = Q2;
        Q2 = Q3;
        Q3 = TmpQ;
       } else if( Q3 > Q1 )
       {
        TmpQ = Q1;
        Q1 = Q3;
        Q3 = TmpQ;
       }

       if( Q2 > Q1 ) 
       {
        TmpQ = Q1;
        Q1 = Q2;
        Q2 = TmpQ;
       }

       G4int NewCode = Q1*1000 + Q2* 100 + Q3*  10 + 2; 
       return NewCode;
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
