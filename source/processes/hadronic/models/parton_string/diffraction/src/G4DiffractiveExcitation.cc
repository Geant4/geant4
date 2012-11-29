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
// $Id$
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

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
#include "G4Neutron.hh"
#include "G4ParticleDefinition.hh"

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
                     G4ElasticHNScattering *theElastic) const  
{
//G4cout<<G4endl<<"ExciteParticipants --------------"<<G4endl;
// -------------------- Projectile parameters -----------------------
     G4LorentzVector Pprojectile=projectile->Get4Momentum();

     if(Pprojectile.z() < 0.)
     {
       target->SetStatus(2);
       return false;
     } 

     G4double ProjectileRapidity = Pprojectile.rapidity();

     G4int    ProjectilePDGcode=projectile->GetDefinition()->GetPDGEncoding();
     G4int    absProjectilePDGcode=std::abs(ProjectilePDGcode);

     G4bool PutOnMassShell(false);
//   G4double M0projectile=projectile->GetDefinition()->GetPDGMass(); // With de-excitation
     G4double M0projectile = Pprojectile.mag();                       // Without de-excitation
//G4cout<<"M0projectile "<<M0projectile<<" "<<ProjectileRapidity<<G4endl;

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
     G4int    absTargetPDGcode=std::abs(TargetPDGcode);
//G4cout<<"Entry to QE or Excit "<<ProjectilePDGcode<<" "<<TargetPDGcode<<G4endl;

     G4LorentzVector Ptarget=target->Get4Momentum();
//G4cout<<"Pproj "<<Pprojectile<<G4endl;
//G4cout<<"Ptarget "<<Ptarget<<G4endl;
     G4double M0target = Ptarget.mag();

//   G4double TargetRapidity = Ptarget.rapidity();
//G4cout<<"M0target "<<M0target<<" "<<TargetRapidity<<G4endl;
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
     G4double ProbLogDistr=theParameters->GetProbLogDistr();         // Uzhi 21.05.2012

//     G4double ProbOfDiffraction=ProbProjectileDiffraction +
//                                ProbTargetDiffraction;

//     G4double SumMasses=M0projectile+M0target+220.*MeV; // 200->220 7 June 2011
     G4double SumMasses=M0projectile+M0target+220.*MeV; // 200->220 7 June 2011

// Kinematical properties of the interactions --------------
     G4LorentzVector Psum;      // 4-momentum in CMS
     Psum=Pprojectile+Ptarget;
     G4double S=Psum.mag2(); 

//Uzhi_SqrtS=std::sqrt(S);

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

/*
G4cout<<"Proj "<<Pprojectile<<G4endl;
G4cout<<"Targ "<<Ptarget<<G4endl;
G4cout<<"SqrtS     "<<SqrtS<<G4endl;
G4cout<<"M0pr M0tr "<<M0projectile<<" "<<M0target<<" "<<SumMasses<<G4endl;
*/
//G4cout<<"SqrtS < 2300*MeV Bary "<<SqrtS<<G4endl;
//   if(absProjectilePDGcode > 1000 && (SqrtS < 2300*MeV || SqrtS < SumMasses)) // 7.06.11
     if(absProjectilePDGcode > 1000 && (SqrtS < SumMasses))
     {target->SetStatus(2);  return false;}  // The model cannot work for
                                             // p+p-interactions
                                             // at Plab < 1.62 GeV/c.

//G4cout<<"SqrtS < 1600*MeV Pion "<<SqrtS<<G4endl;

//   if(( absProjectilePDGcode == 211 || ProjectilePDGcode ==  111) &&        // 7.06.11
//     ((SqrtS < 1600*MeV) || (SqrtS < SumMasses)))
     if(( absProjectilePDGcode == 211 || ProjectilePDGcode ==  111) && 
       (SqrtS < SumMasses))
     {target->SetStatus(2);  return false;}  // The model cannot work for
                                             // Pi+p-interactions
                                             // at Plab < 1. GeV/c.
//G4cout<<"SqrtS < 1600*MeV "<<SqrtS<<G4endl;
//SumMasses=M0projectile+M0target+20.*MeV;
     if(( (absProjectilePDGcode == 321) || (ProjectilePDGcode == -311)   ||
          (absProjectilePDGcode == 311) || (absProjectilePDGcode == 130) ||
          (absProjectilePDGcode == 310)) && (SqrtS < SumMasses)) 
//        (absProjectilePDGcode == 310)) && ((SqrtS < 1600*MeV) || (SqrtS < SumMasses))) 
     {target->SetStatus(2);  return false;}  // The model cannot work for
                                             // K+p-interactions
                                             // at Plab < ??? GeV/c.  ???

     PZcms2=(S*S+M0projectile2*M0projectile2+M0target2*M0target2-
             2*S*M0projectile2 - 2*S*M0target2 - 2*M0projectile2*M0target2)
             /4./S;
//G4cout<<"PZcms2 "<<PZcms2<<G4endl;
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


//G4cout<<"Proj "<<Pprojectile<<G4endl;
//G4cout<<"Targ "<<Ptarget<<G4endl;
     G4double maxPtSquare; // = PZcms2;
/*
Uzhi_targetdiffraction=0;
Uzhi_projectilediffraction=0;
Uzhi_Mx2=1.0;
Uzhi_modT=0.;
G4int Uzhi_QE=0;
*/
/*
G4cout<<"Start --------------------"<<G4endl;
G4cout<<"Proj "<<M0projectile<<" "<<ProjectileDiffStateMinMass<<"  "<<ProjectileNonDiffStateMinMass<<G4endl;
G4cout<<"Targ "<<M0target    <<" "<<TargetDiffStateMinMass    <<" "<<TargetNonDiffStateMinMass<<G4endl;
G4cout<<"SqrtS "<<SqrtS<<G4endl;
G4cout<<"Rapid "<<ProjectileRapidity<<G4endl; //" "<<TargetRapidity<<G4endl;
*/
// Charge exchange can be possible for baryons -----------------

// Getting the values needed for exchange ----------------------
     G4double MagQuarkExchange        =theParameters->GetMagQuarkExchange();//0.045; //
     G4double SlopeQuarkExchange      =theParameters->GetSlopeQuarkExchange();//0.; //
// 777777777777777777777777777777777777777777777777777777777777777777777777777777
     G4double DeltaProbAtQuarkExchange=theParameters->GetDeltaProbAtQuarkExchange();

//     G4double NucleonMass=
//              (G4ParticleTable::GetParticleTable()->FindParticle(2112))->GetPDGMass();     
     G4double DeltaMass=
              (G4ParticleTable::GetParticleTable()->FindParticle(2224))->GetPDGMass();

//G4double TargetRapidity(0.);
//G4cout<<"Prob Q Exch "<<MagQuarkExchange*std::exp(-SlopeQuarkExchange*ProjectileRapidity)<<G4endl;

//G4cout<<"Q exc Mag Slop Wdelta"<<MagQuarkExchange<<" "<<SlopeQuarkExchange<<" "<<DeltaProbAtQuarkExchange<<G4endl;
//G4cout<<"ProjectileRapidity "<<ProjectileRapidity<<G4endl;
//G4cout<<"Prob Exc "<<MagQuarkExchange*std::exp(-SlopeQuarkExchange*(ProjectileRapidity))<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
// Check for possible quark exchange -----------------------------------

     if(G4UniformRand() < MagQuarkExchange*
        std::exp(-SlopeQuarkExchange*ProjectileRapidity))  //TargetRapidity))) 1.45
     {    
//        std::exp(-SlopeQuarkExchange*(ProjectileRapidity - 1.36)))  //TargetRapidity))) 1.45
//G4cout<<"Q exchange"<<G4endl;
//Uzhi_QE=1;
      G4int NewProjCode(0), NewTargCode(0);

      G4int ProjQ1(0), ProjQ2(0), ProjQ3(0);

//  Projectile unpacking --------------------------
      if(absProjectilePDGcode < 1000 )
      {    // projectile is meson ----------------- 
       UnpackMeson(ProjectilePDGcode, ProjQ1, ProjQ2);  
      } else 
      {    // projectile is baryon ----------------
       UnpackBaryon(ProjectilePDGcode, ProjQ1, ProjQ2, ProjQ3);
      } // End of the hadron's unpacking ----------

//  Target unpacking ------------------------------
      G4int TargQ1(0), TargQ2(0), TargQ3(0);
      UnpackBaryon(TargetPDGcode, TargQ1, TargQ2, TargQ3); 

//G4cout<<ProjQ1<<" "<<ProjQ2<<" "<<ProjQ3<<G4endl;
//G4cout<<TargQ1<<" "<<TargQ2<<" "<<TargQ3<<G4endl;
// Sampling of exchanged quarks -------------------
      G4int ProjExchangeQ(0);
      G4int TargExchangeQ(0);

      if(absProjectilePDGcode < 1000 )
      {    // projectile is meson ----------------- 

       if(ProjQ1 > 0 ) // ProjQ1 is quark
       {  
        G4int Navailable=0;
        ProjExchangeQ = ProjQ1;
        if(ProjExchangeQ != TargQ1) Navailable++;
        if(ProjExchangeQ != TargQ2) Navailable++;
        if(ProjExchangeQ != TargQ3) Navailable++;

        G4int Nsampled=CLHEP::RandFlat::shootInt(G4long(Navailable))+1;
//G4cout<<"Navailable Nsampled "<<Navailable<<" "<<Nsampled<<G4endl;
        Navailable=0;
        if(ProjExchangeQ != TargQ1) 
        {
         Navailable++;
         if(Navailable == Nsampled) 
         {TargExchangeQ = TargQ1; TargQ1=ProjExchangeQ; ProjQ1=TargExchangeQ;}
        }

        if(ProjExchangeQ != TargQ2)
        {
         Navailable++;
         if(Navailable == Nsampled)
         {TargExchangeQ = TargQ2; TargQ2=ProjExchangeQ; ProjQ1=TargExchangeQ;}
        }          

        if(ProjExchangeQ != TargQ3)
        {
         Navailable++;
         if(Navailable == Nsampled)
         {TargExchangeQ = TargQ3;  TargQ3=ProjExchangeQ; ProjQ1=TargExchangeQ;}
        }
       } else          // ProjQ2 is quark
       {  
        G4int Navailable=0;
        ProjExchangeQ = ProjQ2;
        if(ProjExchangeQ != TargQ1) Navailable++;
        if(ProjExchangeQ != TargQ2) Navailable++;
        if(ProjExchangeQ != TargQ3) Navailable++;

        G4int Nsampled=CLHEP::RandFlat::shootInt(G4long(Navailable))+1;
//G4cout<<"Navailable Nsampled "<<Navailable<<" "<<Nsampled<<G4endl;
        Navailable=0;
        if(ProjExchangeQ != TargQ1) 
        {
         Navailable++;
         if(Navailable == Nsampled) 
         {TargExchangeQ = TargQ1; TargQ1=ProjExchangeQ; ProjQ2=TargExchangeQ;}
        }

        if(ProjExchangeQ != TargQ2)
        {
         Navailable++;
         if(Navailable == Nsampled)
         {TargExchangeQ = TargQ2; TargQ2=ProjExchangeQ; ProjQ2=TargExchangeQ;}
        }

        if(ProjExchangeQ != TargQ3)
        {
         Navailable++;
         if(Navailable == Nsampled)
         {TargExchangeQ = TargQ3;  TargQ3=ProjExchangeQ; ProjQ2=TargExchangeQ;}
        }
       } // End of if(ProjQ1 > 0 ) // ProjQ1 is quark

//G4cout<<"Exch Pr Tr "<<ProjExchangeQ<<" "<<TargExchangeQ<<G4endl;
//G4cout<<ProjQ1<<" "<<ProjQ2<<" "<<ProjQ3<<G4endl;
//G4cout<<TargQ1<<" "<<TargQ2<<" "<<TargQ3<<G4endl;

       G4int aProjQ1=std::abs(ProjQ1);
       G4int aProjQ2=std::abs(ProjQ2);
       if(aProjQ1 == aProjQ2)          {NewProjCode = 111;} // Pi0-meson 
       else  // |ProjQ1| # |ProjQ2|
       {
        if(aProjQ1 > aProjQ2)          {NewProjCode = aProjQ1*100+aProjQ2*10+1;}
        else                           {NewProjCode = aProjQ2*100+aProjQ1*10+1;}
//      NewProjCode *=(ProjectilePDGcode/absProjectilePDGcode);
       }

G4bool ProjExcited=false;
//G4cout<<"NewProjCode "<<NewProjCode<<G4endl;
        if(G4UniformRand() < 0.5) 
        {
         NewProjCode +=2; // Excited Pi0-meson 
         ProjExcited=true;
        }
        if(aProjQ1 != aProjQ2) NewProjCode *=(ProjectilePDGcode/absProjectilePDGcode); // Uzhi
//G4cout<<"NewProjCode +2 or 0 "<<NewProjCode<<G4endl;

G4ParticleDefinition* TestParticle=0;
TestParticle=G4ParticleTable::GetParticleTable()->FindParticle(NewProjCode);
//G4cout<<"TestParticle ? "<<TestParticle<<G4endl;

if(TestParticle) 
{
 G4double MtestPart=                             // 31.05.2012
 (G4ParticleTable::GetParticleTable()->FindParticle(NewProjCode))->GetPDGMass();
/*
G4cout<<"TestParticle Name "<<NewProjCode<<" "<<TestParticle->GetParticleName()<<G4endl;
G4cout<<"MtestPart M0projectile projectile->GetDefinition()->GetPDGMass() "<<MtestPart<<" "<<M0projectile<<" "<<projectile->GetDefinition()->GetPDGMass()<<G4endl;
G4bool Test =M0projectile <= projectile->GetDefinition()->GetPDGMass(); 
G4cout<<"M0projectile <= projectile->GetDefinition()->GetPDGMass() "<<Test<<G4endl;
*/

  if(MtestPart > M0projectile) 
  {M0projectile = MtestPart;}
  else 
  {
   if(std::abs(M0projectile - projectile->GetDefinition()->GetPDGMass()) < 140.*MeV)
   {M0projectile = MtestPart;}
  }
//G4cout<<"M0projectile After check "<<M0projectile<<G4endl;
  M0projectile2 = M0projectile * M0projectile;

 ProjectileDiffStateMinMass   =M0projectile+210.*MeV; //210 MeV=m_pi+70 MeV 
 ProjectileNonDiffStateMinMass=M0projectile+210.*MeV; //210 MeV=m_pi+70 MeV
} else
{return false;}

//G4cout<<"New TrQ "<<TargQ1<<" "<<TargQ2<<" "<<TargQ3<<G4endl;
       NewTargCode = NewNucleonId(TargQ1, TargQ2, TargQ3);
//G4cout<<"NewTargCode "<<NewTargCode<<G4endl;

//       if( (TargQ1 != TargQ2) && (TargQ1 != TargQ3) && (TargQ2 != TargQ3) // Lambda or Sigma0
//       {if(G4UniformRand() < 0.5) NewTargCode=


       if( (TargQ1 == TargQ2) && (TargQ1 == TargQ3) &&
           (SqrtS > M0projectile+DeltaMass))           {NewTargCode +=2;  //Create Delta isobar
                                                        ProjExcited=true;}
       else if( target->GetDefinition()->GetPDGiIsospin() == 3 )          //Delta was the target
       { if(G4UniformRand() > DeltaProbAtQuarkExchange){NewTargCode +=2;  //Save   Delta isobar
                                                        ProjExcited=true;}
         else                                          {}     // De-excite initial Delta isobar
       }

//       else if((!CreateDelta)                               &&
       else if((!ProjExcited)                               &&
               (G4UniformRand() < DeltaProbAtQuarkExchange) &&         //Nucleon was the target
               (SqrtS > M0projectile+DeltaMass))      {NewTargCode +=2;}  //Create Delta isobar
//       else if( CreateDelta)                          {NewTargCode +=2;}
       else                                           {}                 //Save initial nucleon

//       target->SetDefinition(                                          // Fix 15.12.09
//       G4ParticleTable::GetParticleTable()->FindParticle(NewTargCode));// Fix 15.12.09 

TestParticle=G4ParticleTable::GetParticleTable()->FindParticle(NewTargCode);
//G4cout<<"New targ "<<NewTargCode<<" "<<TestParticle->GetParticleName()<<G4endl;
if(TestParticle) 
{
 G4double MtestPart=                             // 31.05.2012
 (G4ParticleTable::GetParticleTable()->FindParticle(NewTargCode))->GetPDGMass();

 if(MtestPart > M0target)
 {M0target=MtestPart;}
 else
 {
  if(std::abs(M0target - target->GetDefinition()->GetPDGMass()) < 140.*MeV)
  {M0target=MtestPart;}
 }

 TargetDiffStateMinMass   =M0target+220.*MeV;         //220 MeV=m_pi+80 MeV;    
 TargetNonDiffStateMinMass=M0target+220.*MeV;         //220 MeV=m_pi+80 MeV; 
} else
{return false;}
      } else 
      {    // projectile is baryon ----------------
//=========================================================================
       G4double Same=theParameters->GetProbOfSameQuarkExchange(); //0.3; //0.5; 0.
       G4bool ProjDeltaHasCreated(false);
       G4bool TargDeltaHasCreated(false);
 
       G4double Ksi=G4UniformRand();
       if(G4UniformRand() < 0.5)     // Sampling exchange quark from proj. or targ.
       {   // Sampling exchanged quark from the projectile ---
        if( Ksi < 0.333333 ) 
        {ProjExchangeQ = ProjQ1;}
        else if( (0.333333 <= Ksi) && (Ksi < 0.666667))
        {ProjExchangeQ = ProjQ2;}
        else
        {ProjExchangeQ = ProjQ3;}

//G4cout<<"ProjExchangeQ "<<ProjExchangeQ<<G4endl;
        if((ProjExchangeQ != TargQ1)||(G4UniformRand()<Same)) 
        {
         TargExchangeQ = TargQ1; TargQ1=ProjExchangeQ; ProjExchangeQ=TargExchangeQ;
        } else
        if((ProjExchangeQ != TargQ2)||(G4UniformRand()<Same)) 
        {
         TargExchangeQ = TargQ2; TargQ2=ProjExchangeQ; ProjExchangeQ=TargExchangeQ;
        } else 
        {
         TargExchangeQ = TargQ3;  TargQ3=ProjExchangeQ; ProjExchangeQ=TargExchangeQ;
        }

//G4cout<<"ProjExchangeQ "<<ProjExchangeQ<<G4endl;
//G4cout<<"TargExchangeQ "<<TargExchangeQ<<G4endl;
        if( Ksi < 0.333333 ) 
        {ProjQ1=ProjExchangeQ;}
        else if( (0.333333 <= Ksi) && (Ksi < 0.666667))
        {ProjQ2=ProjExchangeQ;}
        else
        {ProjQ3=ProjExchangeQ;}

       } else
       {   // Sampling exchanged quark from the target -------
        if( Ksi < 0.333333 ) 
        {TargExchangeQ = TargQ1;}
        else if( (0.333333 <= Ksi) && (Ksi < 0.666667)) 
        {TargExchangeQ = TargQ2;}
        else
        {TargExchangeQ = TargQ3;}

        if((TargExchangeQ != ProjQ1)||(G4UniformRand()<Same)) 
        {
         ProjExchangeQ = ProjQ1; ProjQ1=TargExchangeQ; TargExchangeQ=ProjExchangeQ;
        } else
        if((TargExchangeQ != ProjQ2)||(G4UniformRand()<Same)) 
        {
         ProjExchangeQ = ProjQ2; ProjQ2=TargExchangeQ; TargExchangeQ=ProjExchangeQ;
        } else 
        {
         ProjExchangeQ = ProjQ3;  ProjQ3=TargExchangeQ; TargExchangeQ=ProjExchangeQ;
        }

        if( Ksi < 0.333333 ) 
        {TargQ1=TargExchangeQ;}
        else if( (0.333333 <= Ksi) && (Ksi < 0.666667)) 
        {TargQ2=TargExchangeQ;}
        else
        {TargQ3=TargExchangeQ;}

       } // End of sampling baryon

       NewProjCode = NewNucleonId(ProjQ1, ProjQ2, ProjQ3); // *****************************

       if((ProjQ1==ProjQ2) && (ProjQ1==ProjQ3)) {NewProjCode +=2; ProjDeltaHasCreated=true;}
       else if(projectile->GetDefinition()->GetPDGiIsospin() == 3)// Projectile was Delta
       { if(G4UniformRand() > DeltaProbAtQuarkExchange)
                                                {NewProjCode +=2; ProjDeltaHasCreated=true;} 
         else                                   {NewProjCode +=0; ProjDeltaHasCreated=false;}
       }
       else                                                       // Projectile was Nucleon
       {
        if((G4UniformRand() < DeltaProbAtQuarkExchange) && (SqrtS > DeltaMass+M0target)) 
                                                {NewProjCode +=2; ProjDeltaHasCreated=true;}
        else                                    {NewProjCode +=0; ProjDeltaHasCreated=false;}
       } 

//G4ParticleDefinition* NewTestParticle=
//                      G4ParticleTable::GetParticleTable()->FindParticle(NewProjCode);
//G4cout<<"TestParticleMass NewTestParticle->GetPDGMass() "<<TestParticleMass<<" "<< NewTestParticle->GetPDGMass()<<G4endl;
//if(TestParticleMass < NewTestParticle->GetPDGMass()) {NewProjCode=TestParticleID;}
 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=

       NewTargCode = NewNucleonId(TargQ1, TargQ2, TargQ3); // *****************************

//G4cout<<"TargQ1, TargQ2, TargQ3 "<<TargQ1<<" "<<TargQ2<<" "<<TargQ3<<" "<<NewTargCode<<G4endl;

//TestParticleID=NewTargCode;
//TestParticleMass=DBL_MAX;

//TestParticle=G4ParticleTable::GetParticleTable()->FindParticle(NewTargCode);
//if(TestParticle) TestParticleMass=TestParticle->GetPDGMass(); 

       if((TargQ1==TargQ2) && (TargQ1==TargQ3)) {NewTargCode +=2; TargDeltaHasCreated=true;}  
       else if(target->GetDefinition()->GetPDGiIsospin() == 3)    // Target was Delta
       { if(G4UniformRand() > DeltaProbAtQuarkExchange)
                                                {NewTargCode +=2; TargDeltaHasCreated=true;} 
         else                                   {NewTargCode +=0; TargDeltaHasCreated=false;}
       }
       else                                                       // Target was Nucleon
       {
        if((G4UniformRand() < DeltaProbAtQuarkExchange) && (SqrtS > M0projectile+DeltaMass)) 
                                                {NewTargCode +=2; TargDeltaHasCreated=true;}
        else                                    {NewTargCode +=0; TargDeltaHasCreated=false;}
       }         

//NewTestParticle=G4ParticleTable::GetParticleTable()->FindParticle(NewTargCode);
//G4cout<<"TestParticleMass NewTestParticle->GetPDGMass() "<<TestParticleMass<<" "<< NewTestParticle->GetPDGMass()<<G4endl;
//if(TestParticleMass < NewTestParticle->GetPDGMass()) {NewTargCode=TestParticleID;}

//G4cout<<"NewProjCode NewTargCode "<<NewProjCode<<" "<<NewTargCode<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;

       if((absProjectilePDGcode == NewProjCode) && (absTargetPDGcode == NewTargCode))
       { // Nothing was changed! It is not right!?
       }
// Forming baryons --------------------------------------------------
if(ProjDeltaHasCreated) {ProbProjectileDiffraction=1.; ProbTargetDiffraction=0.;}
if(TargDeltaHasCreated) {ProbProjectileDiffraction=0.; ProbTargetDiffraction=1.;}
       if(ProjDeltaHasCreated) 
       {
        G4double MtestPart=                             // 31.05.2012
          (G4ParticleTable::GetParticleTable()->FindParticle(NewProjCode))->GetPDGMass();

        if(MtestPart >= M0projectile)                   // 31.05.2012
        {                                               // 31.05.2012
         M0projectile = MtestPart;                      // 31.05.2012
         M0projectile2 = M0projectile * M0projectile;   // 31.05.2012
        }                                               // 31.05.2012

        ProjectileDiffStateMinMass   =M0projectile+210.*MeV; //210 MeV=m_pi+70 MeV
        ProjectileNonDiffStateMinMass=M0projectile+210.*MeV; //210 MeV=m_pi+70 MeV
       }

//      if(M0target < 
//         (G4ParticleTable::GetParticleTable()->FindParticle(NewTargCode))->GetPDGMass())
       if(TargDeltaHasCreated)
       {
        G4double MtestPart=                             // 31.05.2012
          (G4ParticleTable::GetParticleTable()->FindParticle(NewTargCode))->GetPDGMass();

        if(MtestPart >=M0target)                        // 31.05.2012
        {                                               // 31.05.2012
         M0target=MtestPart;                            // 31.05.2012
         M0target2 = M0target * M0target;               // 31.05.2012
        }                                               // 31.05.2012

        TargetDiffStateMinMass   =M0target+210.*MeV;         //210 MeV=m_pi+70 MeV;    
        TargetNonDiffStateMinMass=M0target+210.*MeV;         //210 MeV=m_pi+70 MeV; 
       }
      } // End of if projectile is baryon ---------------------------

//G4cout<<"At end// NewProjCode "<<NewProjCode<<G4endl;
//G4cout<<"At end// NewTargCode "<<NewTargCode<<G4endl;

// If we assume that final state hadrons after the charge exchange will be
// in the ground states, we have to put ----------------------------------
//G4cout<<"M0pr M0tr SqS "<<M0projectile<<" "<<M0target<<" "<<SqrtS<<G4endl;

      PZcms2=(S*S+M0projectile2*M0projectile2+M0target2*M0target2-
             2*S*M0projectile2 - 2*S*M0target2 - 2*M0projectile2*M0target2)
             /4./S;
//G4cout<<"PZcms2 1 "<<PZcms2<<G4endl<<G4endl;
      if(PZcms2 < 0) {return false;}  // It can be if energy is not sufficient for Delta
//----------------------------------------------------------
      projectile->SetDefinition(
                  G4ParticleTable::GetParticleTable()->FindParticle(NewProjCode)); 

      target->SetDefinition(
                  G4ParticleTable::GetParticleTable()->FindParticle(NewTargCode)); 
//----------------------------------------------------------
      PZcms = std::sqrt(PZcms2);

      Pprojectile.setPz( PZcms);
      Pprojectile.setE(std::sqrt(M0projectile2+PZcms2));

      Ptarget.setPz(    -PZcms);
      Ptarget.setE(std::sqrt(M0target2+PZcms2));

// ----------------------------------------------------------

      if(absProjectilePDGcode < 1000)
      { // For projectile meson
       G4double Wexcit=1.-2.256*std::exp(-0.6*ProjectileRapidity);
       Wexcit=0.;
       if(G4UniformRand() > Wexcit)
       {                             // Make elastic scattering
//G4cout<<"Make elastic scattering of new hadrons"<<G4endl;
        Pprojectile.transform(toLab);
        Ptarget.transform(toLab);

        projectile->SetTimeOfCreation(target->GetTimeOfCreation());
        projectile->SetPosition(target->GetPosition());

        projectile->Set4Momentum(Pprojectile);
        target->Set4Momentum(Ptarget);

        G4bool Result= theElastic->ElasticScattering (projectile,target,theParameters);
//G4cout<<"Result of el. scatt "<<Result<<G4endl;
        return Result;
       } // end of if(Make elastic scattering for projectile meson?)
      } else
      { // For projectile baryon
       G4double Wexcit=1.-2.256*std::exp(-0.6*ProjectileRapidity);
       //Wexcit=0.;
       if(G4UniformRand() > Wexcit)
       {                             // Make elastic scattering
//G4cout<<"Make elastic scattering of new hadrons"<<G4endl;
        Pprojectile.transform(toLab);
        Ptarget.transform(toLab);

        projectile->SetTimeOfCreation(target->GetTimeOfCreation());
        projectile->SetPosition(target->GetPosition());

        projectile->Set4Momentum(Pprojectile);
        target->Set4Momentum(Ptarget);

        G4bool Result= theElastic->ElasticScattering (projectile,target,theParameters);
        return Result;
       } // end of if(Make elastic scattering for projectile baryon?)
      }
//G4cout<<"Make excitation of new hadrons"<<G4endl;
     }  // End of charge exchange part ------------------------------

// ------------------------------------------------------------------
     G4double ProbOfDiffraction=ProbProjectileDiffraction + ProbTargetDiffraction;
/*
G4cout<<"Excite --------------------"<<G4endl;
G4cout<<"Proj "<<M0projectile<<" "<<ProjectileDiffStateMinMass<<"  "<<ProjectileNonDiffStateMinMass<<G4endl;
G4cout<<"Targ "<<M0target    <<" "<<TargetDiffStateMinMass    <<" "<<TargetNonDiffStateMinMass<<G4endl;
G4cout<<"SqrtS "<<SqrtS<<G4endl;

G4cout<<"Prob ProjDiff TargDiff "<<ProbProjectileDiffraction<<" "<<ProbTargetDiffraction<<" "<<ProbOfDiffraction<<G4endl;
G4cout<<"Pr Y "<<Pprojectile.rapidity()<<" Tr Y "<<Ptarget.rapidity()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
*/
/*
     if(ProjectileNonDiffStateMinMass + TargetNonDiffStateMinMass > SqrtS) // 24.07.10
     {
      if(ProbOfDiffraction!=0.)
      {
       ProbProjectileDiffraction/=ProbOfDiffraction;
       ProbOfDiffraction=1.;
      } else {return false;}      
     } 

*/

     if(ProbOfDiffraction!=0.)
     {
      ProbProjectileDiffraction/=ProbOfDiffraction;
     }
     else
     {
      ProbProjectileDiffraction=0.;
     }

//G4cout<<"Prob ProjDiff TargDiff "<<ProbProjectileDiffraction<<" "<<ProbTargetDiffraction<<" "<<ProbOfDiffraction<<G4endl;

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
//ProbOfDiffraction=1.;                 // Uzhi Difr
//ProbProjectileDiffraction=1.;         // Uzhi 
     if(G4UniformRand() < ProbOfDiffraction)
       {
        if(G4UniformRand() < ProbProjectileDiffraction)
        { //-------- projectile diffraction ---------------
//G4cout<<"projectile diffraction"<<G4endl;

         do {
/*
Uzhi_projectilediffraction=1;
Uzhi_targetdiffraction=0;
Uzhi_Mx2=1.;
*/
//             Generate pt
//             if (whilecount++ >= 500 && (whilecount%100)==0)
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants possibly looping"
//	   	 << ", loop count/ maxPtSquare : "
//           	 << whilecount << " / " << maxPtSquare << G4endl;

//             whilecount++;
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
//G4cout<<"Masses "<<ProjMassT<<" "<<TargMassT<<" "<<SqrtS<<" "<<ProjMassT+TargMassT<<G4endl;
             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                     2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                    /4./S;

//G4cout<<"PZcms2 PrD"<<PZcms2<<G4endl;
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

//G4cout<<"Masses "<<ProjMassT<<" "<<TargMassT<<" "<<SqrtS<<" "<<ProjMassT+TargMassT<<G4endl;

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                     2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                    /4./S;
//G4cout<<"PZcms2 PrD"<<PZcms2<<G4endl;
             if(PZcms2 < 0 ) continue;
             PZcms =std::sqrt(PZcms2);

             PMinusMin=std::sqrt(ProjMassT2+PZcms2)-PZcms;
             PMinusMax=SqrtS-TargMassT;

             PMinusNew=ChooseP(PMinusMin, PMinusMax);
// PMinusNew=1./sqrt(1./PMinusMin-G4UniformRand()*(1./PMinusMin-1./PMinusMax));
//PMinusNew=1./sqr(1./std::sqrt(PMinusMin)-G4UniformRand()*(1./std::sqrt(PMinusMin)-1./std::sqrt(PMinusMax)));

             TMinusNew=SqrtS-PMinusNew;
             Qminus=Ptarget.minus()-TMinusNew;
             TPlusNew=TargMassT2/TMinusNew;
             Qplus=Ptarget.plus()-TPlusNew;

             Qmomentum.setPz( (Qplus-Qminus)/2 );
             Qmomentum.setE(  (Qplus+Qminus)/2 );

          } while ((Pprojectile+Qmomentum).mag2() <  ProjectileDiffStateMinMass2); //||  
                  //Repeat the sampling because there was not any excitation
//((Ptarget    -Qmomentum).mag2() <  M0target2                  )) );
          projectile->SetStatus(1*projectile->GetStatus());          // VU 10.04.2012
        }
        else
        { // -------------- Target diffraction ----------------

//G4cout<<"Target diffraction"<<G4endl;
         do {
/*
Uzhi_projectilediffraction=0;
Uzhi_targetdiffraction=1;
Uzhi_Mx2=1.;
*/
//             Generate pt
//             if (whilecount++ >= 500 && (whilecount%100)==0)
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants possibly looping"
//	   	 << ", loop count/ maxPtSquare : "
//           	 << whilecount << " / " << maxPtSquare << G4endl;

//             whilecount++;
             if (whilecount > 1000 )
             {
              Qmomentum=G4LorentzVector(0.,0.,0.,0.);
              target->SetStatus(2);  return false;    //  Ignore this interaction
             };
//G4cout<<"Qm while "<<Qmomentum<<" "<<whilecount<<G4endl;
// --------------- Check that the interaction is possible -----------
             ProjMassT2=M0projectile2;
             ProjMassT =M0projectile;

             TargMassT2=TargetDiffStateMinMass2;
             TargMassT =TargetDiffStateMinMass;

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                     2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                    /4./S;

//G4cout<<"PZcms2 TrD <0 "<<PZcms2<<" return"<<G4endl;
             if(PZcms2 < 0 ) 
             {
               target->SetStatus(2);  
               return false;
             }
             maxPtSquare=PZcms2;

             Qmomentum=G4LorentzVector(GaussianPt(AveragePt2,maxPtSquare),0);

//G4cout<<"Qm while "<<Qmomentum<<" "<<whilecount<<G4endl;
             Pt2=G4ThreeVector(Qmomentum.vect()).mag2();

             ProjMassT2=M0projectile2+Pt2;
             ProjMassT =std::sqrt(ProjMassT2);

             TargMassT2=TargetDiffStateMinMass2+Pt2;
             TargMassT =std::sqrt(TargMassT2);

             PZcms2=(S*S + ProjMassT2*ProjMassT2 + TargMassT2*TargMassT2-
                     2.*S*ProjMassT2-2.*S*TargMassT2-2.*ProjMassT2*TargMassT2)
                    /4./S;

//G4cout<<"PZcms2 <0 "<<PZcms2<<" continue"<<G4endl;
             if(PZcms2 < 0 ) continue;
             PZcms =std::sqrt(PZcms2);

             TPlusMin=std::sqrt(TargMassT2+PZcms2)-PZcms;
             TPlusMax=SqrtS-ProjMassT;

             TPlusNew=ChooseP(TPlusMin, TPlusMax);
//TPlusNew=1./sqr(1./std::sqrt(TPlusMin)-G4UniformRand()*(1./std::sqrt(TPlusMin)-1./std::sqrt(TPlusMax)));

//TPlusNew=TPlusMax;

             PPlusNew=SqrtS-TPlusNew;
             Qplus=PPlusNew-Pprojectile.plus();
             PMinusNew=ProjMassT2/PPlusNew;
             Qminus=PMinusNew-Pprojectile.minus();

             Qmomentum.setPz( (Qplus-Qminus)/2 );
             Qmomentum.setE(  (Qplus+Qminus)/2 );

/*
G4cout<<(Pprojectile+Qmomentum).mag()<<" "<<M0projectile<<G4endl;
G4bool First=(Pprojectile+Qmomentum).mag2() <  M0projectile2;
G4cout<<First<<G4endl;

G4cout<<(Ptarget    -Qmomentum).mag()<<" "<<TargetDiffStateMinMass<<" "<<TargetDiffStateMinMass2<<G4endl;
G4bool Seco=(Ptarget    -Qmomentum).mag2() < TargetDiffStateMinMass2;
G4cout<<Seco<<G4endl;
*/

         } while ((Ptarget    -Qmomentum).mag2() <  TargetDiffStateMinMass2);
                 // Repeat the sampling because there was not any excitation
// (((Pprojectile+Qmomentum).mag2() <  M0projectile2          ) ||  //No without excitation
//  ((Ptarget    -Qmomentum).mag2() <  TargetDiffStateMinMass2)) );
//G4cout<<"Go out"<<G4endl;
          target->SetStatus(1*target->GetStatus());     // VU 10.04.2012
         } // End of if(G4UniformRand() < ProbProjectileDiffraction)
        }
        else  //----------- Non-diffraction process ------------
        {

//G4cout<<"Non-diffraction process"<<G4endl;
         do {
/*
Uzhi_projectilediffraction=0;
Uzhi_targetdiffraction=0;
Uzhi_Mx2=1.;
*/
//             Generate pt
//             if (whilecount++ >= 500 && (whilecount%100)==0)
//	   	 G4cout << "G4DiffractiveExcitation::ExciteParticipants possibly looping"
//	   	 << ", loop count/ maxPtSquare : "
//           	 << whilecount << " / " << maxPtSquare << G4endl;

//             whilecount++;
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
//G4cout<<"PZcms2 ND"<<PZcms2<<G4endl;

             if(PZcms2 < 0 ) continue;
             PZcms =std::sqrt(PZcms2);

             PMinusMin=std::sqrt(ProjMassT2+PZcms2)-PZcms;
             PMinusMax=SqrtS-TargMassT;

             if(G4UniformRand() < ProbLogDistr)                       // Uzhi 25.04.2012
             {     PMinusNew=ChooseP(PMinusMin, PMinusMax);}               // 12.06.11
             else {PMinusNew=(PMinusMax-PMinusMin)*G4UniformRand() + PMinusMin;}
             Qminus=PMinusNew-Pprojectile.minus();

             TPlusMin=std::sqrt(TargMassT2+PZcms2)-PZcms;
             TPlusMax=SqrtS-PMinusNew;                   // Why is it closed???
//             TPlusMax=SqrtS-ProjMassT;      

             if(G4UniformRand() < 0.5) //ProbLogDistr)                 // Uzhi 29.05.2012 0.5)
             {     TPlusNew=ChooseP(TPlusMin, TPlusMax);}                   // 12.06.11
             else {TPlusNew=(TPlusMax-TPlusMin)*G4UniformRand() +TPlusMin;} 

             Qplus=-(TPlusNew-Ptarget.plus());

             Qmomentum.setPz( (Qplus-Qminus)/2 );
             Qmomentum.setE(  (Qplus+Qminus)/2 );
/*
G4cout<<(Pprojectile+Qmomentum).mag2()<<" "<<ProjectileNonDiffStateMinMass2<<G4endl;
G4cout<<(Ptarget    -Qmomentum).mag2()<<" "<<TargetNonDiffStateMinMass2<<G4endl;
G4int Uzhi; G4cin>>Uzhi;
*/
       } while ( 
 ((Pprojectile+Qmomentum).mag2() <  ProjectileNonDiffStateMinMass2) || //No double Diffraction
 ((Ptarget    -Qmomentum).mag2() <  TargetNonDiffStateMinMass2    ));

      projectile->SetStatus(0*projectile->GetStatus());     // VU 10.04.2012
      target->SetStatus(0*target->GetStatus());             // VU 10.04.2012
     }

     Pprojectile += Qmomentum;
     Ptarget     -= Qmomentum;

//G4cout<<"Pr Y "<<Pprojectile.rapidity()<<" Tr Y "<<Ptarget.rapidity()<<G4endl;

// Transform back and update SplitableHadron Participant.
     Pprojectile.transform(toLab);
     Ptarget.transform(toLab);

// Calculation of the creation time ---------------------
     projectile->SetTimeOfCreation(target->GetTimeOfCreation());
     projectile->SetPosition(target->GetPosition());
// Creation time and position of target nucleon were determined at
// ReggeonCascade() of G4FTFModel
// ------------------------------------------------------
/*
if(Uzhi_projectilediffraction != 0) 
{Uzhi_Mx2=Pprojectile.mag2(); Uzhi_modT=(target->Get4Momentum()-Ptarget).mag2();}

if(Uzhi_targetdiffraction     != 0) 
{Uzhi_Mx2=Ptarget.mag2(); Uzhi_modT=(projectile->Get4Momentum()-Pprojectile).mag2();}

if(Uzhi_QE!= 0) 
{
 Uzhi_projectilediffraction=0;
 Uzhi_targetdiffraction    =0;
 Uzhi_Mx2                  =1.;
}
*/
//G4cout<<"Mproj "<<Pprojectile.mag()<<G4endl;
//G4cout<<"Mtarg "<<Ptarget.mag()<<G4endl;
     projectile->Set4Momentum(Pprojectile);
     target->Set4Momentum(Ptarget);

     projectile->IncrementCollisionCount(1);
     target->IncrementCollisionCount(1);

     return true;
}

// ---------------------------------------------------------------------
void G4DiffractiveExcitation::CreateStrings(G4VSplitableHadron * hadron, 
                                            G4bool isProjectile,
                                            G4ExcitedString * &FirstString, 
                                            G4ExcitedString * &SecondString,
                                            G4FTFParameters *theParameters) const
{
/*
G4cout<<"Create Strings SplitUp "<<hadron<<G4endl;
G4cout<<"Defin "<<hadron->GetDefinition()<<G4endl;
G4cout<<"Defin "<<hadron->GetDefinition()->GetPDGEncoding()<<G4endl;
*/
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
//G4cout<<start<<" "<<start->GetPDGcode()<<" "<<end<<" "<<end->GetPDGcode()<<G4endl;
//G4cout<<"Create string "<<start->GetPDGcode()<<" "<<end->GetPDGcode()<<G4endl;
        G4LorentzVector Phadron=hadron->Get4Momentum();
//G4cout<<"String mom "<<Phadron<<G4endl;
        G4LorentzVector Pstart(0.,0.,0.,0.);
        G4LorentzVector Pend(0.,0.,0.,0.);
        G4LorentzVector Pkink(0.,0.,0.,0.);
        G4LorentzVector PkinkQ1(0.,0.,0.,0.);
        G4LorentzVector PkinkQ2(0.,0.,0.,0.);

        G4int PDGcode_startQ = std::abs(start->GetDefinition()->GetPDGEncoding());
        G4int PDGcode_endQ   = std::abs(  end->GetDefinition()->GetPDGEncoding());
//G4cout<<"PDGcode_startQ "<<PDGcode_startQ<<" PDGcode_endQ   "<<PDGcode_endQ  <<G4endl;

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
//G4cout<<"Wmin W "<<Wmin<<" "<<W<<G4endl;
        G4double W2=W*W;

        G4double Pt(0.), x1(0.), x3(0.); // x2(0.), 
        G4bool Kink=false;

       if(!(((start->GetDefinition()->GetParticleSubType() == "di_quark") &&
             (  end->GetDefinition()->GetParticleSubType() == "di_quark")   ) ||
            ((start->GetDefinition()->GetParticleSubType() == "quark")    &&
             (  end->GetDefinition()->GetParticleSubType() == "quark")      )))
       {   // Kinky strings are allowed only for qq-q strings
           // Kinky strings are impossible for other systems (qq-qqbar, q-qbar)
           // according to the analysis of Pbar P interactions
//G4cout<<G4endl<<"Check for Kink!##############"<<G4endl<<G4endl;
        if(W > Wmin)
        {                                        // Kink is possible
          if(hadron->GetStatus() == 0)           // VU 10.04.2012
          {
           G4double Pt2kink=theParameters->GetPt2Kink(); // For non-diffractive
           Pt = std::sqrt(Pt2kink*(std::pow(W2/16./Pt2kink+1.,G4UniformRand()) - 1.));
          }
          else
          {Pt=0.;}                                       // For diffractive

          if(Pt > 500.*MeV)
          {
             G4double Ymax = std::log(W/2./Pt + std::sqrt(W2/4./Pt/Pt - 1.));
             G4double Y=Ymax*(1.- 2.*G4UniformRand());

             x1=1.-Pt/W*std::exp( Y);
             x3=1.-Pt/W*std::exp(-Y);
//             x2=2.-x1-x3;

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

             if((P2_1 <= 0.) || (P2_3 <= 0.))
             { Kink=false;}
             else
             {
               G4double P_1=std::sqrt(P2_1);
               G4double P_2=std::sqrt(P2_2);
               G4double P_3=std::sqrt(P2_3);

               G4double CosT12=(P2_3-P2_1-P2_2)/(2.*P_1*P_2);
               G4double CosT13=(P2_2-P2_1-P2_3)/(2.*P_1*P_3);
//             Pt=P_2*std::sqrt(1.-CosT12*CosT12);  // because system was rotated 11.12.09

               if((std::abs(CosT12) >1.) || (std::abs(CosT13) > 1.)) 
               { Kink=false;}
               else
               { 
                 Kink=true; 
                 Pt=P_2*std::sqrt(1.-CosT12*CosT12);  // because system was rotated 11.12.09
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

                 G4double Cos2Psi=(sqr(x1) -sqr(x3)+2.*sqr(x3*CosT13))/
                          std::sqrt(sqr(sqr(x1)-sqr(x3)) + sqr(2.*x1*x3*CosT13));
                 G4double Psi=std::acos(Cos2Psi);

G4LorentzRotation Rotate;
if(isProjectile) {Rotate.rotateY(Psi);}
else             {Rotate.rotateY(pi-Psi);}                   
Rotate.rotateZ(twopi*G4UniformRand());

Pstart*=Rotate;
Pkink*=Rotate;
PkinkQ1*=Rotate;
PkinkQ2*=Rotate;
Pend*=Rotate;

               }
             }      // end of if((P2_1 < 0.) || (P2_3 <0.))
          }         // end of if(Pt > 500.*MeV)
        }           // end of if(W > Wmin) Check for a kink
       }            // end of qq-q string selection
//--------------------------------------------------------------------------------
/*
G4cout<<"Kink "<<Kink<<" "
<<start->GetDefinition()->GetParticleSubType()<<" "
<<  end->GetDefinition()->GetParticleSubType() <<G4endl;

G4cout<<"Kink "<<Kink<<" "
<<start->GetDefinition()->GetPDGEncoding()<<" "
<<  end->GetDefinition()->GetPDGEncoding() <<G4endl;
G4int Uzhi; G4cin>>Uzhi;
*/

        if(Kink)
        {                                        // Kink is possible
//G4cout<<"Kink is sampled!"<<G4endl;
          std::vector<G4double> QuarkProbabilitiesAtGluonSplitUp =
              theParameters->GetQuarkProbabilitiesAtGluonSplitUp();

          G4int QuarkInGluon(1); G4double Ksi=G4UniformRand();
          for(unsigned int Iq=0; Iq <3; Iq++)
          {
//G4cout<<"Iq "<<Iq<<G4endl;

if(Ksi > QuarkProbabilitiesAtGluonSplitUp[Iq]) QuarkInGluon++;}

//G4cout<<"Last Iq "<<QuarkInGluon<<G4endl;
          G4Parton * Gquark = new G4Parton(QuarkInGluon);
          G4Parton * Ganti_quark = new G4Parton(-QuarkInGluon);
//G4cout<<"Lorentz "<<G4endl;

//-------------------------------------------------------------------------------
          G4LorentzRotation toCMS(-1*Phadron.boostVector());

          G4LorentzRotation toLab(toCMS.inverse());
//G4cout<<"Pstart "<<Pstart<<G4endl;
//G4cout<<"Pend   "<<Pend<<G4endl;
          Pstart.transform(toLab);  start->Set4Momentum(Pstart);
          PkinkQ1.transform(toLab);
          PkinkQ2.transform(toLab);
          Pend.transform(toLab);    end->Set4Momentum(Pend);
//G4cout<<"Pstart "<<Pstart<<G4endl;
//G4cout<<"Pend   "<<Pend<<G4endl;
//G4cout<<"Defin "<<hadron->GetDefinition()<<G4endl;
//G4cout<<"Defin "<<hadron->GetDefinition()->GetPDGEncoding()<<G4endl;

//        G4int absPDGcode=std::abs(hadron->GetDefinition()->GetPDGEncoding());
          G4int absPDGcode=1500;  // 23 Dec
if((start->GetDefinition()->GetParticleSubType() == "quark") &&
   (  end->GetDefinition()->GetParticleSubType() == "quark")  )
  absPDGcode=110;

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

              } else
              {                            // Anti_Quark on the end
                FirstString = new G4ExcitedString(end        ,Gquark, +1);
                SecondString= new G4ExcitedString(Ganti_quark,start ,+1);
                Gquark->Set4Momentum(PkinkQ1);
                Ganti_quark->Set4Momentum(PkinkQ2);

              }   // end of if(end->GetPDGcode() > 0)
            } else {                      // Target
              if(end->GetDefinition()->GetPDGEncoding() > 0 )  // A quark on the end
              {                           // Quark on the end
                FirstString = new G4ExcitedString(Ganti_quark,end   ,-1);
                SecondString= new G4ExcitedString(start      ,Gquark,-1);
                Ganti_quark->Set4Momentum(PkinkQ2);
                Gquark->Set4Momentum(PkinkQ1);

              } else
              {                            // Anti_Quark on the end
                FirstString = new G4ExcitedString(Gquark,end        ,-1);
                SecondString= new G4ExcitedString(start ,Ganti_quark,-1);
                Gquark->Set4Momentum(PkinkQ2);
                Ganti_quark->Set4Momentum(PkinkQ1);

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

              } else
              {                            // Anti_DiQuark on the end or quark
                FirstString = new G4ExcitedString(end   ,Ganti_quark, +1);
                SecondString= new G4ExcitedString(Gquark,start      ,+1);
                Ganti_quark->Set4Momentum(PkinkQ1);
                Gquark->Set4Momentum(PkinkQ2);

              }   // end of if(end->GetPDGcode() > 0)
            } else {                      // Target

              if((end->GetDefinition()->GetParticleType() == "diquarks") &&
                 (end->GetDefinition()->GetPDGEncoding() > 0           )   ) 
              {                            // DiQuark on the end
                FirstString = new G4ExcitedString(Gquark,end        ,-1);

                SecondString= new G4ExcitedString(start ,Ganti_quark,-1);
                Gquark->Set4Momentum(PkinkQ1);
                Ganti_quark->Set4Momentum(PkinkQ2);

              } else
              {                            // Anti_DiQuark on the end or Q
                FirstString = new G4ExcitedString(Ganti_quark,end   ,-1);
                SecondString= new G4ExcitedString(start      ,Gquark,-1);
                Gquark->Set4Momentum(PkinkQ2);
                Ganti_quark->Set4Momentum(PkinkQ1);

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
//G4cout<<start<<" "<<start->GetPDGcode()<<" "<<end<<" "<<end->GetPDGcode()<<G4endl;
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

          G4LorentzVector Pstart1(tmp,0.);
          G4LorentzVector   Pend1(tmp,0.);

          if(isProjectile)
          {
           Pstart1*=(-1.)*Minus/2.;
           Pend1  *=(+1.)*Plus /2.;
          } 
          else
          {
           Pstart1*=(+1.)*Plus/2.;
           Pend1  *=(-1.)*Minus/2.;
          }

          Momentum=-Pstart1.mag();
          Pstart1.setT(Momentum);  // It is assumed that quark has m=0.

          Momentum=-Pend1.mag();
          Pend1.setT(Momentum);    // It is assumed that di-quark has m=0.

	  start->Set4Momentum(Pstart1);
	  end->Set4Momentum(Pend1);
          SecondString=0;
        }            // End of kink is impossible 

//G4cout<<"Quarks in the string at creation"<<FirstString->GetRightParton()->GetPDGcode()<<" "<<FirstString->GetLeftParton()->GetPDGcode()<<G4endl;
//G4cout<<FirstString<<" "<<SecondString<<G4endl;

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
//	G4double P=(Pmax-Pmin)*G4UniformRand()+Pmin;
	return P;
}

// ---------------------------------------------------------------------
G4ThreeVector G4DiffractiveExcitation::GaussianPt(G4double AveragePt2, 
                                                  G4double maxPtSquare) const
{            //  @@ this method is used in FTFModel as well. Should go somewhere common!

	G4double Pt2(0.);
        if(AveragePt2 <= 0.) {Pt2=0.;}
        else
        {
         Pt2 = -AveragePt2 * std::log(1. + G4UniformRand() * 
                            (std::exp(-maxPtSquare/AveragePt2)-1.));
        }
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
