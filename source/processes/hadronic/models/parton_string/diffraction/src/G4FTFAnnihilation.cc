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
// $Id: G4FTFAnnihilation.cc,v 1.1 2010-12-07 10:42:40 vuzhinsk Exp $
// ------------------------------------------------------------
//      GEANT 4 class implemetation file
//
//      ---------------- G4FTFAnnihilation --------------
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

#include "G4DiffractiveSplitableHadron.hh"
#include "G4DiffractiveExcitation.hh"
#include "G4FTFParameters.hh"
#include "G4ElasticHNScattering.hh"
#include "G4FTFAnnihilation.hh"

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

G4FTFAnnihilation::G4FTFAnnihilation()
{
}

// ---------------------------------------------------------------------
G4bool G4FTFAnnihilation::
          Annihilate(G4VSplitableHadron *projectile, 
                     G4VSplitableHadron *target,
                     G4VSplitableHadron *&AdditionalString,
                     G4FTFParameters    *theParameters) const  
{
// -------------------- Projectile parameters -----------------------
     G4LorentzVector Pprojectile=projectile->Get4Momentum();

     G4int    ProjectilePDGcode=projectile->GetDefinition()->GetPDGEncoding();
     if(ProjectilePDGcode > 0)
     {
       target->SetStatus(2);
       return false;
     } 

     G4double M0projectile = Pprojectile.mag();  
//     G4double M0projectile2 = M0projectile * M0projectile;

// -------------------- Target parameters -------------------------
     G4int    TargetPDGcode=target->GetDefinition()->GetPDGEncoding();

//G4cout<<"Annihilate "<<ProjectilePDGcode<<" "<<TargetPDGcode<<G4endl;

     G4LorentzVector Ptarget=target->Get4Momentum();
//G4cout<<"Ptarget "<<Ptarget<<G4endl;
     G4double M0target = Ptarget.mag();
//G4cout<<"M0target "<<M0target<<G4endl;
//     G4double M0target2 = M0target * M0target; 

     G4double AveragePt2=theParameters->GetAveragePt2();

//   G4double TargetRapidity = Ptarget.rapidity();

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

     G4double SqrtS=std::sqrt(S);

     G4double maxPtSquare; // = PZcms2;
/*
G4cout<<"Start --------------------"<<G4endl;
G4cout<<"SqrtS "<<SqrtS<<G4endl;
*/
//
//G4cout<<"M0projectile+M0target "<<M0projectile/GeV<<" "<<M0target/GeV<<" "<<(M0projectile+M0target)/GeV<<" "<<SqrtS/GeV<<G4endl;
     G4double FlowF=M0target/M0projectile/std::sqrt(S-
                (M0projectile+M0target)*(M0projectile+M0target))*GeV;
//G4cout<<"Annig FlowF "<<FlowF<<" sqrt "<<SqrtS/GeV<<G4endl;
     G4double X_a=100.*FlowF*(1.-1.12*GeV/SqrtS); // mb 3-shirt diagram
     G4double X_b= 50.*FlowF*std::pow(GeV*(1./(M0projectile+M0target) - 1./SqrtS),0.6); 
                                               // mb anti-quark-quark annihilation
     G4double X_e=10. *GeV*GeV/S+
                  197.*GeV*GeV*GeV/S/SqrtS-
                  181.*GeV*GeV*GeV*GeV/S/S;  
                                               // mb anti-quark-quark string creation
     G4double Xannihilation=X_a+X_b+X_e;

//G4cout<<"X a b e Annih "<<X_a<<" "<<X_b<<" "<<X_e<<" "<<Xannihilation<<G4endl;
// ------------------------------------------------------------------
// ------ Projectile unpacking --------------------------
     G4int AQ[3];
     UnpackBaryon(ProjectilePDGcode, AQ[0], AQ[1], AQ[2]);

// ------ Target unpacking ------------------------------
     G4int Q[3];
     UnpackBaryon(TargetPDGcode, Q[0], Q[1], Q[2]); 

//G4cout<<"Proj "<<AQ[0]<<" "<<AQ[1]<<" "<<AQ[2]<<G4endl;
//G4cout<<"Targ "<<Q[0]<<" "<<Q[1]<<" "<<Q[2]<<G4endl;

//     G4double Pt2;

     G4double Ksi=G4UniformRand();
//============================================================
// Simulation of anti-quark-quark string creation
//
     if(Ksi < X_e/Xannihilation)
     {
//G4cout<<"Process e"<<G4endl;
      G4int CandidatsN(0), CandAQ[9], CandQ[9];
      G4int LeftAQ(0), LeftQ(0);
//------------------------------------------------------------
      for(G4int iAQ1=0; iAQ1<3; iAQ1++)
      {
       for(G4int iAQ2=0; iAQ2<3; iAQ2++)
       {
        if(iAQ1 != iAQ2)
        {
         for(G4int iQ1=0; iQ1<3; iQ1++)
         {
          for(G4int iQ2=0; iQ2<3; iQ2++)
          {
           if(iQ1 != iQ2)
           {
            if((-AQ[iAQ1] == Q[iQ1]) && (-AQ[iAQ2] == Q[iQ2]))
            {
             if((iAQ1 == 0) && (iAQ2 == 1)){CandAQ[CandidatsN]=2;}
             if((iAQ1 == 1) && (iAQ2 == 0)){CandAQ[CandidatsN]=2;}

             if((iAQ1 == 0) && (iAQ2 == 2)){CandAQ[CandidatsN]=1;}
             if((iAQ1 == 2) && (iAQ2 == 0)){CandAQ[CandidatsN]=1;}

             if((iAQ1 == 1) && (iAQ2 == 2)){CandAQ[CandidatsN]=0;}
             if((iAQ1 == 2) && (iAQ2 == 1)){CandAQ[CandidatsN]=0;}
//----------------------------------------------------------------
             if((iQ1 == 0) && (iQ2 == 1)){CandQ[CandidatsN]=2;}
             if((iQ1 == 1) && (iQ2 == 0)){CandQ[CandidatsN]=2;}

             if((iQ1 == 0) && (iQ2 == 2)){CandQ[CandidatsN]=1;}
             if((iQ1 == 2) && (iQ2 == 0)){CandQ[CandidatsN]=1;}

             if((iQ1 == 1) && (iQ2 == 2)){CandQ[CandidatsN]=0;}
             if((iQ1 == 2) && (iQ2 == 1)){CandQ[CandidatsN]=0;}
             CandidatsN++;
            }//--------------------------
           } //end of if(jQ1 != jQ2)
          }  //end of for(G4int jQ2=0; j<3; j++)
         }   //end of for(G4int jQ=0; j<3; j++)
        }    //end of if(iAQ1 != iAQ2)
       }     //end of for(G4int iAQ2=0; i<3; i++)
      }      //end of for(G4int iAQ1=0; i<3; i++)
//------------------------------------------------------------

      if(CandidatsN != 0) 
      {
       CandidatsN--;
       G4int SampledCase=CLHEP::RandFlat::shootInt(G4long(CandidatsN));

       LeftAQ=AQ[CandAQ[SampledCase]];

       LeftQ =Q[CandQ[SampledCase]];

// --------------- Set the string properties ---------------
//G4cout<<"Left Aq Q "<<LeftAQ<<" "<<Q<<G4endl;
       projectile->SplitUp();

       projectile->SetFirstParton(LeftAQ);
       projectile->SetSecondParton(LeftQ);

       projectile->SetStatus(1);
       target->SetStatus(0);

       Pprojectile.setPz(0.);
       Pprojectile.setE(SqrtS);
       Pprojectile.transform(toLab);

// Calculation of the creation time ---------------------
       projectile->SetTimeOfCreation(target->GetTimeOfCreation());
       projectile->SetPosition(target->GetPosition());
// Creation time and position of target nucleon were determined at
// ReggeonCascade() of G4FTFModel
// ------------------------------------------------------

//G4cout<<"Mproj "<<Pprojectile.mag()<<G4endl;
//G4cout<<"Mtarg "<<Ptarget.mag()<<G4endl;
       projectile->Set4Momentum(Pprojectile);

       projectile->IncrementCollisionCount(1);
       return true;
      }  // end of if(CandidatsN != 0)
     }  // if(Ksi < X_e/Xannihilation)
//
//============================================================
// Simulation of anti-diquark-diquark string creation
//
     if(Ksi < (X_e+X_b)/Xannihilation)
     {
//G4cout<<"Process b"<<G4endl;
      G4int CandidatsN(0), CandAQ[9][2], CandQ[9][2];
      G4int LeftAQ1(0), LeftAQ2(0), LeftQ1(0), LeftQ2(0);
//------------------------------------------------------------
      for(G4int iAQ=0; iAQ<3; iAQ++)
      {
       for(G4int iQ=0; iQ<3; iQ++)
       {
        if(-AQ[iAQ] == Q[iQ])
        {
         if(iAQ == 0) {CandAQ[CandidatsN][0]=1; CandAQ[CandidatsN][1]=2;}
         if(iAQ == 1) {CandAQ[CandidatsN][0]=0; CandAQ[CandidatsN][1]=2;}
         if(iAQ == 2) {CandAQ[CandidatsN][0]=0; CandAQ[CandidatsN][1]=1;}
         if(iQ  == 0) {CandQ[CandidatsN][0] =1;  CandQ[CandidatsN][1]=2;}
         if(iQ  == 1) {CandQ[CandidatsN][0] =0;  CandQ[CandidatsN][1]=2;}
         if(iQ  == 2) {CandQ[CandidatsN][0] =0;  CandQ[CandidatsN][1]=1;}
         CandidatsN++;
        }  //end of if(-AQ[i] == Q[j])
       }   //end of cycle on targ. quarks
      }    //end of cycle on proj. anti-quarks
//------------------------------------------------------------

      if(CandidatsN != 0) 
      {
       CandidatsN--;
       G4int SampledCase=CLHEP::RandFlat::shootInt(G4long(CandidatsN));

       LeftAQ1=AQ[CandAQ[SampledCase][0]];
       LeftAQ2=AQ[CandAQ[SampledCase][1]];

       LeftQ1=Q[CandQ[SampledCase][0]];
       LeftQ2=Q[CandQ[SampledCase][1]];

// -------- Build anti-diquark and diquark
       G4int Anti_DQ(0), DQ(0);

       if(std::abs(LeftAQ1) > std::abs(LeftAQ2))
       { 
        Anti_DQ=1000*LeftAQ1+100*LeftAQ2-3;   // 1
       } else
       {
        Anti_DQ=1000*LeftAQ2+100*LeftAQ1-3;   // 1
       }
//       if(G4UniformRand() > 0.5) Anti_DQ-=2;
 
       if(std::abs(LeftQ1) > std::abs(LeftQ2))
       { 
        DQ=1000*LeftQ1+100*LeftQ2+3;  // 1
       } else
       {
        DQ=1000*LeftQ2+100*LeftQ1+3; // 1
       }
//       if(G4UniformRand() > 0.5) DQ+=2;

// --------------- Set the string properties ---------------
//G4cout<<"Left ADiQ DiQ "<<Anti_DQ<<" "<<DQ<<G4endl;
       projectile->SplitUp();

       projectile->SetFirstParton(Anti_DQ);
       projectile->SetSecondParton(DQ);

       projectile->SetStatus(1);
       target->SetStatus(0);

       Pprojectile.setPz(0.);
       Pprojectile.setE(SqrtS);
       Pprojectile.transform(toLab);

// Calculation of the creation time ---------------------
       projectile->SetTimeOfCreation(target->GetTimeOfCreation());
       projectile->SetPosition(target->GetPosition());
// Creation time and position of target nucleon were determined at
// ReggeonCascade() of G4FTFModel
// ------------------------------------------------------

//G4cout<<"Mproj "<<Pprojectile.mag()<<G4endl;
//G4cout<<"Mtarg "<<Ptarget.mag()<<G4endl;
       projectile->Set4Momentum(Pprojectile);

       projectile->IncrementCollisionCount(1);

       return true;
      }  // end of if(CandidatsN != 0)
     }  // if(Ksi < (X_e+X_b)/Xannihilation)
//
//============================================================
// Simulation of 3 anti-quark-quark strings creation
//   Sampling of anti-quark order in projectile

//G4cout<<"Process a"<<G4endl;
     G4int SampledCase=CLHEP::RandFlat::shootInt(G4long(5));

     G4int Tmp1(0), Tmp2(0);
     if(SampledCase == 0) {}
     if(SampledCase == 1) {Tmp1=AQ[1]; AQ[1]=AQ[2]; AQ[2]=Tmp1;}
     if(SampledCase == 2) {Tmp1=AQ[0]; AQ[0]=AQ[1]; AQ[1]=Tmp1;}
     if(SampledCase == 3) {Tmp1=AQ[0]; Tmp2=AQ[1]; AQ[0]=AQ[2]; AQ[1]=Tmp1; AQ[2]=Tmp2;}
     if(SampledCase == 4) {Tmp1=AQ[0]; Tmp2=AQ[1]; AQ[0]=Tmp2; AQ[1]=AQ[2]; AQ[2]=Tmp1;}
     if(SampledCase == 5) {Tmp1=AQ[0]; Tmp2=AQ[1]; AQ[0]=AQ[2]; AQ[1]=Tmp2; AQ[2]=Tmp1;}

// --------------- Set the string properties ---------------
//G4cout<<"String 1 "<<AQ[0]<<" "<<Q[0]<<G4endl;
       projectile->SplitUp();

       projectile->SetFirstParton(AQ[0]);
       projectile->SetSecondParton(Q[0]);
       projectile->SetStatus(1);

//G4cout<<"String 2 "<<Q[1]<<" "<<AQ[1]<<G4endl;
       target->SplitUp();

       target->SetFirstParton(Q[1]);
       target->SetSecondParton(AQ[1]);
       target->SetStatus(1);

//G4cout<<"String 3 "<<AQ[2]<<" "<<Q[2]<<G4endl;
AdditionalString=new G4DiffractiveSplitableHadron();
AdditionalString->SplitUp();
AdditionalString->SetFirstParton(AQ[2]);
AdditionalString->SetSecondParton(Q[2]);
AdditionalString->SetStatus(1);
//G4cout<<G4endl<<"*AdditionalString in Annih"<<AdditionalString<<G4endl;

// Sampling kinematical properties
// 1 string AQ[0]-Q[0]// 2 string AQ[1]-Q[1]// 3 string AQ[2]-Q[2]
     G4ThreeVector Quark_Mom[6];
     G4double         ModMom[6];
//     do   // while(SumMt_anti+SumMt_bary > SqrtS)
//     {
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//      G4double AveragePt2(200.),maxPtSquare(1000.);
      AveragePt2=200.*200.; maxPtSquare=S;
//AveragePt2*=4.; // ######################################################
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      G4double ModPtSum(0.), ProdMod(1.);
      do  // while((ModPtSum >SqrtS)||(ProdMod==0.))
      {
       G4ThreeVector PtSum(0.,0.,0.);
       for(G4int i=0; i<6; i++)
       {
        Quark_Mom[i]=GaussianPt(AveragePt2, maxPtSquare);
        PtSum+=Quark_Mom[i];
       }

       PtSum/=6.;

       ModPtSum=0.; ProdMod=1.;     
       for(G4int i=0; i<6; i++)
       {
        Quark_Mom[i]-=PtSum;
        ModMom[i]=Quark_Mom[i].mag();
        ModPtSum+=ModMom[i]; ProdMod*=ModMom[i];
       }
      }  while((ModPtSum >SqrtS)||(ProdMod==0.));

//------------------------------------
/*
// Sampling X's of anti-baryon -------
      G4bool Succes(true);
      do  // while(!Succes==false)
      {
       G4double SumX(0), Xupper(1.), SampledX(0.);
       for(G4int i=0; i<2; i++)
       {
        SampledX=Xupper*G4UniformRand();
        Quark_Mom[i].setZ(SampledX);
        Xupper=SampledX;
        SumX+=SampledX;
       }

       Quark_Mom[2].setZ(1.-SumX);

       for(G4int i=0; i<3; i++)
       {
        SampledX=Quark_Mom[i].getZ();
        if((SampledX==0.)||(SampledX==1.)) {Succes=false; break;}
       }
      } while(Succes==false);

// Sampling X's of baryon ------------
      if(Succes == true)
      {
       do  // while(Succes==false)
       {
        G4double SumX(0), Xupper(1.), SampledX(0.);
        for(G4int i=3; i<5; i++)
        {
         SampledX=Xupper*G4UniformRand();
         Quark_Mom[i].setZ(SampledX);
         Xupper=SampledX;
         SumX+=SampledX;
        }

        Quark_Mom[5].setZ(1.-SumX);

        for(G4int i=3; i<6; i++)
        {
         SampledX=Quark_Mom[i].getZ();
         if((SampledX==0.)||(SampledX==1.)) {Succes=false; break;}
        }
       } while(Succes==false);
      }  // if(Succes == true)
*/
       G4double SumMod_anti=ModMom[0]+ModMom[1]+ModMom[2];
       Quark_Mom[0].setZ(ModMom[0]/SumMod_anti);
       Quark_Mom[1].setZ(ModMom[1]/SumMod_anti);
       Quark_Mom[2].setZ(ModMom[2]/SumMod_anti);

       G4double SumMod_bary=ModMom[3]+ModMom[4]+ModMom[5];
       Quark_Mom[3].setZ(ModMom[3]/SumMod_bary);
       Quark_Mom[4].setZ(ModMom[4]/SumMod_bary);
       Quark_Mom[5].setZ(ModMom[5]/SumMod_bary);

       G4double M2anti=SumMod_anti*SumMod_anti;
       G4double M2bary=SumMod_bary*SumMod_bary;
       G4double DecayMomentum2= S*S+M2anti*M2anti+M2bary*M2bary
                          -2.*S*M2anti - 2.*S*M2bary -2.*M2anti*M2bary;

       G4double WminusTarget=(S-M2anti+M2bary+std::sqrt(DecayMomentum2))/2./SqrtS;
       G4double WplusProjectile=SqrtS - M2bary/WminusTarget;

       for(G4int i=0; i<3; i++)
       {
        G4double Pz=WplusProjectile*Quark_Mom[i].getZ()/2.-
                    ModMom[i]*ModMom[i]/(2.*WplusProjectile*Quark_Mom[i].getZ()); // SumMod_anti
        Quark_Mom[i].setZ(Pz);
//G4cout<<"Anti_Quarks mom "<<Quark_Mom[i]<<G4endl;
       }

       for(G4int i=3; i<6; i++)
       {
        G4double Pz=-WminusTarget*Quark_Mom[i].getZ()/2.+
                     ModMom[i]*ModMom[i]/(2.*WminusTarget*Quark_Mom[i].getZ()); // SumMod_bary
        Quark_Mom[i].setZ(Pz);
//G4cout<<"Quarks mom "<<Quark_Mom[i]<<G4endl;
       }
//G4cout<<Quark_Mom[0]+Quark_Mom[1]+Quark_Mom[2]<<G4endl;
//G4cout<<Quark_Mom[3]+Quark_Mom[4]+Quark_Mom[5]<<G4endl;
//-------------------------------------

       G4ThreeVector tmp=Quark_Mom[0]+Quark_Mom[3];
       G4LorentzVector Pstring1(tmp,Quark_Mom[0].mag()+Quark_Mom[3].mag());
       G4double        Ystring1=Pstring1.rapidity();
/*
G4cout<<"Mom 1 string "<<G4endl;
G4cout<<Quark_Mom[0]<<G4endl;
G4cout<<Quark_Mom[3]<<G4endl;
G4cout<<tmp<<" "<<tmp.mag()<<G4endl;
G4cout<<Pstring1<<" "<<Pstring1.mag()<<" "<<Ystring1<<G4endl;
*/
                     tmp=Quark_Mom[1]+Quark_Mom[4];
       G4LorentzVector Pstring2(tmp,Quark_Mom[1].mag()+Quark_Mom[4].mag());
       G4double        Ystring2=Pstring2.rapidity();
/*
G4cout<<"Mom 2 string "<<G4endl;
G4cout<<Quark_Mom[1]<<G4endl;
G4cout<<Quark_Mom[4]<<G4endl;
G4cout<<tmp<<" "<<tmp.mag()<<G4endl;
G4cout<<Pstring2<<" "<<Pstring2.mag()<<" "<<Ystring2<<G4endl;
*/

                     tmp=Quark_Mom[2]+Quark_Mom[5];
       G4LorentzVector Pstring3(tmp,Quark_Mom[2].mag()+Quark_Mom[5].mag());
       G4double        Ystring3=Pstring3.rapidity();
/*
G4cout<<"Mom 3 string "<<G4endl;
G4cout<<Quark_Mom[2]<<G4endl;
G4cout<<Quark_Mom[5]<<G4endl;
G4cout<<tmp<<" "<<tmp.mag()<<G4endl;
G4cout<<Pstring3<<" "<<Pstring3.mag()<<" "<<Ystring3<<G4endl;
*/
//--------------------------------
       G4LorentzVector LeftString(0.,0.,0.,0.);
//-----
       if((Ystring1 > Ystring2)&&(Ystring2 > Ystring3))
       {
        Pprojectile=Pstring1;
        LeftString =Pstring2;
        Ptarget    =Pstring3;
       }

       if((Ystring1 > Ystring3)&&(Ystring3 > Ystring2))
       {
        Pprojectile=Pstring1;
        LeftString =Pstring3;
        Ptarget    =Pstring2;
       }
//-----
       if((Ystring2 > Ystring1)&&(Ystring1 > Ystring3))
       {
        Pprojectile=Pstring2;
        LeftString =Pstring1;
        Ptarget    =Pstring3;
       }
    
       if((Ystring2 > Ystring3)&&(Ystring3 > Ystring1))
       {
        Pprojectile=Pstring2;
        LeftString =Pstring3;
        Ptarget    =Pstring1;
       }
//-----
       if((Ystring3 > Ystring1)&&(Ystring1 > Ystring2))
       {
        Pprojectile=Pstring3;
        LeftString =Pstring1;
        Ptarget    =Pstring2;
       }

       if((Ystring3 > Ystring2)&&(Ystring2 > Ystring1))
       {
        Pprojectile=Pstring3;
        LeftString =Pstring2;
        Ptarget    =Pstring1;
       }
//-----
/*
       G4double SumMt_anti(0.), SumMt_baryon(0.)
       G4double SumMt_anti=Quark_Mom[0].perp2()/Quark_Mom[0].getZ()+
                           Quark_Mom[1].perp2()/Quark_Mom[1].getZ()+
                           Quark_Mom[2].perp2()/Quark_Mom[2].getZ();

       G4double SumMt_bary=Quark_Mom[3].perp2()/Quark_Mom[3].getZ()+
                           Quark_Mom[4].perp2()/Quark_Mom[4].getZ()+
                           Quark_Mom[5].perp2()/Quark_Mom[5].getZ();
       if(SumMt_anti+SumMtbary > SqrtS) Succes=false;
*/
//     } while(!Succes);
//-------------------------------------------------------
//G4cout<<"SumP "<<Pprojectile+LeftString+Ptarget<<" "<<SqrtS<<G4endl;

     Pprojectile.transform(toLab);
     LeftString.transform(toLab);
     Ptarget.transform(toLab);
//G4cout<<"SumP "<<Pprojectile+LeftString+Ptarget<<" "<<SqrtS<<G4endl;

// Calculation of the creation time ---------------------
     projectile->SetTimeOfCreation(target->GetTimeOfCreation());
     projectile->SetPosition(target->GetPosition());

     AdditionalString->SetTimeOfCreation(target->GetTimeOfCreation());
     AdditionalString->SetPosition(target->GetPosition());
// Creation time and position of target nucleon were determined at
// ReggeonCascade() of G4FTFModel
// ------------------------------------------------------

//G4cout<<"Mproj "<<Pprojectile.mag()<<G4endl;
//G4cout<<"Mtarg "<<Ptarget.mag()<<G4endl;
     projectile->Set4Momentum(Pprojectile);
     AdditionalString->Set4Momentum(LeftString);
     target->Set4Momentum(Ptarget);

     projectile->IncrementCollisionCount(1);
     AdditionalString->IncrementCollisionCount(1);
     target->IncrementCollisionCount(1);

     return true;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//G4cout<<"Pr Y "<<Pprojectile.rapidity()<<" Tr Y "<<Ptarget.rapidity()<<G4endl;
}

// ---------------------------------------------------------------------
G4double G4FTFAnnihilation::ChooseX(G4double Alpha, G4double Beta) const
{
// choose an x between Xmin and Xmax with P(x) ~ 1/x
//  to be improved...

/*
	G4double range=Pmax-Pmin;                    

	if ( Pmin <= 0. || range <=0. )
	{
		G4cout << " Pmin, range : " << Pmin << " , " << range << G4endl;
		throw G4HadronicException(__FILE__, __LINE__, "G4DiffractiveExcitation::ChooseP : Invalid arguments ");
	}

	G4double P=Pmin * std::pow(Pmax/Pmin,G4UniformRand()); 
	return P;
*/
        G4double tmp=Alpha*Beta;
        tmp*=1.;
        return 0.5;
}

// ---------------------------------------------------------------------
G4ThreeVector G4FTFAnnihilation::GaussianPt(G4double AveragePt2, 
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
void G4FTFAnnihilation::UnpackBaryon(G4int IdPDG, 
                              G4int &Q1, G4int &Q2, G4int &Q3) const // Uzhi 7.09.09
    {
       G4int AbsId=std::abs(IdPDG);

       Q1 =  AbsId           / 1000;
       Q2 = (AbsId % 1000)  / 100;
       Q3 = (AbsId % 100)   / 10;
       
       if(IdPDG < 0 ) {Q1=-Q1; Q2=-Q2; Q3=-Q3;} // Anti-baryon
       
       return;
    }

// ---------------------------------------------------------------------
G4FTFAnnihilation::G4FTFAnnihilation(const G4FTFAnnihilation &)
{
	throw G4HadronicException(__FILE__, __LINE__, "G4FTFAnnihilation copy contructor not meant to be called");
}


G4FTFAnnihilation::~G4FTFAnnihilation()
{
}


const G4FTFAnnihilation & G4FTFAnnihilation::operator=(const G4FTFAnnihilation &)
{
	throw G4HadronicException(__FILE__, __LINE__, "G4FTFAnnihilation = operator meant to be called");
	return *this;
}


int G4FTFAnnihilation::operator==(const G4FTFAnnihilation &) const
{
	throw G4HadronicException(__FILE__, __LINE__, "G4FTFAnnihilation == operator meant to be called");
	return false;
}

int G4FTFAnnihilation::operator!=(const G4FTFAnnihilation &) const
{
	throw G4HadronicException(__FILE__, __LINE__, "G4DiffractiveExcitation != operator meant to be called");
	return true;
}
