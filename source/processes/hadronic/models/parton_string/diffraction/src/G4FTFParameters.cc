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
// $Id: G4FTFParameters.cc,v 1.16 2010/12/07 10:42:40 vuzhinsk Exp $
// GEANT4 tag $Name:  $
//

#include "G4FTFParameters.hh"

#include "G4ios.hh"
#include <utility>                                        

G4FTFParameters::G4FTFParameters()
{;}


G4FTFParameters::~G4FTFParameters()
{;}
//**********************************************************************************************
G4FTFParameters::G4FTFParameters(const G4ParticleDefinition * particle, 
                                                   G4int   theA,
                                                   G4int   theZ,
                                                   G4double   s) 
{
    G4int    ProjectilePDGcode    = particle->GetPDGEncoding();
    G4int    ProjectileabsPDGcode = std::abs(ProjectilePDGcode);
    G4double ProjectileMass       = particle->GetPDGMass();
    G4double ProjectileMass2      =ProjectileMass*ProjectileMass;

    G4int    ProjectileBaryonNumber(0), AbsProjectileBaryonNumber(0);
    G4int                               AbsProjectileCharge(0);
    G4bool   ProjectileIsNucleus=false;

    if(std::abs(particle->GetBaryonNumber()) > 1)
    { // The projectile is a nucleus
     ProjectileIsNucleus      =true;
     ProjectileBaryonNumber   =particle->GetBaryonNumber();
     AbsProjectileBaryonNumber=std::abs(ProjectileBaryonNumber);
     AbsProjectileCharge      =(G4int) particle->GetPDGCharge();

     if(ProjectileBaryonNumber > 1)
     {      ProjectilePDGcode= 2212; ProjectileabsPDGcode=2212;} // Proton
     else { ProjectilePDGcode=-2212; ProjectileabsPDGcode=2212;} // Anti-Proton

     ProjectileMass =G4Proton::Proton()->GetPDGMass();
     ProjectileMass2=sqr(ProjectileMass);
    } 

    G4double TargetMass     = G4Proton::Proton()->GetPDGMass();
    G4double TargetMass2    = TargetMass*TargetMass;

    G4double Elab = (s - ProjectileMass2 - TargetMass2)/(2*TargetMass);
    G4double Plab = std::sqrt(Elab * Elab - ProjectileMass2);

//G4cout<<"Proj S Plab "<<ProjectilePDGcode<<" "<<s/GeV/GeV<<" "<<Plab/GeV<<G4endl;
//G4cout<<" A Z "<<theA<<" "<<theZ<<G4endl;

    G4double Ylab,Xtotal,Xelastic,Xannihilation;
    G4int NumberOfTargetNucleons;

    Ylab=0.5*std::log((Elab+Plab)/(Elab-Plab));

    G4double ECMSsqr=s/GeV/GeV;
    G4double SqrtS  =std::sqrt(s)/GeV;

    TargetMass     /=GeV; TargetMass2     /=(GeV*GeV);
    ProjectileMass /=GeV; ProjectileMass2 /=(GeV*GeV);

    Plab/=GeV;
    G4double LogPlab    = std::log( Plab );
    G4double sqrLogPlab = LogPlab * LogPlab;

    G4int NumberOfTargetProtons  = theZ; 
    G4int NumberOfTargetNeutrons = theA-theZ;

    NumberOfTargetNucleons = NumberOfTargetProtons + NumberOfTargetNeutrons;

    if( ProjectilePDGcode > 1000 )                        //------Projectile is baryon --------
      {        
       G4double XtotPP = 48.0 +  0. *std::pow(Plab, 0.  ) + 0.522*sqrLogPlab - 4.51*LogPlab;
       G4double XtotPN = 47.3 +  0. *std::pow(Plab, 0.  ) + 0.513*sqrLogPlab - 4.27*LogPlab;

       G4double XelPP  = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;
       G4double XelPN  = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;

       if(!ProjectileIsNucleus)
       { // Projectile is hadron
        Xtotal          = ( NumberOfTargetProtons  * XtotPP + 
                            NumberOfTargetNeutrons * XtotPN  ) / NumberOfTargetNucleons;
        Xelastic        = ( NumberOfTargetProtons  * XelPP  + 
                            NumberOfTargetNeutrons * XelPN   ) / NumberOfTargetNucleons;
       } else
       { // Projectile is a nucleus
        Xtotal  = (
                 AbsProjectileCharge                           *NumberOfTargetProtons *XtotPP + 
                (AbsProjectileBaryonNumber-AbsProjectileCharge)*NumberOfTargetNeutrons*XtotPP + 
               ( AbsProjectileCharge                           *NumberOfTargetNeutrons +
                (AbsProjectileBaryonNumber-AbsProjectileCharge)*NumberOfTargetProtons)*XtotPN
                   )/(AbsProjectileBaryonNumber*NumberOfTargetNucleons);

        Xelastic= (
                 AbsProjectileCharge                           *NumberOfTargetProtons *XelPP + 
                (AbsProjectileBaryonNumber-AbsProjectileCharge)*NumberOfTargetNeutrons*XelPP + 
               ( AbsProjectileCharge                           *NumberOfTargetNeutrons +
                (AbsProjectileBaryonNumber-AbsProjectileCharge)*NumberOfTargetProtons)*XelPN
                   )/(AbsProjectileBaryonNumber*NumberOfTargetNucleons);
      }

       Xannihilation   = 0.;
      }
    else if( ProjectilePDGcode < -1000 )         //------Projectile is anti_baryon --------
      {        

// Total and elastic cross section of PbarP interactions a'la Arkhipov
       G4double LogS=std::log(ECMSsqr/33.0625);
       G4double Xasmpt=36.04+0.304*LogS*LogS;    // mb

                LogS=std::log(SqrtS/20.74);
       G4double Basmpt=11.92+0.3036*LogS*LogS;   // GeV^(-2)
       G4double R0=std::sqrt(0.40874044*Xasmpt-Basmpt); // GeV^(-1)

       G4double FlowF=SqrtS/
       std::sqrt(ECMSsqr*ECMSsqr+ProjectileMass2*ProjectileMass2+TargetMass2*TargetMass2-
       2.*ECMSsqr*ProjectileMass2 -2.*ECMSsqr*TargetMass2 -2.*ProjectileMass2*TargetMass2);

       Xtotal=Xasmpt*(1.+13.55*FlowF/R0/R0/R0*
                         (1.-4.47/SqrtS+12.38/ECMSsqr-12.43/SqrtS/ECMSsqr)); // mb

       Xasmpt=4.4+0.101*LogS*LogS;    // mb
       Xelastic=Xasmpt*(1.+59.27*FlowF/R0/R0/R0*
                         (1.-6.95/SqrtS+23.54/ECMSsqr-25.34/SqrtS/ECMSsqr)); // mb
//G4cout<<"Param Xelastic "<<Xelastic<<G4endl;
//G4cout<<"FlowF "<<FlowF<<" SqrtS "<<SqrtS<<G4endl;
//G4cout<<"Param Xelastic-NaN "<<Xelastic<<" "<<1.5*16.654/pow(ECMSsqr/2.176/2.176,2.2)<<" "<<ECMSsqr<<G4endl;
       G4double X_a=25.*FlowF;               // mb, 3-shirts diagram

       G4double X_b(0.);
       G4double MesonProdThreshold=ProjectileMass+TargetMass+
                                     (2.*0.14+0.016); // 2 Mpi +DeltaE
       if(SqrtS < MesonProdThreshold)
       {
        X_b=3.13+140.*pow(MesonProdThreshold-SqrtS,2.5);// mb anti-quark-quark annihilation
        Xelastic-=3.*X_b;  // Xel-X(PbarP->NNbar)
       } else
       {
        X_b=6.8/SqrtS;                                 // mb anti-quark-quark annihilation
        Xelastic-=3.*X_b;  // Xel-X(PbarP->NNbar)
       }

       G4double X_c=2.*FlowF*sqr(ProjectileMass+TargetMass)/ECMSsqr; // mb rearrangement

//G4cout<<"Old new Xa "<<35.*FlowF<<" "<<25.*FlowF<<G4endl;

       G4double X_d=23.3/ECMSsqr;                       // mb anti-quark-quark string creation

//---------------------------------------------------------------
//G4cout<<"Para a b c d "<<X_a<<" "<<5.*X_b<<" "<<5.*X_c<<" "<<6.*X_d<<G4endl;
       G4double Xann_on_P(0.), Xann_on_N(0.);

       if(ProjectilePDGcode == -2212)       // Pbar+P/N
       {Xann_on_P=X_a + X_b*5. + X_c*5. + X_d*6.; Xann_on_N=X_a + X_b*4. + X_c*4. + X_d*4.;} 
       else if(ProjectilePDGcode == -2112) // NeutrBar+P/N
       {Xann_on_P=X_a + X_b*4. + X_c*4. + X_d*4.; Xann_on_N=X_a + X_b*5. + X_c*5. + X_d*6.;} 
       else if(ProjectilePDGcode == -3122) // LambdaBar+P/N
       {Xann_on_P=X_a + X_b*3. + X_c*3. + X_d*2.; Xann_on_N=X_a + X_b*3. + X_c*3. + X_d*2.;} 
       else if(ProjectilePDGcode == -3112) // Sigma-Bar+P/N
       {Xann_on_P=X_a + X_b*2. + X_c*2. + X_d*0.; Xann_on_N=X_a + X_b*4. + X_c*4. + X_d*2.;} 
       else if(ProjectilePDGcode == -3212) // Sigma0Bar+P/N
       {Xann_on_P=X_a + X_b*3. + X_c*3. + X_d*2.; Xann_on_N=X_a + X_b*3. + X_c*3. + X_d*2.;} 
       else if(ProjectilePDGcode == -3222) // Sigma+Bar+P/N
       {Xann_on_P=X_a + X_b*4. + X_c*4. + X_d*2.; Xann_on_N=X_a + X_b*2. + X_c*2. + X_d*0.;} 
       else if(ProjectilePDGcode == -3312) // Xi-Bar+P/N
       {Xann_on_P=X_a + X_b*1. + X_c*1. + X_d*0.; Xann_on_N=X_a + X_b*2. + X_c*2. + X_d*0.;} 
       else if(ProjectilePDGcode == -3322) // Xi0Bar+P/N
       {Xann_on_P=X_a + X_b*2. + X_c*2. + X_d*0.; Xann_on_N=X_a + X_b*1. + X_c*1. + X_d*0.;} 
       else if(ProjectilePDGcode == -3334) // Omega-Bar+P/N
       {Xann_on_P=X_a + X_b*0. + X_c*0. + X_d*0.; Xann_on_N=X_a + X_b*0. + X_c*0. + X_d*0.;} 
       else {G4cout<<"Unknown anti-baryon for FTF annihilation"<<G4endl;}
//---------------------------------------------------------------

//G4cout<<"Sum          "<<Xann_on_P<<G4endl;

       if(!ProjectileIsNucleus)
       { // Projectile is anti-baryon
        Xannihilation   = ( NumberOfTargetProtons  * Xann_on_P  + 
                            NumberOfTargetNeutrons * Xann_on_N   ) / NumberOfTargetNucleons;
       } else
       { // Projectile is a nucleus
        Xannihilation=(
          ( AbsProjectileCharge                           *NumberOfTargetProtons+ 
           (AbsProjectileBaryonNumber-AbsProjectileCharge)*NumberOfTargetNeutrons )*Xann_on_P + 
          ( AbsProjectileCharge                           *NumberOfTargetNeutrons+
           (AbsProjectileBaryonNumber-AbsProjectileCharge)*NumberOfTargetProtons  )*Xann_on_N
                      )/(AbsProjectileBaryonNumber*NumberOfTargetNucleons);
       }
       G4double Xftf=0.;  
       MesonProdThreshold=ProjectileMass+TargetMass+(0.14+0.08); // Mpi +DeltaE
       if(SqrtS > MesonProdThreshold) {Xftf=36.*(1.-MesonProdThreshold/SqrtS);}

       Xtotal = Xelastic + Xannihilation + Xftf;
/*
G4cout<<"Plab Xtotal, Xelastic  Xinel Xftf "<<Plab<<" "<<Xtotal<<" "<<Xelastic<<" "<<Xtotal-Xelastic<<" "<<Xtotal-Xelastic-Xannihilation<<G4endl;
G4cout<<"Plab Xelastic/Xtotal,  Xann/Xin "<<Plab<<" "<<Xelastic/Xtotal<<" "<<Xannihilation/(Xtotal-Xelastic)<<G4endl;
G4int Uzhi; G4cin>>Uzhi;
*/
//---------------------------------------------------------------
      }
    else if( ProjectilePDGcode ==  211 )                     //------Projectile is PionPlus -------
      {
       G4double XtotPiP = 16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 0.0 *LogPlab;
       G4double XtotPiN = 33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab - 4.03*LogPlab;
           
       G4double XelPiP  =  0.0 + 11.4*std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 0.0 *LogPlab;
       G4double XelPiN  = 1.76 + 11.2*std::pow(Plab,-0.64) + 0.043*sqrLogPlab - 0.0 *LogPlab;

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP + 
                            NumberOfTargetNeutrons * XtotPiN  ) / NumberOfTargetNucleons;
       Xelastic         = ( NumberOfTargetProtons  * XelPiP  + 
                            NumberOfTargetNeutrons * XelPiN   ) / NumberOfTargetNucleons; 
       Xannihilation   = 0.;
      }
    else if( ProjectilePDGcode == -211 )            //------Projectile is PionMinus -------
      {
       G4double XtotPiP = 33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab - 4.03*LogPlab;
       G4double XtotPiN = 16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 0.0 *LogPlab;
           
       G4double XelPiP  = 1.76 + 11.2*std::pow(Plab,-0.64) + 0.043*sqrLogPlab - 0.0 *LogPlab;
       G4double XelPiN  =  0.0 + 11.4*std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 0.0 *LogPlab;

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP + 
                            NumberOfTargetNeutrons * XtotPiN  ) / NumberOfTargetNucleons;
       Xelastic         = ( NumberOfTargetProtons  * XelPiP  + 
                            NumberOfTargetNeutrons * XelPiN   ) / NumberOfTargetNucleons;
       Xannihilation   = 0.;
      }

    else if( ProjectilePDGcode ==  111 )          //------Projectile is PionZero  -------
      {
       G4double XtotPiP =(16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab - 
                                                              0.0 *LogPlab +    //Pi+
                          33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab -
                                                              4.03*LogPlab)/2;  //Pi-

       G4double XtotPiN =(33.0 + 14.0 *std::pow(Plab,-1.36) + 0.456*sqrLogPlab -
                                                              4.03*LogPlab +    //Pi+
                          16.4 + 19.3 *std::pow(Plab,-0.42) + 0.19 *sqrLogPlab -
                                                              0.0 *LogPlab)/2;  //Pi-
           
       G4double XelPiP  =( 0.0 + 11.4 *std::pow(Plab,-0.40) + 0.079*sqrLogPlab - 
                                                              0.0 *LogPlab +     //Pi+
                           1.76 +11.2 *std::pow(Plab,-0.64) + 0.043*sqrLogPlab -
                                                              0.0 *LogPlab)/2;  //Pi-
       G4double XelPiN  =( 1.76 +11.2 *std::pow(Plab,-0.64) + 0.043*sqrLogPlab -
                                                              0.0 *LogPlab +    //Pi+
                           0.0  +11.4 *std::pow(Plab,-0.40) + 0.079*sqrLogPlab -
                                                              0.0 *LogPlab)/2;  //Pi-

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP + 
                            NumberOfTargetNeutrons * XtotPiN  ) / NumberOfTargetNucleons;
       Xelastic         = ( NumberOfTargetProtons  * XelPiP  + 
                            NumberOfTargetNeutrons * XelPiN   ) / NumberOfTargetNucleons; 
       Xannihilation   = 0.;
      }
    else if( ProjectilePDGcode == 321 )             //------Projectile is KaonPlus -------
      {
       G4double XtotKP = 18.1 +  0. *std::pow(Plab, 0.  ) + 0.26 *sqrLogPlab - 1.0 *LogPlab;
       G4double XtotKN = 18.7 +  0. *std::pow(Plab, 0.  ) + 0.21 *sqrLogPlab - 0.89*LogPlab;

       G4double XelKP  =  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab;
       G4double XelKN  =  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotKP + 
                           NumberOfTargetNeutrons * XtotKN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelKP  + 
                           NumberOfTargetNeutrons * XelKN   ) / NumberOfTargetNucleons;
       Xannihilation   = 0.;
      }
    else if( ProjectilePDGcode ==-321 )             //------Projectile is KaonMinus ------
      {
       G4double XtotKP = 32.1 +  0. *std::pow(Plab, 0.  ) + 0.66 *sqrLogPlab - 5.6 *LogPlab;
       G4double XtotKN = 25.2 +  0. *std::pow(Plab, 0.  ) + 0.38 *sqrLogPlab - 2.9 *LogPlab;

       G4double XelKP  =  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab - 2.4 *LogPlab;
       G4double XelKN  =  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab - 1.3 *LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotKP + 
                           NumberOfTargetNeutrons * XtotKN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelKP  + 
                           NumberOfTargetNeutrons * XelKN   ) / NumberOfTargetNucleons;
       Xannihilation   = 0.;
      }
    else if((ProjectilePDGcode == 311) || 
            (ProjectilePDGcode == 130) || 
            (ProjectilePDGcode == 310))               //Projectile is KaonZero
      {
       G4double XtotKP =( 18.1 +  0. *std::pow(Plab, 0.  ) + 0.26 *sqrLogPlab -
                                                             1.0 *LogPlab +    //K+
                          32.1 +  0. *std::pow(Plab, 0.  ) + 0.66 *sqrLogPlab -
                                                             5.6 *LogPlab)/2;  //K-
       G4double XtotKN =( 18.7 +  0. *std::pow(Plab, 0.  ) + 0.21 *sqrLogPlab -
                                                             0.89*LogPlab +    //K+
                          25.2 +  0. *std::pow(Plab, 0.  ) + 0.38 *sqrLogPlab -
                                                             2.9 *LogPlab)/2;  //K-

       G4double XelKP  =(  5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab -
                                                             1.3 *LogPlab +    //K+
                           7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab -
                                                             2.4 *LogPlab)/2;  //K-
       G4double XelKN  =(  7.3 +  0. *std::pow(Plab,-0.  ) + 0.29 *sqrLogPlab -
                                                             2.4 *LogPlab +    //K+
                           5.0 +  8.1*std::pow(Plab,-1.8 ) + 0.16 *sqrLogPlab -
                                                             1.3 *LogPlab)/2;  //K-
       Xtotal          = ( NumberOfTargetProtons  * XtotKP + 
                           NumberOfTargetNeutrons * XtotKN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelKP  + 
                           NumberOfTargetNeutrons * XelKN   ) / NumberOfTargetNucleons;
       Xannihilation   = 0.;
      }
    else                 //------Projectile is undefined, Nucleon assumed
      {
       G4double XtotPP = 48.0 +  0. *std::pow(Plab, 0.  ) + 0.522*sqrLogPlab - 4.51*LogPlab;
       G4double XtotPN = 47.3 +  0. *std::pow(Plab, 0.  ) + 0.513*sqrLogPlab - 4.27*LogPlab;

       G4double XelPP  = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;
       G4double XelPN  = 11.9 + 26.9*std::pow(Plab,-1.21) + 0.169*sqrLogPlab - 1.85*LogPlab;

       Xtotal          = ( NumberOfTargetProtons  * XtotPP + 
                           NumberOfTargetNeutrons * XtotPN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelPP  + 
                           NumberOfTargetNeutrons * XelPN   ) / NumberOfTargetNucleons;
       Xannihilation   = 0.;
      };

//----------- Geometrical parameters ------------------------------------------------
      SetTotalCrossSection(Xtotal);
      SetElastisCrossSection(Xelastic);
      SetInelasticCrossSection(Xtotal-Xelastic);

/*
G4cout<<"Plab Xtotal, Xelastic Xinel Xftf "<<Plab<<" "<<Xtotal<<" "<<Xelastic<<" "<<Xtotal-Xelastic<<" "<<Xtotal-Xelastic-Xannihilation<<G4endl;
G4cout<<"Plab Xelastic/Xtotal,  Xann/Xin "<<Plab<<" "<<Xelastic/Xtotal<<" "<<Xannihilation/(Xtotal-Xelastic)<<G4endl;
G4int Uzhi; G4cin>>Uzhi;
*/
//  // Interactions with elastic and inelastic collisions
      SetProbabilityOfElasticScatt(Xtotal, Xelastic);
      SetRadiusOfHNinteractions2(Xtotal/pi/10.);

      SetProbabilityOfAnnihilation(Xannihilation/(Xtotal-Xelastic));
//
/* //==== No elastic scattering ============================
      SetProbabilityOfElasticScatt(Xtotal, 0.);
      SetRadiusOfHNinteractions2((Xtotal-Xelastic)/pi/10.);
      SetProbabilityOfAnnihilation(1.);
//        SetProbabilityOfAnnihilation(0.);
*/ //=======================================================

//-----------------------------------------------------------------------------------  

      SetSlope( Xtotal*Xtotal/16./pi/Xelastic/0.3894 ); // Slope parameter of elastic scattering
                                                        //      (GeV/c)^(-2))
//-----------------------------------------------------------------------------------
      SetGamma0( GetSlope()*Xtotal/10./2./pi );

//----------- Parameters of elastic scattering --------------------------------------
                                                        // Gaussian parametrization of
                                                        // elastic scattering amplitude assumed
      SetAvaragePt2ofElasticScattering(1./(Xtotal*Xtotal/16./pi/Xelastic/0.3894)*GeV*GeV);

//----------- Parameters of excitations ---------------------------------------------

//G4cout<<"Param ProjectilePDGcode "<<ProjectilePDGcode<<G4endl;
           if( ProjectilePDGcode > 1000 )             //------Projectile is baryon --------
             {
              SetMagQuarkExchange(1.84);//(3.63);
              SetSlopeQuarkExchange(0.7);//(1.2);
              SetDeltaProbAtQuarkExchange(0.);
              if(NumberOfTargetNucleons > 26) {SetProbOfSameQuarkExchange(1.);}
              else                            {SetProbOfSameQuarkExchange(0.);}

              SetProjMinDiffMass(1.16);                   // GeV 
              SetProjMinNonDiffMass(1.16);                // GeV 
              SetProbabilityOfProjDiff(0.805*std::exp(-0.35*Ylab));// 0.5

              SetTarMinDiffMass(1.16);                    // GeV
              SetTarMinNonDiffMass(1.16);                 // GeV 
              SetProbabilityOfTarDiff(0.805*std::exp(-0.35*Ylab));// 0.5

              SetAveragePt2(0.15);                        // 0.15 GeV^2
             }
           else if( ProjectilePDGcode < -1000 )  //------Projectile is anti_baryon --------
             {
              SetMagQuarkExchange(0.);
              SetSlopeQuarkExchange(0.);
              SetDeltaProbAtQuarkExchange(0.);
              SetProbOfSameQuarkExchange(0.);

              SetProjMinDiffMass(ProjectileMass+0.22);             // GeV 
              SetProjMinNonDiffMass(ProjectileMass+0.22);          // GeV
              SetProbabilityOfProjDiff(0.805*std::exp(-0.35*Ylab));// 0.5
//SetProbabilityOfProjDiff(0.5);
//G4cout<<"PrDif "<<GetProbabilityOfProjDiff()<<" "<<1.-2.*GetProbabilityOfProjDiff()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
              SetTarMinDiffMass(TargetMass+0.22);                  // GeV
              SetTarMinNonDiffMass(TargetMass+0.22);               // GeV
              SetProbabilityOfTarDiff(0.805*std::exp(-0.35*Ylab)); // 0.5
//SetProbabilityOfTarDiff(0.5);
              SetAveragePt2(0.15);                        // 0.15 GeV^2
             }
           else if( ProjectileabsPDGcode == 211 || 
                    ProjectilePDGcode ==  111)     //------Projectile is Pion -----------
             {
              SetMagQuarkExchange(240.); 
              SetSlopeQuarkExchange(2.);
              SetDeltaProbAtQuarkExchange(0.56); //(0.35);

              SetProjMinDiffMass(0.5);                    // GeV
              SetProjMinNonDiffMass(0.5);                 // GeV 0.3
              SetProbabilityOfProjDiff(0.);//(0.*0.62*std::pow(s/GeV/GeV,-0.51)); // 40/32 X-dif/X-inel

              SetTarMinDiffMass(1.16);                     // GeV
              SetTarMinNonDiffMass(1.16);                  // GeV
//              SetProbabilityOfTarDiff(1.);//(2.*0.62*std::pow(s/GeV/GeV,-0.51));
//              SetProbabilityOfTarDiff(2.6*std::exp(-0.46*Ylab));
              SetProbabilityOfTarDiff(0.8*std::exp(-0.6*(Ylab-3.)));

              SetAveragePt2(0.3);                         // GeV^2
             }
           else if( (ProjectileabsPDGcode == 321) || 
                    (ProjectileabsPDGcode == 311) || 
                    (ProjectilePDGcode == 130)    || 
                    (ProjectilePDGcode == 310))        //Projectile is Kaon
             {
// Must be corrected, taken from PiN
              SetMagQuarkExchange(120.);
              SetSlopeQuarkExchange(2.0);
              SetDeltaProbAtQuarkExchange(0.6);

              SetProjMinDiffMass(0.7);                    // GeV 1.1
              SetProjMinNonDiffMass(0.7);                 // GeV
              SetProbabilityOfProjDiff(0.85*std::pow(s/GeV/GeV,-0.5)); // 40/32 X-dif/X-inel

              SetTarMinDiffMass(1.1);                     // GeV
              SetTarMinNonDiffMass(1.1);                  // GeV
              SetProbabilityOfTarDiff(0.85*std::pow(s/GeV/GeV,-0.5)); // 40/32 X-dif/X-inel

              SetAveragePt2(0.3);                         // GeV^2
             }
           else                                           //------Projectile is undefined,
                                                          //------Nucleon assumed
             {
              SetMagQuarkExchange(3.5);
              SetSlopeQuarkExchange(1.0);
              SetDeltaProbAtQuarkExchange(0.1);

              SetProjMinDiffMass((particle->GetPDGMass()+160.*MeV)/GeV); 
              SetProjMinNonDiffMass((particle->GetPDGMass()+160.*MeV)/GeV);
              SetProbabilityOfProjDiff(0.95*std::pow(s/GeV/GeV,-0.35)); // 40/32 X-dif/X-inel

              SetTarMinDiffMass(1.1);                     // GeV
              SetTarMinNonDiffMass(1.1);                  // GeV
              SetProbabilityOfTarDiff(0.95*std::pow(s/GeV/GeV,-0.35)); // 40/32 X-dif/X-inel

              SetAveragePt2(0.3);                         // GeV^2
             }

// ---------- Set parameters of a string kink -------------------------------
             SetPt2Kink(6.*GeV*GeV);
             G4double Puubar(1./3.), Pddbar(1./3.), Pssbar(1./3.); // SU(3) symmetry
//           G4double Puubar(0.41 ), Pddbar(0.41 ), Pssbar(0.18 ); // Broken SU(3) symmetry
             SetQuarkProbabilitiesAtGluonSplitUp(Puubar, Pddbar, Pssbar);

// --------- Set parameters of nuclear destruction--------------------

    if( ProjectileabsPDGcode < 1000 )               // Meson projectile
    {
      SetMaxNumberOfCollisions(Plab,2.); //3.); ##############################
      SetCofNuclearDestruction(1.*std::exp(4.*(Ylab-2.1))/
                              (1.+std::exp(4.*(Ylab-2.1)))); //0.62 1.0

      SetR2ofNuclearDestruction(1.5*fermi*fermi);

      SetDofNuclearDestruction(0.4);
      SetPt2ofNuclearDestruction((0.035+0.04*std::exp(4.*(Ylab-2.5))/
                                         (1.+std::exp(4.*(Ylab-2.5))))*GeV*GeV); //0.09
      SetMaxPt2ofNuclearDestruction(1.0*GeV*GeV);

      SetExcitationEnergyPerWoundedNucleon(75.*MeV);
    } else if( ProjectilePDGcode < -1000 )             // for anti-baryon projectile
    {
//G4cout<<"Nucl destruct Anti Bar"<<G4endl;

      SetMaxNumberOfCollisions(Plab,2.); //3.); ##############################
      SetCofNuclearDestruction(1.*std::exp(4.*(Ylab-2.1))/
                              (1.+std::exp(4.*(Ylab-2.1)))); //0.62 1.0

      SetR2ofNuclearDestruction(1.5*fermi*fermi);

      SetDofNuclearDestruction(0.4);
      SetPt2ofNuclearDestruction((0.035+0.04*std::exp(4.*(Ylab-2.5))/
                                         (1.+std::exp(4.*(Ylab-2.5))))*GeV*GeV); //0.09
      SetMaxPt2ofNuclearDestruction(1.0*GeV*GeV);

      SetExcitationEnergyPerWoundedNucleon(75.*MeV);
      if(Plab < 2.*GeV)
      { // For slow anti-baryon we have to garanty putting on mass-shell
       SetCofNuclearDestruction(0.);
       SetR2ofNuclearDestruction(0.*fermi*fermi);
       SetDofNuclearDestruction(0.01);
       SetPt2ofNuclearDestruction(0.035*GeV*GeV);
       SetMaxPt2ofNuclearDestruction(0.04*GeV*GeV);

//       SetExcitationEnergyPerWoundedNucleon(0.);   // ?????
      }
    } else                                        // Projectile baryon assumed
    {
      SetMaxNumberOfCollisions(Plab,2.); //3.); ##############################
      SetCofNuclearDestruction(1.*std::exp(4.*(Ylab-2.1))/
                              (1.+std::exp(4.*(Ylab-2.1)))); //0.62 1.0

      SetR2ofNuclearDestruction(1.5*fermi*fermi);

      SetDofNuclearDestruction(0.4);
      SetPt2ofNuclearDestruction((0.035+0.04*std::exp(4.*(Ylab-2.5))/
                                         (1.+std::exp(4.*(Ylab-2.5))))*GeV*GeV); //0.09
      SetMaxPt2ofNuclearDestruction(1.0*GeV*GeV);

      SetExcitationEnergyPerWoundedNucleon(75.*MeV);
    }

//SetCofNuclearDestruction(0.47*std::exp(2.*(Ylab-2.5))/(1.+std::exp(2.*(Ylab-2.5)))); 
//SetPt2ofNuclearDestruction((0.035+0.1*std::exp(4.*(Ylab-3.))/(1.+std::exp(4.*(Ylab-3.))))*GeV*GeV);

//SetMagQuarkExchange(120.); // 210. PipP
//SetSlopeQuarkExchange(2.0);
//SetDeltaProbAtQuarkExchange(0.6);
//SetProjMinDiffMass(0.7);                    // GeV 1.1
//SetProjMinNonDiffMass(0.7);                 // GeV
//SetProbabilityOfProjDiff(0.85*std::pow(s/GeV/GeV,-0.5)); // 40/32 X-dif/X-inel
//SetTarMinDiffMass(1.1);                     // GeV
//SetTarMinNonDiffMass(1.1);                  // GeV
//SetProbabilityOfTarDiff(0.85*std::pow(s/GeV/GeV,-0.5)); // 40/32 X-dif/X-inel
//
//SetAveragePt2(0.3);                         // GeV^2
//------------------------------------
//SetProbabilityOfElasticScatt(1.,1.); //(Xtotal, Xelastic);
//SetProbabilityOfProjDiff(1.*0.62*std::pow(s/GeV/GeV,-0.51)); // 0->1
//SetProbabilityOfTarDiff(4.*0.62*std::pow(s/GeV/GeV,-0.51)); // 2->4
//SetAveragePt2(0.3);                              //(0.15);
//SetAvaragePt2ofElasticScattering(0.);

//SetMaxNumberOfCollisions(4.*(Plab+0.01),Plab); //6.); // ##############################
//SetCofNuclearDestruction(0.2); //(0.4);                  
//SetExcitationEnergyPerWoundedNucleon(0.*MeV); //(75.*MeV); 
//SetDofNuclearDestruction(0.4); //(0.4);                  
//SetPt2ofNuclearDestruction(0.1*GeV*GeV); //(0.168*GeV*GeV); 

} 
//**********************************************************************************************
