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
// GEANT4 tag $Name:  $
//

#include <utility>                                        

#include "G4FTFParameters.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"             // 31 May 2011

#include "G4Proton.hh"                         // 31 May 2011
#include "G4Neutron.hh"                        // 31 May 2011

#include "G4PionPlus.hh"                       // 31 May 2011
#include "G4PionMinus.hh"                      // 31 May 2011
#include "G4KaonPlus.hh"                       // 31 May 2011
#include "G4KaonMinus.hh"                      // 31 May 2011

G4FTFParameters::G4FTFParameters() :
  //A.R. 14-Aug-2012 : Coverity fix.
  FTFhNcmsEnergy(0.0), 
  FTFxsManager(0),
  FTFXtotal(0.0), FTFXelastic(0.0), FTFXinelastic(0.0), FTFXannihilation(0.0),
  ProbabilityOfAnnihilation(0.0), ProbabilityOfElasticScatt(0.0),
  RadiusOfHNinteractions2(0.0), FTFSlope(0.0), 
  AvaragePt2ofElasticScattering(0.0), FTFGamma0(0.0),
  MagQuarkExchange(0.0), SlopeQuarkExchange(0.0), DeltaProbAtQuarkExchange(0.0),
  ProbOfSameQuarkExchange(0.0), ProjMinDiffMass(0.0), ProjMinNonDiffMass(0.0),
  ProbabilityOfProjDiff(0.0), TarMinDiffMass(0.0), TarMinNonDiffMass(0.0),
  ProbabilityOfTarDiff(0.0), AveragePt2(0.0), ProbLogDistr(0.0),
  Pt2kink(0.0),
  MaxNumberOfCollisions(0.0), ProbOfInelInteraction(0.0), CofNuclearDestruction(0.0),
  R2ofNuclearDestruction(0.0), ExcitationEnergyPerWoundedNucleon(0.0),
  DofNuclearDestruction(0.0), Pt2ofNuclearDestruction(0.0), MaxPt2ofNuclearDestruction(0.0) 
{}


G4FTFParameters::~G4FTFParameters()
{}
//**********************************************************************************************
G4FTFParameters::G4FTFParameters(const G4ParticleDefinition * particle, 
                                                   G4int      theA,
                                                   G4int      theZ,
                                                   G4double   PlabPerParticle) 
{

    //A.R. 25-Jul-2012 Coverity fix.
    FTFXannihilation = 0.0;
    FTFhNcmsEnergy = 0.0;
    ProbOfSameQuarkExchange = 0.0;

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

    G4double Plab = PlabPerParticle;
    G4double Elab = std::sqrt(Plab*Plab+ProjectileMass2);
    G4double KineticEnergy = Elab-ProjectileMass;                 // 31 May 2011

    G4double S=ProjectileMass2 + TargetMass2 + 2.*TargetMass*Elab;

//G4cout<<"Proj Plab "<<ProjectilePDGcode<<" "<<Plab<<G4endl;
//G4cout<<"Mass KinE "<<ProjectileMass<<" "<<KineticEnergy<<G4endl;
//G4cout<<" A Z "<<theA<<" "<<theZ<<G4endl;

    G4double Ylab,Xtotal,Xelastic,Xannihilation;
    G4int NumberOfTargetNucleons;

    Ylab=0.5*std::log((Elab+Plab)/(Elab-Plab));

    G4double ECMSsqr=S/GeV/GeV;
    G4double SqrtS  =std::sqrt(S)/GeV;
//G4cout<<"Sqrt(s) "<<SqrtS<<G4endl;

    TargetMass     /=GeV; TargetMass2     /=(GeV*GeV);
    ProjectileMass /=GeV; ProjectileMass2 /=(GeV*GeV);

    static G4ChipsComponentXS* _instance = new G4ChipsComponentXS();  // Witek Pokorski
    FTFxsManager = _instance;

    Plab/=GeV;
//  G4double LogPlab    = std::log( Plab );
//  G4double sqrLogPlab = LogPlab * LogPlab;

    G4int NumberOfTargetProtons  = theZ; 
    G4int NumberOfTargetNeutrons = theA-theZ;

    NumberOfTargetNucleons = NumberOfTargetProtons + NumberOfTargetNeutrons;

    if( (ProjectilePDGcode == 2212) || 
        (ProjectilePDGcode == 2112)   )    //------Projectile is nucleon --------
      {        
       G4double XtotPP = FTFxsManager->
                  GetTotalElementCrossSection(  particle,KineticEnergy,1,0);
       G4ParticleDefinition* Neutron=G4Neutron::Neutron();
       G4double XtotPN = FTFxsManager->
                  GetTotalElementCrossSection(   Neutron,KineticEnergy,1,0);


       G4double XelPP  = FTFxsManager->
                  GetElasticElementCrossSection(particle,KineticEnergy,1,0);
       G4double XelPN  = FTFxsManager->
                  GetElasticElementCrossSection(   Neutron,KineticEnergy,1,0);
//G4cout<<"Xs "<<XtotPP/millibarn<<" "<<XelPP/millibarn<<G4endl;
//G4cout<<"Xs "<<XtotPN/millibarn<<" "<<XelPN/millibarn<<G4endl;
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
       Xtotal/=millibarn;
       Xelastic/=millibarn;
      }
    else if( ProjectilePDGcode < -1000 )         //------Projectile is anti_baryon --------
      {        

       G4double X_a(0.), X_b(0.), X_c(0.), X_d(0.);
       G4double MesonProdThreshold=ProjectileMass+TargetMass+(2.*0.14+0.016); // 2 Mpi +DeltaE;

       if(PlabPerParticle < 40.*MeV)
       { // Low energy limits. Projectile at rest.
        Xtotal=   1512.9;    // mb
        Xelastic=  473.2;    // mb
        X_a=       625.1;    // mb
        X_b=         9.780;  // mb
        X_c=        49.989;  // mb
        X_d=         6.614;  // mb
       }
       else
       { // Total and elastic cross section of PbarP interactions a'la Arkhipov
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
//G4cout<<"Param Xtotal Xelastic "<<Xtotal<<" "<<Xelastic<<G4endl;
//G4cout<<"FlowF "<<FlowF<<" SqrtS "<<SqrtS<<G4endl;
//G4cout<<"Param Xelastic-NaN "<<Xelastic<<" "<<1.5*16.654/pow(ECMSsqr/2.176/2.176,2.2)<<" "<<ECMSsqr<<G4endl;
        X_a=25.*FlowF;               // mb, 3-shirts diagram

        if(SqrtS < MesonProdThreshold)
        {
         X_b=3.13+140.*std::pow(MesonProdThreshold-SqrtS,2.5);// mb anti-quark-quark annihilation
         Xelastic-=3.*X_b;  // Xel-X(PbarP->NNbar)
        } else
        {
         X_b=6.8/SqrtS;                                 // mb anti-quark-quark annihilation
         Xelastic-=3.*X_b;  // Xel-X(PbarP->NNbar)
        }

        X_c=2.*FlowF*sqr(ProjectileMass+TargetMass)/ECMSsqr; // mb rearrangement

//G4cout<<"Old new Xa "<<35.*FlowF<<" "<<25.*FlowF<<G4endl;

        X_d=23.3/ECMSsqr;                       // mb anti-quark-quark string creation
       }
//---------------------------------------------------------------
//G4cout<<"Param Xtotal Xelastic "<<Xtotal<<" "<<Xelastic<<G4endl;
//G4cout<<"Para a b c d "<<X_a<<" "<<X_b<<" "<<X_c<<" "<<X_d<<G4endl;
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
//G4int Uzhi; G4cin>>Uzhi;
*/
//---------------------------------------------------------------
      }
    else if( ProjectilePDGcode ==  211 )    //------Projectile is PionPlus -------
      {
       G4double XtotPiP = FTFxsManager->
                  GetTotalElementCrossSection(  particle,KineticEnergy,1,0); 
       G4ParticleDefinition* PionMinus=G4PionMinus::PionMinus();
       G4double XtotPiN =  FTFxsManager->
                  GetTotalElementCrossSection( PionMinus,KineticEnergy,1,0);
           
       G4double XelPiP  = FTFxsManager->
                  GetElasticElementCrossSection(particle,KineticEnergy,1,0); 
       G4double XelPiN  = FTFxsManager->
                  GetElasticElementCrossSection( PionMinus,KineticEnergy,1,0);

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP + 
                            NumberOfTargetNeutrons * XtotPiN  ) / NumberOfTargetNucleons;
       Xelastic         = ( NumberOfTargetProtons  * XelPiP  + 
                            NumberOfTargetNeutrons * XelPiN   ) / NumberOfTargetNucleons; 
       Xannihilation   = 0.;
       Xtotal/=millibarn;
       Xelastic/=millibarn;
      }
    else if( ProjectilePDGcode == -211 )            //------Projectile is PionMinus -------
      {
       G4double XtotPiP = FTFxsManager->
                  GetTotalElementCrossSection(  particle,KineticEnergy,1,0);
       G4ParticleDefinition* PionPlus=G4PionPlus::PionPlus();
       G4double XtotPiN = FTFxsManager->
                  GetTotalElementCrossSection(  PionPlus,KineticEnergy,1,0);
           
       G4double XelPiP  = FTFxsManager->
                  GetElasticElementCrossSection(particle,KineticEnergy,1,0);
       G4double XelPiN  = FTFxsManager->
                  GetElasticElementCrossSection(  PionPlus,KineticEnergy,1,0);

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP + 
                            NumberOfTargetNeutrons * XtotPiN  ) / NumberOfTargetNucleons;
       Xelastic         = ( NumberOfTargetProtons  * XelPiP  + 
                            NumberOfTargetNeutrons * XelPiN   ) / NumberOfTargetNucleons;
       Xannihilation   = 0.;
       Xtotal/=millibarn;
       Xelastic/=millibarn;
      }

    else if( ProjectilePDGcode ==  111 )          //------Projectile is PionZero  -------
      {
       G4ParticleDefinition* PionPlus=G4PionPlus::PionPlus();
       G4double XtotPipP= FTFxsManager->
                  GetTotalElementCrossSection(   PionPlus,KineticEnergy,1,0);

       G4ParticleDefinition* PionMinus=G4PionMinus::PionMinus();
       G4double XtotPimP= FTFxsManager->
                  GetTotalElementCrossSection(   PionMinus,KineticEnergy,1,0);
           
       G4double XelPipP = FTFxsManager->
                  GetElasticElementCrossSection(  PionPlus,KineticEnergy,1,0);
       G4double XelPimP = FTFxsManager->
                  GetElasticElementCrossSection( PionMinus,KineticEnergy,1,0);

       G4double XtotPiP= (XtotPipP + XtotPimP)/2.;
       G4double XtotPiN=XtotPiP;
       G4double XelPiP = (XelPipP  + XelPimP )/2.;
       G4double XelPiN = XelPiP;

       Xtotal           = ( NumberOfTargetProtons  * XtotPiP + 
                            NumberOfTargetNeutrons * XtotPiN  ) / NumberOfTargetNucleons;
       Xelastic         = ( NumberOfTargetProtons  * XelPiP  + 
                            NumberOfTargetNeutrons * XelPiN   ) / NumberOfTargetNucleons; 
       Xannihilation   = 0.;

       Xtotal/=millibarn;
       Xelastic/=millibarn;
      }
    else if( ProjectilePDGcode == 321 )             //------Projectile is KaonPlus -------
      {
       G4double XtotKP = FTFxsManager->
                  GetTotalElementCrossSection(  particle,KineticEnergy,1,0);

       G4ParticleDefinition* KaonMinus=G4KaonMinus::KaonMinus();
       G4double XtotKN = FTFxsManager->
                  GetTotalElementCrossSection( KaonMinus,KineticEnergy,1,0);

       G4double XelKP  = FTFxsManager->
                  GetElasticElementCrossSection(particle,KineticEnergy,1,0);
       G4double XelKN  = FTFxsManager->
                  GetElasticElementCrossSection( KaonMinus,KineticEnergy,1,0);

       Xtotal          = ( NumberOfTargetProtons  * XtotKP + 
                           NumberOfTargetNeutrons * XtotKN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelKP  + 
                           NumberOfTargetNeutrons * XelKN   ) / NumberOfTargetNucleons;
       Xannihilation   = 0.;

       Xtotal/=millibarn;
       Xelastic/=millibarn;
      }
    else if( ProjectilePDGcode ==-321 )             //------Projectile is KaonMinus ------
      {
       G4double XtotKP = FTFxsManager->
                  GetTotalElementCrossSection(  particle,KineticEnergy,1,0);

       G4ParticleDefinition* KaonPlus=G4KaonPlus::KaonPlus();
       G4double XtotKN = FTFxsManager->
                  GetTotalElementCrossSection(  KaonPlus,KineticEnergy,1,0);

       G4double XelKP  = FTFxsManager->
                  GetElasticElementCrossSection(particle,KineticEnergy,1,0);

       G4double XelKN  = FTFxsManager->
                  GetElasticElementCrossSection(KaonPlus,KineticEnergy,1,0); 

       Xtotal          = ( NumberOfTargetProtons  * XtotKP + 
                           NumberOfTargetNeutrons * XtotKN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelKP  + 
                           NumberOfTargetNeutrons * XelKN   ) / NumberOfTargetNucleons;
       Xannihilation   = 0.;
       
       Xtotal/=millibarn;
       Xelastic/=millibarn;
      }
    else if((ProjectilePDGcode == 311) || 
            (ProjectilePDGcode == 130) || 
            (ProjectilePDGcode == 310))               //Projectile is KaonZero
      {
       G4ParticleDefinition* KaonPlus=G4KaonPlus::KaonPlus();
       G4double XtotKpP= FTFxsManager->
                  GetTotalElementCrossSection(    KaonPlus,KineticEnergy,1,0);

       G4ParticleDefinition* KaonMinus=G4KaonMinus::KaonMinus();
       G4double XtotKmP= FTFxsManager->
                  GetTotalElementCrossSection(   KaonMinus,KineticEnergy,1,0);

       G4double XelKpP = FTFxsManager->
                  GetElasticElementCrossSection(  KaonPlus,KineticEnergy,1,0);
       G4double XelKmP = FTFxsManager->
                  GetElasticElementCrossSection(   KaonMinus,KineticEnergy,1,0);

       G4double XtotKP=(XtotKpP+XtotKmP)/2.;
       G4double XtotKN=XtotKP;
       G4double XelKP =(XelKpP +XelKmP )/2.; 
       G4double XelKN =XelKP;

       Xtotal          = ( NumberOfTargetProtons  * XtotKP + 
                           NumberOfTargetNeutrons * XtotKN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelKP  + 
                           NumberOfTargetNeutrons * XelKN   ) / NumberOfTargetNucleons;
       Xannihilation   = 0.;

       Xtotal/=millibarn;
       Xelastic/=millibarn;
      }
    else                 //------Projectile is undefined, Nucleon assumed
      {
       G4ParticleDefinition* Proton=G4Proton::Proton();
       G4double XtotPP = FTFxsManager->
                  GetTotalElementCrossSection(  Proton,KineticEnergy,1,0);

       G4ParticleDefinition* Neutron=G4Neutron::Neutron();
       G4double XtotPN = FTFxsManager->
                  GetTotalElementCrossSection( Neutron,KineticEnergy,1,0);

       G4double XelPP  = FTFxsManager->
                  GetElasticElementCrossSection(Proton,KineticEnergy,1,0);
       G4double XelPN  = FTFxsManager->
                  GetElasticElementCrossSection( Neutron,KineticEnergy,1,0);

       Xtotal          = ( NumberOfTargetProtons  * XtotPP + 
                           NumberOfTargetNeutrons * XtotPN  ) / NumberOfTargetNucleons;
       Xelastic        = ( NumberOfTargetProtons  * XelPP  + 
                           NumberOfTargetNeutrons * XelPN   ) / NumberOfTargetNucleons;
       Xannihilation   = 0.;

       Xtotal/=millibarn;
       Xelastic/=millibarn;
      };

//----------- Geometrical parameters ------------------------------------------------
      SetTotalCrossSection(Xtotal);
      SetElastisCrossSection(Xelastic);
      SetInelasticCrossSection(Xtotal-Xelastic);

/*
G4cout<<"Plab Xtotal, Xelastic Xinel Xftf "<<Plab<<" "<<Xtotal<<" "<<Xelastic<<" "<<Xtotal-Xelastic<<" "<<Xtotal-Xelastic-Xannihilation<<G4endl;
if(Xtotal-Xelastic != 0.)
{
  G4cout<<"Plab Xelastic/Xtotal,  Xann/Xin "<<Plab<<" "<<Xelastic/Xtotal<<" "<<Xannihilation/
  (Xtotal-Xelastic)<<G4endl;
} else 
{
  G4cout<<"Plab Xelastic/Xtotal,  Xann     "<<Plab<<" "<<Xelastic/Xtotal<<" "<<
  Xannihilation<<G4endl;
}
//G4int Uzhi; G4cin>>Uzhi;
*/
//  // Interactions with elastic and inelastic collisions
      SetProbabilityOfElasticScatt(Xtotal, Xelastic);
      SetRadiusOfHNinteractions2(Xtotal/pi/10.);

      if(Xtotal-Xelastic == 0.)
      {
       SetProbabilityOfAnnihilation(0.);
      } else
      {SetProbabilityOfAnnihilation(Xannihilation/(Xtotal-Xelastic));}
//
//SetProbabilityOfElasticScatt(Xtotal, 0.);
// //==== No elastic scattering ============================
//      SetProbabilityOfElasticScatt(Xtotal, 0.);
//      SetRadiusOfHNinteractions2((Xtotal-Xelastic)/pi/10.);
//      SetProbabilityOfAnnihilation(1.);
//        SetProbabilityOfAnnihilation(0.);
// //=======================================================

//-----------------------------------------------------------------------------------  

      SetSlope( Xtotal*Xtotal/16./pi/Xelastic/0.3894 ); // Slope parameter of elastic scattering
                                                        //      (GeV/c)^(-2))
//G4cout<<"Slope "<<GetSlope()<<G4endl;
//-----------------------------------------------------------------------------------
      SetGamma0( GetSlope()*Xtotal/10./2./pi );

//----------- Parameters of elastic scattering --------------------------------------
                                                        // Gaussian parametrization of
                                                        // elastic scattering amplitude assumed
      SetAvaragePt2ofElasticScattering(1./(Xtotal*Xtotal/16./pi/Xelastic/0.3894)*GeV*GeV);
//G4cout<<"AvaragePt2ofElasticScattering "<<GetAvaragePt2ofElasticScattering()<<G4endl;
//----------- Parameters of excitations ---------------------------------------------

            G4double Xinel=Xtotal-Xelastic;                        // Uzhi 25.04.2012
//G4cout<<"Param ProjectilePDGcode "<<ProjectilePDGcode<<G4endl;
           if( ProjectilePDGcode > 1000 )             //------Projectile is baryon --------
             {
              SetMagQuarkExchange(1.84);//(3.63);
              SetSlopeQuarkExchange(0.7);//(1.2);
              SetDeltaProbAtQuarkExchange(0.);
              if(NumberOfTargetNucleons > 26) {SetProbOfSameQuarkExchange(1.);}
              else                            {SetProbOfSameQuarkExchange(0.);}

              SetProjMinDiffMass(1.16);                              // GeV 
              SetProjMinNonDiffMass(1.16);                           // GeV 
//            SetProbabilityOfProjDiff(0.805*std::exp(-0.35*Ylab));  // Uzhi 21.05.2012
              SetProbabilityOfProjDiff(6./Xinel+1.5/ECMSsqr);        // Uzhi 25.04.2012

              SetTarMinDiffMass(1.16);                               // GeV
              SetTarMinNonDiffMass(1.16);                            // GeV 
//            SetProbabilityOfTarDiff(0.805*std::exp(-0.35*Ylab));   // Uzhi 21.05.2012
              SetProbabilityOfTarDiff(6./Xinel+1.5/ECMSsqr);         // Uzhi 25.04.2012
//            SetAveragePt2(0.15);                                   // 0.15 GeV^2
              SetAveragePt2(0.3);                         // 0.30 GeV^2 Uzhi 21.05.2012

              SetProbLogDistr(0.5);                                  // Uzhi 21.05.2012
             }
           else if( ProjectilePDGcode < -1000 )  //------Projectile is anti_baryon --------
             {
              SetMagQuarkExchange(0.);
              SetSlopeQuarkExchange(0.);
              SetDeltaProbAtQuarkExchange(0.);
              SetProbOfSameQuarkExchange(0.);

              SetProjMinDiffMass(ProjectileMass+0.22);               // GeV 
              SetProjMinNonDiffMass(ProjectileMass+0.22);            // GeV
//            SetProbabilityOfProjDiff(0.805*std::exp(-0.35*Ylab));  // Uzhi 21.05.2012
              SetProbabilityOfProjDiff(6./Xinel+1.5/ECMSsqr);        // Uzhi 25.04.2012

              SetTarMinDiffMass(TargetMass+0.22);                  // GeV
              SetTarMinNonDiffMass(TargetMass+0.22);               // GeV
//            SetProbabilityOfTarDiff(0.805*std::exp(-0.35*Ylab));   // Uzhi 21.05.2012
              SetProbabilityOfTarDiff(6./Xinel+1.5/ECMSsqr);         // Uzhi 25.04.2012

              SetAveragePt2(0.3);                   // 0.15 GeV^2    // Uzhi 21.05.2012

              SetProbLogDistr(0.5);                                  // Uzhi 21.05.2012
             }
           else if( ProjectileabsPDGcode == 211 || 
                    ProjectilePDGcode ==  111)     //------Projectile is Pion -----------
             {
              SetMagQuarkExchange(240.); 
              SetSlopeQuarkExchange(2.);         
              SetDeltaProbAtQuarkExchange(0.56); //(0.35);

              SetProjMinDiffMass(0.5);                               // GeV
              SetProjMinNonDiffMass(0.5);                            // GeV 0.3
//              SetProbabilityOfProjDiff(0.);                        // Uzhi 3.06.2012 
              SetProbabilityOfProjDiff((6.2-3.7*std::exp(-sqr(SqrtS-7.)/16.))/Xinel*0.);

              SetTarMinDiffMass(1.16);                               // GeV
              SetTarMinNonDiffMass(1.16);                            // GeV
//            SetProbabilityOfTarDiff(0.8*std::exp(-0.6*(Ylab-3.))); // Uzhi 3.06.2012
              SetProbabilityOfTarDiff((2.+22./ECMSsqr)/Xinel);

              SetAveragePt2(0.3);                                   // GeV^2 7 June 2011
              SetProbLogDistr(0.);                                   // Uzhi 21.05.2012

              SetProbLogDistr(1.);
             }
           else if( (ProjectileabsPDGcode == 321) || 
                    (ProjectileabsPDGcode == 311) || 
                    (ProjectilePDGcode == 130)    || 
                    (ProjectilePDGcode == 310))        //Projectile is Kaon
             {
              SetMagQuarkExchange(40.);
              SetSlopeQuarkExchange(2.25);
              SetDeltaProbAtQuarkExchange(0.6);

              SetProjMinDiffMass(0.6);                               // GeV 0.7 0.6
              SetProjMinNonDiffMass(0.6);                            // GeV 0.7 0.6
//            SetProbabilityOfProjDiff(0.85*std::pow(s/GeV/GeV,-0.5)); // 40/32 X-dif/X-inel
              SetProbabilityOfProjDiff(0.*4.7/Xinel);                   // Uzhi 5.06.2012

              SetTarMinDiffMass(1.1);                                // GeV
              SetTarMinNonDiffMass(1.1);                             // GeV
//            SetProbabilityOfTarDiff(0.45*std::pow(s/GeV/GeV,-0.5));// 40/32 X-dif/X-inel
              SetProbabilityOfTarDiff(1.5/Xinel);                    // Uzhi 5.06.2012
              SetAveragePt2(0.3);                                    // GeV^2 7 June 2011
              SetProbLogDistr(1.);                                   // Uzhi 5.06.2012
             }
           else                                           //------Projectile is undefined,
                                                          //------Nucleon assumed
             {
/*                 // Uzhi 6.06.2012
              SetMagQuarkExchange(1.85);       // 7 June 2011
              SetSlopeQuarkExchange(0.7);      // 7 June 2011
              SetDeltaProbAtQuarkExchange(0.); // 7 June 2011

              SetProjMinDiffMass((940.+160.*MeV)/GeV);     // particle->GetPDGMass()
              SetProjMinNonDiffMass((940.+160.*MeV)/GeV);  // particle->GetPDGMass()
              SetProbabilityOfProjDiff(0.805*std::pow(s/GeV/GeV,-0.35)); // 40/32 X-dif/X-inel

              SetTarMinDiffMass(1.16);                     // GeV
              SetTarMinNonDiffMass(1.16);                  // GeV
              SetProbabilityOfTarDiff(0.805*std::pow(s/GeV/GeV,-0.35)); // 40/32 X-dif/X-inel
*/

              SetMagQuarkExchange(0.);
              SetSlopeQuarkExchange(0.);
              SetDeltaProbAtQuarkExchange(0.);
              SetProbOfSameQuarkExchange(0.);

              SetProjMinDiffMass(ProjectileMass+0.22);               // GeV 
              SetProjMinNonDiffMass(ProjectileMass+0.22);            // GeV
              SetProbabilityOfProjDiff(6./Xinel+1.5/ECMSsqr);        // Uzhi 25.04.2012

              SetTarMinDiffMass(TargetMass+0.22);                    // GeV
              SetTarMinNonDiffMass(TargetMass+0.22);                 // GeV
              SetProbabilityOfTarDiff(6./Xinel+1.5/ECMSsqr);         // Uzhi 25.04.2012

              SetAveragePt2(0.3);                         // 0.15 GeV^2 Uzhi 21.05.2012
              SetProbLogDistr(0.5);                                  // Uzhi 21.05.2012
             }

//           if(theA > 4) SetProbabilityOfProjDiff(0.);     // Uzhi 6.07.2012 Closed

//G4cout<<"Param Get Min Dif "<<GetProjMinNonDiffMass()<<G4endl;

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

      SetDofNuclearDestruction(0.3);
      SetPt2ofNuclearDestruction((0.035+0.04*std::exp(4.*(Ylab-2.5))/
                                         (1.+std::exp(4.*(Ylab-2.5))))*GeV*GeV); //0.09
//G4cout<<"Parm Pt2 Y "<<(0.035+0.04*std::exp(4.*(Ylab-2.5))/(1.+std::exp(4.*(Ylab-2.5))))<<" "<<Ylab<<G4endl;
      SetMaxPt2ofNuclearDestruction(1.0*GeV*GeV);

      SetExcitationEnergyPerWoundedNucleon(100.*MeV);
    } else if( ProjectilePDGcode < -1000 )             // for anti-baryon projectile
    {
//G4cout<<"Nucl destruct Anti Bar"<<G4endl;

      SetMaxNumberOfCollisions(Plab,2.); //3.); ##############################
      SetCofNuclearDestruction(1.*std::exp(4.*(Ylab-2.1))/
                              (1.+std::exp(4.*(Ylab-2.1)))); //0.62 1.0

      SetR2ofNuclearDestruction(1.5*fermi*fermi);

      SetDofNuclearDestruction(0.3);
      SetPt2ofNuclearDestruction((0.035+0.04*std::exp(4.*(Ylab-2.5))/
                                         (1.+std::exp(4.*(Ylab-2.5))))*GeV*GeV); //0.09
      SetMaxPt2ofNuclearDestruction(1.0*GeV*GeV);

      SetExcitationEnergyPerWoundedNucleon(100.*MeV);

      if(Plab < 2.)   // 2 GeV/c
      { // For slow anti-baryon we have to garanty putting on mass-shell
       SetCofNuclearDestruction(0.);
       SetR2ofNuclearDestruction(1.5*fermi*fermi);
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

      SetDofNuclearDestruction(0.3);
      SetPt2ofNuclearDestruction((0.035+0.04*std::exp(4.*(Ylab-2.5))/
                                         (1.+std::exp(4.*(Ylab-2.5))))*GeV*GeV); //0.09
      SetMaxPt2ofNuclearDestruction(1.0*GeV*GeV);

      SetExcitationEnergyPerWoundedNucleon(100.*MeV);
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

//SetMaxNumberOfCollisions(Plab,6.); //(4.*(Plab+0.01),Plab); //6.); // ##########
//SetAveragePt2(0.15);
//G4cout<<"Cnd "<<GetCofNuclearDestruction()<<G4endl;
//SetCofNuclearDestruction(0.4);// (0.2); //(0.4);                  
//SetExcitationEnergyPerWoundedNucleon(0.*MeV); //(75.*MeV); 
//SetDofNuclearDestruction(0.);                  
//SetPt2ofNuclearDestruction(0.*GeV*GeV); //(0.168*GeV*GeV); 
//G4cout<<"Pt2 "<<GetPt2ofNuclearDestruction()/GeV/GeV<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
} 
//**********************************************************************************************
