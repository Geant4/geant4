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
// $Id: G4Reggeons.cc 99348 2016-09-19 08:39:04Z vuzhinsk $
//

#include "G4Reggeons.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"
#include "G4Log.hh"


G4Reggeons::G4Reggeons(const G4ParticleDefinition * particle) 
{
//                    KP              Orig
 Alpha_pomeron      = 1.12;         //0.9808;
 Alphaprime_pomeron = 0.22/GeV/GeV; //0.25/GeV/GeV; 
 S0_pomeron         = 1.0*GeV*GeV;  //2.7*GeV*GeV; 

 Alpha_pomeronHard  = 1.47;
 Gamma_pomeronHard  = 0.0/GeV/GeV;

 G4int PDGcode = particle->GetPDGEncoding();
 G4int absPDGcode = std::abs(PDGcode);

//-------------------------------------------------------
//                           KP                               Orig
 G4double C_pomeron_NN     = 1.5; 
 G4double C_pomeron_N      = std::sqrt(C_pomeron_NN);

 G4double Gamma_pomeron_NN = 2.14/GeV/GeV;                 //(2.6+3.96)
 G4double Gamma_pomeron_N  = std::sqrt(Gamma_pomeron_NN);
 G4double Gamma_pomeron_Pr(0.), Gamma_pomeron_Tr(0.);

 G4double Rsquare_pomeron_NN =  3.30/GeV/GeV;              //3.56
 G4double Rsquare_pomeron_N  = Rsquare_pomeron_NN/2.;
 G4double Rsquare_pomeron_Pr(0.), Rsquare_pomeron_Tr(0.); 
//-------------------------------------------------------

 if( absPDGcode > 1000 ) { //  Projectile is baryon or anti_baryon --------
	Cpr_pomeron = C_pomeron_N;                              // Shower enhancement coefficient for projectile 
	Ctr_pomeron = C_pomeron_N;                              // Shower enhancement coefficient for target
        C_pomeron   = Cpr_pomeron*Ctr_pomeron;

	Gamma_pomeron_Pr = Gamma_pomeron_N;                     // vertex constant for projectile
	Gamma_pomeron_Tr = Gamma_pomeron_N;                     // vertex constant for target
	Gamma_pomeron    = Gamma_pomeron_Pr * Gamma_pomeron_Tr;

	Rsquare_pomeron_Pr = Rsquare_pomeron_N;                 // R^2 of pomeron-projectile interaction
	Rsquare_pomeron_Tr = Rsquare_pomeron_N;                 // R^2 of pomeron-target interaction
	Rsquare_pomeron    = Rsquare_pomeron_Pr + Rsquare_pomeron_Tr;

        Freggeon_Alpha      =   0.7;                            // Intersept of f-trajectory
        Freggeon_Alphaprime =   0.8/GeV/GeV;                    // Slope of f-trajectory
        Freggeon_Gamma      =   sqr(2.871)/GeV/GeV;             // Vertex constant of f-meson - nucleon interactions
        Freggeon_Rsquare    =   2*0.916/GeV/GeV;                // R^2 of f-meson - nucleon interactions
        Freggeon_C          =   1.0;                            // Shower enhancement coefficient
        FParity             =  +1;                              // Parity of the trajectory

        Wreggeon_Alpha      =  0.4;                             // Intersept of omega-trajectory (w)
        Wreggeon_Alphaprime =  0.9/GeV/GeV;                     // Slope of w-trajectory
        Wreggeon_Gamma      =  sqr(2.241)/GeV/GeV;              // Vertex constant of w-meson - nucleon interactions
        Wreggeon_Rsquare    =  2*0.945/GeV/GeV *0.5;            // R^2 of w-meson - nucleon interactions
        Wreggeon_C          =   1.0;                            // Shower enhancement coefficient
        if(PDGcode > 0) WParity = -1;                           // Parity +1 for Pbar P, and -1 for PP interactions
        if(PDGcode < 0) WParity = +1;
 }
 else if ( absPDGcode == 211  ||  PDGcode ==  111 ) {  // Projectile is Pion
	Cpr_pomeron = 1.352; 
	Ctr_pomeron = C_pomeron_N; 
        C_pomeron   = Cpr_pomeron*Ctr_pomeron;
//                         KP
	Gamma_pomeron_Pr = 0.89/GeV;   // 0.85 -> 0.89    // Uzhi
	Gamma_pomeron_Tr = Gamma_pomeron_N;
	Gamma_pomeron    = Gamma_pomeron_Pr * Gamma_pomeron_Tr;

	Rsquare_pomeron_Pr = 0.5/GeV/GeV;
	Rsquare_pomeron_Tr = Rsquare_pomeron_N;
	Rsquare_pomeron    = Rsquare_pomeron_Pr + Rsquare_pomeron_Tr;

        Freggeon_Alpha      =   0.7;
        Freggeon_Gamma      =   3.524/GeV/GeV;
        Freggeon_Rsquare    =   1.0/GeV/GeV;
        Freggeon_Alphaprime =   0.8/GeV/GeV;
        Freggeon_C          =   1.0;
        FParity             =  +1;

        Wreggeon_Alpha      =   0.5;
        Wreggeon_Gamma      =   0.56/GeV/GeV;   // 1.12 -> 0.56 Uzhi
        Wreggeon_Rsquare    =   9.19/GeV/GeV;
        Wreggeon_Alphaprime =   0.9/GeV/GeV;
        Wreggeon_C          =   1.0;
        if(PDGcode > 0) WParity = -1;
        if(PDGcode < 0) WParity = +1;
 }
 else if ( absPDGcode == 321  ||  absPDGcode == 311  || 
              PDGcode == 130  ||  PDGcode == 310 )      {  // Projectile is Kaon

	Cpr_pomeron = 1.522; 
	Ctr_pomeron = C_pomeron_N; 
        C_pomeron   = Cpr_pomeron*Ctr_pomeron;

	Gamma_pomeron_Pr = 1.312/GeV;
	Gamma_pomeron_Tr = Gamma_pomeron_N;
	Gamma_pomeron    = Gamma_pomeron_Pr * Gamma_pomeron_Tr;

	Rsquare_pomeron_Pr = 0.31/GeV/GeV;
	Rsquare_pomeron_Tr = Rsquare_pomeron_N;
	Rsquare_pomeron    = Rsquare_pomeron_Pr + Rsquare_pomeron_Tr;

        Freggeon_Alpha      =   0.0;
        Freggeon_Gamma      =   0.0/GeV/GeV;
        Freggeon_Rsquare    =   1.0/GeV/GeV;
        Freggeon_Alphaprime =   0.0/GeV/GeV;
        Freggeon_C          =   1.0;
        FParity             =  +1;

        Wreggeon_Alpha      =   0.0;
        Wreggeon_Gamma      =   0.0/GeV/GeV;
        Wreggeon_Rsquare    =   1.0/GeV/GeV;
        Wreggeon_Alphaprime =   0.0/GeV/GeV;
        Wreggeon_C          =   1.0;
        WParity             =  -1;
 }
 else if ( absPDGcode == 22 ) {                 // Projectile is Gamma

	Cpr_pomeron = 1.437; 
	Ctr_pomeron = C_pomeron_N; 
        C_pomeron   = Cpr_pomeron*Ctr_pomeron;

	Gamma_pomeron_Pr = 1.415/GeV/GeV;
	Gamma_pomeron_Tr = Gamma_pomeron_N;
	Gamma_pomeron    = Gamma_pomeron_Pr * Gamma_pomeron_Tr;

	Rsquare_pomeron_Pr = 0.51/GeV/GeV;
	Rsquare_pomeron_Tr = Rsquare_pomeron_N;
	Rsquare_pomeron    = Rsquare_pomeron_Pr + Rsquare_pomeron_Tr;

        Freggeon_Alpha      =   0.0;
        Freggeon_Gamma      =   0.0/GeV/GeV;
        Freggeon_Rsquare    =   1.0/GeV/GeV;
        Freggeon_Alphaprime =   0.0/GeV/GeV;
        Freggeon_C          =   1.0;
        FParity             =  +1;

        Wreggeon_Alpha      =   0.0;
        Wreggeon_Gamma      =   0.0/GeV/GeV;
        Wreggeon_Rsquare    =   1.0/GeV/GeV;
        Wreggeon_Alphaprime =   0.0/GeV/GeV;
        Wreggeon_C          =   1.0;
        WParity             =  -1;
 }
 else {  // Projectile is undefined, Nucleon assumed
	Cpr_pomeron = C_pomeron_N; 
	Ctr_pomeron = C_pomeron_N; 
        C_pomeron   = Cpr_pomeron*Ctr_pomeron;

	Gamma_pomeron_Pr = Gamma_pomeron_N;
	Gamma_pomeron_Tr = Gamma_pomeron_N;
	Gamma_pomeron    = Gamma_pomeron_Pr * Gamma_pomeron_Tr;

	Rsquare_pomeron_Pr = Rsquare_pomeron_N;
	Rsquare_pomeron_Tr = Rsquare_pomeron_N;
	Rsquare_pomeron    = Rsquare_pomeron_Pr + Rsquare_pomeron_Tr;

        Freggeon_Alpha      =   0.723;
        Freggeon_Gamma      =   8.801/GeV/GeV;
        Freggeon_Rsquare    =   0.396/GeV/GeV;
        Freggeon_Alphaprime =   1.324/GeV/GeV;
        Freggeon_C          =   1.0;
        FParity             =  +1;

        Wreggeon_Alpha      =   0.353;
        Wreggeon_Gamma      =   8.516/GeV/GeV;
        Wreggeon_Rsquare    =  24.40/GeV/GeV;
        Wreggeon_Alphaprime =   1.5/GeV/GeV;
        Wreggeon_C          =   1.0;
        WParity             =  -1;
 }
/*
G4cout<<G4endl<<"Reggeon's parameters for Particle "<<particle->GetParticleName()<<" "<<PDGcode<<G4endl<<G4endl;
G4cout<<"Alpha_pomeron "<<Alpha_pomeron;
G4cout<<" Alphaprime_pomeron "<<Alphaprime_pomeron*GeV*GeV;
G4cout<<" S0_pomeron "<<S0_pomeron/GeV/GeV<<G4endl;
G4cout<<"Gamma_pomeron "<<Gamma_pomeron*GeV*GeV;
G4cout<<" Rsquare_pomeron "<<Rsquare_pomeron*GeV*GeV;
G4cout<<" C_pomeron     "<<C_pomeron<<G4endl<<G4endl;
*/ 
}

G4double G4Reggeons::Get_Cprojectile() {return Cpr_pomeron;}

G4double G4Reggeons::Get_Ctarget()     {return Ctr_pomeron;}

G4Reggeons::~G4Reggeons() {}

void G4Reggeons::SetS(G4double S) {Sint = S;}

void G4Reggeons::CalculateXs()
{
	Xtotal  =0.; XtotalP=0.; XtotalR=0.;
	Xelastic=0.; Xpr_Diff=0.; Xtr_Diff=0.; XDDiff=0.; G4double XDiff=0.; 
        Xinel   =0.; Xnd=0.; XndP=0.; XndR=0.;

        G4double AmplitudeP(0.), AmplitudeR(0.);

        G4double B_max = 10.*fermi;
        G4double dB    = B_max/10000.;

        G4double B =-dB/2.;

        G4double chiP(0.), chiR(0.), chiRin(0.);   // chiPin Pomeron inelastic phase is a data member
        chiPin=0.;
	for(G4int i=0; i<10000;i++)
	{
	 B += dB;

         chiP   = Chi_pomeron(1.,B); chiR   = Chi_reggeon(1.,B);
         chiPin = Chi_pomeron(2.,B); chiRin = Chi_reggeon(2.,B);

	 AmplitudeP = (1.0/C_pomeron)*(1.0 - G4Exp(-chiP))*G4Exp(-chiR);
	 AmplitudeR =                 (1.0 - G4Exp(-chiR));

	 Xtotal   += 2 * (AmplitudeP + AmplitudeR) * B * dB;
	 XtotalP  += 2 * (AmplitudeP + 0.        ) * B * dB;
	 XtotalR  += 2 * (0.         + AmplitudeR) * B * dB;

         Xelastic +=                       sqr(AmplitudeP + AmplitudeR) * B * dB;
         Xpr_Diff +=  (Cpr_pomeron - 1.0) * sqr(AmplitudeP)   * B * dB;
         Xtr_Diff +=  (Ctr_pomeron - 1.0) * sqr(AmplitudeP)   * B * dB;
         XDiff    +=  (Cpr_pomeron - 1.0) * (Ctr_pomeron - 1.0) * sqr(AmplitudeP)   * B * dB;

// ----------------------------------
	 AmplitudeP = (1.0/C_pomeron)*(1.0 - G4Exp(-chiPin))*G4Exp(-chiRin);
	 AmplitudeR =                 (1.0 - G4Exp(-chiRin));

	 Xnd      += (AmplitudeP + AmplitudeR) * B * dB;
	 XndP     += (AmplitudeP + 0.        ) * B * dB;
	 XndR     += (0.         + AmplitudeR) * B * dB;
	}

	Xtotal *=twopi; XtotalP *=twopi; XtotalR *=twopi;
	Xelastic *=twopi; Xpr_Diff *=twopi; Xtr_Diff *=twopi; XDiff *=twopi;
	Xinel = Xtotal - Xelastic; 
        (void)Xinel;  // To avoid compiler warning "variable not used"

	Xnd  *=twopi; XndP *=twopi; XndR *=twopi;
	XDDiff = XDiff-Xpr_Diff-Xtr_Diff;

/*
G4cout<<"Total totalP totalR  "<<Xtotal/millibarn  <<" "<<XtotalP/millibarn <<" "<<XtotalR/millibarn<<" mb"<<G4endl;
G4cout<<"Elastic              "<<Xelastic/millibarn                                                        <<G4endl;
G4cout<<"PrDiff TrDiff W_Diff "<<Xpr_Diff/millibarn<<" "<<Xtr_Diff/millibarn<<" "<<(XDiff-Xpr_Diff-Xtr_Diff)/millibarn<<G4endl;
G4cout<<"Inelastic            "<<Xinel/millibarn                                                           <<G4endl;
G4cout<<"NonDiff Pom & Reg    "<<Xnd/millibarn     <<" "<<XndP/millibarn    <<" "<<XndR/millibarn          <<G4endl;
*/
}

G4double G4Reggeons::Chi_pomeron(G4double Mult, G4double B)
{
        G4double R2 = Rsquare_pomeron + Alphaprime_pomeron * G4Log(Sint/S0_pomeron); 
	G4double Eikonal = Mult * C_pomeron * Gamma_pomeron/R2 * 
			   G4Pow::GetInstance()->powA(Sint/S0_pomeron, Alpha_pomeron -1.) *
                           G4Exp(-sqr(B)/4.0/R2/hbarc_squared);
	return Eikonal;
}

G4double G4Reggeons::Chi_reggeon(G4double Mult, G4double B)
{
        G4double R2F = Freggeon_Rsquare + Freggeon_Alphaprime * G4Log(Sint/S0_pomeron); 
        G4double R2W = Wreggeon_Rsquare + Wreggeon_Alphaprime * G4Log(Sint/S0_pomeron); 

	G4double Eikonal = Mult * FParity * Freggeon_C * Freggeon_Gamma/R2F * 
			   G4Pow::GetInstance()->powA(Sint/S0_pomeron, Freggeon_Alpha -1.) *
                           G4Exp(-sqr(B)/4.0/R2F/hbarc_squared);

	         Eikonal+= Mult * WParity * Wreggeon_C * Wreggeon_Gamma/R2W * 
			   G4Pow::GetInstance()->powA(Sint/S0_pomeron, Wreggeon_Alpha -1.) *
                           G4Exp(-sqr(B)/4.0/R2W/hbarc_squared);
	return Eikonal;
}

G4double G4Reggeons::GetTotalX()  {	return Xtotal;             }
G4double G4Reggeons::GetTotalXp() {	return XtotalP;            }
G4double G4Reggeons::GetTotalXr() {	return XtotalR;            }

G4double G4Reggeons::GetElasticX(){	return Xelastic;           }
G4double G4Reggeons::GetPrDiffX() {	return Xpr_Diff;           }
G4double G4Reggeons::GetTrDiffX() {	return Xtr_Diff;           }
G4double G4Reggeons::GetDDiffX()  {	return XDDiff;             }

G4double G4Reggeons::GetInelX()   {	return Xinel;              }

G4double G4Reggeons::GetND_X()    {	return Xnd;                }
G4double G4Reggeons::GetNDp_X()   {	return XndP;               }
G4double G4Reggeons::GetNDr_X()   {	return XndR;               }

//----------------------------------------------------------------------------------------------
void G4Reggeons::GetProbabilities(G4double B, G4int Mode,
				  G4double & Pint, 
                                  G4double & Pprd, G4double & Ptrd, G4double & Pdd, 
                                  G4double & Pnd,  G4double & Pnvr)
{
 // Puprose of the method is a calculation of inelastic interaction probability     (Pint),
 //                                           probability of projectile diffraction (Pprd),
 //                                           probability of target diffraction     (Ptrd),
 //                                           probability of double diffraction     (Pdd ),
 //                                           probability of non-diffractive inter. (Pnd ),
 //                                           probability of quark-exc. inter.      (Pnvr),
 //                                           number of cutted pomerons             (NcutPomerons).
 // The input parameters are B - impact parameter, and Mode = All/WITHOUT_R/NON_DIFF
 //
	if( B > 2.* fermi ) { Pint=0.; Pprd=0.; Ptrd=0.; Pdd=0.; Pnd=0.; Pnvr=0.; return;}
        // At large B for hN collisions it is better to return zero inter. probability
	
        G4double chiP   = Chi_pomeron(1.,B); G4double chiR   = Chi_reggeon(1.,B);
                 chiPin = Chi_pomeron(2.,B); G4double chiRin = Chi_reggeon(2.,B);
	       //chiPin is data member of the class

	G4double Exp_ChiR   = G4Exp(-chiR);

	G4double AmplitudeP = (1.0/C_pomeron)*(1.0 - G4Exp(-chiP))*Exp_ChiR;
	G4double AmplitudeR =                 (1.0 - Exp_ChiR);

	G4double AmplitudeP2, Apr_Diff, Atr_Diff, ADiff;

      //Aelastic = sqr(AmplitudeP + AmplitudeR);  
	AmplitudeP2 = sqr(AmplitudeP);
        Apr_Diff = (Cpr_pomeron - 1.0) * AmplitudeP2;
        Atr_Diff = (Ctr_pomeron - 1.0) * AmplitudeP2;
        ADiff    = (Cpr_pomeron - 1.0) * (Ctr_pomeron - 1.0) * AmplitudeP2;

// ----------------------------------
	Exp_ChiR   = G4Exp(-chiRin);
	AmplitudeP = (1.0/C_pomeron)*(1.0 - G4Exp(-chiPin))*Exp_ChiR;
	AmplitudeR =                 (1.0 - Exp_ChiR);

	G4double And, AndP, AndR;

	And      = (AmplitudeP + AmplitudeR);
	AndP     = (AmplitudeP + 0.        );
	AndR     = (0.         + AmplitudeR);

// ----------------------------------
	if( Mode == ALL)
	{
	  Pint = Apr_Diff + Atr_Diff + ADiff + And;
	  Pprd = Apr_Diff/Pint;                     // Probability of projectile diffraction
	  Ptrd = Atr_Diff/Pint;                     // Probability of target diffraction 
	  Pdd  = ADiff   /Pint;                     // Probability of double diffraction
	  Pnd  = AndP    /Pint;                     // Probability of non-diffractive inelastic 
                                                    // interaction
	  Pnvr = AndR    /Pint;                     // Probability of non-vacuum reggeon (nvr) 
                                                    // inelastic interaction
	}
	else if( Mode == WITHOUT_R)
	{
	  Pint = Apr_Diff + Atr_Diff + ADiff + AndP;
	  Pprd = Apr_Diff/Pint; 
	  Ptrd = Atr_Diff/Pint; 
	  Pdd  = ADiff   /Pint;
	  Pnd  = AndP    /Pint;
	  Pnvr = 0.;
	}
	else
	{// Mode == NON_DIFF (of projectile)
	  Pint =            Atr_Diff         + AndP;
	  Pprd = 0.;
	  Ptrd = Atr_Diff/Pint; 
	  Pdd  = 0.;
	  Pnd  = AndP    /Pint;
	  Pnvr = 0.;
	}

	return;
}

G4int G4Reggeons::ncPomerons()         // Non-complite Poisson distribution
{
	if( chiPin < 0.001 ) return 0; // At small average multiplicity of cutted pomerons
				       // it is better to return 0 to avoid problems with
				       // calculation exactness.
	G4double ksi  = G4UniformRand() * (1.0-G4Exp(-chiPin)) * G4Exp(chiPin);
	G4double Term = chiPin;
	G4double Sum  = Term;
	G4int nCuts   = 1;

	while( Sum < ksi)
	{
	  nCuts++;
	  Term *= chiPin/(G4double) nCuts;
	  Sum += Term;
	}

	return nCuts;
}

