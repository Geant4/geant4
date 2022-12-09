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
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: 
//             Gunter Folger, August/September 2001
//               Create class; algorithm previously in G4VLongitudinalStringDecay.
// -----------------------------------------------------------------------------

#include "G4HadronBuilder.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4HadronicException.hh"
#include "G4ParticleTable.hh"

//#define debug_Hbuilder
//#define debug_heavyHadrons

G4HadronBuilder::G4HadronBuilder(const std::vector<G4double> & mesonMix, const G4double barionMix,
		                 const std::vector<G4double> & scalarMesonMix,
		                 const std::vector<G4double> & vectorMesonMix,
                                 const G4double Eta_cProb, const G4double Eta_bProb)
{
	mesonSpinMix       = mesonMix;
	barionSpinMix      = barionMix;
	scalarMesonMixings = scalarMesonMix;
	vectorMesonMixings = vectorMesonMix;
        ProbEta_c          = Eta_cProb;
        ProbEta_b          = Eta_bProb;
}

//-------------------------------------------------------------------------

G4ParticleDefinition * G4HadronBuilder::Build(G4ParticleDefinition * black, G4ParticleDefinition * white)
{
	if (black->GetParticleSubType()== "di_quark" || white->GetParticleSubType()== "di_quark" ) {
           // Barion
	   Spin spin = (G4UniformRand() < barionSpinMix) ? SpinHalf : SpinThreeHalf;
	   return Barion(black,white,spin);
	} else {
           // Meson
	   G4int StrangeQ = 0;
	   if( std::abs(black->GetPDGEncoding()) >= 3 ) StrangeQ++;
	   if( std::abs(white->GetPDGEncoding()) >= 3 ) StrangeQ++;
           Spin spin = (G4UniformRand() < mesonSpinMix[StrangeQ]) ? SpinZero : SpinOne;
	   return Meson(black,white,spin);
	}
}

//-------------------------------------------------------------------------

G4ParticleDefinition * G4HadronBuilder::BuildLowSpin(G4ParticleDefinition * black, G4ParticleDefinition * white)
{
	if ( black->GetParticleSubType()== "quark" && white->GetParticleSubType()== "quark" ) {
		return Meson(black,white, SpinZero);
	} else {
                // will return a SpinThreeHalf Barion if all quarks the same
		return Barion(black,white, SpinHalf); 
	}	
}

//-------------------------------------------------------------------------

G4ParticleDefinition * G4HadronBuilder::BuildHighSpin(G4ParticleDefinition * black, G4ParticleDefinition * white)
{
	if ( black->GetParticleSubType()== "quark" && white->GetParticleSubType()== "quark" ) {
		return Meson(black,white, SpinOne);
	} else {
		return Barion(black,white,SpinThreeHalf);
	}
}

//-------------------------------------------------------------------------

G4ParticleDefinition * G4HadronBuilder::Meson(G4ParticleDefinition * black, 
					      G4ParticleDefinition * white, Spin theSpin)
{
       #ifdef debug_Hbuilder
       // Verify Input Charge
       G4double charge =  black->GetPDGCharge() + white->GetPDGCharge();	 
       if (std::abs(charge) > 2 || std::abs(3.*charge - 3*G4int(charge*1.001)) > perCent )   // 1.001 to avoid int(.9999) -> 0
       {
	    G4cerr << " G4HadronBuilder::Build()" << G4endl;
	    G4cerr << "    Invalid total charge found for on input: " 
			<< charge<< G4endl;
	    G4cerr << "    PGDcode input quark1/quark2 : " <<
			black->GetPDGEncoding() << " / "<< 
			white->GetPDGEncoding() << G4endl;
	    G4cerr << G4endl;
	} 
        #endif	

        G4int id1 = black->GetPDGEncoding();
	G4int id2 = white->GetPDGEncoding();

        // G4int ifl1= std::max(std::abs(id1), std::abs(id2));
	if ( std::abs(id1) < std::abs(id2) )
	{
	   G4int xchg = id1; 
	   id1 = id2;  
	   id2 = xchg;
	}

        G4int abs_id1 = std::abs(id1); 
	
	if ( abs_id1 > 5 )
	   throw G4HadronicException(__FILE__, __LINE__, "G4HadronBuilder::Meson : Illegal Quark content as input");
	
        G4int PDGEncoding=0;

	if (id1 + id2 == 0) {     
           if ( abs_id1 < 4) {  // light quarks: u, d or s
	      G4double rmix = G4UniformRand();
	      G4int    imix = 2*std::abs(id1) - 1;
	      if (theSpin == SpinZero) {
	         PDGEncoding = 110*(1 + (G4int)(rmix + scalarMesonMixings[imix - 1])
        	                      + (G4int)(rmix + scalarMesonMixings[imix])
		        	    ) +  theSpin;
	      } else {
	         PDGEncoding = 110*(1 + (G4int)(rmix + vectorMesonMixings[imix - 1])
		        	      + (G4int)(rmix + vectorMesonMixings[imix])
				    ) +  theSpin;
	      }
           } else {  // for c and b quarks

              PDGEncoding = abs_id1*100 + abs_id1*10;

              if (PDGEncoding == 440) {
                if ( G4UniformRand() < ProbEta_c ) {
                  PDGEncoding +=1;
                } else {
                  PDGEncoding +=3;
                }
              }
              if (PDGEncoding == 550) {
                if ( G4UniformRand() < ProbEta_b ) {
                  PDGEncoding +=1;
                } else {
                  PDGEncoding +=3;
                }
              }
           }
	} else {
	   PDGEncoding = 100 * std::abs(id1) + 10 * std::abs(id2) +  theSpin;  
	   G4bool IsUp = (std::abs(id1)&1) == 0;	// quark 1 up type quark (u or c)
	   G4bool IsAnti = id1 < 0; 		        // quark 1 is antiquark?
	   if ( (IsUp && IsAnti ) || (!IsUp && !IsAnti ) ) PDGEncoding = - PDGEncoding;
 	}

        // ---------------------------------------------------------------------
        // Special treatment for charmed and bottom mesons : in Geant4 there are
	// no excited charmed or bottom mesons, therefore we need to transform these
	// into existing charmed and bottom mesons in Geant4. Whenever possible,
	// we use the corresponding ground state mesons with the same quantum numbers;
	// else, we prefer to conserve the electric charge rather than other flavor numbers.
        #ifdef debug_heavyHadrons
	G4int initialPDGEncoding = PDGEncoding;
	#endif
        if      ( std::abs( PDGEncoding ) == 10411 )  // D*0(2400)+   ->  D+
	  ( PDGEncoding > 0 ? PDGEncoding = 411 : PDGEncoding = -411 );
        else if ( std::abs( PDGEncoding ) == 10421 )  // D*0(2400)0   ->  D0
	  ( PDGEncoding > 0 ? PDGEncoding = 421 : PDGEncoding = -421 );
        else if ( std::abs( PDGEncoding ) == 413 )    // D*(2010)+    ->  D+
	  ( PDGEncoding > 0 ? PDGEncoding = 411 : PDGEncoding = -411 );
        else if ( std::abs( PDGEncoding ) == 423 )    // D*(2007)0    ->  D0
	  ( PDGEncoding > 0 ? PDGEncoding = 421 : PDGEncoding = -421 );
        else if ( std::abs( PDGEncoding ) == 10413 )  // D1(2420)+    ->  D+
	  ( PDGEncoding > 0 ? PDGEncoding = 411 : PDGEncoding = -411 );
        else if ( std::abs( PDGEncoding ) == 10423 )  // D1(2420)0    ->  D0
	  ( PDGEncoding > 0 ? PDGEncoding = 421 : PDGEncoding = -421 );
        else if ( std::abs( PDGEncoding ) == 20413 )  // D1(H)+       ->  D+
	  ( PDGEncoding > 0 ? PDGEncoding = 411 : PDGEncoding = -411 );
        else if ( std::abs( PDGEncoding ) == 20423 )  // D1(2430)0    ->  D0
	  ( PDGEncoding > 0 ? PDGEncoding = 421 : PDGEncoding = -421 );
        else if ( std::abs( PDGEncoding ) == 415 )    // D2*(2460)+   ->  D+
	  ( PDGEncoding > 0 ? PDGEncoding = 411 : PDGEncoding = -411 );
        else if ( std::abs( PDGEncoding ) == 425 )    // D2*(2460)0   ->  D0
	  ( PDGEncoding > 0 ? PDGEncoding = 421 : PDGEncoding = -421 );
        else if ( std::abs( PDGEncoding ) == 10431 )  // Ds0*(2317)+  ->  Ds+
	  ( PDGEncoding > 0 ? PDGEncoding = 431 : PDGEncoding = -431 );
        else if ( std::abs( PDGEncoding ) == 433 )    // Ds*+         ->  Ds+
	  ( PDGEncoding > 0 ? PDGEncoding = 431 : PDGEncoding = -431 );
        else if ( std::abs( PDGEncoding ) == 10433 )  // Ds1(2536)+   ->  Ds+
	  ( PDGEncoding > 0 ? PDGEncoding = 431 : PDGEncoding = -431 );
        else if ( std::abs( PDGEncoding ) == 20433 )  // Ds1(2460)+   ->  Ds+
	  ( PDGEncoding > 0 ? PDGEncoding = 431 : PDGEncoding = -431 );
        else if ( std::abs( PDGEncoding ) == 435 )    // Ds2*(2573)+  ->  Ds+
	  ( PDGEncoding > 0 ? PDGEncoding = 431 : PDGEncoding = -431 );
        else if ( std::abs( PDGEncoding ) ==   10441 ) PDGEncoding = 441;  // chi_c0(1P) ->  eta_c
        else if ( std::abs( PDGEncoding ) ==  100441 ) PDGEncoding = 441;  // eta_c(2S)  ->  eta_c
        else if ( std::abs( PDGEncoding ) ==   10443 ) PDGEncoding = 443;  // h_c(1P)    ->  J/psi
        else if ( std::abs( PDGEncoding ) ==   20443 ) PDGEncoding = 443;  // chi_c1(1P) ->  J/psi
        else if ( std::abs( PDGEncoding ) ==  100443 ) PDGEncoding = 443;  // psi(2S)    ->  J/psi
        else if ( std::abs( PDGEncoding ) ==   30443 ) PDGEncoding = 443;  // psi(3770)  ->  J/psi
        else if ( std::abs( PDGEncoding ) == 9000443 ) PDGEncoding = 443;  // psi(4040)  ->  J/psi
        else if ( std::abs( PDGEncoding ) == 9010443 ) PDGEncoding = 443;  // psi(4160)  ->  J/psi
        else if ( std::abs( PDGEncoding ) == 9020443 ) PDGEncoding = 443;  // psi(4415)  ->  J/psi
        else if ( std::abs( PDGEncoding ) ==     445 ) PDGEncoding = 443;  // chi_c2(1P) ->  J/psi
        else if ( std::abs( PDGEncoding ) ==  100445 ) PDGEncoding = 443;  // chi_c2(2P) ->  J/psi
        // Bottom mesons
        else if ( std::abs( PDGEncoding ) == 10511 )  // B0*0         ->  B0
	  ( PDGEncoding > 0 ? PDGEncoding = 511 : PDGEncoding = -511 );
        else if ( std::abs( PDGEncoding ) == 10521 )  // B0*+         ->  B+
	  ( PDGEncoding > 0 ? PDGEncoding = 521 : PDGEncoding = -521 );
        else if ( std::abs( PDGEncoding ) == 513 )    // B*0          ->  B0
	  ( PDGEncoding > 0 ? PDGEncoding = 511 : PDGEncoding = -511 );
        else if ( std::abs( PDGEncoding ) == 523 )    // B*+          ->  B+
	  ( PDGEncoding > 0 ? PDGEncoding = 521 : PDGEncoding = -521 );
        else if ( std::abs( PDGEncoding ) == 10513 )  // B1(L)0       ->  B0
	  ( PDGEncoding > 0 ? PDGEncoding = 511 : PDGEncoding = -511 );
        else if ( std::abs( PDGEncoding ) == 10523 )  // B1(L)+       ->  B+
	  ( PDGEncoding > 0 ? PDGEncoding = 521 : PDGEncoding = -521 );
        else if ( std::abs( PDGEncoding ) == 20513 )  // B1(H)0       ->  B0
	  ( PDGEncoding > 0 ? PDGEncoding = 511 : PDGEncoding = -511 );
        else if ( std::abs( PDGEncoding ) == 20523 )  // B1(H)+       ->  B+
	  ( PDGEncoding > 0 ? PDGEncoding = 521 : PDGEncoding = -521 );
        else if ( std::abs( PDGEncoding ) == 515 )    // B2*0         ->  B0
	  ( PDGEncoding > 0 ? PDGEncoding = 511 : PDGEncoding = -511 );
        else if ( std::abs( PDGEncoding ) == 525 )    // B2*+         ->  B+
	  ( PDGEncoding > 0 ? PDGEncoding = 521 : PDGEncoding = -521 );
        else if ( std::abs( PDGEncoding ) == 10531 )  // Bs0*0        ->  Bs0
	  ( PDGEncoding > 0 ? PDGEncoding = 531 : PDGEncoding = -531 );
        else if ( std::abs( PDGEncoding ) == 533 )    // Bs*0         ->  Bs0
	  ( PDGEncoding > 0 ? PDGEncoding = 531 : PDGEncoding = -531 );
        else if ( std::abs( PDGEncoding ) == 10533 )  // Bs1(L)0      ->  Bs0
	  ( PDGEncoding > 0 ? PDGEncoding = 531 : PDGEncoding = -531 );
        else if ( std::abs( PDGEncoding ) == 20533 )  // Bs1(H)0      ->  Bs0
	  ( PDGEncoding > 0 ? PDGEncoding = 531 : PDGEncoding = -531 );
        else if ( std::abs( PDGEncoding ) == 535 )    // Bs2*0        ->  Bs0
	  ( PDGEncoding > 0 ? PDGEncoding = 531 : PDGEncoding = -531 );
        else if ( std::abs( PDGEncoding ) == 10541 )  // Bc0*+        ->  Bc+
	  ( PDGEncoding > 0 ? PDGEncoding = 541 : PDGEncoding = -541 );
        else if ( std::abs( PDGEncoding ) == 543 )    // Bc*+         ->  Bc+
	  ( PDGEncoding > 0 ? PDGEncoding = 541 : PDGEncoding = -541 );
        else if ( std::abs( PDGEncoding ) == 10543 )  // Bc1(L)+      ->  Bc+
	  ( PDGEncoding > 0 ? PDGEncoding = 541 : PDGEncoding = -541 );
	else if ( std::abs( PDGEncoding ) == 20543 )  // Bc1(H)+      ->  Bc+
	  ( PDGEncoding > 0 ? PDGEncoding = 541 : PDGEncoding = -541 );
        else if ( std::abs( PDGEncoding ) == 545 )    // Bc2*+        ->  Bc+
	  ( PDGEncoding > 0 ? PDGEncoding = 541 : PDGEncoding = -541 );
        else if ( std::abs( PDGEncoding ) ==     551 )  PDGEncoding = 553;  // eta_b(1S)       ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==   10551 )  PDGEncoding = 553;  // chi_b0(1P)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  100551 )  PDGEncoding = 553;  // eta_b(2S)       ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  110551 )  PDGEncoding = 553;  // chi_b0(2P)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  200551 )  PDGEncoding = 553;  // eta_b(3S)       ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  210551 )  PDGEncoding = 553;  // chi_b0(3P)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==   10553 )  PDGEncoding = 553;  // h_b(1P)         ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==   20553 )  PDGEncoding = 553;  // chi_b1(1P)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==   30553 )  PDGEncoding = 553;  // Upsilon_1(1D)   ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  100553 )  PDGEncoding = 553;  // Upsilon(2S)     ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  110553 )  PDGEncoding = 553;  // h_b(2P)         ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  120553 )  PDGEncoding = 553;  // chi_b1(2P)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  130553 )  PDGEncoding = 553;  // Upsilon_1(2D)   ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  200553 )  PDGEncoding = 553;  // Upsilon(3S)     ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  210553 )  PDGEncoding = 553;  // h_b(3P)         ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  220553 )  PDGEncoding = 553;  // chi_b1(3P)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  300553 )  PDGEncoding = 553;  // Upsilon(4S)     ->  Upsilon
        else if ( std::abs( PDGEncoding ) == 9000553 )  PDGEncoding = 553;  // Upsilon(10860)  ->  Upsilon
        else if ( std::abs( PDGEncoding ) == 9010553 )  PDGEncoding = 553;  // Upsilon(11020)  ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==     555 )  PDGEncoding = 553;  // chi_b2(1P)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==   10555 )  PDGEncoding = 553;  // eta_b2(1D)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==   20555 )  PDGEncoding = 553;  // Upsilon_2(1D)   ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  100555 )  PDGEncoding = 553;  // chi_b2(2P)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  110555 )  PDGEncoding = 553;  // eta_b2(2D)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  120555 )  PDGEncoding = 553;  // Upsilon_2(2D)   ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  200555 )  PDGEncoding = 553;  // chi_b2(3P)      ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==     557 )  PDGEncoding = 553;  // Upsilon_3(1D)   ->  Upsilon
        else if ( std::abs( PDGEncoding ) ==  100557 )  PDGEncoding = 553;  // Upsilon_3(2D)   ->  Upsilon
        #ifdef debug_heavyHadrons
	if ( initialPDGEncoding != PDGEncoding ) {
	  G4cout << "G4HadronBuilder::Meson : forcing (inexisting in G4) heavy meson with pdgCode="
		 << initialPDGEncoding << " into pdgCode=" << PDGEncoding << G4endl;
	}
	#endif
        // ---------------------------------------------------------------------
	
	G4ParticleDefinition * MesonDef=
		G4ParticleTable::GetParticleTable()->FindParticle(PDGEncoding);

        #ifdef debug_Hbuilder
	if (MesonDef == 0 ) {
		G4cerr << " G4HadronBuilder - Warning: No particle for PDGcode= "
		       << PDGEncoding << G4endl;
	} else if  ( (  black->GetPDGCharge() + white->GetPDGCharge()
	       		- MesonDef->GetPDGCharge() ) > perCent   ) {
	      	G4cerr << " G4HadronBuilder - Warning: Incorrect Charge : "
			<< " Quark1/2 = " 
			<< black->GetParticleName() << " / "
			<< white->GetParticleName() 
			<< " resulting Hadron " << MesonDef->GetParticleName() 
			<< G4endl;
	}
        #endif

	return MesonDef;
}

//-------------------------------------------------------------------------

G4ParticleDefinition * G4HadronBuilder::Barion(G4ParticleDefinition * black, 
					       G4ParticleDefinition * white,Spin theSpin)
{
        #ifdef debug_Hbuilder
        // Verify Input Charge
        G4double charge =  black->GetPDGCharge() + white->GetPDGCharge();	 
        if (std::abs(charge) > 2 || std::abs(3.*charge - 3*G4int(charge*1.001)) > perCent )
   	{
	    G4cerr << " G4HadronBuilder::Build()" << G4endl;
	    G4cerr << "    Invalid total charge found for on input: " 
			<< charge<< G4endl;
	    G4cerr << "    PGDcode input quark1/quark2 : " <<
			black->GetPDGEncoding() << " / "<< 
			white->GetPDGEncoding() << G4endl;
	    G4cerr << G4endl;
	} 
        #endif	

        G4int id1 = black->GetPDGEncoding();
	G4int id2 = white->GetPDGEncoding();

	if ( std::abs(id1) < std::abs(id2) )
	{
	   G4int xchg = id1; 
	   id1 = id2;  
	   id2 = xchg;
	}

	if (std::abs(id1) < 1000 || std::abs(id2) > 5 )
	   throw G4HadronicException(__FILE__, __LINE__, "G4HadronBuilder::Barion: Illegal quark content as input");   

	G4int ifl1= std::abs(id1)/1000;
	G4int ifl2 = (std::abs(id1) - ifl1 * 1000)/100;
	G4int diquarkSpin = std::abs(id1)%10; 
	G4int ifl3 = id2;
	if (id1 < 0)
	{
	   ifl1 = - ifl1;
	   ifl2 = - ifl2;
	}
	//... Construct barion, distinguish Lambda and Sigma barions.
	G4int kfla = std::abs(ifl1);
	G4int kflb = std::abs(ifl2);
	G4int kflc = std::abs(ifl3);

	G4int kfld = std::max(kfla,kflb);
	      kfld = std::max(kfld,kflc);
	G4int kflf = std::min(kfla,kflb);
	      kflf = std::min(kflf,kflc);

	G4int kfle = kfla + kflb + kflc - kfld - kflf;

	//... barion with content uuu or ddd or sss has always spin = 3/2
	theSpin = (kfla == kflb && kflb == kflc)? SpinThreeHalf : theSpin;   

	G4int kfll = 0;
        if (kfld < 6) {
	   if (theSpin == SpinHalf && kfld > kfle && kfle > kflf) { 
              // Spin J=1/2 and all three quarks different
              // Two states exist: (uds -> lambda or sigma0)
              //   -  lambda: s(ud)0 s : 3122; ie. reverse the two lighter quarks
              //   -  sigma0: s(ud)1 s : 3212
	      if (diquarkSpin == 1 ) {
	         if ( kfla == kfld) {   // heaviest quark in diquark
	            kfll = 1;
	         } else {
	            kfll = (G4int)(0.25 + G4UniformRand());
	         }
	      }   
	      if (diquarkSpin == 3 && kfla != kfld)
	          kfll = (G4int)(0.75 + G4UniformRand());
	   }
        }
	
	G4int PDGEncoding;
	if (kfll == 1)
	   PDGEncoding = 1000 * kfld + 100 * kflf + 10 * kfle + theSpin;
	else    
	   PDGEncoding = 1000 * kfld + 100 * kfle + 10 * kflf + theSpin;

	if (id1 < 0)
	   PDGEncoding = -PDGEncoding;

        // ---------------------------------------------------------------------
        // Special treatment for charmed and bottom baryons : in Geant4 there are
	// neither excited charmed or bottom baryons, nor baryons with two or three
	// heavy (c, b) constitutent quarks:
	//   Sigma_c* , Xi_c' , Xi_c* , Omega_c* ,
	//   Xi_cc , Xi_cc* , Omega_cc , Omega_cc* , Omega_ccc ;
	//   Sigma_b* , Xi_b' , Xi_b* , Omega_b*,
	//   Xi_bc , Xi_bc' , Xi_bc* , Omega_bc , Omega_bc' , Omega_bc* ,
	//   Omega_bcc , Omega_bcc* , Xi_bb, Xi_bb* , Omega_bb, Omega_bb* ,
	//   Omega_bbc , Omega_bbc* , Omega_bbb
	// therefore we need to transform these into existing charmed and bottom
	// baryons in Geant4. Whenever possible, we use the corresponding ground state
	// baryons with the same quantum numbers; else, we prefer to conserve the
	// electric charge rather than other flavor numbers.
        #ifdef debug_heavyHadrons
	G4int charmViolation = 0, bottomViolation = 0;  // Only positive
	G4int initialPDGEncoding = PDGEncoding;
	#endif
        if      ( std::abs( PDGEncoding ) == 4224 ) {  // Sigma_c*++   ->  Sigma_c++
	  ( PDGEncoding > 0 ? PDGEncoding = 4222 : PDGEncoding = -4222 );
	} else if ( std::abs( PDGEncoding ) == 4214 ) {  // Sigma_c*+    ->  Sigma_c+
	  ( PDGEncoding > 0 ? PDGEncoding = 4212 : PDGEncoding = -4212 );
	} else if ( std::abs( PDGEncoding ) == 4114 ) {  // Sigma_c*0    ->  Sigma_c0
	  ( PDGEncoding > 0 ? PDGEncoding = 4112 : PDGEncoding = -4112 );
	} else if ( std::abs( PDGEncoding ) == 4322 ) {  // Xi_c'+       ->  Xi_c+
	  ( PDGEncoding > 0 ? PDGEncoding = 4232 : PDGEncoding = -4232 );
	} else if ( std::abs( PDGEncoding ) == 4312 ) {  // Xi_c'0       ->  Xi_c0
	  ( PDGEncoding > 0 ? PDGEncoding = 4132 : PDGEncoding = -4132 );
	} else if ( std::abs( PDGEncoding ) == 4324 ) {  // Xi_c*+       ->  Xi_c+
	  ( PDGEncoding > 0 ? PDGEncoding = 4232 : PDGEncoding = -4232 );
	} else if ( std::abs( PDGEncoding ) == 4314 ) {  // Xi_c*0       ->  Xi_c0
	  ( PDGEncoding > 0 ? PDGEncoding = 4132 : PDGEncoding = -4132 );
	} else if ( std::abs( PDGEncoding ) == 4334 ) {  // Omega_c*0    ->  Omega_c0
	  ( PDGEncoding > 0 ? PDGEncoding = 4332 : PDGEncoding = -4332 );
	} else if ( std::abs( PDGEncoding ) == 4412 ) {  // Xi_cc+       ->  Xi_c+
	  ( PDGEncoding > 0 ? PDGEncoding = 4232 : PDGEncoding = -4232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 4422 ) {  // Xi_cc++      ->  Sigma_c++ (use Sigma to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 4222 : PDGEncoding = -4222 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
	} else if ( std::abs( PDGEncoding ) == 4414 ) {  // Xi_cc*+      ->  Xi_c+
	  ( PDGEncoding > 0 ? PDGEncoding = 4232 : PDGEncoding = -4232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 4424 ) {  // Xi_cc*++     ->  Sigma_c++ (use Sigma to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 4222 : PDGEncoding = -4222 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
	} else if ( std::abs( PDGEncoding ) == 4432 ) {  // Omega_cc+    ->  Xi_c+ (use Xi to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 4232 : PDGEncoding = -4232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
	} else if ( std::abs( PDGEncoding ) == 4434 ) {  // Omega_cc*+   ->  Xi_c+ (use Xi to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 4232 : PDGEncoding = -4232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 4444 ) {  // Omega_ccc++  ->  Sigma_c++ (use Sigma to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 4222 : PDGEncoding = -4222 );
          #ifdef debug_heavyHadrons
	  charmViolation = 2;
	  #endif
        // Bottom baryons
        } else if ( std::abs( PDGEncoding ) == 5114 ) {  // Sigma_b*-    ->  Sigma_b-
	  ( PDGEncoding > 0 ? PDGEncoding = 5112 : PDGEncoding = -5112 );
        } else if ( std::abs( PDGEncoding ) == 5214 ) {  // Sigma_b*0    ->  Sigma_b0
	  ( PDGEncoding > 0 ? PDGEncoding = 5212 : PDGEncoding = -5212 );
        } else if ( std::abs( PDGEncoding ) == 5224 ) {  // Sigma_b*+    ->  Sigma_b+
	  ( PDGEncoding > 0 ? PDGEncoding = 5222 : PDGEncoding = -5222 );
        } else if ( std::abs( PDGEncoding ) == 5312 ) {  // Xi_b'-       ->  Xi_b-
	  ( PDGEncoding > 0 ? PDGEncoding = 5132 : PDGEncoding = -5132 );
        } else if ( std::abs( PDGEncoding ) == 5322 ) {  // Xi_b'0       ->  Xi_b0
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
        } else if ( std::abs( PDGEncoding ) == 5314 ) {  // Xi_b*-       ->  Xi_b-
	  ( PDGEncoding > 0 ? PDGEncoding = 5132 : PDGEncoding = -5132 );
        } else if ( std::abs( PDGEncoding ) == 5324 ) {  // Xi_b*0       ->  Xi_b0
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
        } else if ( std::abs( PDGEncoding ) == 5334 ) {  // Omega_b*-    ->  Omega_b-
	  ( PDGEncoding > 0 ? PDGEncoding = 5332 : PDGEncoding = -5332 );
        } else if ( std::abs( PDGEncoding ) == 5142 ) {  // Xi_bc0       ->  Xi_b0
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5242 ) {  // Xi_bc+       ->  Sigma_b+ (use Sigma to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 5222 : PDGEncoding = -5222 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5412 ) {  // Xi_bc'0      ->  Xi_b0
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5422 ) {  // Xi_bc'+      ->  Sigma_b+ (use Sigma to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 5222 : PDGEncoding = -5222 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5414 ) { // Xi_bc*0      ->  Xi_b0
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5424 ) {  // Xi_bc*+      ->  Sigma_b+ (use Sigma to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 5222 : PDGEncoding = -5222 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5342 ) {  // Omega_bc0    ->  Xi_b0 (use Xi to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5432 ) {  // Omega_bc'0   ->  Xi_b0 (use Xi to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5434 ) {  // Omega_bc*0   ->  Xi_b0 (use Xi to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5442 ) {  // Omega_bcc+   ->  Sigma_b+ (use Sigma to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 5222 : PDGEncoding = -5222 );
          #ifdef debug_heavyHadrons
	  charmViolation = 2;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5444 ) {  // Omega_bcc*+  ->  Sigma_b+ (use Sigma to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 5222 : PDGEncoding = -5222 );
          #ifdef debug_heavyHadrons
	  charmViolation = 2;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5512 ) {  // Xi_bb-       ->  Xi_b-
	  ( PDGEncoding > 0 ? PDGEncoding = 5132 : PDGEncoding = -5132 );
          #ifdef debug_heavyHadrons
	  bottomViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5522 ) {  // Xi_bb0       ->  Xi_b0
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
          #ifdef debug_heavyHadrons
	  bottomViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5514 ) {  // Xi_bb*-      ->  Xi_b-
	  ( PDGEncoding > 0 ? PDGEncoding = 5132 : PDGEncoding = -5132 );
          #ifdef debug_heavyHadrons
	  bottomViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5524 ) {  // Xi_bb*0      ->  Xi_b0
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
          #ifdef debug_heavyHadrons
	  bottomViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5532 ) {  // Omega_bb-    ->  Omega_b-
	  ( PDGEncoding > 0 ? PDGEncoding = 5332 : PDGEncoding = -5332 );
          #ifdef debug_heavyHadrons
	  bottomViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5534 ) {  // Omega_bb*-   ->  Omega_b-
	  ( PDGEncoding > 0 ? PDGEncoding = 5332 : PDGEncoding = -5332 );
          #ifdef debug_heavyHadrons
	  bottomViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5542 ) {  // Omega_bbc0   ->  Xi_b0 (use Xi to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1; bottomViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5544 ) {  // Omega_bbc*0  ->  Xi_b0 (use Xi to conserve charge)
	  ( PDGEncoding > 0 ? PDGEncoding = 5232 : PDGEncoding = -5232 );
          #ifdef debug_heavyHadrons
	  charmViolation = 1; bottomViolation = 1;
	  #endif
        } else if ( std::abs( PDGEncoding ) == 5554 ) {   // Omega_bbb-   ->  Omega_b-
	  ( PDGEncoding > 0 ? PDGEncoding = 5332 : PDGEncoding = -5332 );
          #ifdef debug_heavyHadrons
	  bottomViolation = 2;
	  #endif
	}
        #ifdef debug_heavyHadrons
	if ( initialPDGEncoding != PDGEncoding ) {
	  G4cout << "G4HadronBuilder::Barion : forcing (inexisting in G4) heavy baryon with pdgCode="
		 << initialPDGEncoding << " into pdgCode=" << PDGEncoding << G4endl;
	  if ( charmViolation != 0  ||  bottomViolation != 0 ) {
	    G4cout << "\t --> VIOLATION of " << ( charmViolation != 0 ? " CHARM " : " " )
	           << ( charmViolation != 0  &&  bottomViolation != 0 ? " and " : " " )
	           << ( bottomViolation != 0 ? " BOTTOM " : " " ) << " quantum number ! " << G4endl;
	  }
	}
	#endif
        // ---------------------------------------------------------------------
	
	G4ParticleDefinition * BarionDef=
		G4ParticleTable::GetParticleTable()->FindParticle(PDGEncoding);

        #ifdef debug_Hbuilder
	if (BarionDef == 0 ) {
		G4cerr << " G4HadronBuilder - Warning: No particle for PDGcode= "
		       << PDGEncoding << G4endl;
	} else if  ( (  black->GetPDGCharge() + white->GetPDGCharge()
	       		- BarionDef->GetPDGCharge() ) > perCent   ) {
	      	G4cerr << " G4HadronBuilder - Warning: Incorrect Charge : "
			<< " DiQuark/Quark = " 
			<< black->GetParticleName() << " / "
			<< white->GetParticleName() 
			<< " resulting Hadron " << BarionDef->GetParticleName() 
			<< G4endl;
	}
        #endif

	return BarionDef;
}
