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
// $Id: G4HadronBuilder.cc 102717 2017-02-20 10:37:13Z gcosmo $
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

G4HadronBuilder::G4HadronBuilder(G4double mesonMix, G4double barionMix,
		     std::vector<double> scalarMesonMix,
		     std::vector<double> vectorMesonMix)
{
	mesonSpinMix=mesonMix;	     
	barionSpinMix=barionMix;
	scalarMesonMixings=scalarMesonMix;
	vectorMesonMixings=vectorMesonMix;
}

G4ParticleDefinition * G4HadronBuilder::Build(G4ParticleDefinition * black, G4ParticleDefinition * white)
{

	if (black->GetParticleSubType()== "di_quark" || white->GetParticleSubType()== "di_quark" ) {

//    Barion
	   Spin spin = (G4UniformRand() < barionSpinMix) ? SpinHalf : SpinThreeHalf;
	   return Barion(black,white,spin);

	} else {	

//    Meson
	   Spin spin = (G4UniformRand() < mesonSpinMix) ? SpinZero : SpinOne;
	   return Meson(black,white,spin);
	}
}

//-------------------------------------------------------------------------

G4ParticleDefinition * G4HadronBuilder::BuildLowSpin(G4ParticleDefinition * black, G4ParticleDefinition * white)
{
	if ( black->GetParticleSubType()== "quark" && white->GetParticleSubType()== "quark" ) {
		return Meson(black,white, SpinZero);
	} else {
//		      will return a SpinThreeHalf Barion if all quarks the same
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
#ifdef G4VERBOSE
//  Verify Input Charge
   
   G4double charge =  black->GetPDGCharge() 
                    + white->GetPDGCharge();	 
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
	
	G4int id1= black->GetPDGEncoding();
	G4int id2= white->GetPDGEncoding();
//	G4int ifl1= std::max(std::abs(id1), std::abs(id2));
	if ( std::abs(id1) < std::abs(id2) )
	   {
	   G4int xchg = id1; 
	   id1 = id2;  
	   id2 = xchg;
	   }
	
	if (std::abs(id1) > 3 ) 
	   throw G4HadronicException(__FILE__, __LINE__, "G4HadronBuilder::Meson : Illegal Quark content as input");
	
        G4int PDGEncoding=0;

	if (id1 + id2 == 0) {     
	   G4double rmix = G4UniformRand();
	   G4int    imix = 2*std::abs(id1) - 1;
	   if(theSpin == SpinZero) {
	      PDGEncoding = 110*(1 + (G4int)(rmix + scalarMesonMixings[imix - 1])
        	                   + (G4int)(rmix + scalarMesonMixings[imix])
		 		) +  theSpin;
	   } else {
	      PDGEncoding = 110*(1 + (G4int)(rmix + vectorMesonMixings[imix - 1])
				   + (G4int)(rmix + vectorMesonMixings[imix])
				) +  theSpin;
	   }
	} else {
	   PDGEncoding = 100 * std::abs(id1) + 10 * std::abs(id2) +  theSpin;  
	   G4bool IsUp = (std::abs(id1)&1) == 0;	// quark 1 up type quark (u or c)
	   G4bool IsAnti = id1 < 0; 		// quark 1 is antiquark?
	   if( (IsUp && IsAnti ) || (!IsUp && !IsAnti ) ) 
	      PDGEncoding = - PDGEncoding;
 	}
	   
	   
	G4ParticleDefinition * MesonDef=
		G4ParticleTable::GetParticleTable()->FindParticle(PDGEncoding);
#ifdef G4VERBOSE
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


G4ParticleDefinition * G4HadronBuilder::Barion(G4ParticleDefinition * black, 
					      G4ParticleDefinition * white,Spin theSpin)
{

#ifdef G4VERBOSE
//  Verify Input Charge
   G4double charge =  black->GetPDGCharge() 
                    + white->GetPDGCharge();	 
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
	G4int id1= black->GetPDGEncoding();
	G4int id2= white->GetPDGEncoding();
	if ( std::abs(id1) < std::abs(id2) )
	   {
	   G4int xchg = id1; 
	   id1 = id2;  
	   id2 = xchg;
	   }

	if (std::abs(id1) < 1000 || std::abs(id2) > 3 ) 
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
	if(theSpin == SpinHalf && kfld > kfle && kfle > kflf) { 
// Spin J=1/2 and all three quarks different
// Two states exist: (uds -> lambda or sigma0)
//   -  lambda: s(ud)0 s : 3122; ie. reverse the two lighter quarks
//   -  sigma0: s(ud)1 s : 3212
	   if(diquarkSpin == 1 ) {
	      if ( kfla == kfld) {   // heaviest quark in diquark
	         kfll = 1;
	      } else {
	         kfll = (G4int)(0.25 + G4UniformRand());
	      }
	   }   
	   if(diquarkSpin == 3 && kfla != kfld)
	       kfll = (G4int)(0.75 + G4UniformRand());
	}
	
	G4int PDGEncoding;
	if (kfll == 1)
	   PDGEncoding = 1000 * kfld + 100 * kflf + 10 * kfle + theSpin;
	else    
	   PDGEncoding = 1000 * kfld + 100 * kfle + 10 * kflf + theSpin;

	if (id1 < 0)
	   PDGEncoding = -PDGEncoding;


	G4ParticleDefinition * BarionDef=
		G4ParticleTable::GetParticleTable()->FindParticle(PDGEncoding);
#ifdef G4VERBOSE
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
