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
// $Id: G4QHadronBuilder.cc,v 1.2 2006/12/12 11:02:22 mkossov Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -----------------------------------------------------------------------------
//      GEANT4 class implementation file
//
//      Created by Mikhail Kosov, October 2006 
//      Simple algorithm making a hadron out of two partons
//      For comparison mirror member functions are taken from G4 class:
//      G4HadronBuilder
// -----------------------------------------------------------------------------

//#define debug

#include "G4QHadronBuilder.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"

G4QHadronBuilder::G4QHadronBuilder()
{
	 mesonSpinMix = 0.5;	                    // probability to create vector meson 
	 baryonSpinMix= 0.5;                     // probability to create 3/2 baryon 
  scalarMesonMixings.resize(6);
  scalarMesonMixings[0] = 0.5; 
  scalarMesonMixings[1] = 0.25; 
  scalarMesonMixings[2] = 0.5; 
  scalarMesonMixings[3] = 0.25; 
  scalarMesonMixings[4] = 1.0; 
  scalarMesonMixings[5] = 0.5; 
  vectorMesonMixings.resize(6);
  vectorMesonMixings[0] = 0.5;
  vectorMesonMixings[1] = 0.0;
  vectorMesonMixings[2] = 0.5;
  vectorMesonMixings[3] = 0.0;
  vectorMesonMixings[4] = 1.0;
  vectorMesonMixings[5] = 1.0; 
}

G4QHadron* G4QHadronBuilder::Build(G4QParton* black, G4QParton* white)
{
	 if(black->GetParticleSubType()=="di_quark" || white->GetParticleSubType()=="di_quark")
  {
    // Baryon consists of quark and at least one di-quark
	   Spin spin = (G4UniformRand() < baryonSpinMix) ? SpinHalf : SpinThreeHalf;
	   return Baryon(black,white,spin);
	 }
  else
  {	
    // Meson consists of quark and abnti-quark
	   Spin spin = (G4UniformRand() < mesonSpinMix) ? SpinZero : SpinOne;
	   return Meson(black,white,spin);
	 }
}

G4QHadron* G4QHadronBuilder::BuildLowSpin(G4QParton* black, G4QParton* white)
{
	 if(black->GetParticleSubType() == "quark" && white->GetParticleSubType() == "quark")
		      return Meson(black,white, SpinZero);
  //	returns a SpinThreeHalf Baryon if all quarks are the same
	 else		return Baryon(black,white, SpinHalf);
}

G4QHadron* G4QHadronBuilder::BuildHighSpin(G4QParton* black, G4QParton* white)
{
 	if(black->GetParticleSubType() == "quark" && white->GetParticleSubType() == "quark")
		     return Meson(black,white, SpinOne);
	 else return Baryon(black,white,SpinThreeHalf);
}

G4QHadron* G4QHadronBuilder::Meson(G4QParton* black, G4QParton* white, Spin theSpin)
{
#ifdef debug
  //  Verify Input Charge
  G4double charge =  black->GetPDGCharge() + white->GetPDGCharge();	 
  if (std::abs(charge)>2 || std::abs(3*(charge-G4int(charge*1.001))) > perCent)
	   G4cerr<<"-Warning-G4QHadronBuilder::Meson: invalid TotalCharge="<<charge<<", PDG_q1=" 
          <<black->GetPDGEncoding()<< ", PDG_q2="<<white->GetPDGEncoding()<<G4endl;
#endif	
	
	 G4int id1= black->GetPDGCode();
	 G4int id2= white->GetPDGCode();
	 if (std::abs(id1) < std::abs(id2))     // exchange black and white
	 {
	   G4int xchg = id1; 
	   id1 = id2;  
	   id2 = xchg;
	 }
	 if(std::abs(id1)>3) 
	 {
				G4cerr<<"***G4QHadronBuilder::Meson: q1="<<id1<<", q2="<<id2
          <<" while CHIPS is only SU(3)"<<G4endl;
    G4Exception("G4QHadronBuilder::Meson:","72",FatalException,"HeavyQuarkFound");
  }	
  G4int PDGEncoding=0;
	 if(!(id1+id2))                      // annihilation case (neutral)
  {     
	   G4double rmix = G4UniformRand();
	   G4int    imix = 2*std::abs(id1) - 1;
	   if(theSpin == SpinZero)
	      PDGEncoding = 110*(1 + G4int(rmix + scalarMesonMixings[imix - 1])
        	                   + G4int(rmix + scalarMesonMixings[imix]    ) ) +  theSpin;
		  else
	      PDGEncoding = 110*(1 + G4int(rmix + vectorMesonMixings[imix - 1])
                    				    + G4int(rmix + vectorMesonMixings[imix]    ) ) +  theSpin;
	 }
  else
  {
	   PDGEncoding = 100 * std::abs(id1) + 10 * std::abs(id2) +  theSpin;  
	   G4bool IsUp = (std::abs(id1)&1) == 0;	// quark 1 is up type quark (u or c?)
	   G4bool IsAnti = id1 < 0; 		           // quark 1 is an antiquark?
	   if( (IsUp && IsAnti) || (!IsUp && !IsAnti) )  PDGEncoding = - PDGEncoding;
	 }
	 G4QHadron* Meson= new G4QHadron(PDGEncoding);
#ifdef debug
	 if(std::abs(black->GetPDGCharge() + white->GetPDGCharge() - Meson->GetCharge()) > .001)
	   G4cout<<"-Warning-G4QHadronBuilder::Meson:wrongCharge, q1="<<black->GetPDGCharge()<<"("
          <<black->->GetParticleName()<<"), q2="<<white->GetPDGCharge()<<"("
          <<white->->GetParticleName()<<"), qM="<<Meson->GetCharge()<<"/"<<PDGEncoding
          <<G4endl;
#endif
	return Meson;
}

G4QHadron* G4QHadronBuilder::Baryon(G4QParton* black, G4QParton* white, Spin theSpin)
{
#ifdef debug
  //  Verify Input Charge
  G4double charge =  black->GetPDGCharge() + white->GetPDGCharge();	 
  if(std::abs(charge) > 2 || std::abs(3*(charge - G4int(charge*1.001))) > perCent )
	   G4cerr<<"-Warning-G4QHadronBuilder::Baryon: invalid TotalCharge="<<charge<<", PDG_q1=" 
          <<black->GetPDGEncoding()<< ", PDG_q2="<<white->GetPDGEncoding()<<G4endl;
#endif	
	 G4int id1= black->GetPDGCode();
	 G4int id2= white->GetPDGCode();
	 if(std::abs(id1) < std::abs(id2))
	 {
	   G4int xchg = id1; 
	   id1 = id2;  
	   id2 = xchg;
	 }
	 if(std::abs(id1)<1000 || std::abs(id2)> 3) 
	 {
				G4cerr<<"***G4QHadronBuilder::Baryon: q1="<<id1<<", q2="<<id2
          <<" can't create a Baryon"<<G4endl;
    G4Exception("G4QHadronBuilder::Baryon:","72",FatalException,"WrongQdQSequence");
  }
	 G4int ifl1= std::abs(id1)/1000;
	 G4int ifl2 = (std::abs(id1) - ifl1 * 1000)/100;
	 G4int diquarkSpin = std::abs(id1)%10; 
	 G4int ifl3 = id2;
	 if (id1 < 0) {ifl1 = - ifl1; ifl2 = - ifl2;}
	 //... Construct baryon, distinguish Lambda and Sigma baryons.
	 G4int kfla = std::abs(ifl1);
	 G4int kflb = std::abs(ifl2);
	 G4int kflc = std::abs(ifl3);
	 G4int kfld = std::max(kfla,kflb);
	       kfld = std::max(kfld,kflc);
	 G4int kflf = std::min(kfla,kflb);
	       kflf = std::min(kflf,kflc);
	 G4int kfle = kfla + kflb + kflc - kfld - kflf;
	 //... baryon with content uuu or ddd or sss has always spin = 3/2
	 if(kfla==kflb && kflb==kflc) theSpin=SpinThreeHalf;   

	 G4int kfll = 0;
	 if(theSpin == SpinHalf && kfld > kfle && kfle > kflf)
  { 
    // Spin J=1/2 and all three quarks different
    // Two states exist: (uds -> lambda or sigma0)
    //   -  lambda: s(ud)0 s : 3122; ie. reverse the two lighter quarks
    //   -  sigma0: s(ud)1 s : 3212
	   if(diquarkSpin == 1 )
    {
	      if ( kfla == kfld) kfll = 1; // heaviest quark in diquark
	      else kfll = G4int(0.25 + G4UniformRand());
	   }   
	   if(diquarkSpin==3 && kfla!=kfld) kfll = G4int(0.75+G4UniformRand());
	 }
	 G4int PDGEncoding=0;
	 if (kfll == 1) PDGEncoding = 1000 * kfld + 100 * kflf + 10 * kfle + theSpin;
	 else           PDGEncoding = 1000 * kfld + 100 * kfle + 10 * kflf + theSpin;
 	if (id1 < 0) PDGEncoding = -PDGEncoding;
	 G4QHadron* Baryon= new G4QHadron(PDGEncoding);
#ifdef debug
	 if(std::abs(black->GetPDGCharge() + white->GetPDGCharge()- Baryon->GetCharge()) > .001)
	   G4cout<<"-Warning-G4QHadronBuilder::Baryon:wrongCharge,dq="<<black->GetPDGCharge()<<"("
          <<black->->GetParticleName()<<"), q="<<white->GetPDGCharge()<<"("
          <<white->->GetParticleName()<<"), qB="<<Baryon->GetCharge()<<"/"<<PDGEncoding
          <<G4endl;
#endif
	return Baryon;
}
