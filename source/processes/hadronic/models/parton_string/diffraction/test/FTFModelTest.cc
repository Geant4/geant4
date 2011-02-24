//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: FTFModelTest.cc,v 1.1 2003-10-08 13:48:52 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ------------------------------------------------------------
//      GEANT 4 file
//
//
//             by Gunter Folger, June 1998.
//       class exercising G4FTFModel class.
// ------------------------------------------------------------
	   
//#define USE_FTFmodel
#define USE_DECAY

#if defined(USE_DECAY)
# undef USE_FTFmodel
#endif


#include "G4FTFModel.hh"
#include "G4PionPlus.hh"
#include "G4ReactionProductVector.hh"
#include "G4ExcitedStringDecay.hh"
#include "Randomize.hh"
#include "G4StableIsotopes.hh"
// Solve templates.


#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/HBookHistogram.h"

#include "G4ParticleTable.hh"
#include "G4LeptonConstructor.hh" 
#include "G4BarionConstructor.hh" 
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
//#include "G4QuarkAnnihilator.hh"

//#define __USE_STD_IOSTREAM
#ifdef __DECCXX
#include <strstream>
#endif

void ConstructParticle();
G4double Rapidity(const G4double p, const  G4double E);


// Code in g4templates.hh controlled by the macro G4_SOLVE_TEMPLATES
#include "g4templates.hh"
#ifdef G4_SOLVE_TEMPLATES
template class G4RWTPtrOrderedVector<G4Parton>;
template class G4RWTPtrOrderedVector<G4ExcitedString>;
template class G4RWTPtrOrderedVector<G4KineticTrack>;
template class G4RWTPtrOrderedVector<G4VSplitableHadron>;
#endif

#ifdef __DECCXX
// why is this needed on DEC6 @@
template class G4RWTPtrOrderedVector<G4Parton>;
#endif


int main()
{

	G4cout << " Welcome; How many events shal I do ? " ;
	G4int maxEvents;
	G4cin >> maxEvents;

	G4Nucleus lead(207.2, 82.);	// lead
	G4Nucleus Be(9.,4.); 		// Beryllium
	G4Nucleus H(1.,1.);		// Hydrogen (Proton)
	G4Nucleus nuclei[] = { H, H, Be, lead };  // 0 is dummy
	
	G4StableIsotopes theIso;

	G4int inucleus;
	cout << G4endl;
	for ( G4int achoice=1; achoice< sizeof(nuclei)/sizeof(H); achoice++)
	{   cout  << achoice << ": " << theIso.GetName(nuclei[achoice].GetZ()) 
//	         << "(" << nuclei[achoice].GetN() << "," << nuclei[achoice].GetZ() << ")"
	         << G4endl;
	}
	G4cin >> inucleus;
	
	G4Nucleus theNucleus = nuclei[inucleus];
	
	G4double ptparameter;
	cout << G4endl << " sigma pt in GeV" ;
	G4cin >> ptparameter; 

	G4double minExtraMass;
	cout << G4endl << " minExtraMass " ;
	G4cin >> minExtraMass; 
 
	cout << G4endl << " Projectile Momentum " ;
        G4double proj_momentum;
        G4cin >> proj_momentum;
         		
	G4cout << G4endl;
	
       
	ConstructParticle();
	
	G4FTFModel model(ptparameter*GeV, minExtraMass);
	G4VLongitudinalStringDecay * aStringDecayMethod=
				new G4LundStringFragmentation(ptparameter*GeV);
	G4VStringFragmentation * stringDecay=new G4ExcitedStringDecay(aStringDecayMethod);
	model.SetFragmentationModel(stringDecay);

	
	
#ifdef __DECCXX
	std::ostrstream hbfilename;
	hbfilename << "ftf1-" << theIso.GetName(theNucleus.GetZ()) << "-"
	           << proj_momentum << "-"
	           << ptparameter << "-"
	           << minExtraMass 
#if defined(USE_FTFmodel)
		   << "-String"
#elif defined(USE_DECAY)
		   << "-decay"
#endif
	           << ".hbook";
 
 	HBookFile hbfile(hbfilename.str(), 29);
 
// 	char * 	hbfilename[256];
// 	sprintf(hbfilename,"ftf1-%s-%d-%d.hbook\0", 
// 			theIso.GetName(theNucleus.GetZ()),
// 			proj_momentum,
// 			ptparameter);

#else
	HBookFile hbfile("ftf-test.hbook", 29);
#endif	
		
	proj_momentum *=GeV;

	HBookHistogram numberOfStrings("Number of Strings",100,0.,100.,10);

	HBookHistogram stringMass("String mass(MeV)",500,0.,5000.,20);
	HBookHistogram projMass("projectile string mass(MeV)",500,0.,5000.,21);
	HBookHistogram tgtMass("target string mass(MeV)",500,0.,5000.,22);
	HBookHistogram barionMass("barion mass(MeV)",500,0.,5000.,25);
	HBookHistogram mesonMass("meson mass(MeV)",500,0.,5000.,26);

	HBookHistogram stringEnergy("string energy(GeV)",500,0.,proj_momentum/GeV,30);
	HBookHistogram projEnergy("projectile string energy(GeV)",500,0.,proj_momentum/GeV,31);
	HBookHistogram tgtEnergy("target string energy(GeV)",500,0.,proj_momentum/GeV/10.,32);
	
	HBookHistogram stringEt("string transverse energy(GeV)",500,0.,proj_momentum/GeV,40);
	HBookHistogram projEt("string transverse energy(GeV)",500,0.,proj_momentum/GeV,42);
	HBookHistogram cutstringEt("string transverse energy(GeV),-0.1@<etarap@<2.9",500,0.,proj_momentum/GeV,45);
	HBookHistogram cutstringEta("string transverse energy(GeV),-0.1@<etarap@<5.5",500,0.,proj_momentum/GeV,46);
		
	HBookHistogram stringRap("rapidity (all)",200,-10,10,50);
	HBookHistogram stringRapHE("rapidity (all @> 50MeV)",200,-10,10,51);
	HBookHistogram projRap("projectile rapidity",200,-10,10,53);
	HBookHistogram barionRap("barion rapidity",200,-10,10,55);
	HBookHistogram mesonRap("meson rapidity",200,-10,10,56);
	
	for ( G4int ntimes=0; ntimes< maxEvents; ntimes++) {
	   G4ParticleMomentum direction(0.,0.,1.);
//	   G4ParticleMomentum direction(G4UniformRand(),G4UniformRand(),G4UniformRand());
	   direction.setMag(1.);
	   G4ParticleDefinition * proton = G4Proton::Proton();
	   G4double Ekinetic=
	        sqrt(sqr(proj_momentum) + sqr(proton->GetPDGMass())) - proton->GetPDGMass();
	         
	   G4DynamicParticle primary(proton, 
				     direction,
				     Ekinetic);
#if defined(USE_FTFmodel)
	   G4ExcitedStringVector * result = NULL;
	   G4int attempts = 0, maxAttempts=20;
	   while ( result  == NULL )
	   {
		model.Init(theNucleus,primary);
		result =model.GetStrings();
		if (attempts++ > maxAttempts ) 
		{
			G4cout << "G4VPartonStringModel::Scatter(): fails to generate strings"
		       		<< attempts << G4endl;
			maxAttempts *=2;       
		}
	   }
	   
	   
#else

	   G4KineticTrackVector * result;			  
	   result = model.Scatter(theNucleus,primary);
	   
#if defined(USE_DECAY)    


	   G4KineticTrackVector *result1, *secondaries;
	   result1=result;
	   result=new G4KineticTrackVector();
	   			  
	   for (G4int aResult=0; aResult < result1->entries(); aResult++)
	   {
		G4ParticleDefinition * pdef;
		pdef=result1->at(aResult)->GetDefinition();
		
// 		cout << G4endl<< " Primary " << pdef->GetParticleName() << G4endl;
// 		cout << "    width, lifetime " << pdef->GetPDGWidth();
// 		cout << " , " << pdef->GetPDGLifeTime() << G4endl;
		secondaries=NULL;
		if ( pdef->GetPDGWidth() > 0 || pdef->GetPDGLifeTime() < 1*ns )
		{
		   secondaries = result1->at(aResult)->Decay();
		}
		if ( secondaries == NULL )
		{
		   result->insert(result1->at(aResult));
		   
		   G4String pname=result1->at(aResult)->GetDefinition()->GetParticleName();
		   if ( pname == "pi-" || pname == "pi+" 
		     || pname == "proton" || pname == "neutron"
		     || pname == "anti_proton" || pname == "anti_neutron"
		     || pname == "e-"|| pname == "e+"|| pname == "gamma"
		     || pname == "kaon+" || pname == "kaon-"|| pname == "kaon0L"
	//	     || pname == "" || pname == ""
		     )
		   {;
		   } else {
		   	G4ParticleDefinition * pdef1=result1->at(aResult)->GetDefinition();
		   	G4cout << "found a stable " << pname 
		   	       << "  ( life time / width : " 
		   	       << pdef1->GetPDGLifeTime() << " / "
		   	       << pdef1->GetPDGWidth()  << " )" << G4endl;
		   }
		} else
		{
		  for (G4int aSecondary=0; aSecondary<secondaries->entries(); aSecondary++)
		  {
		      result1->append(secondaries->at(aSecondary));
		      pdef=secondaries->at(aSecondary)->GetDefinition();
// 		      cout << " Secondary " << pdef->GetParticleName() << G4endl;
// 		      cout << "    width, lifetime " << pdef->GetPDGWidth();
// 		      cout << " , " << pdef->GetPDGLifeTime() << G4endl;
		  }
		  delete result1->at(aResult);
		  delete secondaries;
		}
	   }
	   delete result1;
#endif
#endif	   
	   numberOfStrings.accumulate(float(result->entries()));
	   
// 	   G4LorentzVector cms;
// 	   for ( G4int aResult=0; aResult < result->entries() ; aResult++ )
// 	   {
// 	   	cms += result->at(aResult)->Get4Momentum();
// 	   }
// 	   
// 	   cout << " total  Momentum " << cms << G4endl;
// 	    
// 	   G4LorentzRotation toCms(-1*cms.boostVector());
// 
// 	   for ( G4int aResult=0; aResult < result->entries() ; aResult++ )
// 	   {
// 		G4LorentzVector p=result->at(aResult)->Get4Momentum();
// 	   	cout << "string 4 Momentum in cms" 
// 	   	     << p.transform(toCms)
// 	   	     << G4endl 
// 	   	     << " Mass " 
// 	   	     << result->at(aResult)->Get4Momentum().mag() 
// 	   	     << G4endl;;
// 	   }
	   
	   
	
	   G4double Et=0., cutEt=0.,cutEta=0.;
	   G4double Etcurrent;
	   G4double rapidity;
	   G4LorentzVector Sum=0;
	   G4int astring;
	   for (astring=0; astring < result->entries(); astring++)
	   {
	   	if( (*result)[astring] == NULL ) G4cout << "got NULL" << G4endl;
		Etcurrent= (*result)[astring]->Get4Momentum().e() * abs(sin((*result)[astring]->Get4Momentum().theta()));
		Et+= Etcurrent;
		Sum += (*result)[astring]->Get4Momentum();
		
		rapidity=Rapidity((*result)[astring]->Get4Momentum().pz(),(*result)[astring]->Get4Momentum().e());

		cutEt += (rapidity> -0.1 && rapidity<2.9 ) ? Etcurrent : 0.;
		cutEta += (rapidity> -0.1 && rapidity<5.5 ) ? Etcurrent : 0.;

		stringRap.accumulate(rapidity);
		if ( result->at(astring)->Get4Momentum().e() > 50*MeV ) 
		   {   stringRapHE.accumulate(rapidity);  
		   }
		G4double mass=(*result)[astring]->Get4Momentum().mag();
#ifdef USE_DECAY
		if ( mass > 1000) 
		{
		   G4String aname=result->at(astring)->GetDefinition()->GetParticleName();
		   G4cout << " high mass for "<< aname << " : " << mass << G4endl;
		}
#endif
 		stringMass.accumulate((*result)[astring]->Get4Momentum().mag());
		stringEnergy.accumulate((*result)[astring]->Get4Momentum().e()/GeV);
		
		if (astring == 0 )
		{
		    projMass.accumulate((*result)[astring]->Get4Momentum().mag());
		    projEnergy.accumulate((*result)[astring]->Get4Momentum().e()/GeV);
		    projEt.accumulate(Et/GeV);
		    projRap.accumulate(rapidity);
		} else
		{
		    tgtMass.accumulate((*result)[astring]->Get4Momentum().mag());
		    tgtEnergy.accumulate((*result)[astring]->Get4Momentum().e()/GeV);
		}
		
		G4bool isBarion=false;
#ifndef USE_FTFmodel
		isBarion = (*result)[astring]->GetDefinition()->GetPDGEncoding()%10000 >= 1000;
#endif
		if
		( isBarion ) 
		{
		   // Barion
		   
		   barionRap.accumulate(rapidity);
		   barionMass.accumulate((*result)[astring]->Get4Momentum().mag());
		} else
		{
		   // Meson 
		   mesonRap.accumulate(rapidity);
		   mesonMass.accumulate((*result)[astring]->Get4Momentum().mag());
		}
		
	   }
	   
	   stringEt.accumulate(Et/GeV);
	   cutstringEta.accumulate(cutEta/GeV);
	   cutstringEt.accumulate(cutEt/GeV);

//	   cout << "Total 4 Momentum: " << Sum << G4endl;
	   result->clearAndDestroy();
	   delete result;
	   
	   G4V3DNucleus * hitNucleus;
	   
	   hitNucleus= model.GetWoundedNucleus();
	   
	   hitNucleus->StartLoop();
	   G4Nucleon * nucleon;
	   G4int allnucleons=0,hitnucleons=0;
	   
	   while ( (nucleon=hitNucleus->GetNextNucleon()) != NULL )
	   {
	   	allnucleons++;
	   	if (nucleon->AreYouHit()) hitnucleons++;
	   }
//	   cout << "Nucleons, hitNucleons " << allnucleons << ", " << hitnucleons << G4endl;
	   
	}

	hbfile.write();
	return 0;
}

G4double Rapidity(const G4double p, const G4double E)
{
	return 0.5*log((E+p)/(E-p));
} 

//**************************************************
void ConstructParticle()
  {
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4LeptonConstructor Leptons;
  Leptons.ConstructParticle();

  G4BosonConstructor Bosons;
  Bosons.ConstructParticle();

  G4MesonConstructor Mesons;
  Mesons.ConstructParticle();

  G4BarionConstructor Barions;
  Barions.ConstructParticle();

  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();
  }
