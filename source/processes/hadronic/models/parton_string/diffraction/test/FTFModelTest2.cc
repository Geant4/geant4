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
// $Id: FTFModelTest2.cc,v 1.1 2003-10-08 13:48:52 hpw Exp $
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

#if defined(USE_G4RW)
# define STL_size entries
# define STL_push_back insert
#else
# define STL_size size
# define STL_push_back push_back
#endif

#define debout if(false) G4cout
//#define debout if(debprint) G4cout
#include <strstream> 
  
#include "G4FTFModel.hh"
#include "G4PionPlus.hh"
#include "G4ReactionProductVector.hh"
#include "G4ExcitedStringDecay.hh"
#include "Randomize.hh"
#include "G4StableIsotopes.hh"

#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/HBookHistogram.h"

#include "G4Neutron.hh"
#include "G4Proton.hh"

#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ExcitedMesonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ExcitedBaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Quarks.hh"
#include "G4DiQuarks.hh"

#include "G4ParticleTable.hh"

// #include "G4ParticleTable.hh"
// #include "G4LeptonConstructor.hh" 
// #include "G4BarionConstructor.hh" 
// #include "G4MesonConstructor.hh"
// #include "G4BosonConstructor.hh"
// #include "G4ShortLivedConstructor.hh"
//#include "G4QuarkAnnihilator.hh"

//void ConstructParticle();
G4double Rapidity(const G4double p, const  G4double E);


int main()
{

	G4cout << " Welcome; How many events shal I do ? " << std::flush;
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
	{   G4cout  << achoice << ": " << theIso.GetName(nuclei[achoice].GetZ()) 
//	         << "(" << nuclei[achoice].GetN() << "," << nuclei[achoice].GetZ() << ")"
	         << G4endl;
	}
	G4cout  << std::flush;
	G4cin >> inucleus;
	
	G4Nucleus theNucleus = nuclei[inucleus];
	G4cout <<  " Nucleus : " << theIso.GetName(theNucleus.GetZ()) << G4endl;
	
	G4double ptparameter;
	G4cout << G4endl << " sigma pt in GeV "  << std::flush;
	G4cin >> ptparameter; 
	G4cout <<  " ptparameter = " << ptparameter <<  G4endl;

	G4double minExtraMass;
	G4cout << G4endl << " minExtraMass (MeV) "  << std::flush;
	G4cin >> minExtraMass;
	G4cout <<  " minExtraMass = " << minExtraMass << " MeV" << G4endl;
 
	G4double x0Mass;
	G4cout << G4endl << " x0Mass (MeV) " << std::flush;
	G4cin >> x0Mass; 
	G4cout <<  " x0Mass = " << x0Mass << " MeV" << G4endl;

	G4double pspin;
	G4cout << G4endl << " pspsin (vector/scalar mesons) " << std::flush;
	G4cin >> pspin; 
	G4cout <<  " pspin = " << pspin <<  G4endl;

	G4cout << G4endl << " Projectile Momentum (GeV) " << std::flush;
        G4double proj_momentum;
        G4cin >> proj_momentum;
	G4cout <<  " proj_momentum = " << proj_momentum << " GeV" << G4endl;
         		
	G4cout << G4endl;
	
       
	G4FTFModel model;
//	G4FTFModel model(ptparameter*GeV, minExtraMass*MeV,x0Mass*MeV);
	G4VLongitudinalStringDecay * aStringDecayMethod=
//				new G4LundStringFragmentation(ptparameter*GeV);
				new G4LundStringFragmentation();
//	aStringDecayMethod->SetVectorMesonProbability(pspin);
	G4VStringFragmentation * stringDecay=new G4ExcitedStringDecay(aStringDecayMethod);
	model.SetFragmentationModel(stringDecay);

	
	
#ifndef __DECCXX
	std::ostrstream hbfilename;
	hbfilename << "ftf2-" << theIso.GetName(theNucleus.GetZ()) << "-"
	           << proj_momentum << "-"
	           << ptparameter << "-"
	           << minExtraMass <<  "-"
	           << x0Mass <<  "-"
	           << pspin   
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

	G4int Spim = 0;
	G4int Spip = 0;
	G4int Spi0 = 0;

	G4int SKm = 0;
	G4int SKp = 0;
	G4int SK0s = 0;
	G4int SK0l = 0;

	G4int SLambda = 0;
	G4int SLambdaBar = 0;

	G4int SSigma0 = 0;
	G4int SSigma0Bar = 0;

	G4int Sproton = 0;
	G4int SprotonBar = 0;

	G4int Sneutron = 0;
	G4int SneutronBar = 0;
// ----------- now get all particles of interest ---------
//	ConstructParticle();
   
	G4BosonConstructor Bosons;
	Bosons.ConstructParticle();

	G4MesonConstructor Mesons;
	Mesons.ConstructParticle();

	G4LeptonConstructor Leptons;
	Leptons.ConstructParticle();

	G4BaryonConstructor Baryons;
	Baryons.ConstructParticle();

	G4ShortLivedConstructor ShortLived;
	ShortLived.ConstructParticle();
	
	G4int Acoding[]={0,2212,2112,2222,2122,1212,1112,0};
	G4int Encoding;
	
	for ( G4int aCode=0; (Encoding=Acoding[aCode]) !=0; aCode++ )
	{
	   G4ParticleDefinition* ptr = G4ParticleTable::GetParticleTable()->FindParticle(Encoding);
	   if (ptr == NULL)
	   {
             G4cout << "Particle with encoding "<<Encoding<<" does not exist!!!"<<G4endl;
	   } else {
             G4cout << " code / name : "<< Encoding<<" / "<<ptr->GetParticleName()<<G4endl;
	   }
	}
//------- 
	   			  
	HBookTuple Tuple_mom(" momentum",1,"//");
	HBookTuple Tuple_fxp("protons fx vs p",3,"//"); 
//	HBookTuple Tuple_rap(" rapidity",2,"//");

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
		
	HBookHistogram stringRap  ("rapidity (all)",       200,-10,10,50);
	HBookHistogram stringRapHE("rapidity, all@>50MeV)",200,-10,10,51);
	HBookHistogram projRap    ("projectile rapidity",  200,-10,10,53);
	HBookHistogram barionRap  ("barion rapidity",      200,-10,10,54);
	HBookHistogram mesonRap   ("meson rapidity",       200,-10,10,55);
	HBookHistogram chargedRap ("rapidity charged",     200,-10,10,56);
	HBookHistogram chposRap   ("rapidity positives",   200,-10,10,57);
	HBookHistogram chnegRap   ("rapidity negatives",   200,-10,10,58);
	HBookHistogram neutralRap ("rapidity neutrals",    200,-10,10,59);

	HBookHistogram hfeynmanX("feynman x",220,-1.1,1.1,60);
	HBookHistogram hfeynmanX_cplus("feynman x, positve",220,-1.1,1.1,61);
	HBookHistogram hfeynmanX_cminus("feynman x, negative",220,-1.1,1.1,62);
	HBookHistogram hfeynmanX_neutral("feynman x, neutral",220,-1.1,1.1,63);
	HBookHistogram hfeynmanX_p("feynman x, protons",220,-1.1,1.1,65);
	HBookHistogram hfeynmanX_pim("feynman x, pi-",220,-1.1,1.1,66);
	HBookHistogram hfeynmanX_pip("feynman x, pi+",220,-1.1,1.1,67);
	HBookHistogram hfeynmanX_Km("feynman x, K-",220,-1.1,1.1,68);
	HBookHistogram hfeynmanX_Kp("feynman x, K+",220,-1.1,1.1,69);

	HBookHistogram ifeynmanX("feynman x",220,-1.1,1.1,70);
	HBookHistogram ifeynmanX_cplus("feynman x, positve",220,-1.1,1.1,71);
	HBookHistogram ifeynmanX_cminus("feynman x, negative",220,-1.1,1.1,72);
	HBookHistogram ifeynmanX_neutral("feynman x, neutral",220,-1.1,1.1,73);
	HBookHistogram ifeynmanX_p("feynman x, protons",220,-1.1,1.1,75);
	HBookHistogram ifeynmanX_pim("feynman x, pi-",220,-1.1,1.1,76);
	HBookHistogram ifeynmanX_pip("feynman x, pi+",220,-1.1,1.1,77);
	HBookHistogram ifeynmanX_Km("feynman x, K-",220,-1.1,1.1,78);
	HBookHistogram ifeynmanX_Kp("feynman x, K+",220,-1.1,1.1,79);

	HBookHistogram hptSq("ptsquare (GeV**2)",200,0.,10.,80);
	HBookHistogram hptSq_cplus("ptsquare (GeV**2), positve",200,0.,10.,81);
	HBookHistogram hptSq_cminus("ptsquare (GeV**2), negative",200,0.,10.,82);
	HBookHistogram hptSq_neutral("ptsquare (GeV**2), neutral",200,0.,10.,83);
	HBookHistogram hptSq_p("ptsquare (GeV**2), protons",200,0.,10.,85);
	HBookHistogram hptSq_pim("ptsquare (GeV**2), pi-",200,0.,10.,86);
	HBookHistogram hptSq_pip("ptsquare (GeV**2), pi+",200,0.,10.,87);
	HBookHistogram hptSq_Km("ptsquare (GeV**2), K-",200,0.,10.,88);
	HBookHistogram hptSq_Kp("ptsquare (GeV**2), K+",200,0.,10.,89);

	HBookHistogram iptSq("ptsquare (GeV**2)",200,0.,10.,90);
	HBookHistogram iptSq_cplus("ptsquare (GeV**2), positve",200,0.,10.,91);
	HBookHistogram iptSq_cminus("ptsquare (GeV**2), negative",200,0.,10.,92);
	HBookHistogram iptSq_neutral("ptsquare (GeV**2), neutral",200,0.,10.,93);
	HBookHistogram iptSq_p("ptsquare (GeV**2), protons",200,0.,10.,95);
	HBookHistogram iptSq_pim("ptsquare (GeV**2), pi-",200,0.,10.,96);
	HBookHistogram iptSq_pip("ptsquare (GeV**2), pi+",200,0.,10.,97);
	HBookHistogram iptSq_Km("ptsquare (GeV**2), K-",200,0.,10.,98);
	HBookHistogram iptSq_Kp("ptsquare (GeV**2), K+",200,0.,10.,99);

	HBookHistogram hpt("pt (GeV)",200,0.,4.,180);
	HBookHistogram hpt_cplus("pt (GeV), positve",200,0.,4.,181);
	HBookHistogram hpt_cminus("pt (GeV), negative",200,0.,4.,182);
	HBookHistogram hpt_neutral("pt (GeV), neutral",200,0.,4.,183);
	HBookHistogram hpt_p("pt (GeV), protons",200,0.,4.,185);
	HBookHistogram hpt_pim("pt (GeV), pi-",200,0.,4.,186);
	HBookHistogram hpt_pip("pt (GeV), pi+",200,0.,4.,187);
	HBookHistogram hpt_Km("pt (GeV), K-",200,0.,4.,188);
	HBookHistogram hpt_Kp("pt (GeV), K+",200,0.,4.,189);

	HBookHistogram ipt("pt (GeV)",200,0.,4.,190);
	HBookHistogram ipt_cplus("pt (GeV), positve",200,0.,4.,191);
	HBookHistogram ipt_cminus("pt (GeV), negative",200,0.,4.,192);
	HBookHistogram ipt_neutral("pt (GeV), neutral",200,0.,4.,193);
	HBookHistogram ipt_p("pt (GeV), protons",200,0.,4.,195);
	HBookHistogram ipt_pim("pt (GeV), pi-",200,0.,4.,196);
	HBookHistogram ipt_pip("pt (GeV), pi+",200,0.,4.,197);
	HBookHistogram ipt_Km("pt (GeV), K-",200,0.,4.,198);
	HBookHistogram ipt_Kp("pt (GeV), K+",200,0.,4.,199);

	HBookHistogram Mpim("pi-",30,-0.25,14.75,100);
	HBookHistogram Mpip("pi+",30,-0.25,14.75,101);
	HBookHistogram Mpi0("pi0",30,-0.25,14.75,102);
	
	HBookHistogram MKm("K-",30,-0.25,14.75,110);
	HBookHistogram MKp("K+",30,-0.25,14.75,111);
	HBookHistogram MK0s("K0s",30,-0.25,14.75,112);
	HBookHistogram MK0l("K0l",30,-0.25,14.75,113);
	
	HBookHistogram MLambda("Lambda",30,-0.25,14.75,120);
	HBookHistogram MLambdaBar("LambdaBar",30,-0.25,14.75,121);
	
	HBookHistogram MSigma0("Sigma0",30,-0.25,14.75,130);
	HBookHistogram MSigma0Bar("Sigma0Bar",30,-0.25,14.75,131);
	
	HBookHistogram Mproton("proton",30,-0.25,14.75,140);
	HBookHistogram MprotonBar("antiproton",30,-0.25,14.75,141);
	
	HBookHistogram Mneutron("neutron",30,-0.25,14.75,142);
	HBookHistogram MneutronBar("antineutron",30,-0.25,14.75,143);	
	
	for ( G4int ntimes=0; ntimes< maxEvents; ntimes++) {
	
	   G4bool debprint=ntimes<5;

	
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

#endif

	   G4LorentzVector pNucleus(0.,0.,0.,model.GetWoundedNucleus()->GetMass());
	   G4LorentzRotation toCMS(-1*(pNucleus + primary.Get4Momentum()).boostVector());
	   G4double Ecms=(pNucleus + primary.Get4Momentum()).mag();

	   debout << " Ecms " << Ecms<< G4endl;
// 	   debout << " toCMS xx,xy,xz,xt " << toCMS.xx() <<" "<< toCMS.xy() <<" "<< toCMS.xz() <<" "<< toCMS.xt() << G4endl;
// 	   debout << " toCMS yx,yy,yz,yt " << toCMS.yx() <<" "<< toCMS.yy() <<" "<< toCMS.yz() <<" "<< toCMS.yt() << G4endl;
// 	   debout << " toCMS zx,zy,zz,zt " << toCMS.zx() <<" "<< toCMS.zy() <<" "<< toCMS.zz() <<" "<< toCMS.zt() << G4endl;
// 	   debout << " toCMS tx,ty,tz,tt " << toCMS.tx() <<" "<< toCMS.ty() <<" "<< toCMS.tz() <<" "<< toCMS.tt() << G4endl;
// 	   
// 	   debout << " projectile in cms: " << toCMS*primary.Get4Momentum()<< G4endl;



#if defined(USE_DECAY)    

	   G4KineticTrackVector *result1, *secondaries;
	   result1=result;
	   result=new G4KineticTrackVector();

	   G4int Npim = 0;
	   G4int Npip = 0;
	   G4int Npi0 = 0;

	   G4int NKm = 0;
	   G4int NKp = 0;
	   G4int NK0s = 0;
	   G4int NK0l = 0;

	   G4int NLambda = 0;
	   G4int NLambdaBar = 0;

	   G4int NSigma0 = 0;
	   G4int NSigma0Bar = 0;

	   G4int Nproton = 0;
	   G4int NprotonBar = 0;

	   G4int Nneutron = 0;
	   G4int NneutronBar = 0;
	   
	   debout << " primary kin-tracks produced "<<
	   result1->STL_size()<< G4endl;
	   G4int astring;
	   G4KineticTrackVector &Result =*result;
	   G4KineticTrackVector &Result1=*result1;
/*
 * 	   std::vector<G4KineticTrack *>::iterator IResult=Result.begin();
 * 	   G4int ind=0; 
 * 	   while (IResult++ != Result.end())
 * 	   {
 * 	   	G4ParticleDefinition * pdef = IResult->GetDefinition();
 * 		G4cout << " iter, hand" << pdef->GetParticleName()
 * 			<< Result[ind++]->GetDefinition()->GetParticleName() <<
 * 			G4endl;
 * 	   }
 */
	   for (G4int aResult=0; aResult < result1->STL_size(); aResult++)
	   {
		G4ParticleDefinition * pdef;
		pdef=Result1[aResult]->GetDefinition();
		
		debout << G4endl<< " Primary " << pdef->GetParticleName() ;
		debout << "    width, lifetime " << pdef->GetPDGWidth();
		debout << " , " << pdef->GetPDGLifeTime() << G4endl;

		secondaries=NULL;
		G4String pname=pdef->GetParticleName();
 		debout << " Primary momentumn  : " << Result1[aResult]->Get4Momentum() << G4endl;
		if ( pname == "pi-" || pname == "pi+" || pname == "pi0"
		  || pname == "proton" || pname == "neutron"
		  || pname == "anti_proton" || pname == "anti_neutron"
		  || pname == "kaon+" || pname == "kaon-"
		  || pname == "kaon0L" || pname == "kaon0S"
		  || pname == "lambda" || pname == "anti_lambda"
		  || pname == "sigma0" || pname == "anti_sigma0"
//		  || pdef->GetPDGWidth() <= 0 || pdef->GetPDGLifeTime() > 1*ns 
		  )
		{
		  if ( pname == "pi-" ) Npim++;
		  if ( pname == "pi+" ) Npip++;
		  if ( pname == "pi0" ) Npi0++;
		  
		  if ( pname == "kaon-" ) NKm++;
		  if ( pname == "kaon+" ) NKp++;
		  if ( pname == "kaon0S" ) NK0s++;
		  if ( pname == "kaon0L" ) NK0l++;

		  if ( pname == "lambda" ) NLambda++;
		  if ( pname == "anti_lambda" ) NLambdaBar++;
		  
		  if ( pname == "sigma0" ) NSigma0++;
		  if ( pname == "anti_sigma0" ) NSigma0Bar++;
		  
		  if ( pname == "proton" ) Nproton++;
		  if ( pname == "anti_proton" ) NprotonBar++;
		  
		  if ( pname == "neutron" ) Nneutron++;
		  if ( pname == "anti_neutron" ) NneutronBar++;
		} else {
		   secondaries = Result1[aResult]->Decay();
		   G4LorentzVector checkP=0;
		   for ( G4int aSec=0;secondaries!=NULL && aSec<secondaries->STL_size(); aSec++)
		   {
		   	debout << " P secondary " << (*secondaries)[aSec]->Get4Momentum() << G4endl;
		   	checkP += (*secondaries)[aSec]->Get4Momentum();
		   }
		   debout << " Sum Of decay P " << checkP <<  G4endl;
		}
		if ( secondaries == NULL )
		{
		   Result.STL_push_back(Result1[aResult]);
		   		   
		} else{
		  for (G4int aSecondary=0; aSecondary<secondaries->STL_size(); aSecondary++)
		  {
		      Result1.STL_push_back((*secondaries)[aSecondary]);
		  }
		  delete Result1[aResult];
		  delete secondaries;
		}
	   }

	   delete result1;
	   
	   Mpim.accumulate(Npim);
	   Mpip.accumulate(Npip);
	   Mpi0.accumulate(Npi0);

	   MKm.accumulate(NKm);
	   MKp.accumulate(NKp);
	   MK0s.accumulate(NK0s);
	   MK0l.accumulate(NK0l);

	   MLambda.accumulate(NLambda);
	   MLambdaBar.accumulate(NLambdaBar);

	   MSigma0.accumulate(NSigma0);
	   MSigma0Bar.accumulate(NSigma0Bar);

	   Mproton.accumulate(Nproton);
	   MprotonBar.accumulate(NprotonBar);

	   Mneutron.accumulate(Nneutron);
	   MneutronBar.accumulate(NneutronBar);

	   Spim += Npim;
	   Spip += Npip;
	   Spi0 += Npi0;

	   SKm += NKm;
	   SKp += NKp;
	   SK0s += NK0s;
	   SK0l += NK0l;

	   SLambda += NLambda;
	   SLambdaBar += NLambdaBar;

	   SSigma0 += NSigma0;
	   SSigma0Bar += NSigma0Bar;

	   Sproton += Nproton;
	   SprotonBar += NprotonBar;

	   Sneutron += Nneutron;
	   SneutronBar += NneutronBar;

//  Now decay remaining Sigmas Lamda....

	   result1=result;
	   result=new G4KineticTrackVector();

	   for (G4int bResult=0; bResult < result1->STL_size(); bResult++)
	   {
		G4ParticleDefinition * pdef;
		pdef=(*result1)[bResult]->GetDefinition();
		
		debout << G4endl<< " Primary " << pdef->GetParticleName() ;
		debout << "    width, lifetime " << pdef->GetPDGWidth();
		debout << " , " << pdef->GetPDGLifeTime() << G4endl;

		secondaries=NULL;
		G4String pname=pdef->GetParticleName();
 		debout << " Primary momentumn  : " << (*result)[astring]->Get4Momentum() << G4endl;
		if ( pname == "pi-" || pname == "pi+" || pname == "pi0"
		  || pname == "proton" || pname == "neutron"
		  || pname == "anti_proton" || pname == "anti_neutron"
		  || pname == "kaon+" || pname == "kaon-"
		  || pname == "kaon0L" 
// 		  || pname == "kaon0S"
// 		  || pname == "lambda" || pname == "anti_lambda"
// 		  || pname == "sigma0" || pname == "anti_sigma0"
//		  || pdef->GetPDGWidth() <= 0 || pdef->GetPDGLifeTime() > 1*ns 
		  )
		{;} else {
		   secondaries = (*result1)[bResult]->Decay();
		   G4LorentzVector checkP=0;
		   for ( G4int aSec=0;secondaries!=NULL && aSec<secondaries->STL_size(); aSec++)
		   {
		   	debout << " P secondary... " << (*secondaries)[aSec]->Get4Momentum() << G4endl;
		   	checkP += (*secondaries)[aSec]->Get4Momentum();
		   }
		   debout << " Sum Of decay P... " << checkP <<  G4endl;
		}
		if ( secondaries == NULL )
		{
		   (*result).STL_push_back((*result1)[bResult]);
		   		   
		} else{
		  for (G4int aSecondary=0; aSecondary<secondaries->STL_size(); aSecondary++)
		  {
		      (*result).STL_push_back((*secondaries)[aSecondary]);
		  }
		  delete (*result1)[bResult];
		  delete secondaries;
		}
	   }
	   delete result1;		   
#endif
	



	   numberOfStrings.accumulate(float(result->STL_size()));
	   	
	   G4double Et=0., cutEt=0.,cutEta=0.;
	   G4double Etcurrent;
	   G4double rapidity;
	   G4double Epartsum=0;
	   debout << "Final Kinetic tracks  " << result->STL_size() << G4endl;
	   G4LorentzVector Sum=0;
	   for (astring=0; astring < result->STL_size(); astring++)
	   {
	   	if( (*result)[astring] == NULL ) G4cout << "got NULL" << G4endl;
		Etcurrent= (*result)[astring]->Get4Momentum().e() * abs(sin((*result)[astring]->Get4Momentum().theta()));
		Et+= Etcurrent;
		Sum += (*result)[astring]->Get4Momentum();
		
		rapidity=Rapidity((*result)[astring]->Get4Momentum().pz(),(*result)[astring]->Get4Momentum().e());

		cutEt += (rapidity> -0.1 && rapidity<2.9 ) ? Etcurrent : 0.;
		cutEta += (rapidity> -0.1 && rapidity<5.5 ) ? Etcurrent : 0.;

		stringRap.accumulate(rapidity);
		if ( (*result)[astring]->Get4Momentum().e() > 50*MeV ) 
		   {   stringRapHE.accumulate(rapidity);  
		   }


#ifdef USE_DECAY
		G4String pname=(*result)[astring]->GetDefinition()->GetParticleName();

		debout << " particle is " << pname << G4endl;
		
//		Tuple_rap.column(pname,rapidity);
//		Tuple_rap.dumpData();

		G4double mass=(*result)[astring]->Get4Momentum().mag();

		if ( mass > 1000 
		&& !( pname == "pi-" || pname == "pi+" || pname == "pi0"
		  || pname == "proton" || pname == "neutron"
		  || pname == "anti_proton" || pname == "anti_neutron"
		  || pname == "kaon+" || pname == "kaon-"
		  || pname == "kaon0L" || pname == "kaon0S"
		  || pname == "lambda" || pname == "anti_lambda"
		  || pname == "sigma0" || pname == "anti_sigma0"
		   ))
		{
		   G4String aname=(*result)[astring]->GetDefinition()->GetParticleName();
		   G4cout << " high mass for "<< aname << " : " << mass << G4endl;
		}
#endif

#ifdef USE_FTFmodel
		G4double charge=0;
		G4String pname("string");
#else		
		G4double charge=(*result)[astring]->GetDefinition()->GetPDGCharge();
#endif
		G4double feynmanX=2*(toCMS*(*result)[astring]->Get4Momentum()).z()/Ecms;
		G4LorentzVector px2((*result)[astring]->Get4Momentum().vect(),
		         sqrt((*result)[astring]->Get4Momentum().vect().mag2()+sqr(938.27)));
		G4double fx2=2*(toCMS*px2).z()/Ecms;

		G4double ptSquare= (*result)[astring]->Get4Momentum().perp2()/sqr(GeV);
		G4double pt= sqrt(ptSquare);
				
		hfeynmanX.accumulate(feynmanX);
		hptSq.accumulate(ptSquare);
		hpt.accumulate(pt);
		
#ifdef USE_FTFmodel
		if ( pname == "string" ) 
		{
			Tuple_fxp.column("p",(*result)[astring]->Get4Momentum().vect().mag());
			Tuple_fxp.column("pt",pt);
			Tuple_fxp.column("fx",feynmanX);
			Tuple_fxp.column("fx2",fx2);			
			Tuple_fxp.column("mass",(*result)[astring]->Get4Momentum().mag());
			Tuple_fxp.dumpData();
		}
#endif
				
		if ( pname == "proton" ) 
		{
// 			Tuple_fxp.column("p",(*result)[astring]->Get4Momentum().vect().mag());
// 			Tuple_fxp.column("pt",pt);
// 			Tuple_fxp.column("fx",feynmanX);
// 			Tuple_fxp.dumpData();
			
			hfeynmanX_p.accumulate(feynmanX);
			if ( feynmanX > 1. ) 
			{
			    G4LorentzVector pmom=(*result)[astring]->Get4Momentum();
			    G4cout <<" fx >1 " << pmom.x() << ", "
			    		     << pmom.y() << ", "
					     << pmom.z() << ", "
					     << pmom.e() << G4endl;

			}
			hptSq_p.accumulate(ptSquare);
			hpt_p.accumulate(pt);
		}
		if ( pname == "pi-" )
		{
			hfeynmanX_pim.accumulate(feynmanX);
			hptSq_pim.accumulate(ptSquare);
			hpt_pim.accumulate(pt);
		}

		if ( pname == "pi+" )
		{
			hfeynmanX_pip.accumulate(feynmanX);
			hptSq_pip.accumulate(ptSquare);
			hpt_pip.accumulate(pt);
		}

		if ( pname == "kaon-" )
		{
			hfeynmanX_Km.accumulate(feynmanX);
			hptSq_Km.accumulate(ptSquare);
			hpt_Km.accumulate(pt);
		}

		if ( pname == "kaon+" )
		{
			hfeynmanX_Kp.accumulate(feynmanX);
			hptSq_Kp.accumulate(ptSquare);
			hpt_Kp.accumulate(pt);
		}

		if      ( charge >  0.5 )
		{
			hfeynmanX_cplus.accumulate(feynmanX);
			hptSq_cplus.accumulate(ptSquare);
			hpt_cplus.accumulate(pt);
		
		} else if ( charge < -0.5 )
		{
			hfeynmanX_cminus.accumulate(feynmanX);
			hptSq_cminus.accumulate(ptSquare);
			hpt_cminus.accumulate(pt);
		} else {
 			hfeynmanX_neutral.accumulate(feynmanX);
			hptSq_neutral.accumulate(ptSquare);
			hpt_neutral.accumulate(pt);
 		}
		
		if ( pname != "proton" || (*result)[astring]->Get4Momentum().vect().mag() > 1200. )
		{
			ifeynmanX.accumulate(feynmanX);
			iptSq.accumulate(ptSquare);
			ipt.accumulate(pt);
			
			if ( pname == "proton" ) 
			{
			
				ifeynmanX_p.accumulate(feynmanX);
				iptSq_p.accumulate(ptSquare);
				ipt_p.accumulate(pt);
			}

			if ( pname == "pi-" )
			{
				ifeynmanX_pim.accumulate(feynmanX);
				iptSq_pim.accumulate(ptSquare);
				ipt_pim.accumulate(pt);
			}
			if ( pname == "pi+" )
			{
				ifeynmanX_pip.accumulate(feynmanX);
				iptSq_pip.accumulate(ptSquare);
				ipt_pip.accumulate(pt);
			}

			if ( pname == "kaon-" )
			{
				ifeynmanX_Km.accumulate(feynmanX);
				iptSq_Km.accumulate(ptSquare);
				ipt_Km.accumulate(pt);
			}
			if ( pname == "kaon+" )
			{
				ifeynmanX_Kp.accumulate(feynmanX);
				iptSq_Kp.accumulate(ptSquare);
				ipt_Kp.accumulate(pt);
			}

			if      ( charge >  0.5 )
			{
				ifeynmanX_cplus.accumulate(feynmanX);
				iptSq_cplus.accumulate(ptSquare);
				iptSq_cplus.accumulate(pt);

			} else if ( charge < -0.5 )
			{
				ifeynmanX_cminus.accumulate(feynmanX);
				iptSq_cminus.accumulate(ptSquare);
				iptSq_cminus.accumulate(pt);
			} else {
 				ifeynmanX_neutral.accumulate(feynmanX);
				iptSq_neutral.accumulate(ptSquare);
				iptSq_neutral.accumulate(pt);
 			}
		}
				
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
//		isBarion = (*result)[astring]->GetDefinition()->GetPDGEncoding()%10000 >= 1000;
		isBarion = (*result)[astring]->GetDefinition()->GetParticleType() == G4String("baryon");
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
		if ( charge > 0.5 ) 
		{
			chposRap.accumulate(rapidity);
			chargedRap.accumulate(rapidity);
		} else if ( charge < -0.5 ) 
		{
			chnegRap.accumulate(rapidity);
			chargedRap.accumulate(rapidity);
		} else
		{
			neutralRap.accumulate(rapidity);
		}
		
	   }
	   
	   stringEt.accumulate(Et/GeV);
	   cutstringEta.accumulate(cutEta/GeV);
	   cutstringEt.accumulate(cutEt/GeV);

//	   cout << "Total 4 Momentum: " << Sum << G4endl;
	
	   Tuple_mom.column("px1",Sum.x());
	   Tuple_mom.column("py1",Sum.y());
	   Tuple_mom.column("pz1",Sum.z());
	   Tuple_mom.column("E1",Sum.e());
	   
// 	   std::for_each(result->begin(),result->end(),DeleteKineticTrack());
	   for (G4int i1=0;i1<result->STL_size();i1++)
	   {
	   	delete (*result)[i1];
	   }
	   delete result;
	   
	   G4V3DNucleus * hitNucleus;
	   
	   hitNucleus= model.GetWoundedNucleus();
	   
	   hitNucleus->StartLoop();
	   G4Nucleon * nucleon;
	   G4int allnucleons=0,hitnucleons=0;
	   
	   while ( (nucleon=hitNucleus->GetNextNucleon()) != NULL )
	   {
	   	allnucleons++;
	   	if (nucleon->AreYouHit())
	   	{
//	   		cout << "Nucleon HIT " << nucleon->Get4Momentum() << G4endl;
	   	   hitnucleons++;
	   	} else {
//	   		cout << "Nucleon     " << nucleon->Get4Momentum() << G4endl;
	   	   Sum += nucleon->Get4Momentum();
	   	}
	   }
//	   cout << "Nucleons, hitNucleons " << allnucleons << ", " << hitnucleons << G4endl;
//	   cout << "Total 4 Momentum: " << Sum << G4endl;

	   Tuple_mom.column("px",Sum.x());
	   Tuple_mom.column("py",Sum.y());
	   Tuple_mom.column("pz",Sum.z());
	   Tuple_mom.column("E",Sum.e());

	   Tuple_mom.dumpData();
	   
	}
	
	cout << "  pim : " <<        G4double(Spim) /maxEvents 
			 << " (" << sqrt(G4double(Spim))/maxEvents << ")" << G4endl;
	cout << "  pip : " <<        G4double(Spip) /maxEvents 		
			 << " (" << sqrt(G4double(Spip))/maxEvents << ")" << G4endl;
	cout << "  pi0 : " <<        G4double(Spi0) /maxEvents 		
			 << " (" << sqrt(G4double(Spi0))/maxEvents << ")" << G4endl;

	cout << "  Kp : " <<  	     G4double(SKp) /maxEvents 		
			 << " (" << sqrt(G4double(SKp))/maxEvents << ")" << G4endl;
	cout << "  Km : " <<  	     G4double(SKm) /maxEvents 		
			 << " (" << sqrt(G4double(SKm))/maxEvents << ")" << G4endl;
	cout << "  K0s : " <<        G4double(SK0s) /maxEvents 		
			 << " (" << sqrt(G4double(SK0s))/maxEvents << ")" << G4endl;
	cout << "  K0l : " << 	     G4double(SK0l) /maxEvents 		
			 << " (" << sqrt(G4double(SK0l))/maxEvents << ")" << G4endl;
	
	cout << "  Lambda : " <<     G4double(SLambda) /maxEvents 	
			 << " (" << sqrt(G4double(SLambda))/maxEvents << ")" << G4endl;
	cout << "  LambdaBar : " <<  G4double(SLambdaBar) /maxEvents 	
			 << " (" << sqrt(G4double(SLambdaBar))/maxEvents << ")" << G4endl;
	
	cout << "  Sigma0 : " <<     G4double(SSigma0) /maxEvents 	
			 << " (" << sqrt(G4double(SSigma0))/maxEvents << ")" << G4endl;
	cout << "  Sigma0Bar : " <<  G4double(SSigma0Bar) /maxEvents 	
			 << " (" << sqrt(G4double(SSigma0Bar))/maxEvents << ")" << G4endl;

	cout << "  proton : " <<     G4double(Sproton) /maxEvents 			
			 << " (" << sqrt(G4double(Sproton))/maxEvents << ")" << G4endl;
	cout << "  protonBar : " <<  G4double(SprotonBar) /maxEvents 			
			 << " (" << sqrt(G4double(SprotonBar))/maxEvents << ")" << G4endl;

	cout << "  neutron : " <<    G4double(Sneutron) /maxEvents 			
			 << " (" << sqrt(G4double(Sneutron))/maxEvents << ")" << G4endl;
	cout << "  neutronBar : " << G4double(SneutronBar) /maxEvents 			
			 << " (" << sqrt(G4double(SneutronBar))/maxEvents << ")" << G4endl;

	cout << " "  << G4double(Spim) /maxEvents
	     << " "  << G4double(Spip) /maxEvents
	     << " "  << G4double(Spi0) /maxEvents

	     << " "  << G4double(SKp) /maxEvents
	     << " "  << G4double(SKm) /maxEvents
	     << " "  << G4double(SK0s) /maxEvents
	
	     << " "  << G4double(SLambda + SSigma0) /maxEvents
	     << " "  << G4double(SLambdaBar + SSigma0Bar) /maxEvents
	
	     << " "  <<  G4double(Sproton) /maxEvents
	     << " "  <<  G4double(SprotonBar) /maxEvents
	     << "  | "
	     << proj_momentum/GeV << "-"
	     << ptparameter << "-"
	     << minExtraMass <<  "-"
	     << x0Mass <<  "-"
	     << pspin	
	     << G4endl;


	hbfile.write();
	return 0;
}

G4double Rapidity(const G4double p, const G4double E)
{
	return 0.5*log((E+p)/(E-p));
} 

//**************************************************
// void ConstructParticle()
//   {
//   // In this method, static member functions should be called
//   // for all particles which you want to use.
//   // This ensures that objects of these particle types will be
//   // created in the program. 
// 
//   G4LeptonConstructor Leptons;
//   Leptons.ConstructParticle();
// 
//   G4BosonConstructor Bosons;
//   Bosons.ConstructParticle();
// 
//   G4MesonConstructor Mesons;
//   Mesons.ConstructParticle();
// 
//   G4BarionConstructor Barions;
//   Barions.ConstructParticle();
// 
//   G4ShortLivedConstructor ShortLived;
//   ShortLived.ConstructParticle();
//   }
