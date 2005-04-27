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
// $Id: FTFModelTest2.cc,v 1.3 2005-04-27 15:26:54 gunter Exp $
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
#define USE_QGSM
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
#include "G4QGSModel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4QGSParticipants.hh"
#include "G4PionPlus.hh"
#include "G4ReactionProductVector.hh"
#include "G4ExcitedStringDecay.hh"
#include "Randomize.hh"
#include "G4StableIsotopes.hh"

//#include "CLHEP/Hist/HBookFile.h"
//#include "CLHEP/Hist/HBookHistogram.h"

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

#include "G4Material.hh"
#include "G4HadronCrossSections.hh"
#include "G4VCrossSectionDataSet.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4HadronInelasticDataSet.hh"
#include "G4IonTable.hh"
#include "G4HadTmpUtil.hh"   // G4lrint

#include <memory> // for the auto_ptr(T>
#include "AIDA/AIDA.h"

//void ConstructParticle();
G4double Rapidity(const G4double p, const  G4double E);


int main()
{

	G4cout << " Welcome; How many events shal I do ? " << std::flush;
	G4int maxEvents;
	G4cin >> maxEvents;

	G4Nucleus lead(207.2, 82.);	// lead
	G4Material* m_lead = new G4Material("Pb",    82., 207.19*g/mole, 11.35*g/cm3);
	G4Nucleus Mg(24.3, 12.0);	// Magnesium
	G4Material* m_Mg = new G4Material("Mg",    12., 24.305*g/mole, 1.738*g/cm3);
	G4Nucleus Be(9.,4.); 		// Beryllium
	G4Material* m_Be = new G4Material("Be",    4.,  9.01*g/mole, 1.848*g/cm3);
	G4Nucleus H(1.,1.);		// Hydrogen (Proton)
	G4Material* m_H = new G4Material("H",     1.,  1.0*g/mole, 1.*g/cm3);
	G4Nucleus nuclei[] = { H, H, Be, Mg, lead };  // 0 is dummy
	
	G4Material * materials[] = {m_H, m_H, m_Be, m_Mg, m_lead };
		
	G4StableIsotopes theIso;

	G4int inucleus;
	G4cout << G4endl;
	for ( unsigned int achoice=1; achoice< sizeof(nuclei)/sizeof(H); achoice++)
	{   G4cout  << achoice << ": " << theIso.GetName(int(nuclei[achoice].GetZ())) 
//	         << "(" << nuclei[achoice].GetN() << "," << nuclei[achoice].GetZ() << ")"
	         << G4endl;
	}
	G4cout  << std::flush;
	G4cin >> inucleus;
	
	G4Nucleus theNucleus = nuclei[inucleus];
	G4Material * material = materials[inucleus];
	
	G4cout <<  " Nucleus : " << theIso.GetName(int(theNucleus.GetZ())) << G4endl;
	
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

	std::ostrstream hbfilename;
	G4VPartonStringModel * model;

#ifdef USE_QGSM
	hbfilename << "qgsm-";
	model= new G4QGSModel<G4QGSParticipants>; //ptparameter*GeV, minExtraMass*MeV,x0Mass*MeV);
	G4VLongitudinalStringDecay * aStringDecayMethod=
				new G4QGSMFragmentation();
//	aStringDecayMethod->Setpspin(pspin);
	G4VStringFragmentation * stringDecay=new G4ExcitedStringDecay(aStringDecayMethod);
	model->SetFragmentationModel(stringDecay);
#else
	hbfilename << "ftf2-"; 
	model= new G4FTFModel ;
//	G4FTFModel model(ptparameter*GeV, minExtraMass*MeV,x0Mass*MeV);
	G4VLongitudinalStringDecay * aStringDecayMethod=
//				new G4LundStringFragmentation(ptparameter*GeV);
				new G4LundStringFragmentation();
//	aStringDecayMethod->SetVectorMesonProbability(pspin);
	G4VStringFragmentation * stringDecay=new G4ExcitedStringDecay(aStringDecayMethod);
	model->SetFragmentationModel(stringDecay);

#endif
	
	
	hbfilename << theIso.GetName(int(theNucleus.GetZ())) << "-"
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
 
// 	HBookFile hbfile(hbfilename.str(), 29);

     // Creating the analysis factory
    std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );

    // Creating the tree factory
    std::auto_ptr< AIDA::ITreeFactory > tf( af->createTreeFactory() );

    // Creating a tree mapped to a new hbook file.
    std::auto_ptr< AIDA::ITree > tree( tf->create( hbfilename.str(),  "hbook", false, true ));
    std::cout << "Tree store : " << tree->storeName() << std::endl;

    // Creating a tuple factory, whose tuples will be handled by the tree
    //   std::auto_ptr< AIDA::ITupleFactory > tpf( af->createTupleFactory( *tree ) );

    // Creating a histogram factory, whose histograms will be handled by the tree
      std::auto_ptr< AIDA::IHistogramFactory > hf( af->createHistogramFactory( *tree ) );

    const G4int nhisto = 500;
//    AIDA::IHistogram1D* h[nhisto];
    //    AIDA::IHistogram2D* h2;
    //AIDA::ITuple* ntuple1 = 0;


		
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
	   			  
//	HBookTuple Tuple_mom = hf->createHistogram1D(" momentum",1,"//");
//	HBookTuple Tuple_fxp = hf->createHistogram1D("protons fx vs p",3,"//"); 
//	HBookTuple Tuple_rap = hf->createHistogram1D(" rapidity",2,"//");

 	AIDA::IHistogram1D* numberOfStrings = hf->createHistogram1D("10","Number of Strings",100,0.,100.);

	AIDA::IHistogram1D* stringMass = hf->createHistogram1D("20","String mass(MeV)",500,0.,5000.);
	AIDA::IHistogram1D* projMass   = hf->createHistogram1D("21","projectile string mass(MeV)",500,0.,5000.);
	AIDA::IHistogram1D* tgtMass    = hf->createHistogram1D("22","target string mass(MeV)",500,0.,5000.);
	AIDA::IHistogram1D* barionMass = hf->createHistogram1D("25","barion mass(MeV)",500,0.,5000.);
	AIDA::IHistogram1D* mesonMass  = hf->createHistogram1D("26","meson mass(MeV)",500,0.,5000.);

	AIDA::IHistogram1D* stringEnergy = hf->createHistogram1D("30","string energy(GeV)",500,0.,proj_momentum/GeV);
	AIDA::IHistogram1D* projEnergy   = hf->createHistogram1D("31","projectile string energy(GeV)",500,0.,proj_momentum/GeV);
	AIDA::IHistogram1D* tgtEnergy    = hf->createHistogram1D("32","target string energy(GeV)",500,0.,proj_momentum/GeV/10.);
	
	AIDA::IHistogram1D* stringEt    = hf->createHistogram1D("40","string transverse energy(GeV)",500,0.,proj_momentum/GeV);
	AIDA::IHistogram1D* projEt      = hf->createHistogram1D("42","string transverse energy(GeV)",500,0.,proj_momentum/GeV);
	AIDA::IHistogram1D* cutstringEt = hf->createHistogram1D("45","string transverse energy(GeV),-0.1@<etarap@<2.9",500,0.,proj_momentum/GeV);
	AIDA::IHistogram1D* cutstringEta= hf->createHistogram1D("46","string transverse energy(GeV),-0.1@<etarap@<5.5",500,0.,proj_momentum/GeV);
		
	AIDA::IHistogram1D* stringRap   = hf->createHistogram1D("50","rapidity (all)",	200,-10,10);
	AIDA::IHistogram1D* stringRapHE = hf->createHistogram1D("51","rapidity, all@>50MeV)",200,-10,10);
	AIDA::IHistogram1D* projRap     = hf->createHistogram1D("53","projectile rapidity",  200,-10,10);
	AIDA::IHistogram1D* barionRap   = hf->createHistogram1D("54","barion rapidity",	200,-10,10);
	AIDA::IHistogram1D* mesonRap    = hf->createHistogram1D("55","meson rapidity",	200,-10,10);
	AIDA::IHistogram1D* chargedRap  = hf->createHistogram1D("56","rapidity charged",	200,-10,10);
	AIDA::IHistogram1D* chposRap    = hf->createHistogram1D("57","rapidity positives",	200,-10,10);
	AIDA::IHistogram1D* chnegRap    = hf->createHistogram1D("58","rapidity negatives",	200,-10,10);
	AIDA::IHistogram1D* neutralRap  = hf->createHistogram1D("59","rapidity neutrals",	200,-10,10);

	AIDA::IHistogram1D* hfeynmanX         = hf->createHistogram1D("60","feynman x",220,-1.1,1.1);
	AIDA::IHistogram1D* hfeynmanX_cplus   = hf->createHistogram1D("61","feynman x, positve",220,-1.1,1.1);
	AIDA::IHistogram1D* hfeynmanX_cminus  = hf->createHistogram1D("63","feynman x, negative",220,-1.1,1.1);
	AIDA::IHistogram1D* hfeynmanX_neutral = hf->createHistogram1D("64","feynman x, neutral",220,-1.1,1.1);
	AIDA::IHistogram1D* hfeynmanX_p       = hf->createHistogram1D("65","feynman x, protons",220,-1.1,1.1);
	AIDA::IHistogram1D* hfeynmanX_pim	 = hf->createHistogram1D("66","feynman x, pi-",220,-1.1,1.1);
	AIDA::IHistogram1D* hfeynmanX_pip	 = hf->createHistogram1D("67","feynman x, pi+",220,-1.1,1.1);
	AIDA::IHistogram1D* hfeynmanX_Km	 = hf->createHistogram1D("68","feynman x, K-",220,-1.1,1.1);
	AIDA::IHistogram1D* hfeynmanX_Kp	 = hf->createHistogram1D("69","feynman x, K+",220,-1.1,1.1);

	AIDA::IHistogram1D* ifeynmanX	 = hf->createHistogram1D("70","feynman x",220,-1.1,1.1);
	AIDA::IHistogram1D* ifeynmanX_cplus	 = hf->createHistogram1D("71","feynman x, positve",220,-1.1,1.1);
	AIDA::IHistogram1D* ifeynmanX_cminus	 = hf->createHistogram1D("73","feynman x, negative",220,-1.1,1.1);
	AIDA::IHistogram1D* ifeynmanX_neutral = hf->createHistogram1D("74","feynman x, neutral",220,-1.1,1.1);
	AIDA::IHistogram1D* ifeynmanX_p	 = hf->createHistogram1D("75","feynman x, protons",220,-1.1,1.1);
	AIDA::IHistogram1D* ifeynmanX_pim	 = hf->createHistogram1D("76","feynman x, pi-",220,-1.1,1.1);
	AIDA::IHistogram1D* ifeynmanX_pip	 = hf->createHistogram1D("77","feynman x, pi+",220,-1.1,1.1);
	AIDA::IHistogram1D* ifeynmanX_Km	 = hf->createHistogram1D("78","feynman x, K-",220,-1.1,1.1);
	AIDA::IHistogram1D* ifeynmanX_Kp	 = hf->createHistogram1D("79","feynman x, K+",220,-1.1,1.1);

	AIDA::IHistogram1D* hptSq	     = hf->createHistogram1D("80","ptsquare (GeV**2)",200,0.,10.);
	AIDA::IHistogram1D* hptSq_cplus   = hf->createHistogram1D("81","ptsquare (GeV**2), positve",200,0.,10.);
	AIDA::IHistogram1D* hptSq_cminus  = hf->createHistogram1D("83","ptsquare (GeV**2), negative",200,0.,10.);
	AIDA::IHistogram1D* hptSq_neutral = hf->createHistogram1D("84","ptsquare (GeV**2), neutral",200,0.,10.);
	AIDA::IHistogram1D* hptSq_p	     = hf->createHistogram1D("85","ptsquare (GeV**2), protons",200,0.,10.);
	AIDA::IHistogram1D* hptSq_pim     = hf->createHistogram1D("86","ptsquare (GeV**2), pi-",200,0.,10.);
	AIDA::IHistogram1D* hptSq_pip     = hf->createHistogram1D("87","ptsquare (GeV**2), pi+",200,0.,10.);
	AIDA::IHistogram1D* hptSq_Km      = hf->createHistogram1D("88","ptsquare (GeV**2), K-",200,0.,10.);
	AIDA::IHistogram1D* hptSq_Kp      = hf->createHistogram1D("89","ptsquare (GeV**2), K+",200,0.,10.);

	AIDA::IHistogram1D* iptSq	     = hf->createHistogram1D("90","ptsquare (GeV**2)",200,0.,10.);
	AIDA::IHistogram1D* iptSq_cplus   = hf->createHistogram1D("91","ptsquare (GeV**2), positve",200,0.,10.);
	AIDA::IHistogram1D* iptSq_cminus  = hf->createHistogram1D("93","ptsquare (GeV**2), negative",200,0.,10.);
	AIDA::IHistogram1D* iptSq_neutral = hf->createHistogram1D("94","ptsquare (GeV**2), neutral",200,0.,10.);
	AIDA::IHistogram1D* iptSq_p	     = hf->createHistogram1D("95","ptsquare (GeV**2), protons",200,0.,10.);
	AIDA::IHistogram1D* iptSq_pim     = hf->createHistogram1D("96","ptsquare (GeV**2), pi-",200,0.,10.);
	AIDA::IHistogram1D* iptSq_pip     = hf->createHistogram1D("97","ptsquare (GeV**2), pi+",200,0.,10.);
	AIDA::IHistogram1D* iptSq_Km	     = hf->createHistogram1D("98","ptsquare (GeV**2), K-",200,0.,10.);
	AIDA::IHistogram1D* iptSq_Kp	     = hf->createHistogram1D("99","ptsquare (GeV**2), K+",200,0.,10.);

	AIDA::IHistogram1D* hpt	   = hf->createHistogram1D("180","pt (GeV)",200,0.,4.);
	AIDA::IHistogram1D* hpt_cplus   = hf->createHistogram1D("181","pt (GeV), positve",200,0.,4.);
	AIDA::IHistogram1D* hpt_cminus  = hf->createHistogram1D("183","pt (GeV), negative",200,0.,4.);
	AIDA::IHistogram1D* hpt_neutral = hf->createHistogram1D("184","pt (GeV), neutral",200,0.,4.);
	AIDA::IHistogram1D* hpt_p	   = hf->createHistogram1D("185","pt (GeV), protons",200,0.,4.);
	AIDA::IHistogram1D* hpt_pim	   = hf->createHistogram1D("186","pt (GeV), pi-",200,0.,4.);
	AIDA::IHistogram1D* hpt_pip	   = hf->createHistogram1D("187","pt (GeV), pi+",200,0.,4.);
	AIDA::IHistogram1D* hpt_Km	   = hf->createHistogram1D("188","pt (GeV), K-",200,0.,4.);
	AIDA::IHistogram1D* hpt_Kp	   = hf->createHistogram1D("189","pt (GeV), K+",200,0.,4.);

	AIDA::IHistogram1D* ipt	   = hf->createHistogram1D("190","pt (GeV)",200,0.,4.);
	AIDA::IHistogram1D* ipt_cplus   = hf->createHistogram1D("191","pt (GeV), positve",200,0.,4.);
	AIDA::IHistogram1D* ipt_cminus  = hf->createHistogram1D("193","pt (GeV), negative",200,0.,4.);
	AIDA::IHistogram1D* ipt_neutral = hf->createHistogram1D("194","pt (GeV), neutral",200,0.,4.);
	AIDA::IHistogram1D* ipt_p	   = hf->createHistogram1D("195","pt (GeV), protons",200,0.,4.);
	AIDA::IHistogram1D* ipt_pim	   = hf->createHistogram1D("196","pt (GeV), pi-",200,0.,4.);
	AIDA::IHistogram1D* ipt_pip	   = hf->createHistogram1D("197","pt (GeV), pi+",200,0.,4.);
	AIDA::IHistogram1D* ipt_Km	   = hf->createHistogram1D("198","pt (GeV), K-",200,0.,4.);
	AIDA::IHistogram1D* ipt_Kp	   = hf->createHistogram1D("199","pt (GeV), K+",200,0.,4.);

	AIDA::IHistogram1D* Mpim = hf->createHistogram1D("100","pi-",30,-0.25,14.75);
	AIDA::IHistogram1D* Mpip = hf->createHistogram1D("101","pi+",30,-0.25,14.75);
	AIDA::IHistogram1D* Mpi0 = hf->createHistogram1D("102","pi0",30,-0.25,14.75);
	
	AIDA::IHistogram1D* MKm  = hf->createHistogram1D("110","K-",30,-0.25,14.75);
	AIDA::IHistogram1D* MKp  = hf->createHistogram1D("112","K+",30,-0.25,14.75);
	AIDA::IHistogram1D* MK0s = hf->createHistogram1D("113","K0s",30,-0.25,14.75);
	AIDA::IHistogram1D* MK0l = hf->createHistogram1D("114","K0l",30,-0.25,14.75);
	
	AIDA::IHistogram1D* MLambda    = hf->createHistogram1D("120","Lambda",30,-0.25,14.75);
	AIDA::IHistogram1D* MLambdaBar = hf->createHistogram1D("121","LambdaBar",30,-0.25,14.75);
	
	AIDA::IHistogram1D* MSigma0    = hf->createHistogram1D("130","Sigma0",30,-0.25,14.75);
	AIDA::IHistogram1D* MSigma0Bar = hf->createHistogram1D("131","Sigma0Bar",30,-0.25,14.75);
	
	AIDA::IHistogram1D* Mproton    = hf->createHistogram1D("140","proton",30,-0.25,14.75);
	AIDA::IHistogram1D* MprotonBar = hf->createHistogram1D("141","antiproton",30,-0.25,14.75);
	
	AIDA::IHistogram1D* Mneutron    = hf->createHistogram1D("142","neutron",30,-0.25,14.75);
	AIDA::IHistogram1D* MneutronBar = hf->createHistogram1D("143","antineutron",30,-0.25,14.75);	


	AIDA::IHistogram1D* Info = hf->createHistogram1D("1","information",10,.5,10.5);


	G4ParticleDefinition * proton = G4Proton::Proton();
	G4double Eprojectile = sqrt(sqr(proj_momentum) + sqr(proton->GetPDGMass()));
	G4double Ekinetic= Eprojectile - proton->GetPDGMass();
	G4double beta = proj_momentum / ( proton->GetPDGMass() + Eprojectile );
	G4double delta_eta = atanh ( beta );        
	G4cout << " offset for rapidity = " << delta_eta << G4endl;
	Info->fill(3.,delta_eta);
        
	G4VCrossSectionDataSet* cs = 0;
	if(proton == G4Proton::Proton() && material->GetElement(0)->GetZ() > 1.5) {
	  cs = new G4ProtonInelasticCrossSection();
	} else if(proton == G4Neutron::Neutron() && material->GetElement(0)->GetZ() > 1.5) {
	  cs = new G4NeutronInelasticCrossSection();
	} else {
	  cs = new G4HadronInelasticDataSet();
	}

	G4DynamicParticle dParticle(proton,G4ThreeVector(1., 0.,0.),Ekinetic);
	G4double cross_sec = 0.0;

	if(cs) {
	  cs->BuildPhysicsTable(*proton);
	  cross_sec = cs->GetCrossSection(&dParticle, material->GetElement(0));
	} else {
	  cross_sec = (G4HadronCrossSections::Instance())->
            GetInelasticCrossSection(&dParticle, material->GetElement(0));
	}

//    G4double factor = cross_sec*MeV*1000.0*(G4double)nbinse/(energy*barn*(G4double)nevt);
//    G4double factora= cross_sec*MeV*1000.0*(G4double)nbinsa/(twopi*2.0*barn*(G4double)nevt);
//    G4double factorb= cross_sec*1000.0/(barn*(G4double)nevt);
//    G4cout << "### factor  = " << factor
//           << "### factora = " << factor
      G4cout << "    cross(b)= " << cross_sec/barn << G4endl;

      Info->fill(1.,double(maxEvents));
      Info->fill(2.,cross_sec);
 
// =========================event loop ================================================================ 
 	for ( G4int ntimes=0; ntimes< maxEvents; ntimes++) {
	
	   G4bool debprint=ntimes<5;

	G4ParticleMomentum direction(0.,0.,1.);
//	G4ParticleMomentum direction(G4UniformRand(),G4UniformRand(),G4UniformRand());
	direction.setMag(1.);

	   G4DynamicParticle primary(proton, 
				     direction,
				     Ekinetic);
      

				     
#if defined(USE_FTFmodel)
	   G4ExcitedStringVector * result = NULL;
	   G4int attempts = 0, maxAttempts=20;
	   while ( result  == NULL )
	   {
		model->Init(theNucleus,primary);
		result =model->GetStrings();
		if (attempts++ > maxAttempts ) 
		{
			G4cout << "G4VPartonStringModel::Scatter(): fails to generate strings"
		       		<< attempts << G4endl;
			maxAttempts *=2;       
		}
	   }
	   
	   
#else

	   G4KineticTrackVector * result;			  
	   result = model->Scatter(theNucleus,primary);

#endif

	   G4LorentzVector pNucleus(0.,0.,0.,model->GetWoundedNucleus()->GetMass());
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
	   for (unsigned int aResult=0; aResult < result1->STL_size(); aResult++)
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
		   for ( unsigned int aSec=0;secondaries!=NULL && aSec<secondaries->STL_size(); aSec++)
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
		  for (unsigned int aSecondary=0; aSecondary<secondaries->STL_size(); aSecondary++)
		  {
		      Result1.STL_push_back((*secondaries)[aSecondary]);
		  }
		  delete Result1[aResult];
		  delete secondaries;
		}
	   }

	   delete result1;
	   
	   Mpim->fill(Npim);
	   Mpip->fill(Npip);
	   Mpi0->fill(Npi0);

	   MKm->fill(NKm);
	   MKp->fill(NKp);
	   MK0s->fill(NK0s);
	   MK0l->fill(NK0l);

	   MLambda->fill(NLambda);
	   MLambdaBar->fill(NLambdaBar);

	   MSigma0->fill(NSigma0);
	   MSigma0Bar->fill(NSigma0Bar);

	   Mproton->fill(Nproton);
	   MprotonBar->fill(NprotonBar);

	   Mneutron->fill(Nneutron);
	   MneutronBar->fill(NneutronBar);

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

	   for (unsigned int bResult=0; bResult < result1->STL_size(); bResult++)
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
		   for (unsigned int aSec=0;secondaries!=NULL && aSec<secondaries->STL_size(); aSec++)
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
		  for (unsigned int aSecondary=0; aSecondary<secondaries->STL_size(); aSecondary++)
		  {
		      (*result).STL_push_back((*secondaries)[aSecondary]);
		  }
		  delete (*result1)[bResult];
		  delete secondaries;
		}
	   }
	   delete result1;		   
#endif
	



	   numberOfStrings->fill(float(result->STL_size()));
	   	
	   G4double Et=0., cutEt=0.,cutEta=0.;
	   G4double Etcurrent;
	   G4double rapidity;
	   G4double Epartsum=0;
	   debout << "Final Kinetic tracks  " << result->STL_size() << G4endl;
	   G4LorentzVector Sum=0;
	   for (astring=0; astring < int(result->STL_size()); astring++)
	   {
	   	if( (*result)[astring] == NULL ) G4cout << "got NULL" << G4endl;
		Etcurrent= (*result)[astring]->Get4Momentum().e() * abs(sin((*result)[astring]->Get4Momentum().theta()));
		Et+= Etcurrent;
		Sum += (*result)[astring]->Get4Momentum();
		
		rapidity=Rapidity((*result)[astring]->Get4Momentum().pz(),(*result)[astring]->Get4Momentum().e());
		
		cutEt += (rapidity> -0.1 && rapidity<2.9 ) ? Etcurrent : 0.;
		cutEta += (rapidity> -0.1 && rapidity<5.5 ) ? Etcurrent : 0.;

		stringRap->fill(rapidity);
		if ( (*result)[astring]->Get4Momentum().e() > 50*MeV ) 
		   {   stringRapHE->fill(rapidity);  
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
//		G4double fx2=2*(toCMS*px2).z()/Ecms;

		G4double ptSquare= (*result)[astring]->Get4Momentum().perp2()/sqr(GeV);
		G4double pt= sqrt(ptSquare);
				
		hfeynmanX->fill(feynmanX);
		hptSq->fill(ptSquare);
		hpt->fill(pt);
		
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
			
			hfeynmanX_p->fill(feynmanX);
			if ( feynmanX > 1. ) 
			{
			    G4LorentzVector pmom=(*result)[astring]->Get4Momentum();
			    G4cout <<" fx >1 " << pmom.x() << ", "
			    		     << pmom.y() << ", "
					     << pmom.z() << ", "
					     << pmom.e() << G4endl;

			}
			hptSq_p->fill(ptSquare);
			hpt_p->fill(pt);
		}
		if ( pname == "pi-" )
		{
			hfeynmanX_pim->fill(feynmanX);
			hptSq_pim->fill(ptSquare);
			hpt_pim->fill(pt);
		}

		if ( pname == "pi+" )
		{
			hfeynmanX_pip->fill(feynmanX);
			hptSq_pip->fill(ptSquare);
			hpt_pip->fill(pt);
		}

		if ( pname == "kaon-" )
		{
			hfeynmanX_Km->fill(feynmanX);
			hptSq_Km->fill(ptSquare);
			hpt_Km->fill(pt);
		}

		if ( pname == "kaon+" )
		{
			hfeynmanX_Kp->fill(feynmanX);
			hptSq_Kp->fill(ptSquare);
			hpt_Kp->fill(pt);
		}

		if      ( charge >  0.5 )
		{
			hfeynmanX_cplus->fill(feynmanX);
			hptSq_cplus->fill(ptSquare);
			hpt_cplus->fill(pt);
		
		} else if ( charge < -0.5 )
		{
			hfeynmanX_cminus->fill(feynmanX);
			hptSq_cminus->fill(ptSquare);
			hpt_cminus->fill(pt);
		} else {
 			hfeynmanX_neutral->fill(feynmanX);
			hptSq_neutral->fill(ptSquare);
			hpt_neutral->fill(pt);
 		}
		
		if ( pname != "proton" || (*result)[astring]->Get4Momentum().vect().mag() > 1200. )
		{
			ifeynmanX->fill(feynmanX);
			iptSq->fill(ptSquare);
			ipt->fill(pt);
			
			if ( pname == "proton" ) 
			{
			
				ifeynmanX_p->fill(feynmanX);
				iptSq_p->fill(ptSquare);
				ipt_p->fill(pt);
			}

			if ( pname == "pi-" )
			{
				ifeynmanX_pim->fill(feynmanX);
				iptSq_pim->fill(ptSquare);
				ipt_pim->fill(pt);
			}
			if ( pname == "pi+" )
			{
				ifeynmanX_pip->fill(feynmanX);
				iptSq_pip->fill(ptSquare);
				ipt_pip->fill(pt);
			}

			if ( pname == "kaon-" )
			{
				ifeynmanX_Km->fill(feynmanX);
				iptSq_Km->fill(ptSquare);
				ipt_Km->fill(pt);
			}
			if ( pname == "kaon+" )
			{
				ifeynmanX_Kp->fill(feynmanX);
				iptSq_Kp->fill(ptSquare);
				ipt_Kp->fill(pt);
			}

			if      ( charge >  0.5 )
			{
				ifeynmanX_cplus->fill(feynmanX);
				iptSq_cplus->fill(ptSquare);
				iptSq_cplus->fill(pt);

			} else if ( charge < -0.5 )
			{
				ifeynmanX_cminus->fill(feynmanX);
				iptSq_cminus->fill(ptSquare);
				iptSq_cminus->fill(pt);
			} else {
 				ifeynmanX_neutral->fill(feynmanX);
				iptSq_neutral->fill(ptSquare);
				iptSq_neutral->fill(pt);
 			}
		}
				
 		stringMass->fill((*result)[astring]->Get4Momentum().mag());
		stringEnergy->fill((*result)[astring]->Get4Momentum().e()/GeV);
		
		if (astring == 0 )
		{
		    projMass->fill((*result)[astring]->Get4Momentum().mag());
		    projEnergy->fill((*result)[astring]->Get4Momentum().e()/GeV);
		    projEt->fill(Et/GeV);
		    projRap->fill(rapidity);
		} else
		{
		    tgtMass->fill((*result)[astring]->Get4Momentum().mag());
		    tgtEnergy->fill((*result)[astring]->Get4Momentum().e()/GeV);
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
		   
		   barionRap->fill(rapidity);
		   barionMass->fill((*result)[astring]->Get4Momentum().mag());
		} else
		{
		   // Meson 
		   mesonRap->fill(rapidity);
		   mesonMass->fill((*result)[astring]->Get4Momentum().mag());
		}
		if ( charge > 0.5 ) 
		{
			chposRap->fill(rapidity);
			chargedRap->fill(rapidity);
		} else if ( charge < -0.5 ) 
		{
			chnegRap->fill(rapidity);
			chargedRap->fill(rapidity);
		} else
		{
			neutralRap->fill(rapidity);
		}
		
	   }
	   
	   stringEt->fill(Et/GeV);
	   cutstringEta->fill(cutEta/GeV);
	   cutstringEt->fill(cutEt/GeV);

//	   std::cout << "Total 4 Momentum: " << Sum << G4endl;
	
//	   Tuple_mom.column("px1",Sum.x());
//	   Tuple_mom.column("py1",Sum.y());
//	   Tuple_mom.column("pz1",Sum.z());
//	   Tuple_mom.column("E1",Sum.e());
	   
// 	   std::for_each(result->begin(),result->end(),DeleteKineticTrack());
	   for (unsigned int i1=0;i1<result->STL_size();i1++)
	   {
	   	delete (*result)[i1];
	   }
	   delete result;
	   
	   G4V3DNucleus * hitNucleus;
	   
	   hitNucleus= model->GetWoundedNucleus();
	   
	   hitNucleus->StartLoop();
	   G4Nucleon * nucleon;
	   G4int allnucleons=0,hitnucleons=0;
	   
	   while ( (nucleon=hitNucleus->GetNextNucleon()) != NULL )
	   {
	   	allnucleons++;
	   	if (nucleon->AreYouHit())
	   	{
//	   		std::cout << "Nucleon HIT " << nucleon->Get4Momentum() << G4endl;
	   	   hitnucleons++;
	   	} else {
//	   		std::cout << "Nucleon     " << nucleon->Get4Momentum() << G4endl;
	   	   Sum += nucleon->Get4Momentum();
	   	}
	   }
//	   std::cout << "Nucleons, hitNucleons " << allnucleons << ", " << hitnucleons << G4endl;
//	   std::cout << "Total 4 Momentum: " << Sum << G4endl;

//	   Tuple_mom.column("px",Sum.x());
//	   Tuple_mom.column("py",Sum.y());
//	   Tuple_mom.column("pz",Sum.z());
//	   Tuple_mom.column("E",Sum.e());

//	   Tuple_mom.dumpData();
	   
	}
	
	std::cout << "  pim : " <<        G4double(Spim) /maxEvents 
			 << " (" << sqrt(G4double(Spim))/maxEvents << ")" << G4endl;
	std::cout << "  pip : " <<        G4double(Spip) /maxEvents 		
			 << " (" << sqrt(G4double(Spip))/maxEvents << ")" << G4endl;
	std::cout << "  pi0 : " <<        G4double(Spi0) /maxEvents 		
			 << " (" << sqrt(G4double(Spi0))/maxEvents << ")" << G4endl;

	std::cout << "  Kp : " <<  	     G4double(SKp) /maxEvents 		
			 << " (" << sqrt(G4double(SKp))/maxEvents << ")" << G4endl;
	std::cout << "  Km : " <<  	     G4double(SKm) /maxEvents 		
			 << " (" << sqrt(G4double(SKm))/maxEvents << ")" << G4endl;
	std::cout << "  K0s : " <<        G4double(SK0s) /maxEvents 		
			 << " (" << sqrt(G4double(SK0s))/maxEvents << ")" << G4endl;
	std::cout << "  K0l : " << 	     G4double(SK0l) /maxEvents 		
			 << " (" << sqrt(G4double(SK0l))/maxEvents << ")" << G4endl;
	
	std::cout << "  Lambda : " <<     G4double(SLambda) /maxEvents 	
			 << " (" << sqrt(G4double(SLambda))/maxEvents << ")" << G4endl;
	std::cout << "  LambdaBar : " <<  G4double(SLambdaBar) /maxEvents 	
			 << " (" << sqrt(G4double(SLambdaBar))/maxEvents << ")" << G4endl;
	
	std::cout << "  Sigma0 : " <<     G4double(SSigma0) /maxEvents 	
			 << " (" << sqrt(G4double(SSigma0))/maxEvents << ")" << G4endl;
	std::cout << "  Sigma0Bar : " <<  G4double(SSigma0Bar) /maxEvents 	
			 << " (" << sqrt(G4double(SSigma0Bar))/maxEvents << ")" << G4endl;

	std::cout << "  proton : " <<     G4double(Sproton) /maxEvents 			
			 << " (" << sqrt(G4double(Sproton))/maxEvents << ")" << G4endl;
	std::cout << "  protonBar : " <<  G4double(SprotonBar) /maxEvents 			
			 << " (" << sqrt(G4double(SprotonBar))/maxEvents << ")" << G4endl;

	std::cout << "  neutron : " <<    G4double(Sneutron) /maxEvents 			
			 << " (" << sqrt(G4double(Sneutron))/maxEvents << ")" << G4endl;
	std::cout << "  neutronBar : " << G4double(SneutronBar) /maxEvents 			
			 << " (" << sqrt(G4double(SneutronBar))/maxEvents << ")" << G4endl;

	std::cout << " "  << G4double(Spim) /maxEvents
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


      std::cout << "Committing..." << std::endl;
      tree->commit();
      std::cout << "Closing the tree..." << std::endl;
      tree->close();
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
