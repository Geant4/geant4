// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QEnvironment.cc,v 1.10 2000-09-18 07:47:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QEnvironment ----------------
//             by Mikhail Kossov, August 2000.
//  class for Multy Quasmon Environment used by the CHIPS Model
// ------------------------------------------------------------
 
//#define debug
//#define pdebug
//#define ppdebug
//#define sdebug

#include "G4QEnvironment.hh"

G4QEnvironment::G4QEnvironment(const G4QHadronVector& projHadrons, const G4int targPDG) :
  theEnvironment(90000000)
{
  G4int nHadrons=projHadrons.entries();       // A#of hadrons in the input Vector
  if(nHadrons<1)
  {
    G4cerr<<"***G4QEnvironment: a#of INPUT QHadrons="<<nHadrons<<"<= 0"<<G4endl;
    G4Exception("***G4QEnvironment: There is no one projectile");
  }
#ifdef pdebug
  // ===== Print out of the input information at Creation time =========
  G4cout<<"G4QEnvironment START: targPDG="<<targPDG<<", nProj="<<nHadrons<<G4endl;
  for(G4int ipr=0; ipr<nHadrons; ipr++)
  {
    G4int           hPDG  = projHadrons[ipr]->GetPDGCode();
    G4LorentzVector h4Mom = projHadrons[ipr]->Get4Momentum();
    G4int           hNFrag= projHadrons[ipr]->GetNFragments();
    G4QContent      hQC   = projHadrons[ipr]->GetQC();
    G4cout<<ipr<<". PDG="<<hPDG<<hQC<<", 4M="<<h4Mom<<", hNFrag="<<hNFrag<<G4endl;
  }
#endif
  G4int nP=theWorld.GetQPEntries();           // A#of init'ed particles in CHIPS World
  G4int nCl=nP-73;       // A#of init'ed clusters in CHIPS World
#ifdef pdebug
  G4cout<<"G4QEnvironment: before Nuclear clusters initialization: nP="<<nP<<G4endl;
#endif
  InitClustersVector(nCl);                    // Initialize clusters as Particles
#ifdef pdebug
  G4cout<<"G4QEnvironment: Nuclear clusters are initialized: nCl="<<nCl<<G4endl;
#endif
  if(targPDG>=90000000)                       // ==> Interaction with a nuclear target
  {
    theEnvironment.InitByPDG(targPDG);        // Create nuclear environment
    for(G4int ih=0; ih<nHadrons; ih++)        // ==> The main LOOP over projQHadrons
    {
      G4QHadron* curHadr=projHadrons[ih];     // Pointer to current projectile Hadron
      G4int hNFrag = curHadr->GetNFragments();// #0 means intermediate (skip)
      if(!hNFrag)                             // => "Final hadron" case
	  {
        if(theEnvironment.GetPDG()==90000000) // Vacuum case
        {
          G4int hPDG  = curHadr->GetPDGCode();// A PDG Code of the projQHadron
          if(!hPDG||hPDG==10)                 // Check for the validity of the QHadron
          {
            G4cout<<"***G4QEnvironment::Constructor: wrong PDG("<<ih<<")="<<hPDG
                  <<", HQC="<<curHadr->GetQC()<<", HM="<<curHadr->GetMass()<<G4endl;
            G4Exception("***G4QEnvironment::Constructor: Can not construct QEnvironment");
          }
          else
          {
            G4int hQ = curHadr->GetQCode();  // One more check for valid of the QHadron
            if(hQ<0)
	        {
              G4cout<<"***G4QEnvironment::Constructor: Q<0, projPDG="<<hPDG<<G4endl;
              G4Exception("***G4QEnvironment::Constructor: Can not construct QEnvironment");
	        }
            else
            {
              G4QHadron* newHadr = new G4QHadron(curHadr);
              theQHadrons.insert(newHadr);  // This is a real hadron and can be filled
            } // End of Q-Code check
          } // End of proper PDG for i-th Hadron
        }
        else                                // Nuclear Environment still exists
		{
          G4LorentzVector h4Mom = curHadr->Get4Momentum();
          G4QContent      hQC   = curHadr->GetQC();
#ifdef pdebug
          G4cout<<"G4QEnvironment: (1) CreateQuasmon h4M="<<h4Mom<<",hQC="<<hQC<<G4endl;
#endif
          CreateQuasmon(hQC, h4Mom);
		} // End of Existing Nuclear Environment case
	  } // End of final hadron case
    } // End of the LOOP over input hadrons
  } // End of interaction with nucleus
  else
  {
    G4cerr<<"***G4QEnvironment: not nuclear targPDG="<<targPDG<<G4endl;
    G4Exception("***G4QEnvironment::HadrQEnvironment: Target is not a Nucleus");
  }
} // End of the G4QEnvironment constructor

G4QEnvironment::G4QEnvironment(const G4QEnvironment &right)
{
  // theQHadrons (Vector)
  G4int nQH             = right.theQHadrons.entries();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
    theQHadrons.insert(curQH);
  }

  theWorld              = right.theWorld;
  nBarClust             = right.nBarClust;
  f2all                 = right.f2all;

  // theQuasmons (Vector)
  G4int nQ              = right.theQuasmons.entries();
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ     = new G4Quasmon(right.theQuasmons[iq]);
    theQuasmons.insert(curQ);
  }

  // theQCandidates (Vector)
  G4int nQC             = right.theQCandidates.entries();
  if(nQC) for(G4int ic=0; ic<nQC; ic++)
  {
    G4QCandidate* curQC = new G4QCandidate(right.theQCandidates[ic]);
    theQCandidates.insert(curQC);
  }

  theEnvironment        = right.theEnvironment;
}

G4QEnvironment::G4QEnvironment(G4QEnvironment* right)
{
  // theQHadrons (Vector)
  G4int nQH             = right->theQHadrons.entries();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right->theQHadrons[ih]);
    theQHadrons.insert(curQH);
  }

  theWorld              = right->theWorld;
  nBarClust             = right->nBarClust;
  f2all                 = right->f2all;

  // theQuasmons (Vector)
  G4int nQ              = right->theQuasmons.entries();
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ     = new G4Quasmon(right->theQuasmons[iq]);
    theQuasmons.insert(curQ);
  }

  // theQCandidates (Vector)
  G4int nQC             = right->theQCandidates.entries();
  if(nQC) for(G4int ic=0; ic<nQC; ic++)
  {
    G4QCandidate* curQC = new G4QCandidate(right->theQCandidates[ic]);
    theQCandidates.insert(curQC);
  }

  theEnvironment        = right->theEnvironment;
}

G4QEnvironment::~G4QEnvironment()
{
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQCandidates"<<G4endl;
#endif
  theQCandidates.clearAndDestroy();
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQuasmons nQ="<<theQuasmons.entries()<<G4endl;
#endif
  theQuasmons.clearAndDestroy();
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQHadrons nH="<<theQHadrons.entries()<<G4endl;
#endif
  theQHadrons.clearAndDestroy();
#ifdef debug
  G4cout<<"~G4QEnvironment: === DONE ==="<<G4endl;
#endif
}

G4double G4QEnvironment::SolidAngle=0.8;    // Part of Solid Angle to capture (@@A-dep.)
G4bool   G4QEnvironment::EnergyFlux=false;  // Flag for Energy Flux use (not Multy Quasmon)
G4double G4QEnvironment::PiPrThresh=141.4;  // Pion Production Threshold for gammas
G4double G4QEnvironment::M2ShiftVir=20000.; // Shift for M2=-Q2=m_pi^2 of the virtual gamma
G4double G4QEnvironment::DiNuclMass=1880.;  // Double Nucleon Mass for virtual normalization
// Fill the private static parameters
void G4QEnvironment::SetParameters(G4double solAn, G4bool efFlag, G4double piThresh,
                                   G4double mpisq, G4double dinum)
{//  =======================================================================================
  EnergyFlux=efFlag;       // Flag for Energy Flux use instead of Multy Quasmon
  SolidAngle=solAn;        // Part of Solid Angle to capture secondaries (@@A-dep)
  PiPrThresh=piThresh;     // Pion Production Threshold for gammas
  M2ShiftVir=mpisq;        // Shift for M2=-Q2=m_pi^2 of the virtual gamma
  DiNuclMass=dinum;        // Double Nucleon Mass for virtual normalization
}

const G4QEnvironment& G4QEnvironment::operator=(const G4QEnvironment &right)
{// ========================================================================
  // theQHadrons (Vector)
  G4int nQH             = right.theQHadrons.entries();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
    theQHadrons.insert(curQH);
  }

  theWorld              = right.theWorld;
  nBarClust             = right.nBarClust;
  f2all                 = right.f2all;

  // theQuasmons (Vector)
  G4int nQ              = right.theQuasmons.entries();
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ     = new G4Quasmon(right.theQuasmons[iq]);
    theQuasmons.insert(curQ);
  }

  // theQCandidates (Vector)
  G4int nQC             = right.theQCandidates.entries();
  if(nQC) for(G4int ic=0; ic<nQC; ic++)
  {
    G4QCandidate* curQC = new G4QCandidate(right.theQCandidates[ic]);
    theQCandidates.insert(curQC);
  }

  theEnvironment        = right.theEnvironment;
}

// Member function for Quasmon Creation & Environment nucleus modification
void G4QEnvironment::CreateQuasmon(const G4QContent& projQC, const G4LorentzVector& proj4M)
{//========================================================================================
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QNucleus vacuum(90000000);
  G4QContent valQ(0.,0.,0.,0.,0.,0.);            // Prototype of the Quasmon's Quark Content
  G4LorentzVector q4Mom(0.,0.,0.,0.);            // Prototype of the Quasmon's 4-momentum
  nBarClust = 1;                                 // By default only quasi-free nucleons
  G4double  projE=proj4M.e();                    // energy of the projectile
  if(projE<=0.)
  {
    G4cout<<"***G4QEnvironment::CreateQuasmon: projE="<<projE<<"<=0, QC="<<projQC<<G4endl;
    G4Exception("***G4QEnvironment::CreateQuasmon: Negative or 0 energy of the projectile");
  }
  G4double  projM2=proj4M.m2();                  // squared mass of the projectile (print & v.gamma)
  G4int     targPDG=theEnvironment.GetPDG();     // PDG Code of the target nucleus
  if(targPDG>90000000)                           // Interaction with a nuclear target
  {
    G4double  tgMass=theEnvironment.GetMass();   // mass of the target (QEnvironment) nucleus
#ifdef pdebug
    G4cout<<"G4QEnvironment::CreateQ:Interact "<<projQC<<proj4M<<"(m2="
          <<projM2<<") + A="<<targPDG<<",M="<<tgMass<<G4endl;
#endif
    G4int envZ=theEnvironment.GetZ();            // A#of protons in the nucleus
    G4int envN=theEnvironment.GetN();            // A#of neutrons in the nucleus
    G4int envS=theEnvironment.GetS();            // A#of lambdas in the nucleus
    G4int nP  =theWorld.GetQPEntries();          // A#of initialized particles in CHIPS World
    G4int nCl =nP-73;                            // A#of initialized clusters in CHIPS World
    if(nCl<0)G4cout<<"***G4QEnv::CreateQ: nP="<<nP<<" for NuclTarg="<<targPDG<<G4endl;
    if     (nCl<3) nBarClust=1;                  // Fix the maximum Baryon Number for clusters
    else if(nCl<9) nBarClust=2;
    else
	{
      G4int v=nCl-9;
      G4int d=v/15;
      G4int r=v%15;
      if(r<7) nBarClust=3+d+d;
      else    nBarClust=4+d+d;
    }
#ifdef pdebug
	G4cout<<"G4QEnvironment::CreateQ:TargNuc Z="<<envZ<<",N="<<envN<<",S="<<envS<<G4endl;
#endif
    theEnvironment.UpdateClusters(nBarClust);    // Clusters are calculated up to BN=nBarCl
#ifdef pdebug
	G4cout<<"G4QEnvironment::CreateQ: Nucleus("<<targPDG<<") is created ("<<nBarClust<<" clast's)";
    for(G4int ic=0;ic<nBarClust;ic++)G4cout<<" #"<<ic<<"("<<theEnvironment.GetProbability(ic)<<")";
    cout<<G4endl;
#endif
    PrepareClusters();                           // Calculate the cluster population
#ifdef pdebug
	G4cout<<"G4QEnvironment::CreateQ: Cluster probab is calculated."<<G4endl;
#endif
    G4bool efFlag=false;                         // Flag of Energy Flow case FALSE(@@==DEFOLT==@@)
    //     ************ Change if necessary to compare Energy Flux & Multy Quasmon **************
    G4int efCounter=0;                           // Counter of Energy Flux particles
    G4QContent EnFlQC(0,0,0,0,0,0);              // Quark Content of Energy Flux
    G4LorentzVector ef4Mom(0.,0.,0.,0.);         // Summed 4-momentum of Energy Flux
    G4int projPDG=projQC.GetSPDGCode();          // Minimum hadron for the projectile QC
    //     ---   Pbar     ---    Nbar  ---  LAMBDAbar  ---  SIGMA-bar  ---  SIGMA0bar  ---  SIGMA+bar
    if((projPDG==-2212||projPDG==-2112||projPDG==-3122||projPDG==-3112||projPDG==-3212||projPDG==-3222)
       && abs(proj4M.rho())<10.)                 // Only for at rest interactions (move to interface)
	{
      // @@ Annihilation on only one baryon is implemented (no annihilation on clusters! @@??) @@
      G4double   zpn=envZ+envN;                  // a#of nucleons in the nucleus
      G4double   rnd=(zpn+envS)*G4UniformRand(); // Random number to find a baryon
      G4int      targBPDG = 0;                   // Bary-Prototype of PDG of Periferal Target
      G4int      targNPDG = 90000000;            // Nucl-Prototype of PDG of Periferal Target
      G4QContent targQC(0,0,0,0,0,0);            // Quark Content of Periferal Target
      if     (rnd<envN)                          // Neutron is a Periferal Target
      {
        targBPDG = 2112;
        targNPDG = 90000001;
        targQC   = neutQC;
	  }
      else if(rnd<zpn)                           // Proton is a Periferal Target
      {
        targBPDG = 2212;
        targNPDG = 90001000;
        targQC   = protQC;
	  }
      else                                       // Lambda is a Periferal Target
      {
        targBPDG = 3122;
        targNPDG = 91000000;
        targQC   = lambQC;
	  }
      theEnvironment.Reduce(targNPDG);             // Subtract periferal baryon from Nucleus
      G4double resMass=theEnvironment.GetGSMass(); // Nuclear mass after baryon subtraction
      G4double barMass=tgMass-resMass;             // Mass of the bound baryon for annihilation
      tgMass=resMass;                              // New mass of theEnvironment
      q4Mom=G4LorentzVector(0,0,0,barMass)+proj4M; // 4-momentum of the intermediate B-Bbar Quasmon
      valQ=targQC+projQC;                          // Quark Content of the intermediate B-Bbar Quasmon
      G4Quasmon* pan = new G4Quasmon(valQ,q4Mom);  // The intermediate B-Bbar Quasmon creation
      G4QNucleus vE = vacuum;                      // The annihilation as in vacuum (@@in NuclMatter?)
      G4QHadronVector* output = pan->Fragment(vE); // Output of "in Vacuum" Annihilation **!DESTROY!**
      delete pan;
      G4QHadronVector input;                       // Input for MultyQuasmon **!!DESTROY!!**
      G4int trgPDG = theEnvironment.GetPDG();      // New PDG Code for the Residual Nucleus
      G4LorentzVector trg4M(0.,0.,0.,resMass);     // New 4-momentum for the Residual Nucleus
      G4int tNH = output->entries();               // For the selection LOOP
      G4ThreeVector dir = RndmDir();               // For the selection in the LOOP
      for (G4int ind=0; ind<tNH; ind++)            // Loop over projectile QHadrons
      {
        G4QHadron*   curHadr = output->at(ind);    // Pointer to the current hadron
        G4int           shDFL= curHadr->GetNFragments();
        G4LorentzVector sh4m = curHadr->Get4Momentum();
        G4ThreeVector   shDIR= sh4m.vect().unit(); // Projection to the random direction
#ifdef pdebug
		G4cout<<"G4QE::CrQ:##"<<ind<<",d="<<shDFL<<",PDG="<<curHadr->GetPDGCode()<<G4endl;
#endif
        if(!shDFL)                                 // Final (not decayed) QHadron (d==0)
		{
#ifdef pdebug
		  G4cout<<"G4QEnv::CrQ: efF="<<efFlag<<", dot="<<dir.dot(shDIR)<<", SA="<<SolidAngle<<G4endl;
#endif
		  if(dir.dot(shDIR)>SolidAngle)            // Sum up these hadrons and make Energy Flow
		  {
            if(efFlag)                             // => Case of Energy Flux approach
		    {
              G4QContent shQC = curHadr->GetQC();  // Quark Content of the Current Hadron
              ef4Mom+=sh4m;
              EnFlQC+=shQC;
              efCounter++;
#ifdef pdebug
              G4int hPDG=curHadr->GetPDGCode();
              G4LorentzVector h4M = curHadr->Get4Momentum();
		      G4cout<<"G4QE::CrQ:#"<<efCounter<<", PDG="<<hPDG<<", h4M="<<h4M<<G4endl;
#endif
		    }
		    else                                   // => Case of MultyQuasmon fragmentation
		    {
              G4QHadron* mqHadron = new G4QHadron(curHadr);
              input.insert(mqHadron);              // "input" is filled by new hadron-copies
#ifdef pdebug
              G4int hPDG=curHadr->GetPDGCode();
              G4LorentzVector h4M = curHadr->Get4Momentum();
		      G4cout<<"G4QE::CrQ: Fill #"<<ind<<", PDG="<<hPDG<<", h4M="<<h4M<<G4endl;
#endif
		    }
		  }
		  else                                     // Direct filling of the output vector
		  {
#ifdef pdebug
            G4int hPDG=curHadr->GetPDGCode();
            G4LorentzVector h4M = curHadr->Get4Momentum();
		    G4cout<<"G4QE::CrQ: Fill OUT #"<<ind<<",PDG="<<hPDG<<",h4M="<<h4M<<G4endl;
#endif
            G4QHadron* curHadron = new G4QHadron(curHadr);
            theQHadrons.insert(curHadron);         // TheQHadrons are filled by new hadron-copies
          }
		}
        delete output->at(ind);                    // DESTROY hadrons of "output"
	  } // End of LOOP over "output" of annihilation
      delete output;                               // DELETE "output"
      if(!efFlag)                                  // ==> Not Energy Flux case: MultyQuasmon case
	  {
        if(!(input.entries())) return;             // *** RETURN *** Without Quasmon creation
        G4Quasmon fakeQ;                           // fake Quasmon to get parameters
        G4double QTemper=fakeQ.GetTemper();        // Temperature defined by user for Quasmons
        G4double QSOverU=fakeQ.GetSOverU();        // S/U defined by user for Quasmons
        G4double QEtaSup=fakeQ.GetEtaSup();        // Eta Suppresion defined by user for Quasmons
        G4Quasmon::SetParameters(180.,.1,.3);      //@@Hardwired parameters for N-barN annihilation
        G4QEnvironment* muq = new G4QEnvironment(input,theEnvironment.GetPDG());
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: before input.clearAndDestroy()"<<G4endl;
#endif
        input.clearAndDestroy();                   // Here we are DESTROING input
        theEnvironment =muq->GetEnvironment();     // Get residual Environment after interaction
        G4QuasmonVector* outQ = muq->GetQuasmons();// Copy of quasmons **!!DESTROY!!**
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: before delete muq"<<G4endl;
#endif
        delete muq;
        G4Quasmon::SetParameters(QTemper,QSOverU,QEtaSup); // Recover user's parameters for Quasmons
	    G4int nMQ = outQ->entries();               // A#of Quasmons in MultyQuasmon output
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: after GetQuasmon nMQ="<<nMQ<<G4endl;
#endif
        if(nMQ) for(G4int mh=0; mh<nMQ; mh++)      // @@ One can escape creation/distruction but...
        {
          G4Quasmon* curQ = new G4Quasmon(outQ->at(mh));
          theQuasmons.insert(curQ);                // Here we fill in theQuasmons new Quasmon-copies 
          delete outQ->at(mh);
		}
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: delete outQ"<<G4endl;
#endif
        delete outQ;
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: befor return"<<G4endl;
#endif
        return;                                    // *** RETURN *** 
      }
      else                                         // ==> Energy Flux case
	  {
        if (!efCounter) return;                    // *** RETURN *** Without Quasmon creation
	  }
	} // End of Hyperon annihilation case
    else EnFlQC=projQC;                            // If it's not anti-baryon don't use Energy Flux
    PrepareInteractionProbabilities(EnFlQC);       // Interaction probabilities for clusters
#ifdef pdebug
	G4cout<<"G4QEnvironment::CreateQ: Interaction Probabilities are calculated"<<G4endl;
#endif
    G4int nCandid = theQCandidates.entries();
    G4double maxP = theQCandidates[nCandid-1]->GetIntegProbability();
    G4double totP = maxP * G4UniformRand();
#ifdef pdebug
	G4cout<<"G4QEnvironment::CreateQ: nC="<<nCandid<<", maxP="<<maxP<<", totP="<<totP<<G4endl;
#endif
    if (!totP)
    {
	  G4cerr<<"***G4QEnv::CreateQ: nC="<<nCandid<<", maxP="<<maxP<<", QE="<<theEnvironment<<G4endl;
      G4Exception("***G4QEnvironment::CreateQ: Can not select a cluster");
	}
    G4int i=0;
    while(theQCandidates[i]->GetIntegProbability()<totP) i++;
    G4QCandidate* curCand = theQCandidates[i];     // Pointer to selected cluster to interact
    G4QContent    curQC   = curCand->GetQC();      // Get Quark Cont. of the selected cluster
    G4QNucleus targClust(curQC.GetP(),curQC.GetN(),curQC.GetL());// Define Cluster as a QNucleus
#ifdef pdebug
	G4cout<<"G4QEnv::CrQ:Clust#"<<i<<" is selected("<<targClust<<") from "<<theEnvironment<<G4endl;
#endif
    theEnvironment.Reduce(targClust.GetPDG());     // Subtract selected cluster from Nucleus
    G4double envMass=theEnvironment.GetGSMass();   // Mass of residual nuclear environment
    if(projPDG==22&&projE<PiPrThresh+(M2ShiftVir+projM2)/DiNuclMass)// Gamma+quark Interaction
	{
      q4Mom=G4LorentzVector(0.,0.,0.,tgMass-envMass);// Photon interacts with BoundedCluster
      valQ=targClust.GetQCZNS();
#ifdef pdebug
      G4cout<<"G4QEnv::CreateQ:Q="<<q4Mom<<valQ<<"+vg="<<proj4M<<", QEnv="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom, proj4M); // Interaction gamma+quark inside
      theQuasmons.insert(curQuasmon);              // New Quasmon is inserted without gamma in Q     
	}
    else
	{
      q4Mom=proj4M+G4LorentzVector(0.,0.,0.,tgMass-envMass);//Projectile + BoundCluster
      valQ=EnFlQC+targClust.GetQCZNS();
#ifdef pdebug
      G4cout<<"G4QEnv::CreateQ: Q="<<q4Mom<<valQ<<", QEnv="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom);
      theQuasmons.insert(curQuasmon);              // New Quasmon is inserted with hadron/gamma in Q     
	}
  }
  else
  {
    G4cerr<<"***G4QEnvironment::CreateQuasmon: Strange targPDG="<<targPDG<<G4endl;
    G4Exception("***G4QEnvironment::HadrQEnvironment: Impossible environment");
  }
}

// Calculate a probability to interact with clusters for the givven PDG of the projectile
void G4QEnvironment::PrepareInteractionProbabilities(const G4QContent& projQC)
//   =========================================================================
{
  static const G4int NUCPDG=90000000;
  G4double sum    = 0.;                             // Sum of probabilities of interaction
  G4double probab = 0.;                             // Interaction probability
  G4double denseB  = 0.;                            // A#of*prob baryons in dense part
  G4double allB = 0.;                               // A#of*prob baryons in the nucleus
  for (G4int index=0; index<theQCandidates.entries(); index++)
  {
    G4QCandidate* curCand=theQCandidates[index];    // Intermediate pointer
    G4int cPDG  = curCand->GetPDGCode();
    if(cPDG>NUCPDG)                                 // ===> Cluster case
	{
      G4int zns= cPDG-NUCPDG;
      G4int nc = zns%1000;                          // N of the cluster
      G4int sz = zns/1000;                          // S*1000+Z
      G4int zc = sz %1000;                          // Z of the cluster
      G4int sc = sz /1000;                          // S of the cluster
      G4int ac = zc+nc+sc;                          // A of the cluster
      G4double nOfCl=curCand->GetPreProbability();  // A number of clusters of the type
      G4double dOfCl=curCand->GetDenseProbability();// A number of clusters in dense region
      if(cPDG==91000000||cPDG==90001000||cPDG==90000001)
	  {
        allB+=nOfCl;
        denseB+=dOfCl;
	  }
      G4QContent pQC=curCand->GetQC();              // Quark Content of the candidate
      G4int pC   = projQC.GetCharge();              // Charge of the projectile
      G4QContent qQC=pQC+projQC;                    // Total Quark content of the Compound
      G4QPDGCode qQPDG(qQC);
      G4int qC   = qQPDG.GetQCode();
      G4int rPDG = qQC.GetSPDGCode();
      G4double d = abs(zc-nc);
      G4double baryn = qQC.GetBaryonNumber();
      G4double charge= qQC.GetCharge();
      G4double dq= abs(baryn-charge-charge);
      G4double fact=1./(1+d);
      //if     (!rPDG||qC<0||pC<0&&ac==1&&cPDG!=2212||d>1) probab=0.;     // negative on any
      //else if(!rPDG||qC<0||pC>0&&ac==1&&cPDG!=2112||d>1) probab=0.;     // positive on any
      //if     (!rPDG||qC<-1||pC<0&&ac==1&&cPDG!=2212||d>1) probab=0.;    // negative on any
      //else if(!rPDG||qC<-1||pC>0&&ac==1&&cPDG!=2112||d>1) probab=0.;    // positive on any
      //OLD//if     (!rPDG||qC<-1||pC<0&&ac==1&&cPDG!=2212||dq>2) probab=0.;    // negative on any
      //OLD//else if(pC>0&&ac==1&&cPDG!=2112) probab=0.;    // positive on any
      //if     (!rPDG||qC<-1||pC<0&&ac==1&&cPDG!=2212||d>1||ac>1) probab=0.;  // negative on any
      //else if(!rPDG||qC<-1||pC>0&&ac==1&&cPDG!=2112||d>1) probab=0.;    // positive on any
      //if(!rPDG||ac<2||qC<0) probab=0.;       // A>1
      //if(!rPDG||ac<2||qC<-1) probab=0.;       // A>1 Chipolino
      //else      probab=nOfCl*ac; // Isotopic focusing (use together with d>1 above)
      //OLD//else      probab=nOfCl*ac*fact; // Isotopic focusing (use together with d>1 above)
      //else      probab=nOfCl*fact; // Isotopic focusing (use together with d>1 above)
      //else      probab=nOfCl; // Isotopic focusing (use together with d>1 above)
      //else      probab=dOfCl; // Isotopic focusing (use together with d>1 above)
      //else      probab=dOfCl*fact; // Isotopic focusing (use together with d>1 above)
      //else      probab=dOfCl*ac*fact; // Isotopic focusing (use together with d>1 above)
      probab=nOfCl*ac*fact; //@@Always@@ Isotopic focusing (use together with d>1 above)
#ifdef sdebug
	  G4cout<<"G4QEnvironment::PrepareInteractionProbabilities:C="<<cPDG<<",P="<<probab<<",ac="
            <<ac<<",dq="<<dq<<",qC="<<qC<<",rPDG="<<rPDG<<",b="<<baryn<<",c="<<charge<<G4endl;
#endif
	}
    else probab=0.;
    sum+=probab;
    curCand->SetIntegProbability(sum);
  }
  if(allB>0.)f2all=(allB-denseB)/allB;
  else       f2all=0.;
} // End of PrepareInteractionProbabilities

//Initialize a Clusters Vector for the Nucleus of the QEnvironment
void G4QEnvironment::InitClustersVector(G4int maxClust)
//   ==================================================
{
#ifdef pdebug
  G4cout<<"G4QEnvironment::InitClustersVector called with nC="<<maxClust<<G4endl;
#endif
  if(maxClust>=0) for (G4int i=0; i<maxClust; i++) 
  {
    G4int clustQCode = i+72; // Q-code of the cluster in the CHIPS World
    G4QPDGCode clustQPDG;
#ifdef sdebug
  G4cout<<"G4QEnvironment::InitClustersVector: Before Init Q ="<<clustQCode<<G4endl;
#endif
    clustQPDG.InitByQCode(clustQCode);
    G4int clusterPDG=clustQPDG.GetPDGCode();
#ifdef sdebug
  G4cout<<"G4QEnvironment::InitClustersVector: Before insert ="<<clusterPDG<<G4endl;
#endif
    theQCandidates.insert(new G4QCandidate(clusterPDG));
#ifdef sdebug
    G4cout<<"G4QEnvironment::InitClustersVector: Cluster # "<<i<<" with code = "
          <<clusterPDG<<", QC="<<clustQPDG.GetQuarkContent()<<G4endl;
#endif
  }
} // End of InitClastersVector

// Fragmentation of the QEnvironment with MultyQuasmon (the main member function)
G4QHadronVector G4QEnvironment::HadronizeQEnvironment()
//              ======================================= ************************
{
  static const G4int  NUCPDG = 90000000;
  static const G4QNucleus vacuum(90000000);
  static const G4QContent PiQC(0,1,0,1,0,0);
  static const G4QContent K0QC(1,0,0,0,0,1);
  static const G4QContent KpQC(0,1,0,0,0,1);
  //static const G4QContent alQC(6,6,0,0,0,0);
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  //static const G4double mAlpha = G4QPDGCode(2112).GetNuclMass(2,2,0);
#ifdef pdebug
  G4cout<<"G4QEnv::HadrQEnv:***===>>> START HADRONIZATION OF Q-ENVIRONMENT<<<===***"<<G4endl;
#endif
  G4int nQuasmons = theQuasmons.entries();
  if(nQuasmons<1)                                // "No Quasmons" case -> Fill QEnviron
  {
    G4int nPDG = theEnvironment.GetPDG();        // PDG code of the residual Nucl.Environ.
    if(nPDG>NUCPDG)
	{
      G4QHadron* rNucleus = new G4QHadron(nPDG); // Create a Hadron for the Environment
      G4LorentzVector e4Mom=theEnvironment.Get4Momentum();
      rNucleus->Set4Momentum(e4Mom);             // Copy 4-momentum of Environment to Hadron
      theQHadrons.insert(rNucleus);              // Fill "new rNucleus", GS - no further decay 
	}
    return theQHadrons;
  }
  if(theEnvironment.GetPDG()==NUCPDG)            // ==> "Environment is Vacuum" case
  {
#ifdef pdebug
    G4cout<<"G4QEnv::HadrQEnv:Vacuum #ofQ="<<nQuasmons<<G4endl;
#endif
    G4QNucleus vE = vacuum;
	if (nQuasmons) for (G4int iq=0; iq<nQuasmons; iq++)
	{
	  G4QHadronVector* output=theQuasmons[iq]->Fragment(vE); // **!!DESTROY!!**
      G4int nHadrons = output->entries();
#ifdef pdebug
      G4cout<<"G4QEnv::HadrQEnv:Vacuum Q#"<<iq<<", nHadr="<<nHadrons<<G4endl;
#endif
      if(nHadrons>0)                             // Transfer QHadrons from Quasmon to Output
	  {
    	for (G4int ih=0; ih<nHadrons; ih++)
        {
          G4QHadron* curH = new G4QHadron(output->at(ih));
          theQHadrons.insert(curH);              // Fill in theQHadrons a new hadron-copies
          delete output->at(ih);
        }
        delete output;
	  }
      else
	  {
        // @@ Temporary for heavy quasmons
        G4QContent QQC=theQuasmons[iq]->GetQC();
        G4int tQBN=QQC.GetBaryonNumber();        // Baryon Number of not decayed Quasmon
        G4QNucleus tqN(QQC);                     // Define the quasmon as a nucleus
        G4double   tqM=tqN.GetMZNS();            // GS Mass
        G4double   ttM=theQuasmons[iq]->Get4Momentum().m(); // Real Mass
        if(tQBN>1&&ttM>tqM)                      // => "Quasmon evaporation" case
		{
          G4QHadron* nuclQ = new G4QHadron(QQC,theQuasmons[iq]->Get4Momentum());
          EvaporateResidual(nuclQ);              // Try to evaporate Quasmon
          theQuasmons[iq]->KillQuasmon();        // Kill evaporated Quasmon
		}
        else if(iq+1<nQuasmons)
		{
          theQuasmons[iq+1]->IncreaseBy(theQuasmons[iq]); // Merge with the second Quasmon
          theQuasmons[iq]->KillQuasmon();        // Kill the week Quasmon
		}
		else
		{
		  G4cerr<<"***G4QEnv::HadrQE: nH="<<nHadrons<<", QQC="<<theQuasmons[iq]->GetQC()
                <<",QM="<<theQuasmons[iq]->Get4Momentum().m()<<G4endl;
          G4Exception("G4QEnvironment::HadronizeQEnvironment: Vacuum Quasmon in NOTHING");
		}
	  }
	}
    return theQHadrons;
  }
  else                                           // ==> "Nuclear environment" case
  {
#ifdef pdebug
    G4cout<<"G4QEnv::HadrQEnv:FRAGMENTATION IN NUCLEAR ENVIRONMENT nQ="<<nQuasmons<<G4endl;
#endif
    G4int  sumstat= 2;                           // Sum of statuses of all Quasmons
    G4bool force  = false;                       // Prototype of the Force Major Flag
    while (sumstat)                              // ===***=== The MAIN FOREVER LOOP ===***===
	{
      G4QContent envQC=theEnvironment.GetQCZNS();// QuarkContent of current NuclearEnvironment
      G4QContent totQC=envQC;                    // Total QuarkContent in the system
      G4double   envM =theEnvironment.GetMass(); // mass of Nuclear Environment (@@GetMZNS())
      G4double   sumM =envM;                     // Sum of all residual masses in the System @@??
      G4LorentzVector env4M=theEnvironment.Get4Momentum();      
      G4LorentzVector tot4M=env4M;               // 4-momentum of the Total System
      sumstat         =0;
      G4int     fCount=0;                        // Counter of successful but not final fragm's

	  for (G4int iq=0; iq<nQuasmons; iq++)       // Sum parameters of Quasmons to make a decision
	  {
	    G4Quasmon*      pQ = theQuasmons[iq];
        G4QContent      QQC= pQ->GetQC();
        totQC             += QQC;
        G4LorentzVector Q4M= pQ->Get4Momentum();
        tot4M             += Q4M;
        G4double        QM = Q4M.m();
        sumM              += QM;
        G4int           Qst= pQ->GetStatus();
        sumstat           += Qst;
#ifdef pdebug
	    G4cout<<"G4QEnv::HadrQEnv: #"<<iq<<", Qst="<<Qst<<", QM="<<Q4M.m()<<", QQC="<<QQC
              <<G4endl;
#endif
        if(Qst==1||Qst==3) fCount++;             // Increment a counter of fragmentations 
	  } // End of summation LOOP over Quasmons
      G4QContent    tQC =totQC;                  // Not subtracted copy for error prints
      G4int      totChg =totQC.GetCharge();      // Total Electric Charge of the Total System
      G4int      totS   =totQC.GetStrangeness(); // Total Strangeness of the Total System
      G4int      totBN  =totQC.GetBaryonNumber();// Total Baryon Number of the Total System
      G4int      NaK    =0;                      // a#of additional anti-Kaons
      G4int      aKPDG  =0;                      // PDG of additional anti-Kaons
      G4double   MaK    =0.;                     // Total Mass of additional anti-Kaons
      G4int      NPi    =0;                      // a#of additional anti-Kaons
      G4int      PiPDG  =0;                      // PDG of additional pions
      G4double   MPi    =0.;                     // Total Mass of additional pions
      if    (totBN>1&&totS<0&&totChg+totChg>=totBN)// => "additional K+" case
	  {
        totChg+=totS;                            // Charge reduction (totS<0!)
        aKPDG=321;
        NaK=-totS;
        MaK=mK*NaK;
        totQC+=totS*KpQC;
	  }
      else if (totBN>1&&totS<0)                  // => "additional K0" case
	  {
        aKPDG=311;
        NaK=-totS;
        MaK=mK0*NaK;
        totQC+=totS*K0QC;
	  }
      if      (totBN>1&&totChg>totBN-totS)       // => "additional DELTA++" case
	  {
        PiPDG=211;
        NPi=totChg-totBN+totS;
        MPi=mPi*NPi;
        totQC-=NPi*PiQC;
        totChg=totBN-totS;
	  }
      else if (totBN>1&&totChg<0)                // => "additional DELTA-" case
	  {
        PiPDG=-211;
        NPi=-totChg;
        MPi=mPi*NPi;
        totQC-=totChg*PiQC;                       // Now anti-Pions must be subtracted
        totChg=0;
	  }
      G4QNucleus     totN(totQC,tot4M);          // Excited nucleus for the Total System
      G4double      totRM=totN.GetMZNS();        // min (GroundSt) Mass of the Subtracted System
      G4double       totM=totRM+MPi+MaK;         // min (GroundSt) Mass of the Total System
      G4int        totPDG=totN.GetPDG();         // Total PDG Code for the Current compound @@??
      G4int           nOH=theQHadrons.entries(); // A#of output hadrons
      G4LorentzVector s4M=tot4M;                 // Total 4-momentum (@@ only for checking)
      if(nOH) for(G4int ih=0; ih<nOH; ih++) s4M+=theQHadrons[ih]->Get4Momentum();     
#ifdef pdebug
	  G4cout<<"G4QEnv::HadrQEnv: BN="<<totBN<<", sumstat="<<sumstat<<", fCount="<<fCount
            <<", QEnv="<<theEnvironment<<", nOUT="<<nOH<<", tLV="<<s4M<<G4endl;
#endif
      if(totBN<2)                                // ==> "Baryons or Mesons" case (@@ anti Baryon)
      {
        totPDG =totQC.GetSPDGCode();
        if(totPDG&&totPDG!=10&&totPDG!=1114&&totPDG!=2224) totM=G4QPDGCode(totPDG).GetMass();
        else if(totPDG==1114) totM=mNeut+mPi;
        else if(totPDG==2224) totM=mProt+mPi;
        else if(totPDG==10)
		{
          G4QChipolino totChip(totQC);         // define the residual as a Chipolino
          totM  =totChip.GetQPDG1().GetMass()+totChip.GetQPDG2().GetMass();
		}
        else
		{
          G4cerr<<"***G4QEnv::HadrQE: totPDG="<<totPDG<<", totQC="<<totQC<<G4endl;
          G4Exception("G4QEnvironment::HadronizeQEnvironment: Total System impossible in CHIPS");
		}
	  }
      G4double   totMass= tot4M.m();             // Total effective Mass
      if(sumstat&&fCount&&!force)                // ==> "Still try to decay Quasmons" case
	  {
	    for (G4int jq=0; jq<nQuasmons; jq++)     // Fragmentation LOOP over Quasmons
	    {
	      G4Quasmon* pQ     = theQuasmons[jq];   // Pointer to the current Quasmon
          G4int      status = pQ->GetStatus();   // Old status of the Quasmon
          if(status)                             // Skip dead Quasmons
		  {
 	        G4QHadronVector* output=pQ->Fragment(theEnvironment); // **!!DESTROY!!**
            status = pQ->GetStatus();            // New status after fragmentation attempt
#ifdef pdebug
	        G4cout<<"G4QEnv::HadrQEnv: after FragmAttempt jq="<<jq<<", status="<<status<<G4endl;
#endif
            G4int nHadrons = output->entries();
            if(!status||status==1||nHadrons)     // Something was filled
			{
              if(nHadrons>0)                     // Transfer QHadrons from Quasmon to Output
	          {
    	        for (G4int ih=0; ih<nHadrons; ih++) theQHadrons.insert(output->at(ih)); //@@??
                pQ->ClearOutput();//Hadr's filled// Clear Frag-output for further fragmentation
	          }
			}
            else if(status<0||status==2)         // => "PANIC or NOTHING was done" case
			{
              if(status<0&&nHadrons)
			  {
		        G4cerr<<"***G4QEnv::HadrQE: nH="<<nHadrons<<"< status="<<status<<G4endl;
                G4Exception("G4QEnvironment::HadronizeQEnvironment: Strange PANIC");
			  }
              else if(status==2)                 // Check PANIC conditions for status=2
			  {
                if(theEnvironment!=vacuum)       // "Nuclear Environment" case
				{
                  G4LorentzVector t4M=pQ->Get4Momentum()+theEnvironment.Get4Momentum();
                  G4double      tM=t4M.m();      // Real total (with environment) mass
                  G4QContent   qQC= pQ->GetQC(); // QuarkContent of the Quasmon
                  G4QContent envQC=theEnvironment.GetQCZNS(); // QuarkCont of NucEnviron
                  G4QContent curQC=envQC+qQC;    // Total Quark Content
                  G4QNucleus curE(curQC);        // Pseudo nucleus for the Total System
                  G4double   curM=curE.GetMZNS();// min mass of the Total System
#ifdef pdebug
    		      G4cout<<"G4QEnv::HadrQEnv:#"<<jq<<", tM="<<tM<<" > gstM="<<curM<<curE<<G4endl;
#endif
                  if(tM<curM) status=-1;         // Q+E is below the Mass Shell - PANIC
                }
                else                             // "Vacuum" case
				{
                  G4LorentzVector t4M=pQ->Get4Momentum();
                  G4QPDGCode QPDGQ=pQ->GetQPDG();// QPDG Code for the Quasmon
                  G4int PDGQ=QPDGQ.GetPDGCode(); // PDG Code of the QUASMON
#ifdef pdebug
				  G4cout<<"G4QEnv::HadrQEnv: vacuum PDGQ="<<PDGQ<<G4endl;
#endif
                  if(!PDGQ) status=-1;           // Unknown Quasmon in Vaquum - PANIC
                  else if (PDGQ!=10)             // @@ Chipolino can wait @@
				  {
                    G4double qM =t4M.m();        // Real mass of the Quasmon
                    G4double gsM=QPDGQ.GetMass();// GSmass of the Quasmon
#ifdef pdebug
    		        G4cout<<"G4QEnv::HadrQEnv:#"<<jq<<", qM="<<qM<<" > gsM="<<gsM<<G4endl;
#endif
					if(abs(qM-gsM)<0.0001)       // "Fill & Kill" Case
					{
                      G4QHadron* resQ = new G4QHadron(PDGQ,t4M);
                      theQHadrons.insert(resQ);  // @@ Check Dibarions @@
                      pQ->KillQuasmon();
					}
                    else if(qM<gsM) status=-1;   // Below Mass Shell - PANIC
				  }
				}
			  }
              if(status<0)                       // Panic: Quasmon is below the mass shell
			  {
                G4int    ppm=jq;                 // Initialized by PANIC Quasmon pointer
                G4int    nRQ=0;                  // Prototype of a#of additional real Quasmons
#ifdef pdebug
    		    G4cout<<"G4QEnv::HadrQEnv: ***PANIC*** for jq="<<jq<<G4endl;
#endif
                G4ThreeVector vp= pQ->Get4Momentum().vect(); // momentum of the PANIC Quasmon
                G4double dpm=1.e+30;             // Just a big number (dot product of momenta)
	            for(G4int ir=0; ir<nQuasmons; ir++)
	            {
                  if(ir!=jq)                     // Skip the current (PANIC) Quasmon
				  {
	                G4Quasmon* rQ = theQuasmons[ir];
                    G4int Qst = rQ->GetStatus(); // Status of a Quasmon
                    if(Qst>0)                    // Skip dead Quasmons
				    {
					  nRQ++;                     // Increment real-Quasmon-counter
                      G4double dp=vp.dot(rQ->Get4Momentum().vect());
                      if(dp<dpm)
					  {
                        ppm=ir;
                        dpm=dp;
					  }
				    }
				  }
                }// End of the partner-search-for-the-PANIC-Quasmon LOOP
                if(nRQ)                            // Merge with the best candidate
    		    {
	              G4Quasmon*      rQ = theQuasmons[ppm];
                  G4QContent      rQC= rQ->GetQC();
                  G4LorentzVector r4M= rQ->Get4Momentum();
                  rQC               += pQ->GetQC();
                  r4M               += pQ->Get4Momentum();
                  rQ->InitQuasmon(rQC, r4M);
                  pQ->KillQuasmon();
			    }
                else
			    {
#ifdef pdebug
		          G4cout<<"***G4QEnv::HadrQEnv: Cann't resolve PANIC, try to Evaporate"<<G4endl;
#endif
                  force=true;
                  break;
			    }
			  }
			}
		  }
	    } // End of fragmentation LOOP over Quasmons (jq)
      }
      else if(totMass>totM+.001)                 // ==> "Try Evaporation or decay" case
	  {
#ifdef pdebug
		G4cout<<"G4QEnv::HadrQEnv: Evap. M="<<totMass<<",GSM="<<totM<<",dM="<<totMass-totM
              <<",NaK="<<NaK<<",KPDG="<<aKPDG<<",NPi="<<NPi<<",PiPDG="<<PiPDG<<G4endl;
#endif
        if(totBN<2)                              // ==> "Baryon/Meson residual Quasmon" case
		{
          if(totPDG==1114||totPDG==2224)         // ==> "DELTA- or DELTA++" case (@@anti DELTA)
		  {
            G4double   mBar=mProt;
			G4int      bPDG=2212;
            G4double   mMes=mPi;
			G4int      mPDG=211;
            if(totPDG==1114)                     // "DELTA-" case
			{
              mBar=mNeut;
              bPDG=2112;
              mPDG=-211;
		    }
            if(totMass<mBar+mMes)
			{
		      G4cerr<<"***G4QEnv::HadrQE: tM="<<totMass<<"< GSM+mPi0="<<totM+mPi0<<G4endl;
              G4Exception("G4QEnvironment::HadronizeQEnvironment: (1) Cann't decay QEnv");
			}
            else
			{
              G4QHadron* delta = new G4QHadron(totPDG,tot4M);
              delta->SetNFragments(2);           // Put a#of Fragments=2
              theQHadrons.insert(delta);         // Fill the residual DELTA to output
              G4LorentzVector b4Mom(0.,0.,0.,mBar);
              G4LorentzVector m4Mom(0.,0.,0.,mMes);
              if(!G4QHadron(tot4M).DecayIn2(b4Mom, m4Mom))
              {
                G4cerr<<"***G4QEnv::HadronizeQEnv: B="<<bPDG<<"(m="<<mBar<<") + M="
					  <<mPDG<<"(m="<<mMes<<") >(?) mD="<<totMass<<G4endl;
    	        G4Exception("G4QEnvironment::HadronizeQEnvironment: D->B+M decay failed");
              }
#ifdef pdebug
	          G4cout<<"G4Quasmon::HadronizeQuasmon: DELTA="<<totPDG<<tot4M<<" -> Bar="
                    <<bPDG<<m4Mom<<" + Mes="<<mPDG<<m4Mom<<G4endl;
#endif
              G4QHadron* curBar = new G4QHadron(bPDG,b4Mom);
              theQHadrons.insert(curBar);        // Fill the baryon to output
              G4QHadron* curMes = new G4QHadron(mPDG,m4Mom);
              theQHadrons.insert(curMes);        // Fill the meson to output
              return theQHadrons;
			}
		  }
          else if(totPDG==10)                    // ==> "Chipolino" case
		  {
            G4QChipolino resChip(totQC);         // define the residual as a Chipolino @@ DONE (752?)
            G4QPDGCode h1QPDG=resChip.GetQPDG1();// QPDG of the first hadron
            G4int      h1PDG=h1QPDG.GetPDGCode();// PDG code of the first hadron
            G4double   h1M  =h1QPDG.GetMass();   // Mass of the first hadron
            G4QPDGCode h2QPDG=resChip.GetQPDG2();// QPDG of the second hadron
            G4int      h2PDG=h2QPDG.GetPDGCode();// PDG code of the second hadron
            G4double   h2M  =h2QPDG.GetMass();   // Mass of the second hadron
            G4LorentzVector h14Mom(0.,0.,0.,h1M);
            G4LorentzVector h24Mom(0.,0.,0.,h2M);
            if(!G4QHadron(tot4M).DecayIn2(h14Mom, h14Mom))
            {
              G4cerr<<"***G4QEnv::HadronizeQEnv: h1="<<h1PDG<<"(m="<<h1M<<") + h2="
				    <<h2PDG<<"(m="<<h2M<<") >(?) mChipo="<<totMass<<G4endl;
    	      G4Exception("G4QEnvironment::HadronizeQEnv: Chipo->1+2 decay failed");
            }
#ifdef pdebug
	        G4cout<<"G4Quasmon::HadronizeQuasmon: Chipo="<<tot4M<<" -> h1="
                  <<h1PDG<<h14Mom<<" + Mes="<<h2PDG<<h24Mom<<G4endl;
#endif
            G4QHadron* curH1 = new G4QHadron(h1PDG,h14Mom);
            theQHadrons.insert(curH1);           // Fill the curH1 to output
            G4QHadron* curH2 = new G4QHadron(h2PDG,h24Mom);
            theQHadrons.insert(curH2);           // Fill the curH2 to output
            return theQHadrons;
		  }
          else if(totBN==1&&totPDG&&totMass<totM+mPi0+.001)// ==> "Baryon+gamma" case
		  {
            G4LorentzVector h4Mom(0.,0.,0.,totM);
            G4LorentzVector g4Mom(0.,0.,0.,0.);
            if(!G4QHadron(tot4M).DecayIn2(h4Mom, g4Mom))
            {
              G4cerr<<"***G4QEnv::HadronizeQEnv: h="<<totPDG<<"(m="<<totM
                    <<") + gamma >(?) mTot="<<totMass<<G4endl;
    	      G4Exception("G4QEnvironment::HadronizeQEnv: Decay in gamma failed");
            }
#ifdef pdebug
	        G4cout<<"G4Quasmon::HadronizeQuasmon: "<<tot4M<<" -> h="
                  <<totPDG<<h4Mom<<" + gamma="<<g4Mom<<G4endl;
#endif
            G4QHadron* curH = new G4QHadron(totPDG,h4Mom);
            theQHadrons.insert(curH);            // Fill the baryon to output
            G4QHadron* curG = new G4QHadron(22,g4Mom);
            theQHadrons.insert(curG);            // Fill the gamma to output
            return theQHadrons;
		  }
          else                                   // ==> "|B|<2 new Quasmon" case
		  {
            G4Quasmon* resid = new G4Quasmon(totQC,tot4M);
            G4QNucleus vacuum(90000000);
 	        G4QHadronVector* output=resid->Fragment(vacuum); // **!!DESTROY!!**
            G4int nHadrons = output->entries();
            if(nHadrons>0)                       // Transfer QHadrons to Output
	        {
    	      for (G4int ih=0; ih<nHadrons; ih++)
              {
                G4QHadron* curH = new G4QHadron(output->at(ih));
                theQHadrons.insert(curH);        // Insert new hadron-copy
                delete output->at(ih);
			  }
              delete output;
	        }
			else
			{
              G4cerr<<"***G4QEnv::HadronizeQEnv: MQ="<<tot4M.m()<<",QC="<<totQC<<G4endl;
			  G4Exception("G4QEnvironment::HadronizeQEnv: Quasmon doesn't decay");
			}
            delete resid;
            return theQHadrons;
		  }
		}
        else if(NaK)                             // ==> "Decay in K0 or K+ + NPi" case
	    {//@@ Can be moved to EvaporateResidual ??
          if(NaK==1&&!NPi)                  // ==> "One anti-strange K" case
		  {
            G4LorentzVector m4Mom(0.,0.,0.,MaK);
            G4LorentzVector n4Mom(0.,0.,0.,totRM);
            if(!G4QHadron(tot4M).DecayIn2(m4Mom, n4Mom))
            {
              G4cerr<<"***G4QEnv::HadronizeQEnv: M="<<aKPDG<<"(m="<<MaK<<") + N="
				    <<totPDG<<"(m="<<totRM<<") >(?) mSN="<<totMass<<G4endl;
    	      G4Exception("G4QEnvironment::HadronizeQEnv: Anti-s-Nucleus decay failed");
            }
#ifdef pdebug
	        G4cout<<"G4Quasmon::HadronizeQuasmon: SN="<<tot4M<<" -> M="
                  <<aKPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<totQC<<G4endl;
#endif
            G4QHadron* curK = new G4QHadron(aKPDG,m4Mom);
            theQHadrons.insert(curK);            // Fill the curK to output
            G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
            EvaporateResidual(curN);             // Try to evaporate residual
		  }
          else if(NaK&&NPi)                      // ==> "Anti-strange K's + DELTA's" case
		  {
            G4LorentzVector m4Mom(0.,0.,0.,MPi);
            G4LorentzVector k4Mom(0.,0.,0.,MaK);
            G4LorentzVector n4Mom(0.,0.,0.,totRM);
            if(!G4QHadron(tot4M).DecayIn3(m4Mom, k4Mom, n4Mom))
            {
              G4cerr<<"***G4QEnv::HadronizeQEnv: K="<<aKPDG<<"(m="<<MaK<<") + PI="<<PiPDG
                    <<"(m="<<MPi<<" + N="<<totPDG<<"(m="<<totRM<<") >(?)SN="<<totMass<<G4endl;
    	      G4Exception("G4QEnvironment::HadronizeQEnv: 2anti-S-Nucleus decay failed");
            }
#ifdef pdebug
	        G4cout<<"G4Quasmon::HadronizeQuasmon: SN="<<tot4M<<" -> nK="<<aKPDG<<k4Mom
                  <<" + nPi="<<PiPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<G4endl;
#endif
            G4LorentzVector onePi=(1./NPi)*m4Mom;// 4-mom of one pion  
            for (G4int ip=0; ip<NPi; ip++)
			{
              G4QHadron* curP = new G4QHadron(PiPDG,onePi);
              theQHadrons.insert(curP);          // Fill the curM to output
	 		}
            G4LorentzVector oneK=(1./NaK)*k4Mom; // 4-mom of one kaon  
            for (G4int jp=0; jp<NaK; jp++)
			{
              G4QHadron* curP = new G4QHadron(aKPDG,oneK);
              theQHadrons.insert(curP);          // Fill the curM to output
	 		}
            G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
            EvaporateResidual(curN);             // Try to evaporate residual
		  }
		  else                                   // ==> "Two anti-strange Kaons" case
		  {
            G4int N1K = NaK/2;                   // First kaon cluster
            G4int N2K = NaK-N1K;                 // Second kaon cluster
            G4double mM  = MaK/NaK;              // Mass of Pi
            G4double m1M = mM*N1K;               // Mass of the first Pi-cluster
            G4double m2M = mM*N2K;               // Mass of the second Pi-cluster
            G4LorentzVector m4Mom(0.,0.,0.,m1M);
            G4LorentzVector k4Mom(0.,0.,0.,m2M);
            G4LorentzVector n4Mom(0.,0.,0.,totRM);
            if(!G4QHadron(tot4M).DecayIn3(m4Mom, k4Mom, n4Mom))
            {
              G4cerr<<"***G4QEnv::HadronizeQEnv: N * K="<<aKPDG<<"(m="<<mM<<") + N="
                    <<totPDG<<"(m="<<totRM<<") >(?)SN="<<totMass<<G4endl;
    	      G4Exception("G4QEnvironment::HadronizeQEnv: 2anti-S-Nucleus decay failed");
            }
#ifdef pdebug
	        G4cout<<"G4Quasmon::HadronizeQuasmon: SN="<<tot4M<<" -> N*K="<<aKPDG
                  <<" (4M1="<<m4Mom<<" + 4M2="<<k4Mom<<") + N="<<totPDG<<n4Mom<<G4endl;
#endif
            G4LorentzVector one1=(1./N1K)*m4Mom;  // 4-mom of one kaon  
            for (G4int ip=0; ip<N1K; ip++)
			{
              G4QHadron* curP = new G4QHadron(aKPDG,one1);
              theQHadrons.insert(curP);          // Fill the curP to output
	 		}
            G4LorentzVector one2=(1./N2K)*k4Mom; // 4-mom of one kaon  
            for (G4int jp=0; jp<N2K; jp++)
			{
              G4QHadron* curP = new G4QHadron(aKPDG,one2);
              theQHadrons.insert(curP);          // Fill the curP to output
	 		}
            G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
            EvaporateResidual(curN);             // Try to evaporate residual
		  }
          return theQHadrons;
		}
        else if(NPi)                             // ==> "Decay in Pi+ or Pi-" case
	    {//@@ Can be moved to EvaporateResidual ??
          if(NPi==1)                        // ==> "One isobar" case
		  {
            G4LorentzVector m4Mom(0.,0.,0.,MPi);
            G4LorentzVector n4Mom(0.,0.,0.,totRM);
            if(!G4QHadron(tot4M).DecayIn2(m4Mom, n4Mom))
            {
              G4cerr<<"***G4QEnv::HadronizeQEnv: M="<<PiPDG<<"(m="<<MPi<<") + N="
				    <<totPDG<<"(m="<<totRM<<") >(?) mSN="<<totMass<<G4endl;
    	      G4Exception("G4QEnvironment::HadronizeQEnv: Isobar-Nucleus decay failed");
            }
#ifdef pdebug
	        G4cout<<"G4Quasmon::HadronizeQuasmon: SN="<<tot4M<<" -> M="
                  <<PiPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<totQC<<G4endl;
#endif
            G4QHadron* curK = new G4QHadron(PiPDG,m4Mom);
            theQHadrons.insert(curK);            // Fill the curK to output
            G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
            EvaporateResidual(curN);             // Try to evaporate residual
		  }
		  else                                   // ==> "Many Isobars" case
		  {
            G4int N1Pi = NPi/2;                  // First pion cluster
            G4int N2Pi = NPi-N1Pi;               // Second pion cluster
            G4double mM  = MPi/NPi;              // Mass of Pi
            G4double m1M = mM*N1Pi;              // Mass of the first Pi-cluster
            G4double m2M = mM*N2Pi;              // Mass of the second Pi-cluster
            G4LorentzVector m4Mom(0.,0.,0.,m1M);
            G4LorentzVector k4Mom(0.,0.,0.,m2M);
            G4LorentzVector n4Mom(0.,0.,0.,totRM);
            if(!G4QHadron(tot4M).DecayIn3(m4Mom, k4Mom, n4Mom))
            {
              G4cerr<<"***G4QEnv::HadronizeQEnv: N * Pi="<<PiPDG<<"(m="<<mM<<") + N="
                    <<totPDG<<"(m="<<totRM<<") >(?)SN="<<totMass<<G4endl;
    	      G4Exception("G4QEnvironment::HadronizeQEnv: ManyIsobars-Nucleus decay failed");
            }
#ifdef pdebug
	        G4cout<<"G4Quasmon::HadronizeQuasmon: SN="<<tot4M<<" -> N*PI="<<PiPDG
                  <<" (4M1="<<m4Mom<<" + 4M2="<<k4Mom<<") + N="<<totPDG<<n4Mom<<G4endl;
#endif
            G4LorentzVector one1=(1./N1Pi)*m4Mom;  // 4-mom of one pion  
            for (G4int ip=0; ip<N1Pi; ip++)
			{
              G4QHadron* curP = new G4QHadron(PiPDG,one1);
              theQHadrons.insert(curP);          // Fill the curP to output
	 		}
            G4LorentzVector one2=(1./N2Pi)*k4Mom; // 4-mom of one pion  
            for (G4int jp=0; jp<N2Pi; jp++)
			{
              G4QHadron* curP = new G4QHadron(PiPDG,one2);
              theQHadrons.insert(curP);          // Fill the curP to output
	 		}
            G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
            EvaporateResidual(curN);             // Try to evaporate residual
		  }
          return theQHadrons;
		}
        theEnvironment.InitByPDG(NUCPDG);        // Cancele the Environment 
        G4QHadron* evH = new G4QHadron(totPDG,tot4M);// Create a Hadron for ResidualNucl
        EvaporateResidual(evH);                  // Try to evaporate residual
        return theQHadrons;
	  }
      else                                       // ==> "Only GSEnvironment exists" case
      { 
        G4double dM=totMass-totM;
#ifdef pdebug
		G4cout<<"G4QEnv::HadrQEnv: Ground State tM-GSM="<<dM<<", GSM="<<totM<<G4endl;
#endif
        if(dM>-0.001)
		{
          G4QHadron* evH = new G4QHadron(totPDG,tot4M);// Create a Hadron for ResidualNucl
          if(dM>0.001||totBN>1) EvaporateResidual(evH);// Try to evaporate residual
          else           theQHadrons.insert(evH);// Fill to OUTPUT as it is
		}
        else if(nQuasmons==1)                    // => "Decay in Quasmon + QEnviron" case
		{
	      G4Quasmon*       pQ = theQuasmons[0];  // Pointer to the only Quasmon          
          G4QPDGCode    QQPDG = pQ->GetQPDG();   // QPDG of the Quasmon
          G4int          QPDG = QQPDG.GetPDGCode();
          if(!QPDG)
		  {
			G4cerr<<"***G4QEnv::HadrQE: Quasmon is unknown QHadron: PDG="<<QPDG<<G4endl;
			G4Exception("G4QEnvironment::HadronizeQEnvironment: (2) Cann't decay QEnv");
		  }
          else if(QPDG==10)                      // => "Quasmon-Chipolino" case
		  {
            G4QContent QQC = pQ->GetQC();        // Quark Content of the Quasmon
            G4QChipolino QChip(QQC);             // define the Quasmon as a Chipolino
            G4QPDGCode h1QPDG=QChip.GetQPDG1();  // QPDG of the first hadron
            G4double   h1M  =h1QPDG.GetMass();   // Mass of the first hadron
            G4QPDGCode h2QPDG=QChip.GetQPDG2();  // QPDG of the second hadron
            G4double   h2M  =h2QPDG.GetMass();   // Mass of the second hadron
            if(h1M+h2M+envM<totMass)             // => "Three parts decay" case
			{
              G4LorentzVector h14M(0.,0.,0.,h1M);
              G4LorentzVector h24M(0.,0.,0.,h2M);
              G4LorentzVector e4M(0.,0.,0.,envM);
              if(!G4QHadron(tot4M).DecayIn3(h14M,h24M,e4M))
			  {
                 G4cerr<<"***G4QEnv::HadQEnv: tM="<<tot4M.m()<<"-> h1="<<h1QPDG<<"(M="
					   <<h1M<<") + h2="<<h1QPDG<<"(M="<<h2M<<") + eM="<<envM<<G4endl;
				 G4Exception("G4QEnv::HadrQEnv:QChipo+Env DecayIn3 did not succeed");
			  }
              G4int h1PDG =h1QPDG.GetPDGCode();  // PDG code of the first hadron            
              G4int h2PDG =h2QPDG.GetPDGCode();  // PDG code of the second hadron
              G4int envPDG=theEnvironment.GetPDG();// PDGCode of the NuclQEnvironment
              G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
              theQHadrons.insert(h1H);
              G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
              theQHadrons.insert(h2H);
              G4QHadron* qeH = new G4QHadron(theEnvironment.GetPDG(),e4M);
              theQHadrons.insert(qeH);
			}
		  }
          else                                   // => "Two particles decay" case
		  {
            G4double QM = pQ->Get4Momentum().m();// Real Mass of the Quasmon
            G4double GSM= QQPDG.GetMass();       // GS Mass of the Quasmon
            if(QM<GSM||GSM+envM>totMass)
			{
		      G4cerr<<"***G4QEnv::HadrQE: QM="<<QM<<" < QGSM="<<GSM<<G4endl;
              G4Exception("G4QEnvironment::HadronizeQEnvironment: (3) Cann't decay");
			}
            else                                 // => "Possible decay" case
			{
              G4LorentzVector fq4M(0.,0.,0.,GSM);
              G4LorentzVector qe4M(0.,0.,0.,envM);
              if(!G4QHadron(tot4M).DecayIn2(fq4M,qe4M))
			  {
                 G4cerr<<"***G4QEnv::HadQEnv: tM="<<tot4M.m()<<"-> QPDG="<<QPDG<<"(M="
					   <<GSM<<") + envM="<<envM<<")"<<G4endl;
				 G4Exception("G4QEnv::HadrQEnv: Quasmon+QEnv DecayIn2 did not succeed");
			  }
              G4QHadron* qH = new G4QHadron(QPDG,fq4M);
              theQHadrons.insert(qH);
              G4QHadron* qeH = new G4QHadron(theEnvironment.GetPDG(),qe4M);
              theQHadrons.insert(qeH);
			}
		  }
		}
        else                                     // "Last decay was fatal" case
		{
#ifdef pdebug
		  G4cout<<"***G4QEnv::HadrQE: M="<<totMass<<",dM="<<dM<<",nQ="<<nQuasmons<<G4endl;
#endif
          G4int          nOfOUT  = theQHadrons.entries();
          while(nOfOUT)
          {
            G4QHadron*     theLast = theQHadrons[nOfOUT-1];
            G4LorentzVector last4M = theLast->Get4Momentum();
            G4QContent      lastQC = theLast->GetQC();
            G4int           lastS  = lastQC.GetStrangeness();
            G4int           totS   = totQC.GetStrangeness();
			if(lastS<0&&lastS+totS<0&&nOfOUT>1) // => "Skip K-meson" case @@ Skip a few @@
			{
              G4QHadron* thePrev = theQHadrons[nOfOUT-2];
              theQHadrons.removeLast();         // the last QHadron is excluded from OUTPUT
              theQHadrons.removeLast();         // the prev QHadron is excluded from OUTPUT
              theQHadrons.insert(theLast);      // the Last becomes the prev
              theLast = thePrev;                // Update parameters (Prev instead of Last)
              last4M = theLast->Get4Momentum();
              lastQC = theLast->GetQC();
			}
            else theQHadrons.removeLast();      // the last QHadron is excluded from OUTPUT 
            delete          theLast;            // the last QHadron is deleated as instance
            totQC+=lastQC;                      // Update (increase) the total QC
            tot4M+=last4M;                      // Update (increase) the total 4-momentum
            totMass=tot4M.m();                  // Calculate new real total mass
            G4int bn=totQC.GetBaryonNumber();   // The BaryNum after addition
            totPDG=totQC.GetSPDGCode();
            if(bn>1)
			{
              totS  =totQC.GetStrangeness();    // Total Strangeness of this System
              if(totS>=0)                       // => "This is a normal nucleus" case
			  {
                G4QNucleus newN(totQC,tot4M);
                totPDG=newN.GetPDG();
                totM  =newN.GetMZNS();          // Calculate new minimum (GS) mass
			  }
              else if(totS==-1)                 // => "Try to decay in K+/aK0 and finish" case
			  {
                G4double m1=mK;         
                G4int  PDG1=321;
                G4QNucleus  newNp(totQC-KpQC);
                G4int  PDG2=newNp.GetPDG();
                G4double m2=newNp.GetMZNS();
                G4QNucleus  newN0(totQC-KpQC);
                G4double m3=newN0.GetMZNS();
                if (m3+mK0<m2+mK)                // => "aK0+ResA is better" case
	            {
                  m1  =mK0;
                  PDG1=311;
                  m2  =m3;
                  PDG2=newN0.GetPDG();
	            }
                if(totMass>m1+m2)               // => "can decay" case
                {
                  G4LorentzVector fq4M(0.,0.,0.,m1);
                  G4LorentzVector qe4M(0.,0.,0.,m2);
                  if(!G4QHadron(tot4M).DecayIn2(fq4M,qe4M))
			      {
                    G4cerr<<"***G4QEnv::HadQEnv: tM="<<tot4M.m()<<"-> aK="<<PDG1<<"(M="
					      <<m1<<") + ResA="<<PDG2<<"(M="<<m2<<")"<<G4endl;
				    G4Exception("G4QEnv::HadrQEnv: aK+ResA DecayIn2 did not succeed");
			      }
                  G4QHadron* H1 = new G4QHadron(PDG1,fq4M);
                  theQHadrons.insert(H1);
                  G4QHadron* H2 = new G4QHadron(PDG2,qe4M);
                  theQHadrons.insert(H2);
                  break;
			    }
                else totM=100000.;              // => "Continue reversion" case
			  }
              else totM=200000.;                // => "Continue reversion" case
			}
            else
			{
              if     (totPDG==1114||totPDG==2224||totPDG==10) // Decay right now and finish
			  {
                G4double m1=mNeut;
                G4int  PDG1=2112;
                G4double m2=mPi;
                G4int  PDG2=-211;
				if(totPDG==2224)
				{
                  m1=mProt;
                  PDG1=2212;
                  m2=mPi;
                  PDG2=211;
				}
                else if(totPDG==10)              // "Chipolino" case
				{
                  G4QChipolino resChip(totQC);   // define the residual as a Chipolino
                  G4QPDGCode h1=resChip.GetQPDG1();
                  PDG1=h1.GetPDGCode();          // PDG code of the first hadron
                  m1  =h1.GetMass();             // Mass of the first hadron
                  G4QPDGCode h2=resChip.GetQPDG2();
                  PDG2=h2.GetPDGCode();          // PDG code of the second hadron
                  m2  =h2.GetMass();             // Mass of the second hadron
				}
                if(totMass>m1+m2)
				{
                  G4LorentzVector fq4M(0.,0.,0.,m1);
                  G4LorentzVector qe4M(0.,0.,0.,m2);
                  if(!G4QHadron(tot4M).DecayIn2(fq4M,qe4M))
			      {
                    G4cerr<<"***G4QEnv::HadQEnv: tM="<<tot4M.m()<<"-> h1="<<PDG1<<"(M="
					      <<m1<<") + h2="<<PDG2<<"(M="<<m2<<")"<<G4endl;
				    G4Exception("G4QEnv::HadrQEnv: h1+h2 DecayIn2 did not succeed");
			      }
                  G4QHadron* H1 = new G4QHadron(PDG1,fq4M);
                  theQHadrons.insert(H1);
                  G4QHadron* H2 = new G4QHadron(PDG2,qe4M);
                  theQHadrons.insert(H2);
                  break;
				}
                else totM=300000.;
			  }
			  else if(totPDG) totM=G4QPDGCode(totPDG).GetMass();
              else totM=400000.;
			}
            G4double dM=totMass-totM;
#ifdef pdebug
		    G4cout<<"G4QEnv::HadrQE: Add H="<<last4M<<lastQC<<",tM="<<totM<<",dM="<<dM<<G4endl;
#endif
            if(dM>-0.001&&totPDG)
		    {
              G4QHadron* evH = new G4QHadron(totPDG,tot4M);// Create a Hadron for ResNucl
              if(dM>0.001)EvaporateResidual(evH);// Try to evaporate residual
              else      theQHadrons.insert(evH);// Fill to OUTPUT as it is
              break;
		    }
            nOfOUT  = theQHadrons.entries();    // Update the value of OUTPUT entries
		  }
		  if(!nOfOUT)
		  {
		    G4cerr<<"***G4QEnv::HadrQE: M="<<totMass<<" < gsM="<<totM<<", dM="<<dM<<G4endl;
            G4Exception("G4QEnvironment::HadronizeQEnvironment: (4) Cann't decay QEnv");
	      }
		}
        CleanUp();
        return theQHadrons;
	  }
	}
  }
} // End of the main member function HadronizeQEnvironment

// Clean up the QEnvironment to Zero
void G4QEnvironment::CleanUp()
//   =========================
{
  static const G4QNucleus vacuum(90000000);
  theEnvironment=vacuum;
  G4int nQuasmons = theQuasmons.entries();
  if (nQuasmons) for (G4int iq=0; iq<nQuasmons; iq++)theQuasmons[iq]->KillQuasmon();
}

// Calculate a#of clusters in the nucleus
void G4QEnvironment::PrepareClusters()
//   =================================
{
  static const G4int NUCPDG=90000000;
  G4double ze  = theEnvironment.GetZ();          // For bn<3 (in all nucleus)
  G4double ne  = theEnvironment.GetN();
  G4double se  = theEnvironment.GetS();
  G4double ae  = ze + ne + se;
  G4double ze1 = theEnvironment.GetDZ() + 1;     // For bn>2 (nly in a dense region)
  G4double ne1 = theEnvironment.GetDN() + 1;
  G4double se1 = theEnvironment.GetDS() + 1;
  G4double pos = theEnvironment.GetProbability();// Vacuum value
  for (G4int index=0; index<theQCandidates.entries(); index++)
  {
    G4QCandidate* curCand=theQCandidates[index];
    G4int cPDG  = curCand->GetPDGCode();
    if(cPDG>NUCPDG)                              // ===> Cluster case
	{
      G4int zns= cPDG-NUCPDG;
      G4int nc = zns%1000;                       // N of the cluster
      G4int sz = zns/1000;
      G4int zc = sz %1000;                       // Z of the cluster
      G4int sc = sz /1000;                       // S of the cluster
      G4int ac = zc+nc+sc;                       // A of the cluster
      pos      = theEnvironment.GetProbability(ac);// Get a cluster prob. normalization factor
      G4double dense=1.;
      if(ac==1)dense=theEnvironment.GetProbability(254)/pos;
      if(ac==2)dense=theEnvironment.GetProbability(255)/pos;
#ifdef sdebug
	  G4cout<<"G4QEnvironment::PrepareClusters: cPDG="<<cPDG<<",norm="<<pos<<",zc="<<zc
            <<",nc="<<nc<<",sc="<<sc<<",ze1="<<ze1<<",ne1="<<ne1<<",se1="<<se1<<G4endl;
#endif
      if     (ac==1)
	  {
        if(zc>0) pos*=ze;
        if(nc>0) pos*=ne;
        if(sc>0) pos*=se;
      }
      else if(ac==2)
	  {
        if(ze<zc||ne<nc||se<sc) pos=0.;
        else
		{
          if     (zc==2) pos*=ze*(ze-1)/2;
          else if(nc==2) pos*=ne*(ne-1)/2;
          else if(sc==2) pos*=se*(se-1)/2;
          else if(zc==1&&nc==1) pos*=ze*ne;
          else if(zc==1&&sc==1) pos*=ze*se;
          else if(sc==1&&nc==1) pos*=se*ne;
          else G4cout<<"***G4QEnvironment::PrepClust:z="<<zc<<",n="<<nc<<",s="<<sc<<G4endl;
		}
      }
      else
	  {
        if(ze1<=zc||ne1<=nc||se1<=sc) pos=0.;
        else
		{
          if(zc>0) for(int iz=1; iz<=zc; iz++) pos*=(ze1-iz)/iz;
          if(nc>0) for(int in=1; in<=nc; in++) pos*=(ne1-in)/in;
          if(sc>0) for(int is=1; is<=sc; is++) pos*=(se1-is)/is;
		}
	  }
      curCand->SetPreProbability(pos);
      curCand->SetDenseProbability(pos*dense);
#ifdef sdebug
	  G4cout<<"G4QEnv::PrepClust:ClastPDG="<<cPDG<<",preProb="<<pos<<",d="<<dense<<G4endl;
#endif
	}
	else
    {
      curCand->SetPreProbability(pos);           // ===> Hadronic case in Vacuum     
      curCand->SetDenseProbability(0.);          // ===> Hadronic case in Vacuum
    }
	curCand->SetPossibility(true);               // All candidates are possible at this point
  }
}// End of PrepareClusters

//Evaporate Residual Nucleus
void G4QEnvironment::EvaporateResidual(G4QHadron* qH)
{//  ==================================================
  static const G4double mAlpha = G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mDeutr = G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mNeut  = G4QPDGCode(2112).GetMass();
  static const G4double mProt  = G4QPDGCode(2212).GetMass();
  static const G4double mLamb  = G4QPDGCode(3122).GetMass();
  static const G4double mPi    = G4QPDGCode(211).GetMass();
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4int NUCPDG=90000000;
  G4int        thePDG = qH->GetPDGCode();      // Get PDG code of the Hadron
  G4LorentzVector q4M = qH->Get4Momentum();    // Get 4-momentum of the Hadron
  G4double    totMass = qH->GetMass();         // Real Mass of the nuclear fragment
#ifdef pdebug
  G4cout<<"G4QEnv::EvaporateResidual:Hadron's PDG="<<thePDG<<",4Mom="<<q4M<<G4endl;
#endif
  if     (thePDG==NUCPDG)                      // ==> "Nothing in the INPUT Hadron" case
  {
    delete qH;
    return;
  }
  else if(thePDG>NUCPDG&&thePDG!=90002999)     // ==> "Decay-Evaporation" case
  {
    G4QNucleus qNuc(q4M,thePDG);               // Make a Nucleus out of the Hadron
    G4double GSMass =qNuc.GetGSMass();         // GrState Mass of the nuclear fragment
    G4int    bA     =qNuc.GetA();              // A#of baryons in the Nucleus
#ifdef pdebug
	G4cout<<"G4QEnv::EvapResid:"<<qNuc<<",PDG="<<thePDG<<",M="<<q4M.m()<<",GSM="<<GSMass<<G4endl;
#endif
    if(abs(totMass-GSMass)<.1)                 // ==> Case of possible decay
    {
      G4QContent totQC=qNuc.GetQCZNS();        // Total Quark Content of Residual Nucleus
      G4int    nN     =qNuc.GetN();            // A#of neutrons in the Nucleus
      G4int    nZ     =qNuc.GetZ();            // A#of protons in the Nucleus
      G4int    nS     =qNuc.GetS();            // A#of lambdas in the Nucleus
      G4double nResM  =1000000.;               // Prototype of mass of residual for a neutron
      G4int    nResPDG=0;                      // Prototype of PDGCode of residual for a neutron
      if(nN>0&&bA>1)                           // It's nucleus and there is a neutron
	  {
        G4QContent resQC=totQC-neutQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        nResPDG=resN.GetPDG();                 // PDG of the Residual Nucleus
        if     (nResPDG==90000001) nResM=mNeut;
        else if(nResPDG==90001000) nResM=mProt;
        else if(nResPDG==91000000) nResM=mLamb;
        else nResM=resN.GetMZNS();             // min mass of the Residual Nucleus
	  }
      G4double pResM  =1000000.;               // Prototype of mass of residual for a proton
      G4int    pResPDG=0;                      // Prototype of PDGCode of residual for a proton
      if(nZ>0&&bA>1)                           // It's nucleus and there is a proton
	  {
        G4QContent resQC=totQC-protQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        pResPDG=resN.GetPDG();                 // PDG of the Residual Nucleus
        if     (pResPDG==90000001) pResM=mNeut;
        else if(pResPDG==90001000) pResM=mProt;
        else if(pResPDG==91000000) pResM=mLamb;
        else pResM  =resN.GetMZNS();           // min mass of the Residual Nucleus
	  }
      G4double lResM  =1000000.;               // Prototype of mass of residual for a Lambda
      G4int    lResPDG=0;                      // Prototype of PDGCode of residual for a Lambda
      if(nS>0&&bA>1)                           // It's nucleus and there is a Lambda
	  {
        G4QContent resQC=totQC-lambQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        lResPDG=resN.GetPDG();                 // PDG of the Residual Nucleus
        if     (lResPDG==90000001) lResM=mNeut;
        else if(lResPDG==90001000) lResM=mProt;
        else if(lResPDG==91000000) lResM=mLamb;
        else lResM  =resN.GetMZNS();           // min mass of the Residual Nucleus
	  }
#ifdef ppdebug
      G4cout<<"G4QEnv::EvapResid: rP="<<pResPDG<<",rN="<<nResPDG<<",rL="<<lResPDG<<",nN="
            <<nN<<",nZ="<<nZ<<",nL="<<nS<<",totM="<<totMass<<",n="<<totMass-nResM-mNeut
            <<",p="<<totMass-pResM-mProt<<",l="<<totMass-lResM-mLamb<<G4endl;
#endif
      if(thePDG==90004004||bA>1&&(nN>0&&totMass>nResM+mNeut||nZ>0&&totMass>pResM+mProt
                                                           ||nS>0&&totMass>lResM+mLamb))
	  {
        G4int barPDG = 90002002;
        G4int resPDG = 90002002;
        G4double barM= mAlpha;
        G4double resM= mAlpha;
		if     (totMass>nResM+mNeut)           // Can radiate a neutron (priority 1)
		{
          barPDG=90000001;
          resPDG=nResPDG;
          barM  =mNeut;
          resM  =nResM;
		}
		else if(totMass>pResM+mProt)           // Can radiate a proton (priority 2)
		{
          barPDG=90001000;
          resPDG=pResPDG;
          barM  =mProt;
          resM  =pResM;
		}
		else if(totMass>lResM+mLamb)           // Can radiate a Lambda (priority 3)
		{
          barPDG=91000000;
          resPDG=lResPDG;
          barM  =mLamb;
          resM  =lResM;
		}
        else if(thePDG!=90004004&&totMass>GSMass)// If it's not Be8 decay in gamma
		{
          barPDG=22;
          resPDG=thePDG;
          barM  =0.;
          resM  =pResM;
		}
        else if(thePDG!=90004004)
		{
          G4cerr<<"***G4QEnv::EvapResid:PDG="<<thePDG<<",M="<<totMass<<"< GSM="<<GSMass<<G4endl;
          G4Exception("***G4QEnvironment::EvaporateResidual: M<GSM & cann't decay in p,n,l");
		}
        G4LorentzVector a4Mom(0.,0.,0.,barM);
        G4LorentzVector b4Mom(0.,0.,0.,resM);
        if(!qH->DecayIn2(a4Mom,b4Mom))
        {
          theQHadrons.insert(qH);              // @@ No decay
          G4cout<<"G4QEnv::EvapResid: rP="<<pResPDG<<",rN="<<nResPDG<<",rL="<<lResPDG<<",nN="
                <<nN<<",nZ="<<nZ<<",nL="<<nS<<",totM="<<totMass<<",n="<<totMass-nResM-mNeut
                <<",p="<<totMass-pResM-mProt<<",l="<<totMass-lResM-mLamb<<G4endl;
          G4cout<<"***G4QEnv::EvapResid: Decay failed bPDG="<<barPDG<<", rPDG="<<resPDG<<G4endl;
	    }
        else
        {
          qH->SetNFragments(2);                // Fill a#of fragments to decaying Hadron
          theQHadrons.insert(qH);              // Fill hadron in the HadronVector with nf=2
          G4QHadron* HadrB = new G4QHadron(barPDG,a4Mom);
          theQHadrons.insert(HadrB);           // Fill the baryon to Output Hadr. Vector
          G4QHadron* HadrR = new G4QHadron(resPDG,b4Mom);
          if(HadrR->GetBaryonNumber()>1) EvaporateResidual(HadrR); // Continue decay of residNucl
          else theQHadrons.insert(HadrR);      // Fill ResidNucleus=Baryon to Output HadronVector
        }
	  }
      else theQHadrons.insert(qH);             // No decay (Alpha decays are not implemented)
    }
    else if (totMass<GSMass)
	{
      G4cerr<<"***G4QEnv::EvapRes:M="<<totMass<<" < GSM="<<GSMass<<", d="<<totMass-GSMass<<G4endl;
      G4Exception("***G4QEnv::EvaporateResidual: Nuc.Mass is below the Ground State value");
	}
    else                                       // ===> Evaporation of excited system
	{
#ifdef pdebug
      G4cout<<"G4QEnv::EvapResid:Evaporate "<<thePDG<<",tM="<<totMass<<" > GS="<<GSMass
          <<qNuc.Get4Momentum()<<", m="<<qNuc.Get4Momentum().m()<<G4endl;
#endif
      G4QHadron* bHadron = new G4QHadron;
      G4QHadron* rHadron = new G4QHadron;
      if(!qNuc.EvaporateBaryon(bHadron,rHadron))
	  {
        G4cerr<<"***G4QEnv::EvapResid:Evaporation, PDG="<<thePDG<<",tM="<<totMass<<G4endl;
        delete bHadron;
        delete rHadron;
        theQHadrons.insert(qH);                // Fill hadron in the HadronVector as it is
	  }
      G4int bPDG=bHadron->GetPDGCode();
      G4int rPDG=rHadron->GetPDGCode();
#ifdef pdebug
      G4cout<<"G4QEnv::EvapResid:Done evapFragmPDG="<<bPDG<<", resNuclPDG="<<rPDG<<G4endl;
#endif
      qH->SetNFragments(2);                    // Fill a#of fragments to decaying Hadron
      theQHadrons.insert(qH);                  // Fill hadron in the HadronVector with nf=2
      G4int fB=bHadron->GetBaryonNumber();     // Baryon number of evaporated fragment
      if(fB<2)theQHadrons.insert(bHadron);     // Fill Evapor. Baryon to Output HadronVector
      else if(fB==2)                           // "Dibaryon" case needs decay
	  {
        G4double fGSM = bHadron->GetQPDG().GetMass(); // Ground State mass of the dibaryon
        G4double fM   = bHadron->GetMass();   // Real mass of the dibaryon
        if(fM<fGSM-0.003)
		{
          G4cerr<<"***G4QEnv::EvapRes: M="<<fM<<" < GSM="<<fGSM<<", d="<<fGSM-fM<<G4endl;
          G4Exception("***G4QEnv::EvaporateResidual: Evaporation below mass shell");
		}
        else if(abs(fM-fGSM)<=0.001&&bPDG==90001001) theQHadrons.insert(bHadron);
        else DecayDibaryon(bHadron);           // "Decay" case
	  }
      else
	  {
        G4cerr<<"***G4QEnv::EvapRes:fB="<<fB<<" > 2 - unexpected evaporated fragment"<<G4endl;
        G4Exception("***G4QEnv::EvaporateResidual: Unexpected evaporation act");
	  }
      G4int rB=rHadron->GetBaryonNumber();     // Baryon number of the residual nucleus
      if(rB>2) EvaporateResidual(rHadron);     // Continue evaporation
      else if(rB==2)                           // "Dibaryon" case needs decay
	  {
        G4double fGSM = rHadron->GetQPDG().GetMass(); // Ground State mass of the dibaryon
        G4double fM   = rHadron->GetMass();    // Real mass of the dibaryon
        if(fM<fGSM-0.001)
		{
          G4cerr<<"***G4QEnv::EvapRes: <residual> M="<<fM<<" < GSM="<<fGSM<<G4endl;
          G4Exception("***G4QEnv::EvaporateResidual: Evaporation below mass shell");
		}
        else if(abs(fM-fGSM)<=0.001&&bPDG==90001001)theQHadrons.insert(rHadron);
        else DecayDibaryon(rHadron);           // "Decay" case
	  }
      else theQHadrons.insert(rHadron);        // Fill ResidNucleus=Baryon to Output HadronVector
	}
  }
  else
  {
    G4cout<<"***G4QEnv::EvapResid: inputHadron="<<qH<<", PDG="<<thePDG<<G4endl;
    if(thePDG)
    {
      if(thePDG==10)                           // "Chipolino decay" case 
	  {
        G4QContent totQC = qH->GetMass();      // Quark content of the hadron
        G4QChipolino resChip(totQC);           // define the Residual as a Chipolino
        G4QPDGCode h1=resChip.GetQPDG1();
        G4double m1  =h1.GetMass();            // Mass of the first hadron
        G4QPDGCode h2=resChip.GetQPDG2();
        G4double m2  =h2.GetMass();            // Mass of the second hadron
        if(totMass>=m1+m2)
        {
          delete qH;                           // Chipolino should not be in a sequence
          G4LorentzVector fq4M(0.,0.,0.,m1);
          G4LorentzVector qe4M(0.,0.,0.,m2);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
		  {
            G4cerr<<"***G4QEnv::HadQEnv: tM="<<totMass<<"-> h1M="<<m1<<" + h2M="<<m2<<G4endl;
		    G4Exception("G4QEnv::HadrQEnv: Chip->h1+h2 DecayIn2 did not succeed");
	      }
          G4QHadron* H1 = new G4QHadron(h1.GetPDGCode(),fq4M);
          theQHadrons.insert(H1);
          G4QHadron* H2 = new G4QHadron(h1.GetPDGCode(),qe4M);
          theQHadrons.insert(H2);
		}
        else
	    {
          G4cerr<<"***G4QEnvironment::EvaporateResidual: M="<<totMass<<" < m1="<<m1
                <<" + m2="<<m2<<", d="<<m1+m2-totMass<<G4endl;
          G4Exception("***G4QEnv::EvaporateResidual: Chipolino is under a Mass Shell");
	    }
	  }
      else if(thePDG==89999003||thePDG==90002999)//=> "ISO-dibarion" case
	  {
        G4int  nucPDG = 2112;
        G4double nucM = mNeut;
        G4int   piPDG = -211;
        if(thePDG==90002999)
		{
          nucPDG = 2212;
          nucM   = mProt;
          piPDG  = 211;
		}
        if(totMass>mPi+nucM+nucM)
		{
          G4LorentzVector n14M(0.,0.,0.,nucM);
          G4LorentzVector n24M(0.,0.,0.,nucM);
          G4LorentzVector pi4M(0.,0.,0.,mPi);
          if(!G4QHadron(q4M).DecayIn3(n14M,n24M,pi4M))
		  {
            G4cerr<<"***G4QEnv::HadQEnv: tM="<<totMass<<"-> 2N="<<nucPDG<<"(M="
			      <<nucM<<") + pi="<<piPDG<<"(M="<<mPi<<")"<<G4endl;
		    G4Exception("G4QEnv::HadrQEnv:ISO-dibaryon DecayIn3 did not succeed");
		  }
          qH->SetNFragments(3);              // Put a#of Fragments=3 for input Hadron
          theQHadrons.insert(qH);            // Fill the residual Hadron to output
          G4QHadron* h1H = new G4QHadron(nucPDG,n14M);
          theQHadrons.insert(h1H);
          G4QHadron* h2H = new G4QHadron(nucPDG,n24M);
          theQHadrons.insert(h2H);
          G4QHadron* piH = new G4QHadron(piPDG,pi4M);
          theQHadrons.insert(piH);
		}
		else
	    {
          G4cerr<<"***G4QEnvironment::EvaporateResidual: IdPDG="<<thePDG<<", q4M="<<q4M
                <<", M="<<totMass<<" < M_2N+Pi, d="<<totMass-2*nucM-mPi<<G4endl;
          G4Exception("***G4QEnv::EvaporateResidual: ISO-dibaryon is under a Mass Shell");
	    }
	  }
      else                                     // "Hadron" case
	  {
        G4double totM=G4QPDGCode(thePDG).GetMass();
        if(abs(totMass-totM)<0.001||abs(thePDG)-10*(abs(thePDG)/10)>2)theQHadrons.insert(qH);
        else if (totMass>totM)                 // "Radiative Hadron decay" case
	    {
          qH->SetNFragments(2);                // Put a#of Fragments=2 for input Hadron
          theQHadrons.insert(qH);              // Fill the residual Hadron to output
          G4LorentzVector fq4M(0.,0.,0.,0.);
          G4LorentzVector qe4M(0.,0.,0.,totM);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
		  {
            G4cerr<<"***G4QEnv::HadQEnv: tM="<<totMass<<"-> h1M="<<totM<<" + gamma"<<G4endl;
		    G4Exception("G4QEnv::HadrQEnv: H*->H+gamma DecayIn2 did not succeed");
	      }
          G4QHadron* H1 = new G4QHadron(22,fq4M);
          theQHadrons.insert(H1);
          G4QHadron* H2 = new G4QHadron(thePDG,qe4M);
          theQHadrons.insert(H2);
	    }
        else
	    {
          G4cerr<<"***G4QEnvironment::EvaporateResidual: residualPDG="<<thePDG
                <<", q4M="<<q4M<<", M="<<totMass<<" < GSM="<<totM<<G4endl;
          G4Exception("***G4QEnv::EvaporateResidual: Hadron is under a Mass Shell");
	    }
	  }
	}
    else
    {
      G4cerr<<"***G4QEnvironment::EvaporateResidual: residualPDG="<<thePDG
            <<", q4M="<<q4M<<", qM="<<totMass<<G4endl;
      G4Exception("***G4QEnv::EvaporateResidual: This is not a nucleus nor a hadron");
    }
  }
} // End of EvaporateResidual

//Make Random Unit 3D-Vector
G4ThreeVector G4QEnvironment::RndmDir()
{//  ==================================
  G4double x = G4UniformRand();
  G4double y = G4UniformRand();
  G4double z = G4UniformRand();
  G4double r2= x*x+y*y+z*z;
  while(r2>1.||r2<.01)
  {
    x = G4UniformRand();
    y = G4UniformRand();
    z = G4UniformRand();
    r2=x*x+y*y+z*z;
  }
  G4double r=sqrt(r2);
  return G4ThreeVector(x/r,y/r,z/r);
} // End of RndmDir

//The public Hadronisation routine with delete responsibility of User (!)
G4QHadronVector* G4QEnvironment::Fragment()
{//              ==========================
#ifdef pdebug
  G4cout<<"G4QEnvironment::Fragment is called"<<G4endl;
#endif
  HadronizeQEnvironment();
  G4int nHadr=theQHadrons.entries();
#ifdef pdebug
  G4cout<<"G4QEnvironment::Fragment after HadronizeQEnvironment nH="<<nHadr<<G4endl;
#endif
  G4QHadronVector* theFragments = new G4QHadronVector; // Intermediate
  if(nHadr) for(G4int hadron=0; hadron<nHadr; hadron++)
  {
    G4QHadron* curHadr = new G4QHadron(theQHadrons[hadron]);
    theFragments->insert(curHadr);
  }
#ifdef pdebug
  G4cout<<"G4QEnvironment::Fragment ===OUT==="<<G4endl;
#endif
  return theFragments;
} // End of "Fragment"

//The public Quasmons duplication with delete responsibility of User (!)
G4QuasmonVector* G4QEnvironment::GetQuasmons()
{//              =============================
  G4int nQ=theQuasmons.entries();
#ifdef pdebug
  G4cout<<"G4QEnvironment::GetQuasmons is called nQ="<<nQ<<G4endl;
#endif
  G4QuasmonVector* quasmons = new G4QuasmonVector; // Intermediate
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ = new G4Quasmon(theQuasmons[iq]);
    quasmons->insert(curQ);
  }
#ifdef pdebug
  G4cout<<"G4QEnvironment::GetQuasmons ===OUT==="<<G4endl;
#endif
  return quasmons;
} // End of GetQuasmons

//Decay of the output dibayon
void G4QEnvironment::DecayDibaryon(G4QHadron* qH)
{//  ============================================
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  G4LorentzVector q4M = qH->Get4Momentum();  // Get 4-momentum of the Dibaryon
  G4int          qPDG = qH->GetPDGCode();    // PDG Code of the decayin dybaryon
  G4int          fPDG = 2212;                // Prototype for pp case
  G4int          sPDG = 2212;
  G4double       fMass= mProt;
  G4double       sMass= mProt;
  if     (qPDG==90000002)                    // "dineutron" case
  {
    fPDG = 2112;
    sPDG = 2112;
    fMass= mNeut;
    sMass= mNeut;    
  }
  else if(qPDG==90001001)                    // "exdeutron" case
  {
    if(abs(q4M.m()-mDeut)<0.001)
	{
      theQHadrons.insert(qH);
      return;
	}
    fPDG = 2112;
    sPDG = 2212;
    fMass= mProt;
    sMass= mNeut;    
  }
  else if(qPDG==91000001)                    // "Lambda-neutron" case
  {
    fPDG = 2112;
    sPDG = 3122;
    fMass= mNeut;
    sMass= mLamb;    
  }
  else if(qPDG==91001000)                    // "Lambda-proton" case
  {
    fPDG = 2212;
    sPDG = 3122;
    fMass= mProt;
    sMass= mLamb;    
  }
  else if(qPDG==92000000)                    // "diLambda" case
  {
    fPDG = 3122;
    sPDG = 3122;
    fMass= mLamb;
    sMass= mLamb;    
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,sMass);     // Mass is random since probab. time
  if(!G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
  {
    G4cerr<<"***G4QEnv::DecayDibaryon:fPDG="<<fPDG<<"(fM="<<fMass<<") + sPDG="<<sPDG
          <<"(sM="<<sMass<<")"<<" >? TotM="<<q4M.m()<<G4endl;
    G4Exception("***G4QEnv::DecayDibaryon: DecayIn2 didn't succeed for dibaryon");
  }
#ifdef pdebug
  G4cout<<"G4QEnv::DecayDibaryon: *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
        <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
  qH->SetNFragments(2);                      // Fill a#of fragments to decaying Dibaryon
  theQHadrons.insert(qH);                    // Fill hadron in the HadronVector with nf=2
  G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
  theQHadrons.insert(H1);                    // Fill "new H1"
  G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
  theQHadrons.insert(H2);                    // Fill "new H2"
} // End of DecayDibaryon
