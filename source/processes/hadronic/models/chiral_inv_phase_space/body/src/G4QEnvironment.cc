// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QEnvironment.cc,v 1.30 2001-10-04 20:00:22 hpw Exp $
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
  theEnvironment(90000000)                    // User is responsible for projHadrons(Vector)
{
  G4int nHadrons=projHadrons.size();       // A#of hadrons in the input Vector
  if(!nHadrons||targPDG==90000000)
  {
    G4cerr<<"***G4QEnvironment: a#of INPUT QHadrons="<<nHadrons<<",tPDG="<<targPDG<<G4endl;
    //G4Exception("***G4QEnvironment: There is no one projectile or vacuum target");
    if(nHadrons)                              // The environment is empti (nothing to interact with)
	{
      for(G4int ih=0; ih<nHadrons; ih++)
      {
        G4QHadron* curQH    = new G4QHadron(projHadrons[ih]);
        theQHadrons.push_back(curQH);            // (delete equivalent)
      }
	}
    else if(targPDG!=90000000)                // No hadrons, only the Environ exists
    {
      G4QHadron* curQH    = new G4QHadron(targPDG);
      theQHadrons.push_back(curQH);              // (delete equivalent)
    }
    return;
  }
  G4int    targA=G4QPDGCode(targPDG).GetBaryNum();
  G4double targM=G4QPDGCode(targPDG).GetMass();
  tot4Mom=G4LorentzVector(0.,0.,0.,targM);
  // ===== Print out of the input information at Creation time =========
#ifdef pdebug
  G4cout<<"G4QEnvironment START: targPDG="<<targPDG<<",targM="<<targM<<", nProj="<<nHadrons<<G4endl;
#endif
  for(G4int ipr=0; ipr<nHadrons; ipr++)       // This LOOP is used for the tot4Mom & for the print
  {
    G4LorentzVector h4Mom = projHadrons[ipr]->Get4Momentum();
    tot4Mom      += h4Mom;
#ifdef pdebug
    G4int           hPDG  = projHadrons[ipr]->GetPDGCode();
    G4int           hNFrag= projHadrons[ipr]->GetNFragments();
    G4QContent      hQC   = projHadrons[ipr]->GetQC();
    G4cout<<ipr<<". PDG="<<hPDG<<hQC<<", 4M="<<h4Mom<<", hNFrag="<<hNFrag<<G4endl;
#endif
  }
#ifdef pdebug
  G4cout<<"G4QEnv: nF="<<projHadrons[0]->GetNFragments()<<", tot4Mom="<<tot4Mom<<G4endl;
#endif
  G4int nP=theWorld.GetQPEntries();           // A#of init'ed particles in CHIPS World
  G4int nCl=nP-72;                            // A#of init'ed clusters in CHIPS World
#ifdef pdebug
  G4cout<<"G4QEnv: before NucClustInit: nP="<<nP<<",nF="<<projHadrons[0]->GetNFragments()<<G4endl;
#endif
  InitClustersVector(nCl,targA);              // Initialize clusters as Particles (for interaction!)
#ifdef pdebug
  G4cout<<"G4QEnv: NucClusts are done: nCl="<<nCl<<",nF="<<projHadrons[0]->GetNFragments()<<G4endl;
#endif
  if(targPDG>80000000)                        // ==> Interaction with a nuclear target (inc NUCPDG)
  {
    theEnvironment.InitByPDG(targPDG);        // Create nuclear environment
    for(G4int ih=0; ih<nHadrons; ih++)        // ==> The main LOOP over projQHadrons
    {
      G4QHadron* curHadr=projHadrons[ih];     // Pointer to current projectile Hadron
      G4int hNFrag = curHadr->GetNFragments();// #0 means intermediate (skip)
#ifdef pdebug
      G4cout<<"G4QEnv: h#"<<ih<<",nF="<<hNFrag<<",nF0="<<projHadrons[0]->GetNFragments()<<G4endl;
#endif
      if(!hNFrag)                             // => "Final hadron" case
	  {
        G4int envPDG=theEnvironment.GetPDG();
        if(envPDG==90000000)                  // ==> "Interaction with vacuum" case
        {
          G4int hPDG  = curHadr->GetPDGCode();// A PDG Code of the projQHadron
          if(!hPDG||hPDG==10)                 // Check for the validity of the QHadron (@@ 10 OK?)
          {
            G4cerr<<"***G4QEnvironment::Constructor: wrong PDG("<<ih<<")="<<hPDG
                  <<", HQC="<<curHadr->GetQC()<<", HM="<<curHadr->GetMass()<<G4endl;
            G4Exception("***G4QEnvironment::Constructor: One of input Hadrons is wrong");
          }
          else
          {
            G4int hQ = curHadr->GetQCode();  // One more check for valid of the QHadron
            if(hQ<0)
	        {
              G4cerr<<"***G4QEnvironment::Constructor: Q<0, PDG=("<<ih<<")"<<hPDG<<G4endl;
              G4Exception("***G4QEnvironment::Constructor: One of input Hadrons is wrong");
	        }
            else
            {
              G4QHadron* newHadr = new G4QHadron(curHadr);
              theQHadrons.push_back(newHadr);   // Fill existing hadron (delete equivalent)
#ifdef pdebug
              G4cout<<"G4QEnviron::CreateQuasmon: Fill h="<<hPDG<<curHadr->Get4Momentum()<<G4endl;
              for(G4int ipo=0; ipo<theQHadrons.size(); ipo++) // This LOOP is just for the print
              {
                G4int           hPDG  = theQHadrons[ipo]->GetPDGCode();
                G4LorentzVector h4Mom = theQHadrons[ipo]->Get4Momentum();
                G4int           hNFrag= theQHadrons[ipo]->GetNFragments();
                G4QContent      hQC   = theQHadrons[ipo]->GetQC();
                G4cout<<"h#"<<ipo<<":hPDG="<<hPDG<<hQC<<",h4M="<<h4Mom<<",hNFrag="<<hNFrag<<G4endl;
              }
#endif
            } // End of Q-Code check
          } // End of proper PDG for i-th Hadron
        }
        else                                 // Nuclear Environment still exists
		{
          G4LorentzVector h4Mom = curHadr->Get4Momentum();
          G4QContent      hQC   = curHadr->GetQC();
#ifdef pdebug
          G4cout<<"G4QEnv: Call CreateQuasm h4M="<<h4Mom<<",hQC="<<hQC<<", EnvPDG="<<envPDG<<G4endl;
#endif
          CreateQuasmon(hQC, h4Mom);
		} // End of Existing Nuclear Environment case
	  } // End of final hadron case
    } // End of the LOOP over input hadrons
  } // End of nuclear target case (including neutron=90000001 & proton=90001000)
  else                                      // => "Unique hadron" case
  {
    // the nuclear environment is already initialized as vacuum + get the first hadron
    G4QHadron* curHadr=projHadrons[0];     // Pointer to the first projectile Hadron (checked)
    G4int hPDG  = curHadr->GetPDGCode();   // A PDG Code of the projQHadron
    if(!hPDG||hPDG==10)                    // Check for the validity of the QHadron
    {
      G4cerr<<"***G4QEnvironment::Constructor:Vacuum, 1st Hadron wrong PDG="<<hPDG
            <<", HQC="<<curHadr->GetQC()<<", HM="<<curHadr->GetMass()<<G4endl;
      G4Exception("***G4QEnvironment::Constructor: Fiest input Hadron is wrong");
    }
    else
    {
      G4int hQ = curHadr->GetQCode();     // One more check for valid of the QHadron
      if(hQ<0)
	  {
        G4cerr<<"***G4QEnvironment::Constructor:Vacuum, Q<0, 1st HPDG="<<hPDG<<G4endl;
        G4Exception("***G4QEnvironment::Constructor: First input Hadron is wrong");
	  }
      else                                // Now we can get 4Mom &  QC of incedent particle
      {
        G4LorentzVector h4Mom = curHadr->Get4Momentum();
        G4QContent      hQC   = curHadr->GetQC();
        if(!targPDG||targPDG==10) G4cout<<"G4QEnv::CreateQ; (1) PDG="<<targPDG<<G4endl;
        G4QPDGCode      tQPDG(targPDG);
        G4int           tQ    = tQPDG.GetQCode();
        if(tQ<0||targPDG==10)
		{
          G4cerr<<"***G4QEnvironment::Constructor:Target Q<0 or Chipolino, PDG="<<targPDG<<G4endl;
          G4Exception("***G4QEnvironment::Constructor: Target is wrong");
		}
        else                              // Now we can create a unique Quasmon
		{
          h4Mom+=G4LorentzVector(0.,0.,0.,tQPDG.GetMass());//Projectile + TargetHadron
          hQC+=tQPDG.GetQuarkContent();
#ifdef pdebug
          G4cout<<"G4QEnv::CreateQ:VacuumHadronTarg Q="<<h4Mom<<hQC<<",QE="<<theEnvironment<<G4endl;
#endif
          G4Quasmon* curQuasmon = new G4Quasmon(hQC, h4Mom);
          theQuasmons.push_back(curQuasmon); // Insert Quasmon or even hadron/gamma (delete equivalent)
		}
      } // End of Q-Code check
    } // End of proper PDG for i-th Hadron
    if(nHadrons>1) for(G4int ih=0; ih<nHadrons; ih++) // fill other Hadrons to Output
	{
      G4QHadron* newHadr = new G4QHadron(curHadr);
      theQHadrons.push_back(newHadr);        // Fill existing hadron (delete equivalent)
	}
  } // End of Unique Hadron target treatment
} // End of the G4QEnvironment constructor

G4QEnvironment::G4QEnvironment(const G4QEnvironment &right)
{
  // theQHadrons (Vector)
  G4int nQH             = right.theQHadrons.size();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
    theQHadrons.push_back(curQH);           // (delete equivalent)
  }

  theWorld              = right.theWorld;
  nBarClust             = right.nBarClust;
  f2all                 = right.f2all;
  tot4Mom               = right.tot4Mom;

  // theQuasmons (Vector)
  G4int nQ              = right.theQuasmons.size();
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ     = new G4Quasmon(right.theQuasmons[iq]);
    theQuasmons.push_back(curQ);            // (delete equivalent)
  }

  // theQCandidates (Vector)
  G4int nQC             = right.theQCandidates.size();
  if(nQC) for(G4int ic=0; ic<nQC; ic++)
  {
    G4QCandidate* curQC = new G4QCandidate(right.theQCandidates[ic]);
    theQCandidates.push_back(curQC);        // (delete equivalent)
  }

  theEnvironment        = right.theEnvironment;
}

G4QEnvironment::G4QEnvironment(G4QEnvironment* right)
{
  // theQHadrons (Vector)
  G4int nQH             = right->theQHadrons.size();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right->theQHadrons[ih]);
    theQHadrons.push_back(curQH);           // (delete equivalent)
  }

  theWorld              = right->theWorld;
  nBarClust             = right->nBarClust;
  f2all                 = right->f2all;
  tot4Mom               = right->tot4Mom;

  // theQuasmons (Vector)
  G4int nQ              = right->theQuasmons.size();
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ     = new G4Quasmon(right->theQuasmons[iq]);
    theQuasmons.push_back(curQ);            // (delete equivalent)
  }

  // theQCandidates (Vector)
  G4int nQC             = right->theQCandidates.size();
  if(nQC) for(G4int ic=0; ic<nQC; ic++)
  {
    G4QCandidate* curQC = new G4QCandidate(right->theQCandidates[ic]);
    theQCandidates.push_back(curQC);        // (delete equivalent)
  }

  theEnvironment        = right->theEnvironment;
}

G4QEnvironment::~G4QEnvironment()
{
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQCandidates nC="<<theQCandidates.size()<<G4endl;
#endif
  G4std::for_each(theQCandidates.begin(), theQCandidates.end(), DeleteQCandidate());
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQuasmons nQ="<<theQuasmons.size()<<G4endl;
#endif
  G4std::for_each(theQuasmons.begin(), theQuasmons.end(), DeleteQuasmon());
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQHadrons nH="<<theQHadrons.size()<<G4endl;
#endif
  G4std::for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron());
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
  G4int nQH             = right.theQHadrons.size();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
    theQHadrons.push_back(curQH);           // (delete equivalent)
  }

  theWorld              = right.theWorld;
  nBarClust             = right.nBarClust;
  f2all                 = right.f2all;

  // theQuasmons (Vector)
  G4int nQ              = right.theQuasmons.size();
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ     = new G4Quasmon(right.theQuasmons[iq]);
    theQuasmons.push_back(curQ);            // (delete equivalent)
  }

  // theQCandidates (Vector)
  G4int nQC             = right.theQCandidates.size();
  if(nQC) for(G4int ic=0; ic<nQC; ic++)
  {
    G4QCandidate* curQC = new G4QCandidate(right.theQCandidates[ic]);
    theQCandidates.push_back(curQC);        // (delete equivalent)
  }

  theEnvironment        = right.theEnvironment;
  return *this;
}

// Member function for Quasmon Creation & Environment nucleus modification
void G4QEnvironment::CreateQuasmon(const G4QContent& projQC, const G4LorentzVector& proj4M)
{//========================================================================================
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mPi2 = mPi*mPi;
  static const G4QContent gamQC(0,0,0,0,0,0);
  static const G4QContent pimQC(1,0,0,0,1,0);
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QNucleus vacuum(90000000);
  G4QContent valQ(0,0,0,0,0,0);                  // Prototype of the Quasmon's Quark Content
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
  if(targPDG>80000000&&targPDG!=90000000)        // Interaction with a nuclear target
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
    G4int nCl =nP-72;                            // A#of initialized clusters in CHIPS World
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
	G4cout<<"G4QEnv::CreateQ:TargNuc Z="<<envZ<<",N="<<envN<<",S="<<envS<<",nC="<<nBarClust<<G4endl;
#endif
    G4int projPDG=projQC.GetSPDGCode();          // Minimum hadron for the projectile QC
    G4bool pbpt=projE<PiPrThresh+(M2ShiftVir+projM2)/DiNuclMass; // PhotonBelowPionThrethold
    G4bool din=false;
    G4bool piF=false;
    G4bool gaF=false;
    //if(abs(projM2-mPi2)<.00001&&projE-mPi<0.1&&projPDG==-211) din=true;// InCaseOf ProjPiMinAtRest
    if(abs(projM2-mPi2)<.00001&&projE-mPi<0.1&&projPDG==-211) piF=true; // InCaseOf ProjPiMinAtRest
    //if(pbpt&&projPDG==22) din=true; // InCaseOf GammaBelowPiThresh needs DiNucl (?)
    if(pbpt&&projPDG==22) gaF=true; // InCaseOf GammaBelowPiThresh needs DiNucl (?)
    theEnvironment.SetMaxClust(nBarClust);
    nBarClust=theEnvironment.UpdateClusters(din);// Clusters are calculated up to maxClust
#ifdef pdebug
	G4cout<<"G4QEnvironment::CreateQ: Nucleus("<<targPDG<<") is created ("<<nBarClust<<" clast's)";
    for(G4int ic=0;ic<nBarClust;ic++)G4cout<<" #"<<ic<<"("<<theEnvironment.GetProbability(ic)<<")";
    G4cout<<G4endl;
#endif
    theEnvironment.PrepareCandidates(theQCandidates,piF,gaF); // Calculate the cluster's population
#ifdef pdebug
	G4cout<<"G4QEnvironment::CreateQ: Cluster probab is calculated."<<G4endl;
#endif
    G4bool efFlag=false;                         // Flag of Energy Flow case FALSE(@@==DEFOLT==@@)
    //     ************ Change if necessary to compare Energy Flux & Multy Quasmon **************
    G4int efCounter=0;                           // Counter of Energy Flux particles
    G4QContent EnFlQC(0,0,0,0,0,0);              // Quark Content of Energy Flux
    G4LorentzVector ef4Mom(0.,0.,0.,0.);         // Summed 4-momentum of Energy Flux
    G4double proj3M=proj4M.rho();
    //   ---   Pbar     ---    Nbar  ---  LAMBDAbar  ---  SIGMA-bar  ---  SIGMA0bar  ---  SIGMA+bar
    if((projPDG==-2212||projPDG==-2112||projPDG==-3122||projPDG==-3112||projPDG==-3212||
        projPDG==-3222) && proj3M<10.)           // Only for AtRest interactions (move to interface)
	{
      // @@ Annihilation on only one baryon is implemented (no annihilation on clusters! @@??) @@
#ifdef pdebug
      G4cout<<"G4QEnviron::CreateQ:Annihilation on a perif. nucleon, Z="<<envZ<<",N="<<envN<<G4endl;
#endif
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
#ifdef pdebug
      G4cout<<"G4QEnvironment::CreateQ:"<<targNPDG<<" is selected Env="<<theEnvironment<<G4endl;
#endif
      G4double resMass=theEnvironment.GetGSMass(); // Nuclear mass after baryon subtraction
      G4double barMass=tgMass-resMass;             // Mass of the bound baryon for annihilation
      tgMass=resMass;                              // New mass of theEnvironment
      q4Mom=G4LorentzVector(0,0,0,barMass)+proj4M; // 4-momentum of the intermediate B-Bbar Quasmon
      valQ=targQC+projQC;                          // Quark Content of intermediate B-Bbar Quasmon
      G4Quasmon* pan = new G4Quasmon(valQ,q4Mom);  // N-Nbar Quasmon creation (deleted 3 lines below
      G4QNucleus vE = vacuum;                      // The annihilation as in vacuum (in NuclMatter?)
#ifdef pdebug
      G4cout<<"G4QEnviron::CreateQ: before Fragment vE="<<vE<<",QQC="<<valQ<<",Q4M="<<q4Mom<<G4endl;
#endif
      G4QHadronVector* output=pan->Fragment(vE,1); // Output of "inVac" Annihilation *!DESTROY!*<--+
#ifdef pdebug
	  G4cout<<"G4QEnvironment::CreateQ: Before delet pan."<<G4endl;  //                            ^
#endif
      delete pan;                                  // The N-Nbar tmp Quasmon is deleted A.S.A.P.   ^
      G4QHadronVector input;                       // Input for MultyQuasmon **!!DESTROY!!** <---+ ^
      G4int trgPDG = theEnvironment.GetPDG();      // New PDG Code for the Residual Nucleus      ^ ^
      G4LorentzVector trg4M(0.,0.,0.,resMass);     // New 4-momentum for the Residual Nucleus    ^ ^
      G4int tNH = output->size();               // For the selection LOOP                     ^ ^
      G4ThreeVector dir = RndmDir();               // For the selection in the LOOP              ^ ^
#ifdef pdebug
	  G4cout<<"G4QEnvironment::CreateQ: Loop over "<<tNH<<" hadrons."<<G4endl;  //               ^ ^
#endif
      for (G4int ind=0; ind<tNH; ind++)            // Loop over projectile QHadrons              ^ ^
      {
        G4QHadron*   curHadr = output->operator[](ind);    // Pointer to the current hadron              ^ ^
        G4int           shDFL= curHadr->GetNFragments(); // A#of decay fragments for proj.       ^ ^
        G4LorentzVector sh4m = curHadr->Get4Momentum();  // 4Mom for the projectile              ^ ^
        G4ThreeVector   shDIR= sh4m.vect().unit(); // Projection to the random direction         ^ ^
#ifdef pdebug
		G4cout<<"G4QE::CrQ:##"<<ind<<",d="<<shDFL<<",PDG="<<curHadr->GetPDGCode()<<G4endl; //    ^ ^
#endif
        if(!shDFL)                                 // Final (not decayed) QHadron (d==0)         ^ ^
		{
#ifdef pdebug
		  G4cout<<"G4QE::CrQ:efF="<<efFlag<<",d="<<dir.dot(shDIR)<<",SA="<<SolidAngle<<G4endl;// ^ ^
#endif
		  if(dir.dot(shDIR)>SolidAngle)            // Sum up these hadrons and make Energy Flow  ^ ^
		  {
            if(efFlag)                             // => Case of Energy Flux approach            ^ ^
		    {
              G4QContent shQC = curHadr->GetQC();  // Quark Content of the Current Hadron        ^ ^
              ef4Mom+=sh4m;
              EnFlQC+=shQC;
              efCounter++;
#ifdef pdebug
              G4int hPDG=curHadr->GetPDGCode();    // Only for gebug printing                    ^ ^
              G4LorentzVector h4M = curHadr->Get4Momentum();  // Only for gebug printing         ^ ^
			  G4cout<<"G4QE::CrQ:#"<<efCounter<<", PDG="<<hPDG<<", h4M="<<h4M<<G4endl; //        ^ ^
#endif
		    }
		    else                                   // => "MultyQuasmon fragmentation" (!efFlag)  ^ ^
		    {
              G4QHadron* mqHadron = new G4QHadron(curHadr);
              input.push_back(mqHadron);              // Fill hadron-copy (delete equivalent)   <...^ ^
#ifdef pdebug
              G4int hPDG=curHadr->GetPDGCode();    // Only for debug printing                    ^ ^
              G4LorentzVector h4M = curHadr->Get4Momentum(); // Only for gebug printing          ^ ^
		      G4cout<<"G4QE::CrQ: Fill #"<<ind<<", PDG="<<hPDG<<", h4M="<<h4M<<G4endl; //        ^ ^
#endif
		    }
		  }
		  else                                     // "Direct filling of the output vector" case ^ ^
		  {
#ifdef pdebug
            G4int hPDG=curHadr->GetPDGCode();      // Only for gebug printing                    ^ ^
            G4LorentzVector h4M = curHadr->Get4Momentum(); // Only for gebug printing            ^ ^
		    G4cout<<"G4QE::CrQ: Fill OUT #"<<ind<<",PDG="<<hPDG<<",h4M="<<h4M<<G4endl; //        ^ ^
#endif
            G4QHadron* curHadron = new G4QHadron(curHadr);
            theQHadrons.push_back(curHadron);         // TheQHadrons are filled by new hadr-copies  ^ ^
          }
		} // End of the LOOP over projectiles                                                    ^ ^
	  } // End of LOOP over "output" of annihilation                                             ^ ^
      G4std::for_each(output->begin(), output->end(), DeleteQHadron());// Here we are DESTROING output >---------------^
      delete output;                               // =============================================^
      if(!efFlag)                                  // => Not Energy Flux case: MultyQuasmon case ^
	  {
        if(!(input.size())) return;             // *** RETURN *** Without Quasmon creation----^
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: Creation of a fake Quasmon to restore parameters"<<G4endl;
#endif
        G4Quasmon fakeQ;                           // fake Quasmon to get and restore parameters ^
        G4double QTemper=fakeQ.GetTemper();        // Temperature defined by user for Quasmons   ^
        G4double QSOverU=fakeQ.GetSOverU();        // S/U defined by user for Quasmons           ^
        G4double QEtaSup=fakeQ.GetEtaSup();        // Eta Suppresion defined by user in Quasmons ^
        G4Quasmon::SetParameters(180.,.1,.3);      //  Fixed Parameters for N-barN Annihilation  ^
        G4QEnvironment* muq = new G4QEnvironment(input,theEnvironment.GetPDG());   //------+
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: before input.clearAndDestroy()"<<G4endl; //      ^     ^
#endif
        G4std::for_each(input.begin(), input.end(), DeleteQHadron());// Here we are DESTROING input >--------^-----^
        theEnvironment =muq->GetEnvironment();     // Get residual Environment after interaction
        G4QuasmonVector* outQ = muq->GetQuasmons();// Copy of quasmons **!!DESTROY!!** <---^-----+
        G4QHadronVector* outH = muq->GetQHadrons();// Copy of hadrons **!!DESTROY!!** <----^--+  ^
        G4int noh = outH->size();               // a#oh hadrons in TmpEnviron UpToNow   ^  ^  ^
        if(noh) for(G4int kh=0; kh<noh; kh++)      // One can escape it but...             ^  ^  ^
        {
          G4QHadron* curH = new G4QHadron(outH->operator[](kh)); // Make a copy to destroy tmp     ^  ^  ^
          theQHadrons.push_back(curH);                // Fill new Hadrons in theQHadrons      ^  ^  ^
		}
        G4std::for_each(outH->begin(), outH->end(), DeleteQHadron());// >------------------------------------^==^  ^
        delete outH;                               // >====================================^==^  ^
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: before delete muq"<<G4endl; //                   ^     ^
#endif
        delete muq;                                //======================================^
        G4Quasmon::SetParameters(QTemper,QSOverU,QEtaSup); // Recover parameters for Quasmons    ^
	    G4int nMQ = outQ->size();               // A#of Quasmons in MultyQuasmon output       ^
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: after GetQuasmon nMQ="<<nMQ<<G4endl; //                ^
#endif
        if(nMQ) for(G4int mh=0; mh<nMQ; mh++)      // One can escape creation/distruction but... ^
        {
          G4Quasmon* curQ = new G4Quasmon(outQ->operator[](mh)); // Make a copy to destroy temporary (?) ^
          theQuasmons.push_back(curQ);                // Fill new Quasmon-copies in theQuasmons     ^
	}
        G4std::for_each(outQ->begin(), outQ->end(), DeleteQuasmon());                   // >------------------------------------------^
        delete outQ;                               // >==========================================^
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
    G4double EnFlP=ef4Mom.rho();
    PrepareInteractionProbabilities(EnFlQC,EnFlP); // Interaction probabilities for clusters
#ifdef pdebug
	G4cout<<"G4QEnvironment::CreateQ: Interaction Probabilities are calculated"<<G4endl;
#endif
    G4int nCandid = theQCandidates.size();
    G4double maxP = theQCandidates[nCandid-1]->GetIntegProbability();
    G4int i=0;
    G4QContent    curQC;                           // Quark Content of the selected cluster
    if(nCandid<=0)
	{
	  G4cerr<<"***G4QEnv::CreateQ: nC="<<nCandid<<", maxP="<<maxP<<", QE="<<theEnvironment<<G4endl;
      G4Exception("***G4QEnvironment::CreateQ: Can not select a cluster");
	}
    else if(nCandid==1||maxP==0.)
	{
      curQC=theEnvironment.GetQCZNS();
      theEnvironment=vacuum;
    }
    else
    {
      G4double totP = maxP * G4UniformRand();
#ifdef pdebug
	  G4cout<<"G4QEnvironment::CreateQ: nC="<<nCandid<<", maxP="<<maxP<<", totP="<<totP<<G4endl;
#endif
      while(theQCandidates[i]->GetIntegProbability()<totP) i++;
      G4QCandidate* curCand = theQCandidates[i];     // Pointer to selected cluster to interact
      curQC   = curCand->GetQC();                    // Get Quark Content of the selected cluster
      G4QNucleus targClust(curQC.GetP(),curQC.GetN(),curQC.GetL());// Define Cluster as a QNucleus
#ifdef pdebug
	  G4cout<<"G4QEnv::CQ:Clust#"<<i<<" is selected("<<targClust<<") from "<<theEnvironment<<G4endl;
#endif
      theEnvironment.Reduce(targClust.GetPDG());     // Subtract selected cluster from Nucleus
	}
    G4double envMass=theEnvironment.GetGSMass();     // Mass of residual nuclear environment
    if(projPDG==22&&projE<PiPrThresh+(M2ShiftVir+projM2)/DiNuclMass)// Gamma+quark Interaction
    //if(2>3)                                        //@@ ***TMP*** PhotonAbsorbtion by q is closed
	{
      q4Mom=G4LorentzVector(0.,0.,0.,tgMass-envMass);// Photon interacts with BoundedCluster
      valQ=curQC;
#ifdef pdebug
      G4cout<<"G4QEnv::CreateQ:Q="<<q4Mom<<valQ<<"+vg="<<proj4M<<", QEnv="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom, proj4M); // Interaction gamma+quark inside
      theQuasmons.push_back(curQuasmon);  // Insert Quasmon without incident gamma (delete equivalent)
	}
    else if(abs(projM2-mPi2)<.00001&&projE-mPi<0.1&&projPDG==-211)
    //if(2>3)                                        //@@ ***TMP*** PionAbsorbAtRest by q is closed
	{
      q4Mom=proj4M+G4LorentzVector(0.,0.,0.,tgMass-envMass);// PION + BoundCluster
      valQ=EnFlQC+curQC;
      if(projE<mPi)G4cout<<"***INPUT ERROR***G4QE::CrQ: pi- Energy="<<projE<<" < mPi="<<mPi<<G4endl;
#ifdef pdebug
      G4cout<<"G4QEnv::CreateQ:Q="<<q4Mom<<valQ<<"+pi="<<proj4M<<", QEnv="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom, -proj4M); // Interaction gamma+quark inside
      theQuasmons.push_back(curQuasmon);  // Insert Quasmon without incident gamma (delete equivalent)
	}
    else
	{
      q4Mom=proj4M+G4LorentzVector(0.,0.,0.,tgMass-envMass);//Projectile + BoundCluster
      valQ=EnFlQC+curQC;
#ifdef pdebug
      G4cout<<"G4QEnv::CreateQ: Q="<<q4Mom<<valQ<<", QEnv="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom);
      theQuasmons.push_back(curQuasmon); // Insert Quasmon including hadron/gamma (delete equivalent)
	}
  }
  else
  {
    G4cerr<<"***G4QEnvironment::CreateQuasmon: Strange targPDG="<<targPDG<<G4endl;
    G4Exception("***G4QEnvironment::HadrQEnvironment: Impossible environment");
  }
}

// Calculate a probability to interact with clusters for the givven PDG of the projectile
void G4QEnvironment::PrepareInteractionProbabilities(const G4QContent& projQC, G4double AP)
//   ===========================================================================Proj.3Mom.=
{
  G4double sum    = 0.;                             // Sum of probabilities of interaction
  G4double probab = 0.;                             // Interaction probability
  G4double denseB = 0.;                             // A#of*prob baryons in dense part
  G4double allB   = 0.;                             // A#of*prob baryons in the nucleus
  G4int pPDG      = projQC.GetSPDGCode();           // PDG code of the projectile particle
  for (G4int index=0; index<theQCandidates.size(); index++)
  {
    G4QCandidate* curCand=theQCandidates[index];    // Intermediate pointer
    G4int cPDG  = curCand->GetPDGCode();
    if(cPDG>80000000&&cPDG!=90000000)               // ===> Cluster case
	{
      G4QNucleus cN(cPDG);
      G4int zc = cN.GetZ();                         // "Z" of the cluster
      G4int nc = cN.GetN();                         // "N" of the cluster
      G4int sc = cN.GetS();                         // "S" of the cluster
      G4int ac = cN.GetA();                         // "A" of the cluster
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
      G4double fact=1./pow(2.,d);
      if (qC<-1) probab=0.;     
      //else if((pPDG==-211&&AP<10.||pPDG==22&&AP<150.)&&ac<2) probab=0.; //PiCapAtRest/GamUndrPi(D)
      //else if(pPDG==22&&AP<152.&&ac<2) probab=nOfCl*ac*fact*.5; //GamUnderPi (only quark capture)
      ///////////////////////////////else if((pPDG==-211&&AP<10.)&&ac<2) probab=0;//PiCapAtRest(D)
      //else if(pPDG==-211&&AP<10.)               probab=nOfCl*fact;     // special PiCaptureAtRest
      //else if(pPDG==-211&&AP<10.)               probab=nOfCl*ac*(ac-1)*fact;
      else                                      probab=nOfCl*ac*fact;
      //else                                      probab=dOfCl*ac*fact;
      //if(ac>1) probab=0.;                       // Suppress clusters
      //if(ac>2) probab=0.;                       // Suppress heavy clusters
#ifdef sdebug
	  G4cout<<"G4QE::PrepIntProb:C="<<cPDG<<",P="<<probab<<",ac="<<ac<<",dq="<<dq<<",qC="<<qC
            <<",rPDG="<<rPDG<<",pPDG="<<pPDG<<",nCP="<<nOfCl<<",dCP="<<dOfCl<<G4endl;
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
void G4QEnvironment::InitClustersVector(G4int maxClust, G4int maxA)
//   ==============================================================
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
    G4int clustB=clustQPDG.GetBaryNum();
#ifdef sdebug
	G4cout<<"G4QEnvironment::InitClustersVector: Before insert ="<<clusterPDG<<G4endl;
#endif
	//theQCandidates.push_back(new G4QCandidate(clusterPDG)); // (delete equivalent)
	if(clustB<=maxA)theQCandidates.push_back(new G4QCandidate(clusterPDG)); // (delete equivalent)
#ifdef sdebug
    G4cout<<"G4QEnvironment::InitClustersVector: Cluster # "<<i<<" with code = "
          <<clusterPDG<<", QC="<<clustQPDG.GetQuarkContent()<<G4endl;
#endif
  }
} // End of InitClastersVector

// Fragmentation of the QEnvironment with MultyQuasmon (the main member function)
G4QHadronVector G4QEnvironment::HadronizeQEnvironment()
//              ================***********************
{
  static const G4int  NUCPDG = 90000000;
  static const G4QNucleus vacuum(NUCPDG);
  static const G4QContent PiQC(0,1,0,1,0,0);
  static const G4QContent K0QC(1,0,0,0,0,1);
  static const G4QContent KpQC(0,1,0,0,0,1);
  //static const G4QContent alQC(6,6,0,0,0,0);
  static const G4QPDGCode nQPDG(2112);
  static const G4QPDGCode pQPDG(2212);
  static const G4QPDGCode lQPDG(3122);
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  //static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0);
  G4int nQuasmons = theQuasmons.size();
#ifdef pdebug
  G4cout<<"G4QEnv::HadrQE:***> HADRONIZE Q-ENVIRONMENT="<<theEnvironment<<",nQ="<<nQuasmons<<G4endl;
#endif
  if(nQuasmons<1)                                // "No Quasmons" case -> Fill QEnviron
  {
    G4int nPDG = theEnvironment.GetPDG();        // PDG code of the residual Nucl.Environ.
#ifdef pdebug
	G4cout<<"G4QEnv::HadrQE: ***NO QUASMONS*** Env="<<nPDG<<theEnvironment.Get4Momentum()<<G4endl;
#endif
    if(nPDG==90000000)
    {
      G4std::for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron());
      return theQHadrons;
    }
    if(nPDG>80000000)
	{
      G4QHadron* rNucleus = new G4QHadron(theEnvironment); // Create a Hadron for the Environment
      theQHadrons.push_back(rNucleus);              // Fill GS - no further decay (delete equivalent) 
#ifdef pdebug
	  G4cout<<"G4QEnv::HadrQE: >>>> Fill Environment"<<G4endl;
#endif
	}
    return theQHadrons;
  }
  if(theEnvironment.GetPDG()==NUCPDG)            // ==> "Environment is Vacuum" case
  {
#ifdef pdebug
    G4cout<<"G4QEnv::HadrQE: ***Vacuum*** #ofQ="<<nQuasmons<<G4endl;
#endif
    G4QNucleus vE = vacuum;
    G4int     nlq = 0;                           // Prototype of a#of Living Quasmons
	if(nQuasmons) for(G4int lq=0; lq<nQuasmons; lq++) if(theQuasmons[lq]->GetStatus())nlq++;
	if(nQuasmons) for(G4int iq=0; iq<nQuasmons; iq++)
	{
      G4int ist=theQuasmons[iq]->GetStatus();    // Status of the Quasmon before fragmentation
      if(ist)
	  {
        G4QHadronVector* output=theQuasmons[iq]->Fragment(vE,1);//!!!DESTROY!!! <----------------+
        G4int ast=theQuasmons[iq]->GetStatus();  // Status of the Quasmon after fragmentation    ^
        if(!ast) nlq--;                          // Reduce nlq is Quasmon decayed                ^
        G4int nHadrons = output->size();      // A#of output Hadrons in the Quasmon           ^
#ifdef pdebug
        G4cout<<"G4QEnv::HadrQE: ***Vacuum*** Q#"<<iq<<", nHadr="<<nHadrons<<G4endl; //               ^
#endif
        if(nHadrons>0)                           // Transfer QHadrons from Quasmon to Output     ^
	    {
    	  for (G4int ih=0; ih<nHadrons; ih++)    // LOOP over Hadrons produced by the Quasmon    ^
          {
            G4QHadron* curH = new G4QHadron(output->operator[](ih)); //                                  ^
            //if(curH->GetQPDG().GetPDGCode()==90002002)G4cout<<"G4QE::HQE: Get Alpha"<<G4endl;//^
            //if(curH->GetQPDG().GetPDGCode()==90001001)G4cout<<"G4QE::HQE: Get Deute"<<G4endl;//^
#ifdef pdebug
            G4cout<<"G4QEnv::HadrQE:Vacuum, H#"<<ih<<", QPDG="<<curH->GetQPDG()
                  <<",4M="<<curH->Get4Momentum()<<G4endl; //                                     ^
#endif
            theQHadrons.push_back(curH);            // Fill hadron-copy (delete equivalent)         ^
          }
	    }
        else                                     // => "Quasmon can't decay" case                ^
	    {
          G4QContent      totQC=theQuasmons[iq]->GetQC();
          G4int     tQBN=totQC.GetBaryonNumber();// Baryon Number of not decayed Quasmon         ^
          G4QNucleus     tqN(totQC);             // Define the quasmon as a nucleus              ^
          G4double   gsM=tqN.GetMZNS();          // GS Mass                                      ^
          G4LorentzVector tot4M=theQuasmons[iq]->Get4Momentum();
          G4double totQM=tot4M.m();              // Real Mass of Quasmon                         ^
          if(tQBN>0&&totQM>gsM)                  // => "Try Quasmon evaporation" case            ^
		  {
            G4QHadron* nuclQ = new G4QHadron(totQC,tot4M);
            EvaporateResidual(nuclQ);            // Try to evaporate Quasmon (delete equivalent) ^
            theQuasmons[iq]->KillQuasmon();      // Kill evaporated Quasmon                      ^
            nlq--;
		  }
          else if(iq+1<nQuasmons&&nlq>1)         // => "Try to merge with next" case             ^
		  {
            G4int s=theQuasmons[iq+1]->GetStatus();//Status of the next Quasmon                  ^
            theQuasmons[iq+1]->IncreaseBy(theQuasmons[iq]);// Merge with the next Quasmon        ^
            theQuasmons[iq]->KillQuasmon();      // Kill the week Quasmon                        ^
            if(s) nlq--;                         // Reduce a number of "living Quasmons"         ^
		  }
          else if(iq+1==nQuasmons&&iq&&nlq>1)    // => "Quasmon stack is exhosted" case          ^
		  {
            G4int s=theQuasmons[0]->GetStatus(); // Status of the first Quasmon                  ^
            theQuasmons[0]->IncreaseBy(theQuasmons[iq]);// Merge with the first Quasmon          ^
            theQuasmons[iq]->KillQuasmon();      // Kill the week Quasmon                        ^
            if(s) nlq--;                         // Reduce a number of "living Quasmons"         ^
		  }
		  else                                   // "Have a chance to recover" case              ^
		  {
#ifdef pdebug
		    G4cout<<"***G4QEnv::HadrQE:"<<iq<<",nH="<<nHadrons<<",QC="<<totQC<<",M="<<totQM<<G4endl;
		    for (G4int kq=0; kq<nQuasmons; kq++) // LOOP over Quasmons JUST for DEBUG PRINTING   ^
			  G4cout<<kq<<". Stat="<<theQuasmons[kq]->GetStatus()<<",QC="<<theQuasmons[kq]->GetQC()
                    <<",M="<<theQuasmons[kq]->Get4Momentum().m()<<G4endl; //                     ^
#endif
            G4int nOfOUT = theQHadrons.size();// Total #of QHadrons at this point             ^
            G4double  dM = totQM-gsM;            // Excitation of the Quasmon                    ^
            while(nOfOUT)                        // LOOP over all existing QHadrons              ^
            {
              G4QHadron*     theLast = theQHadrons[nOfOUT-1];     //                             ^
              G4LorentzVector last4M = theLast->Get4Momentum();   //                             ^
              G4QContent      lastQC = theLast->GetQC();          //                             ^
              G4int           lastS  = lastQC.GetStrangeness();   //                             ^
              G4int           totS   = totQC.GetStrangeness();    //                             ^
              G4int           nFr    = theLast->GetNFragments();  //                             ^
              G4int           gam    = theLast->GetPDGCode();     //                             ^
			  if(gam!=22&&!nFr&&lastS<0&&lastS+totS<0&&nOfOUT>1)  // => "Skip K-mes, gam & decayed" 
			  {
                G4QHadron* thePrev = theQHadrons[nOfOUT-2]; // Come back to the previous         ^
                theQHadrons.pop_back();         // the last QHadron is excluded from OUTPUT    ^
                theQHadrons.pop_back();         // the prev QHadron is excluded from OUTPUT    ^
                theQHadrons.push_back(theLast);      // the Last becomes the previouse              ^
                G4QHadron* destrP=theLast;        // destruction Pointer for the QHadron         ^
                delete     destrP;                // the Last QHadron is destructed              ^
                theLast = thePrev;                // Update parameters (Prev becomes the  Last)  ^
                last4M = theLast->Get4Momentum(); // 4Mom of the previouse Quasmon               ^
                lastQC = theLast->GetQC();        // Quark Content of the previouse Quasmon      ^
			  }
              else                                // Just Clear and destroy the Last             ^
              {
                theQHadrons.pop_back();         // the last QHadron is excluded from OUTPUT    ^
                delete         theLast;           // the last QHadron is deleated as instance    ^
			  }
              totQC+=lastQC;                      // Update (increase) the total QC              ^
              tot4M+=last4M;                      // Update (increase) the total 4-momentum      ^
              totQM=tot4M.m();                    // Calculate new real total mass               ^
              G4QNucleus nN(totQC);               // Define the Quasmon as a nucleus             ^
              gsM=nN.GetMZNS();                   // Calculate the new GS Mass                   ^
              dM = totQM-gsM;                     // Escitation energy for the Quasmon           ^
              if(dM>0)                            // "Mass of Q is big enough" case              ^
			  {
                theQuasmons[iq]->InitQuasmon(totQC,tot4M);// Update the week Quasmon             ^
                G4QHadronVector* curout=theQuasmons[iq]->Fragment(vE,1);//!!!DESTROY!!! <----+   ^
                G4int ast=theQuasmons[iq]->GetStatus();  // Status of the Quasmon            ^   ^
                if(!ast) nlq--;                   // Reduce nlq is Quasmon decayed           ^   ^
                G4int nHadrons=curout->size(); // A#of output Hadrons in the week Quasmon ^   ^
#ifdef pdebug
                G4cout<<"G4QEnv::HadrQEnv:VacuumRecoverQ#"<<iq<<",nH="<<nHadrons<<G4endl;//  ^   ^ 
#endif
                if(nHadrons>0)                    // => "QHadrons from Quasmon to Output"    ^   ^
	            {
    	          for (G4int ih=0; ih<nHadrons; ih++) // LOOP over Hadrons of the Quasmon    ^   ^
                  {
                    G4QHadron* curH = new G4QHadron(curout->operator[](ih)); //                      ^   ^
#ifdef pdebug
                    G4cout<<"G4QEnv::HadrQE:Recovered, H#"<<ih<<", QPDG="<<curH->GetQPDG()
                          <<",4M="<<curH->Get4Momentum()<<G4endl;    //                      ^   ^
#endif
                    theQHadrons.push_back(curH);     // Fill hadron-copy (delete equivalent)    ^   ^
                    delete curout->operator[](ih);        // >-** Necessary to delete instances **>--^   ^
                  } // End of LOOP over Hadrons of the Quasmon                               ^   ^
                  G4std::for_each(curout->begin(), curout->end(), DeleteQHadron()); // >---------------------------------------^   ^
                  delete curout;                  // >Necessary to delete Vector structure>==^   ^
                  break;                          // @@ ??                                   ^   ^
	            } // End of check for existing output Hadrons in the Quasmon                 ^   ^
                else                              //                                         ^   ^
                {
                  G4std::for_each(curout->begin(), curout->end(), DeleteQHadron());      // >---------------------------------------^   ^
                  delete curout;                  // >Necessary to delete Vector structure>--^   ^
				}
			  }
              nOfOUT  = theQHadrons.size();    // Update the value of OUTPUT entries          ^
		    }
		    //if(!nOfOUT&&totQM==0.&&nQuasmons==1)  // TRY TO EVAPORATE THE TOTAL SYSTEM           ^
			//{
            //  G4QHadron* evH = new G4QHadron(totPDG,tot4M);// Create a Hadron for ResidualNucl   ^
            //  EvaporateResidual(evH);             // Try to evaporate residual (delete equiv.)   ^
            //  output->clearAndDestroy();          // >-------------------------------------------^
            //  delete output;                      // >===========================================^
            //  CleanUp();                          //                                             ^
            //  return theQHadrons;                 //                                             ^
            //}
		    //else if(!nOfOUT)                      // Still remain not used Quasmons              ^
		    if(!nOfOUT)                           // Still remain not used Quasmons              ^
		    {
		      G4cerr<<"***G4QE::HQE:4M="<<tot4M<<",M="<<totQM<<" < gsM="<<gsM<<",dM="<<dM
                    <<",Env="<<theEnvironment<<G4endl;
              G4Exception("G4QEnvironment::HadronizeQEnvironment: Can't decay the Quasmon");  // ^
	        }
		  } // End of PANIC treatment                                                            ^
		} // End of trouble handling with Quasmon decay in Vacuum                                ^
        G4std::for_each(output->begin(), output->end(), DeleteQHadron());// >--------------------------------------------^
        delete output;                           // >============================================^
	  } // End of check for the already decayed Quasmon
	} // End of the LOOP over Quasmons
    return theQHadrons;
  }
  else                                           // ==> "Nuclear environment" case
  {
#ifdef pdebug
    G4cout<<"G4QEnv::HadrQE:FRAGMENTATION IN NUCLEAR ENVIRONMENT nQ="<<nQuasmons<<G4endl;
#endif
    G4int   c3Max = 27;
    //G4int   c3Max = 3;
    //G4int   premC = 27;
    //G4int   premC = 3;
    G4int   premC = 1;
    G4int   envA=theEnvironment.GetA();
    //if(envA>1&&envA<13) premC = 24/envA;
    if(envA>1&&envA<19) premC = 36/envA;
    //if(envA>1&&envA<25) premC = 48/envA;
    //if(envA>1&&envA<31) premC = 60/envA;
    G4int  sumstat= 2;                           // Sum of statuses of all Quasmons
    G4bool force  = false;                       // Prototype of the Force Major Flag
    G4int totC    = 0;                           // Counter for the "infinit" loop
    G4int totCM   = 227;                         // Limit for this counter
    while (sumstat||totC<totCM)                  // ===***=== The MAIN "FOREVER" LOOP ===***===
	{
      totC++;
      if(totC==totCM)
      {
        G4cerr<<"G4QEnv::HadrQE:***>>"<<totCM<<"<<***:sumstat="<<sumstat<<",nQ="<<nQuasmons<<G4endl;
        //G4Exception("G4QEnv::HadrQEnv: *** TEMPORARY EXCEPTION *** : Too long SUMSTAT loop");
	  }
      G4QContent envQC=theEnvironment.GetQCZNS();// QuarkContent of current NuclearEnvironment
      G4QContent totQC=envQC;                    // Total QuarkContent in the system
      G4double   envM =theEnvironment.GetMass(); // mass of Nuclear Environment (@@GetMZNS())
      G4double   sumM =envM;                     // Sum of all residual masses in the System @@??
      G4LorentzVector env4M=theEnvironment.Get4Momentum();      
      G4LorentzVector tot4M=env4M;               // 4-momentum of the Total System
      sumstat         =0;
      G4int     fCount=0;                        // Counter of successful but not final fragm's
	  G4int     eCount=0;                        // Counter of not yet decayed Quasmons
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
	    G4cout<<"G4QEnv::HadrQE:#"<<iq<<", Qst="<<Qst<<", Q="<<Q4M<<Q4M.m()<<QQC<<", Env="
              <<theEnvironment<<G4endl;
#endif
        if(Qst==1||Qst==3||Qst==4) fCount++;     // Increment a counter of successful fragmentations
        if(Qst>0)                  eCount++;     // Increment a counter of existing Quasmons 
	  } // End of summation LOOP over Quasmons
      // === Now we should be prepared for evaporation ===
      G4int      totChg=totQC.GetCharge();       // Total Electric Charge of the Total System
      G4int      totS  =totQC.GetStrangeness();  // Total Strangeness of the Total System
      G4int      totBN =totQC.GetBaryonNumber(); // Total Baryon Number of the Total System
      G4double   totM  =0.;                      // min (GroundSt) Mass of the Residual System
      G4int      totPDG=0;                       // Total PDG Code for the Current compound
      if(totBN<2)  
	  {
        totPDG=totQC.GetSPDGCode();              // Minimal total PDG Code for the Current compound
        if(totPDG) totM=G4QPDGCode(totPDG).GetMass(); // min Mass of the Residual System
        else G4Exception("G4QEnv::HadrQEnv: Impossible PDG for B=1");
      }
      else
	  {
        G4QNucleus totN(totQC,tot4M);            // Excited nucleus for the Residual System
        totM=totN.GetMZNS();                     // min (GroundSt) Mass of the Residual System
        totPDG=totN.GetPDG();                    // Total PDG Code for the Current compound
	  }
#ifdef pdebug
      if(totPDG==90999999||totPDG==90999000||totPDG==90000999||totPDG==89999001)
		 G4cout<<"***G4QEnv::HadrQEnv: Meson (1) PDG="<<totPDG<<", M="<<tot4M.m()<<G4endl;
#endif
      G4int           nOH=theQHadrons.size(); // A#of output hadrons
      G4LorentzVector s4M=tot4M;                 // Total 4-momentum (@@ only for checking)
      if(nOH) for(G4int ih=0; ih<nOH; ih++) s4M+=theQHadrons[ih]->Get4Momentum();     
#ifdef pdebug
	  G4cout<<"G4QEnv::HadrQE:tBN="<<totBN<<",s="<<sumstat<<",fC="<<fCount<<",eC="<<eCount<<",En="
            <<theEnvironment<<",nH="<<nOH<<",tLV="<<s4M<<totQC<<G4endl;
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
      G4double totMass = tot4M.m();              // Total effective Mass
      G4bool   Premium = eCount&&premC&&envM;    // Premium condition
      G4int    count3  =0;
      if(sumstat&&(fCount||Premium)&&!force&&count3<c3Max)// ==> "Still try to decay Quasmons" case
	  {
        if(!fCount)premC--;                      // Reduce premium efforts counter
	    for (G4int jq=0; jq<nQuasmons; jq++)     // Fragmentation LOOP over Quasmons
	    {
	      G4Quasmon* pQ     = theQuasmons[jq];   // Pointer to the current Quasmon
          G4int      status = pQ->GetStatus();   // Old status of the Quasmon
          if(status)                             // Skip dead Quasmons
		  {
 	        G4QHadronVector* output=pQ->Fragment(theEnvironment,eCount);//**!!DESTROY!!** <------+
            status = pQ->GetStatus();            // New status after fragmentation attempt       ^
#ifdef pdebug
	        G4cout<<"G4QEnv::HadrQE: **FragmAttempt** jq="<<jq<<", status="<<status          //  ^
                  <<", Environment="<<theEnvironment<<theEnvironment.Get4Momentum()<<G4endl; //  ^
#endif
            G4int nHadrons = output->size();
            if(!status||status==1||nHadrons)     // Something was filled                         ^
			{
              if(nHadrons>0)                     // Transfer QHadrons from Quasmon to Output     ^
	          {
    	        for (G4int ih=0; ih<nHadrons; ih++)      // LOOP over Q-output QHadrons          ^
                {
				  G4QHadron* inpH =output->operator[](ih);
                  G4int hC=inpH->GetCharge();   // Charge of the Hadron                          ^
                  G4int hF=inpH->GetNFragments();
                  G4double hCB=0.;              // Coulomb Barrier                               ^
                  G4double hKE=0.;              // Kinetic Energy of the Hadron                  ^
                  G4LorentzVector hLV=inpH->Get4Momentum();
                  G4bool can=hC&&!hF;
                  if(can)
                  {
                    G4int hB=inpH->GetBaryonNumber();
                    hCB=theEnvironment.CoulombBarrier(hC,hB);
                    hKE=hLV.e()-hLV.m();
				  }
                  if(can&&hKE<hCB)              // => "Suck the Hadron in Quasm or Env" case     ^
				  {
                    if(status)                  // => "Suck in the existing Quasmon" case        ^
                    {
                      G4QContent tQC=inpH->GetQC()+pQ->GetQC(); //                               ^
                      G4LorentzVector tLV=hLV+pQ->Get4Momentum();//                              ^
                      pQ->InitQuasmon(tQC,tLV); // Reinitialize the current Quasmon              ^
#ifdef pdebug
	                  G4cout<<"G4QEnv::HadrQE:Medium, H#"<<ih<<", QPDG="<<inpH->GetQPDG()
                            <<",4M="<<inpH->Get4Momentum()<<" is sucked in Quasmon"<<G4endl;
#endif
				    }
                    else                        // => "Suck in the existing Quasmon" case        ^
                    {
                      G4QContent tQC=inpH->GetQC()+theEnvironment.GetQCZNS();//                  ^
                      G4LorentzVector tLV=hLV+theEnvironment.Get4Momentum(); //                  ^
                      theEnvironment=G4QNucleus(tQC,tLV); // Reinit the current Environment      ^
#ifdef pdebug
	                  G4cout<<"G4QEnv::HadrQE:Medium, H#"<<ih<<", QPDG="<<inpH->GetQPDG()<<",4M="
                            <<inpH->Get4Momentum()<<" is sucked in Environment"<<G4endl;
#endif
				    }
				  }
                  else if(!hF)                   // => "Hadron can go out" case                  ^
                  { 
                    G4QHadron* curH = new G4QHadron(inpH); //                                    ^
#ifdef pdebug
	                G4cout<<"G4QEnv::HadrQE:Medium, H#"<<ih<<", QPDG="<<curH->GetQPDG() //       ^
                          <<",4M="<<curH->Get4Momentum()<<G4endl;
#endif
                    theQHadrons.push_back(curH);    // Fill hadron-copy (delete equivalent)         ^
				  }
				}
                pQ->ClearOutput();               // Hadron is filled, Clear Frag-output          ^
                count3=0;                        // Reset counter of empty hadronizations        ^
	          }
              else count3++;                     // Increment counter of empty hadronizations    ^
			}
            else if(status<0||status==2)         // => "PANIC or NOTHING was done" case          ^
			{
              //if(eCount==1 && CheckGroundState(pQ,true))  //                                   ^
              if(eCount==1 && CheckGroundState(pQ))         //                                   ^
              {
                pQ->KillQuasmon();               // If BackFusion succeeded, kill the Quasmon    ^
                return theQHadrons;
              }
              if(status<0&&nHadrons)
			  {
		        G4cerr<<"***G4QEnv::HadrQE: nH="<<nHadrons<<"< status="<<status<<G4endl; //      ^
                G4Exception("G4QEnvironment::HadronizeQEnvironment: Strange PANIC");
			  }
              else if(status==2)                 // Check PANIC conditions for status=2          ^
			  {
                if(theEnvironment!=vacuum)       // "Nuclear Environment" case                   ^
				{
                  G4LorentzVector t4M=pQ->Get4Momentum()+theEnvironment.Get4Momentum(); //       ^
                  G4double      tM=t4M.m();      // Real total (with environment) mass           ^
                  G4QContent   qQC= pQ->GetQC(); // QuarkContent of the Quasmon                  ^
                  G4QContent envQC=theEnvironment.GetQCZNS(); // QuarkCont of NucEnviron         ^
                  G4QContent curQC=envQC+qQC;    // Total Quark Content                          ^
                  G4QNucleus curE(curQC);        // Pseudo nucleus for the Total System          ^
                  G4double   curM=curE.GetMZNS();// min mass of the Total System                 ^
#ifdef pdebug
    		      G4cout<<"G4QEnv::HadrQEnv:#"<<jq<<",tM="<<tM<<" > gstM="<<curM<<curE<<G4endl;
#endif
                  if(tM<curM) status=-1;         // Q+E is below the Mass Shell - PANIC          ^
                }
                else                             // "Vacuum" case                                ^
				{
                  G4LorentzVector t4M=pQ->Get4Momentum(); // 4Mom for the Quasmon                ^
                  G4QPDGCode QPDGQ=pQ->GetQPDG();// QPDG Code for the Quasmon                    ^
                  G4int PDGQ=QPDGQ.GetPDGCode(); // PDG Code of the QUASMON                      ^
#ifdef pdebug
				  G4cout<<"G4QEnv::HadrQEnv: vacuum PDGQ="<<PDGQ<<G4endl; //                     ^
#endif
                  if(!PDGQ) status=-1;           // Unknown Quasmon in Vaquum - PANIC            ^
                  else if (PDGQ!=10)             // @@ Chipolino can wait @@                     ^
				  {
                    G4double qM =t4M.m();        // Real mass of the Quasmon                     ^
                    G4double gsM=QPDGQ.GetMass();// GSmass of the Quasmon                        ^
#ifdef pdebug
    		        G4cout<<"G4QEnv::HadrQEnv:#"<<jq<<", qM="<<qM<<" > gsM="<<gsM<<G4endl; //    ^
#endif
					if(abs(qM-gsM)<0.0001)       // "Fill & Kill" Case                           ^
					{
                      G4QHadron* resQ = new G4QHadron(PDGQ,t4M);
                      theQHadrons.push_back(resQ);  // @@ Check Dibarions @@ (delete equivalent)    ^
                      pQ->KillQuasmon();         // Make done the current Quasmon                ^
					}
                    else if(qM<gsM) status=-1;   // Below Mass Shell - PANIC                     ^
				  }
				}
			  }
              else if(status==3) count3++;
              if(status<0)                       // Panic: Quasmon is below the mass shell       ^
			  {
                //if(eCount==1 && DecayInEnvQ(pQ))
                //{
                //  pQ->KillQuasmon();
                //  return theQHadrons;
                //}
                G4int    ppm=jq;                 // Initialized by PANIC Quasmon pointer         ^
                G4int    nRQ=0;                  // Prototype of a#of additional real Quasmons   ^
#ifdef pdebug
    		    G4cout<<"G4QEnv::HadrQEnv: ***PANIC*** for jq="<<jq<<G4endl;
#endif
                G4ThreeVector vp= pQ->Get4Momentum().vect(); // momentum of the PANIC Quasmon    ^
                G4double dpm=1.e+30;             // Just a big number (dot product of momenta)   ^
	            for(G4int ir=0; ir<nQuasmons; ir++)
	            {
                  if(ir!=jq)                     // Skip the current (PANIC) Quasmon             ^
				  {
	                G4Quasmon* rQ = theQuasmons[ir];
                    G4int Qst = rQ->GetStatus(); // Status of a Quasmon                          ^
                    if(Qst>0)                    // Skip the dead Quasmon                        ^
				    {
					  nRQ++;                     // Increment real-Quasmon-counter               ^
                      G4double dp=vp.dot(rQ->Get4Momentum().vect());
                      if(dp<dpm)
					  {
                        ppm=ir;                  // Remember the index of MinProj Quasmon        ^
                        dpm=dp;                  // Remember the value of Minimum Projection     ^
					  }
				    }
				  }
                }// End of the partner-search-for-the-PANIC-Quasmon LOOP                         ^
                if(nRQ)                          // Merge with the best candidate                ^
    		    {
	              G4Quasmon*      rQ = theQuasmons[ppm];
                  G4QContent      rQC= rQ->GetQC();
                  G4LorentzVector r4M= rQ->Get4Momentum();
                  rQC               += pQ->GetQC();
                  r4M               += pQ->Get4Momentum();
                  rQ->InitQuasmon(rQC, r4M);     // Make new Quasmon                             ^
                  pQ->KillQuasmon();             // Delete old Quasmon                           ^
			    }
                else // No candidate to resolve PANIC was found                                  ^
			    {
                  //if(eCount==1 && CheckGroundState(pQ, true)) //  BackFusion attempt           ^
                  if(eCount==1 && CheckGroundState(pQ)) //  BackFusion attempt                   ^
                  {
                    pQ->KillQuasmon();           //                                              ^
					return theQHadrons;          //                                              ^
                  }
#ifdef pdebug
		          G4cout<<"G4QEnv::HadrQEnv: Cann't resolve PANIC, try to Evaporate"<<G4endl;//  ^
#endif
                  force=true;                    // Make the force decision                      ^
                  break;                         //                                              ^
			    }
			  }
			}
            G4std::for_each(output->begin(), output->end(), DeleteQHadron());         // >--------------------------------------------^
            delete output;                       // >============================================^
		  }
	    } // End of fragmentation LOOP over Quasmons (jq)
      }
      else if(totMass>totM+.001)                 // ==> "Try Evaporation or decay" case
	  {
#ifdef pdebug
		G4cout<<"G4QEnv::HadrQE: M="<<totMass<<",PDG="<<totPDG<<",B="<<totBN<<",GSM="<<totM<<",dM="
              <<totMass-totM<<G4endl;
        //    <<totMass-totM<<",NaK="<<NaK<<",KPDG="<<aKPDG<<",NPi="<<NPi<<",PiPDG="<<PiPDG<<G4endl;
#endif
        if(totBN<2)                              // ==> "Baryon/Meson residual Quasmon" case
		{
          if(totPDG==90999999||totPDG==90999000||totPDG==90000999||totPDG==89999001)//=>"M"case
		  {
		    G4cerr<<"***G4QEnv::HadrQEnv: Meson (2) PDG="<<totPDG<<", M="<<totMass<<G4endl;
		  }
          else if(totPDG==1114||totPDG==2224)    // ==> "DELTA- or DELTA++" case (@@anti DELTA)
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
              //delta->SetNFragments(2);           // Put a#of Fragments=2
              //theQHadrons.push_back(delta);         // Fill the residual DELTA (delete equivalent)
              // Instead
              delete delta;
              //
              G4LorentzVector b4Mom(0.,0.,0.,mBar);
              G4LorentzVector m4Mom(0.,0.,0.,mMes);
              if(!G4QHadron(tot4M).DecayIn2(b4Mom, m4Mom))
              {
                G4cerr<<"***G4QEnv::HadronizeQEnv: B="<<bPDG<<"(m="<<mBar<<") + M="
					  <<mPDG<<"(m="<<mMes<<") >(?) mD="<<totMass<<G4endl;
    	        G4Exception("G4QEnvironment::HadronizeQEnvironment: D->B+M decay failed");
              }
#ifdef pdebug
	          G4cout<<"G4QEnv::HadronizeQEnv: DELTA="<<totPDG<<tot4M<<" -> Bar="
                    <<bPDG<<m4Mom<<" + Mes="<<mPDG<<m4Mom<<G4endl;
#endif
              G4QHadron* curBar = new G4QHadron(bPDG,b4Mom);
              theQHadrons.push_back(curBar);        // Fill the baryon (delete equivalent)
              G4QHadron* curMes = new G4QHadron(mPDG,m4Mom);
              theQHadrons.push_back(curMes);        // Fill the meson (delete equivalent)
              return theQHadrons;
			}
		  }
          else if(totPDG==10)                    // ==> "Chipolino" case
		  {
            G4QChipolino resChip(totQC);         // define the residual as a Chipolino @@DONE (752?)
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
	        G4cout<<"G4QEnv::HadronizeQEnv: Chipo="<<tot4M<<" -> h1="
                  <<h1PDG<<h14Mom<<" + Mes="<<h2PDG<<h24Mom<<G4endl;
#endif
            G4QHadron* curH1 = new G4QHadron(h1PDG,h14Mom);
            theQHadrons.push_back(curH1);           // Fill the curH1 (delete equivalent)
            G4QHadron* curH2 = new G4QHadron(h2PDG,h24Mom);
            theQHadrons.push_back(curH2);           // Fill the curH2 (delete equivalent)
            return theQHadrons;
		  }
          else if(totBN<2&&totPDG&&totMass<totM+mPi0+.001)// ==> "Meson/Baryon+gamma" case
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
	        G4cout<<"G4QEnv::HadrQEnv:"<<tot4M<<" -> h="<<totPDG<<h4Mom<<" + gamma="<<g4Mom<<G4endl;
#endif
            G4QHadron* curG = new G4QHadron(22,g4Mom);
            theQHadrons.push_back(curG);            // Fill the gamma (delete equivalent)
            G4QHadron* curH = new G4QHadron(totPDG,h4Mom);
            if(totPDG==92000000||totPDG==90002000||totPDG==90000002) DecayDibaryon(curH);// (DelEqu)
            else theQHadrons.push_back(curH);       // Fill the baryon (delete equivalent)
            return theQHadrons;
		  }
          else if(totBN<2&&totPDG)                        // ==> "Meson/Baryon+pi" case
		  {
            G4int piPDG=111;
            G4double mpi=mPi0;
            G4int mbPDG=totPDG;
            G4double mbm=totM;
            if(totPDG==1114)
			{
			  piPDG=-211;
              mpi=mPi;
              mbPDG=2112;
              mbm=mNeut;
            }
            else if(totPDG==2224)
			{
			  piPDG=211;
              mpi=mPi;
              mbPDG=2212;
              mbm=mProt;
            }
            else if(totPDG==113)
			{
			  piPDG=-211;
              mpi=mPi;
              mbPDG=211;
              mbm=mPi;
            }
            G4LorentzVector h4Mom(0.,0.,0.,mbm);
            G4LorentzVector g4Mom(0.,0.,0.,mpi);
            if(!G4QHadron(tot4M).DecayIn2(h4Mom, g4Mom))
            {
              G4cerr<<"***G4QEnv::HadronizeQEnv: h="<<mbPDG<<"(m="<<mbm
                    <<") + pi(m="<<mpi<<") >(?) mTot="<<totMass<<G4endl;
    	      G4Exception("G4QEnvironment::HadronizeQEnv: Decay in pi-meson failed");
            }
#ifdef pdebug
	        G4cout<<"G4QEnv::HadQE:"<<tot4M<<" -> h="<<mbPDG<<h4Mom<<" + pi="<<piPDG<<g4Mom<<G4endl;
#endif
            G4QHadron* curH = new G4QHadron(mbPDG,h4Mom);
            if(totPDG==92000000||totPDG==90002000||totPDG==90000002) DecayDibaryon(curH);// (DelEqu)
            else theQHadrons.push_back(curH);       // Fill the baryon (delete equivalent)
            G4QHadron* curG = new G4QHadron(piPDG,g4Mom);
            theQHadrons.push_back(curG);            // Fill the pi0 (delete equivalent)
            return theQHadrons;
		  }
          else                                   // ==> "|B|<2 new Quasmon" case
		  {
            G4Quasmon* resid = new G4Quasmon(totQC,tot4M); // delete is 3 lines below <-+
            G4QNucleus vacuum(90000000);         //                                     ^
 	        G4QHadronVector* curout=resid->Fragment(vacuum,1);//**!!DESTROY!!** <-+     ^
            delete resid;                        //_______________________________^_____^
            G4int nHadrons = curout->size();  // a#of Hadrons in the outHV     ^
            if(nHadrons>0)                       // Transfer QHadrons to Output   ^
	        {
    	      for (G4int ih=0; ih<nHadrons; ih++)// LOOP over output QHadrons     ^
              {
                G4QHadron* curH = new G4QHadron(curout->operator[](ih)); //               ^
#ifdef pdebug
				G4cout<<"G4QEnv::HadrQE:NewB<2, H#"<<ih<<", QPDG="<<curH->GetQPDG()
                      <<",4M="<<curH->Get4Momentum()<<G4endl;
#endif
                theQHadrons.push_back(curH);        // Insert (delete equivalent)    ^
			  }
	        }
			else                                 //                               ^
			{
              G4cerr<<"***G4QEnv::HadronizeQEnv: MQ="<<tot4M.m()<<",QC="<<totQC<<G4endl;
			  G4Exception("G4QEnvironment::HadronizeQEnv: Quasmon doesn't decay");
			}
            G4std::for_each(curout->begin(), curout->end(), DeleteQHadron()); // >-----------------------------^
            delete curout;                       // >=============================^
            return theQHadrons;
		  }
		}
        else
		{
          G4QContent    tQC =totQC;                  // Not subtracted copy for error prints
          G4int      NaK    =0;                      // a#of additional Kaons/anti-Kaons
          G4int      aKPDG  =0;                      // PDG of additional Kaons/anti-Kaons
          G4double   MaK    =0.;                     // Total Mass of additional Kaons/anti-Kaons
          G4int      NPi    =0;                      // a#of additional pions
          G4int      PiPDG  =0;                      // PDG of additional pions
          G4double   MPi    =0.;                     // Total Mass of additional pions
          if    (totBN>0&&totS<0&&totChg+totChg>=totBN)// => "additional K+" case
	      {
            aKPDG=321;
            NaK=-totS;
            MaK=mK*NaK;
            totQC+=totS*KpQC;
            totChg+=totS;                            // Charge reduction (totS<0!)
            totS=0;                                  // Anti-strangness goes to anti-Kaons
	      }
          else if (totBN>0&&totS<0)                 // => "additional aK0" case
	      {
            aKPDG=311;
            NaK=-totS;
            MaK=mK0*NaK;
            totQC+=totS*K0QC;
            totS=0;                                  // Anti-strangness goes to anti-Kaons
	      }
          else if (totBN>0&&totS>totBN&&totBN<totS+totChg)// => "additional K0" case
	      {// @@ Here Ksi0 check should be added totS=2>totBN=1&&totBN=1<totS=2+totChg=0
            aKPDG=-311;
            NaK=totS-totBN;
            MaK=mK0*NaK;
            totQC+=NaK*K0QC;
            totS-=NaK;                               // Reduce residualstrangeness
	      }
          else if (totBN>0&&totS>totBN&&totChg<0)    // => "additional K-" case
	      {// @@ Here Ksi- check should be added totS=2>totBN=1&&totChg=-1<0
            aKPDG=-321;
            NaK=totS-totBN;
            MaK=mK0*NaK;
            totQC+=NaK*KpQC;
            totChg+=NaK;                             // Increase residual charge
            totS-=NaK;                               // Reduce residual strangeness
	      }
          // === Now residual DELTAS should be subtracted === 
          if      (totBN>0&&totChg>totBN-totS)       // => "additional PI+" case
	      {// @@ Here Sigma+ check should be added totChg=1>totBn=1-totS=1
            PiPDG=211;
            NPi=totChg-totBN+totS;
            MPi=mPi*NPi;
            totQC-=NPi*PiQC;
            totChg-=NPi;
	      }
          else if (totBN>0&&totChg<0)                // => "additional PI-" case
	      {// @@ Here Sigma- check should be added totChg<0
            PiPDG=-211;
            NPi=-totChg;
            MPi=mPi*NPi;
            totQC+=NPi*PiQC;                         // Now anti-Pions must be subtracted
            totChg+=NPi;
	      }
          else if (!totBN&&totChg>1-totS)            // => "additional PI+" case
	      {// @@ Here Sigma+ check should be added totChg=1>totBn=1-totS=1
            PiPDG=211;
            NPi=totChg+totS-1;
            MPi=mPi*NPi;
            totQC-=NPi*PiQC;
            totChg-=NPi;
	      }
          else if (!totBN&&totChg<-1-totS)           // => "additional PI-" case
	      {// @@ Here Sigma- check should be added totChg<0
            PiPDG=-211;
            NPi-=totChg+totS+1;
            MPi=mPi*NPi;
            totQC+=NPi*PiQC;                         // Now anti-Pions must be subtracted
            totChg+=NPi;
	      }
          G4double      totRM=0.;                    // min (GroundSt) Mass of the Residual System
          if(totBN<2)
	      {
            totPDG=totQC.GetSPDGCode();              // Minimal PDG Code for the Residual compound
            if(totPDG) totRM=G4QPDGCode(totPDG).GetMass(); // min Mass of the Residual System
            else G4Exception("G4QEnv::HadrQEnv: Impossible PDG for B=1");
          }
          else
	      {
            G4QNucleus totN(totQC,tot4M);            // Excited nucleus for the Residual System
            totRM=totN.GetMZNS();                    // min (GroundSt) Mass of the Residual System
            totPDG=totN.GetPDG();                    // Total PDG Code for the Current compound
	      }
          if(NaK)                                    // ==> "Decay in K0 or K+ + NPi" case
	      {//@@ Can (must) be moved to EvaporateResidual ??
            if(NaK==1&&!NPi)                         // ==> "One anti-strange K" case
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
	          G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> M="
                    <<aKPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<totQC<<G4endl;
#endif
              G4QHadron* curK = new G4QHadron(aKPDG,m4Mom);
              theQHadrons.push_back(curK);            // Fill the curK (delete equivalent)
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
              EvaporateResidual(curN);             // Try to evaporate residual (delete equivalent)
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
	          G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> nK="<<aKPDG<<k4Mom
                    <<" + nPi="<<PiPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<G4endl;
#endif
              G4LorentzVector onePi=(1./NPi)*m4Mom;// 4-mom of one pion  
              for (G4int ip=0; ip<NPi; ip++)
			  {
                G4QHadron* curP = new G4QHadron(PiPDG,onePi);
                theQHadrons.push_back(curP);          // Fill the curM (delete equivalent)
	 		  }
              G4LorentzVector oneK=(1./NaK)*k4Mom; // 4-mom of one kaon  
              for (G4int jp=0; jp<NaK; jp++)
			  {
                G4QHadron* curP = new G4QHadron(aKPDG,oneK);
                theQHadrons.push_back(curP);          // Fill the curM (delete equivalent)
	 		  }
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
              EvaporateResidual(curN);             // Try to evaporate residual (delete equivalent)
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
	          G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> N*K="<<aKPDG
                    <<" (4M1="<<m4Mom<<" + 4M2="<<k4Mom<<") + N="<<totPDG<<n4Mom<<G4endl;
#endif
              G4LorentzVector one1=(1./N1K)*m4Mom;  // 4-mom of one kaon  
              for (G4int ip=0; ip<N1K; ip++)
			  {
                G4QHadron* curP = new G4QHadron(aKPDG,one1);
                theQHadrons.push_back(curP);          // Fill the curP (delete equivalent)
	 		  }
              G4LorentzVector one2=(1./N2K)*k4Mom; // 4-mom of one kaon  
              for (G4int jp=0; jp<N2K; jp++)
			  {
                G4QHadron* curP = new G4QHadron(aKPDG,one2);
                theQHadrons.push_back(curP);          // Fill the curP (delete equivalent)
	 		  }
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
              EvaporateResidual(curN);             // Try to evaporate residual (delete equivalent)
		    }
            return theQHadrons;
		  }
          else if(NPi)                             // ==> "Decay in Pi+ or Pi-" case
	      {
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
	          G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> M="
                    <<PiPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<totQC<<G4endl;
#endif
              G4QHadron* curK = new G4QHadron(PiPDG,m4Mom);
              theQHadrons.push_back(curK);            // Fill the curK (delete equivalent)
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
              EvaporateResidual(curN);             // Evaporate residual (delete equivalent)
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
	          G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> N*PI="<<PiPDG
                    <<" (4M1="<<m4Mom<<" + 4M2="<<k4Mom<<") + N="<<totPDG<<n4Mom<<G4endl;
#endif
              G4LorentzVector one1=(1./N1Pi)*m4Mom;  // 4-mom of one pion  
              for (G4int ip=0; ip<N1Pi; ip++)
			  {
                G4QHadron* curP = new G4QHadron(PiPDG,one1);
                theQHadrons.push_back(curP);          // Fill the curP (delete equivalent)
	 		  }
              G4LorentzVector one2=(1./N2Pi)*k4Mom; // 4-mom of one pion  
              for (G4int jp=0; jp<N2Pi; jp++)
			  {
                G4QHadron* curP = new G4QHadron(PiPDG,one2);
                theQHadrons.push_back(curP);          // Fill the curP (delete equivalent)
	 		  }
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
              EvaporateResidual(curN);             // Try to evaporate residual (delete equivalent)
		    }
            return theQHadrons;
		  }
		}
        theEnvironment.InitByPDG(NUCPDG);        // Cancele the Environment 
        G4QHadron* evH = new G4QHadron(totPDG,tot4M);// Create a Hadron for ResidualNucl
        EvaporateResidual(evH);                  // Try to evaporate residual (delete equivalent)
        return theQHadrons;
	  }
      else                                       // ==> "Only GSEnvironment exists" case
      { 
        if(totPDG==90000000||abs(totMass)<0.000001)
		{
          CleanUp();
          return theQHadrons;
		}
        G4double dM=totMass-totM;
#ifdef pdebug
		G4cout<<"G4QEnv::HadrQEnv:GroundState tM-GSM="<<dM<<",GSM="<<totM<<",tPDG="<<totPDG<<",nQ="
              <<nQuasmons<<G4endl;
#endif
	    G4Quasmon*       pQ = theQuasmons[0];  // Pointer to the only Quasmon          
        G4QPDGCode    QQPDG = pQ->GetQPDG();   // QPDG of the Quasmon
        G4int          QPDG = QQPDG.GetPDGCode();
        if(dM>-0.001)
		{
#ifdef pdebug
		  G4cout<<"G4QEnv::HadrQEnv:ExcitedNucleus, dM="<<dM<<">0, tBN="<<totBN<<",nQ="<<G4endl;
#endif
          G4QHadron* evH = new G4QHadron(totPDG,tot4M);// Create a Hadron for ResidualNucl
          if(dM>0.001&&totBN>0)EvaporateResidual(evH); // Try to evaporate residual (delete equiv.)
          else if(totBN==2)    DecayDibaryon(evH);     // Decay dibaryon (delete equivalent)
          else if(totBN==5)    DecayAlphaBar(evH); //DelEqu
          else if(totPDG==90004004) DecayAlphaAlpha(evH); //DelEqu
          else                 theQHadrons.push_back(evH);// Fill to OUTPUT as it is (delete equiv.)
		}
        else if(nQuasmons==1&&QPDG!=22)          // => "Decay in Quasmon + QEnviron" case
		{
          G4int envPDG = theEnvironment.GetPDG();// PDGCode of the NuclQEnvironment
#ifdef pdebug
		  G4cout<<"G4QEnv::HadrQEnv: nQ=1, QPDG=="<<QPDG<<G4endl;
#endif
          if(!QPDG)
		  {
			G4cerr<<"***G4QEnv::HadrQE: Quasmon is unknown QHadron: PDG="<<QPDG<<G4endl;
			G4Exception("G4QEnvironment::HadronizeQEnvironment: (2) Cann't decay QEnv");
		  }
          // => "Quasmon-Chipolino or Environment-Dibaryon" case
          else if(QPDG==10||QPDG==92000000||QPDG==90002000||QPDG==90000002)
		  {
            G4QContent QQC = pQ->GetQC();        // Quark Content of the Quasmon
            G4QPDGCode h1QPDG=nQPDG;             // QPDG of the first hadron
            G4double   h1M   =mNeut;             // Mass of the first hadron
            G4QPDGCode h2QPDG=h1QPDG;            // QPDG of the second hadron
            G4double   h2M   =mNeut;             // Mass of the second hadron
            if(QPDG==10)
            {
              G4QChipolino QChip(QQC);           // define the Quasmon as a Chipolino
              h1QPDG=QChip.GetQPDG1();           // QPDG of the first hadron
              h1M   =h1QPDG.GetMass();           // Mass of the first hadron
              h2QPDG=QChip.GetQPDG2();           // QPDG of the second hadron
              h2M   =h2QPDG.GetMass();           // Mass of the second hadron
			}
			else if(QPDG==90002000)
			{
              h1QPDG=pQPDG;                      // QPDG of the first hadron
              h1M   =mProt;                      // Mass of the first hadron
              h2QPDG=h1QPDG;                     // QPDG of the second hadron
              h2M   =mProt;                      // Mass of the second hadron
			}
			else if(QPDG==92000000)
			{
              h1QPDG=lQPDG;                      // QPDG of the first hadron
              h1M   =mLamb;                      // Mass of the first hadron
              h2QPDG=h1QPDG;                     // QPDG of the second hadron
              h2M   =mLamb;                      // Mass of the second hadron
			}
            if(h1M+h2M+envM<totMass)             // => "Three parts decay" case
			{
              G4LorentzVector h14M(0.,0.,0.,h1M);
              G4LorentzVector h24M(0.,0.,0.,h2M);
              G4LorentzVector e4M(0.,0.,0.,envM);
              if(!G4QHadron(tot4M).DecayIn3(h14M,h24M,e4M))
			  {
				G4cerr<<"***G4QEnv::HadQEnv: tM="<<tot4M.m()<<","<<totMass<<"-> h1="<<h1QPDG<<"("
					  <<h1M<<")+h2="<<h1QPDG<<"("<<h2M<<")+envM="<<envM<<"=="<<h1M+h2M+envM<<G4endl;
				G4Exception("G4QEnv::HadrQEnv:QChipo+Environment DecayIn3 did not succeed");
			  }
              G4int h1PDG =h1QPDG.GetPDGCode();  // PDG code of the first hadron            
              G4int h2PDG =h2QPDG.GetPDGCode();  // PDG code of the second hadron
              G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
              theQHadrons.push_back(h1H);           // (delete equivalent)
              G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
              theQHadrons.push_back(h2H);           // (delete equivalent)
              G4QHadron* qeH = new G4QHadron(envPDG,e4M);
              theQHadrons.push_back(qeH);           // (delete equivalent)
			}
			else
			{
              //if(eCount==1&&CheckGroundState(pQ,true))// BackFusion attempt
              if(eCount==1&&CheckGroundState(pQ))// BackFusion attempt
              {
                pQ->KillQuasmon();
                return theQHadrons;
              }
#ifdef pdebug
              G4cerr<<"***G4QEnv::HadQEnv: tM="<<tot4M.m()<<"< h1="<<h1QPDG<<"(M="<<h1M<<") + h2="
                    <<h1QPDG<<"(M="<<h2M<<") + envM="<<envM<<"="<<h1M+h2M+envM<<G4endl;
			  //G4Exception("G4QEnv::HadrQEnv:QChipo+Env mass is more than decaying mass");
#endif
              theEnvironment.InitByPDG(NUCPDG);  // Cancele the Environment 
              G4QHadron* evH = new G4QHadron(totQC,tot4M);// Create a Hadron for ResidualNucl
              //if(totBN==2) DecayDibaryon(evH);   // Try to decay unstable Dibaryon
              //else if(totBN==5) DecayAlphaBar(evH);   // Try to decay unstable A5 system
              //else EvaporateResidual(evH);       // Try to evaporate residual (delete equivalent)
              EvaporateResidual(evH);       // Try to evaporate residual (delete equivalent)
              return theQHadrons;
			}
		  }
          else                                   // => "Two particles decay" case
		  {
            G4double QM = pQ->Get4Momentum().m();// Real Mass of the Quasmon
            G4double GSM= QQPDG.GetMass();       // GS Mass of the Quasmon
            if(QM<GSM||GSM+envM>totMass)
			{
#ifdef pdebug
		      G4cerr<<"***G4QEnv::HadrQE: QM="<<QM<<" < QGSM="<<GSM<<" or GSM+envM="<<GSM+envM
                    <<" > totM="<<totMass<<", envM="<<envM<<", Q="<<QQPDG<<",tQC="<<totQC<<G4endl;
              //G4Exception("G4QEnvironment::HadronizeQEnvironment: (3) Cann't decay");
#endif
              theEnvironment.InitByPDG(NUCPDG);  // Cancele the Environment 
              G4QHadron* evH = new G4QHadron(totQC,tot4M);// Create a Hadron for ResidualNucl
              //if(totBN==2) DecayDibaryon(evH);   // Try to decay unstable Dibaryon
              //else if(totBN==5) DecayAlphaBar(evH);   // Try to decay unstable A5 system
              //else EvaporateResidual(evH);       // Try to evaporate residual (delete equivalent)
              EvaporateResidual(evH);       // Try to evaporate residual (delete equivalent)
              return theQHadrons;
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
              theQHadrons.push_back(qH);            // (delete equivalent)
              G4QHadron* qeH = new G4QHadron(envPDG,qe4M);
              if(envPDG==92000000||envPDG==90002000||envPDG==90000002) DecayDibaryon(qeH); //(DelEq)
              else theQHadrons.push_back(qeH);      // (delete equivalent)
			}
		  }
		}
        else                                     // "Last decay was fatal" case
		{
#ifdef pdebug
		  G4cout<<"***G4QEnv::HadrQE: M="<<totMass<<",dM="<<dM<<",nQ="<<nQuasmons<<G4endl;
#endif
          G4int          nOfOUT  = theQHadrons.size();
          while(nOfOUT)
          {
            G4QHadron*     theLast = theQHadrons[nOfOUT-1];
            G4LorentzVector last4M = theLast->Get4Momentum();
            G4QContent      lastQC = theLast->GetQC();
            G4int           lastS  = lastQC.GetStrangeness();
            G4int           totS   = totQC.GetStrangeness();
            G4int           nFr    = theLast->GetNFragments();
            G4int           gam    = theLast->GetPDGCode();
			if(gam!=22&&!nFr&&lastS<0&&lastS+totS<0&&nOfOUT>1)  // => "Skip K-mes, gamma & decayed" 
			{
              G4QHadron* thePrev = theQHadrons[nOfOUT-2];
              theQHadrons.pop_back();         // the last QHadron is excluded from OUTPUT
              theQHadrons.pop_back();         // the prev QHadron is excluded from OUTPUT
              theQHadrons.push_back(theLast);      // the Last becomes the Prev
              theLast = thePrev;                // Update parameters (Prev instead of Last)
              last4M = theLast->Get4Momentum();
              lastQC = theLast->GetQC();
			}
            else theQHadrons.pop_back();      // the last QHadron is excluded from OUTPUT 
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
                G4QNucleus  newN0(totQC-K0QC);
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
                  theQHadrons.push_back(H1);       // (delete equivalent)
                  G4QHadron* H2 = new G4QHadron(PDG2,qe4M);
                  theQHadrons.push_back(H2);       // (delete equivalent)
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
                  theQHadrons.push_back(H1);        // (delete equivalent)
                  G4QHadron* H2 = new G4QHadron(PDG2,qe4M);
                  theQHadrons.push_back(H2);        // (delete equivalent)
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
              G4QHadron* evH = new G4QHadron(totPDG,tot4M);// Create a Hadron for Residual Nucleus
              if(dM>0.001&&totBN>2)EvaporateResidual(evH); // Try to evaporate ResidNucl (del.equiv)
              else if(totBN==2) DecayDibaryon(evH);        // Decay dibaryon (delete equivalent)
              else if(totBN==5) DecayAlphaBar(evH);        // "Alpha+Baryon Decay" case (delete eq.)
              else if(totPDG==90004004) DecayAlphaAlpha(evH); // "Alpha+Alpha Decay" case (del eq.)
              else      theQHadrons.push_back(evH);// Fill to OUTPUT as it is (delete equivalent)
              break;
		    }
            nOfOUT  = theQHadrons.size();    // Update the value of OUTPUT entries
		  }
		  if(!nOfOUT)
		  {
		    G4cerr<<"***G4QEnv::HadrQE: M="<<totMass<<" < gsM="<<totM<<", dM="<<dM<<G4endl;
            G4Exception("G4QEnvironment::HadronizeQEnvironment: Exhosted but can't decay QEnv");
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
  G4int nQuasmons = theQuasmons.size();
  if (nQuasmons) for (G4int iq=0; iq<nQuasmons; iq++)theQuasmons[iq]->KillQuasmon();
}

//Evaporate Residual Nucleus
void G4QEnvironment::EvaporateResidual(G4QHadron* qH, G4bool corFlag)
{//  ================================================================
  static const G4double mHel6 = G4QPDGCode(2112).GetNuclMass(2,4,0);
  static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mDeut = G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mNeut = G4QPDGCode(2112).GetMass();
  static const G4double mProt = G4QPDGCode(2212).GetMass();
  static const G4double mLamb = G4QPDGCode(3122).GetMass();
  static const G4double mPi   = G4QPDGCode(211).GetMass();
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QContent deutQC(3,3,0,0,0,0);
  static const G4QContent alphQC(6,6,0,0,0,0);
  G4int       thePDG = qH->GetPDGCode();       // Get PDG code of the Residual Nucleus
  /// @@@@@@@ *** TEMPORARY TO AVOID HYPERMUCLEI FOR GEANT4 *** @@@@@@@
  if(thePDG!=91000000 && thePDG>90999999)
  {
    G4int S=(thePDG-90000000)/1000000;
    thePDG-=S*999999;
    qH->SetQPDG(G4QPDGCode(thePDG));
  }
  /// @@@ *** ^^^ END OF TEMPORARY ^^^ *** @@@
  G4QContent  theQC  = qH->GetQC();            // Quark Content of the hadron
  if(!thePDG) thePDG = theQC.GetSPDGCode();    // If PDG code's not transferred, get it from QC
  G4LorentzVector q4M = qH->Get4Momentum();    // Get 4-momentum of the Residual Nucleus
  G4double    totMass = qH->GetMass();         // Get Real Mass of the Residual Nucleus
#ifdef pdebug
  G4cout<<"G4QEnvironment::EvaporateResidual(EvaRes): ===IN==> PDG="<<thePDG<<",4Mom="<<q4M<<G4endl;
#endif
  if     (thePDG==90000000)                    // ==> "Nothing in the INPUT Hadron" case KEEP IT!
  {
    delete qH;
#ifdef debug
    G4cerr<<"***G4QEnv::EvaRes: Residual Nucleus is jast a vacuum PDG=90000000, 4Mom="<<q4M<<G4endl;
#endif
    return;
  }
  else if((thePDG>80000000&&thePDG!=90000000||thePDG==2112||thePDG==2212||thePDG==3122) &&
          thePDG!=90002999 && thePDG!=89999003)
  { // @@ Should be improved for Sigma+, Sigma-, Ksi0 & Ksi- content in the Total Residual Nucleus
    if(thePDG<80000000)                        // => "Switch from QHadron coding to QNuclear coding"
    {
      if     (thePDG==2112) thePDG=90000001;   // n
      else if(thePDG==2212) thePDG=90001000;   // p
      else if(thePDG==3122) thePDG=91000000;   // lambda
	}
    G4QNucleus qNuc(q4M,thePDG);               // Make a Nucleus out of the Total Residual Nucleus
    G4double GSMass =qNuc.GetGSMass();         // Ground State Mass of the Total Residual Nucleus
    G4QContent totQC=qNuc.GetQCZNS();          // QuarkContent of the TotalResidualNucleus (=theQC?)
    G4int    bA     =qNuc.GetA();              // A#of baryons in the Total Residual Nucleus
    G4int    bZ     =qNuc.GetZ();              // A#of protons in the Total Residual Nucleus
    G4int    bN     =qNuc.GetN();              // A#of neutrons in the Total Residual Nucleus
    G4int    bS     =qNuc.GetS();              // A#of lambdas in the Total Residual Nucleus
#ifdef debug
    if(bZ==2&&bN==5)G4cout<<"G4QE::EvR: GSM="<<GSMass<<" > "
						  <<G4QPDGCode(2112).GetNuclMass(2,4,0)+mNeut<<G4endl;
    G4double dM=totMass-GSMass;
	G4cout<<"G4QEnv::EvaRes:"<<qNuc<<",PDG="<<thePDG<<",M="<<totMass<<",dM="<<dM<<G4endl;
    ////////if(dM>7.) G4Exception("G4QEnvironment::EvaporateResidual: CALLED");
#endif
    G4bool   bsCond =qNuc.SplitBaryon();       // (Bary/Deut/Alph)SeparetionCondition for TotResNucl
    G4bool   dbsCond=qNuc.Split2Baryons();     // (Two Baryons)SeparetionCondition for TotResidNucl
    if((bA==1||!bsCond&&!dbsCond)&&totMass>GSMass+.003)//==>"Fuse&DecayTechnology to avoid gammaDec"
	//if(2>3)                                    // Close "Fuse&Decay Technology" ***@@@***
	{
#ifdef debug
	  G4cout<<"G4QEnv::EvaR: Can't SplitBar s="<<bsCond<<",M="<<totMass<<" > GSM="<<GSMass<<G4endl;
#endif
      G4int nOfOUT = theQHadrons.size();    // Total #of QHadrons in Vector at this point
      while(nOfOUT && corFlag)                 // Try BackFusionDecays untill something is in Vector
	  {
        G4QHadron*     theLast = theQHadrons[nOfOUT-1];
        G4int          lastBN = theLast->GetBaryonNumber();
        G4int          nFragm = theLast->GetNFragments();
        G4int          gam    = theLast->GetPDGCode();
#ifdef debug
		G4cout<<"G4QEnv::EvaRes:*BackFusion* lBN="<<lastBN<<",lnF="<<nFragm<<",nH="<<nOfOUT<<G4endl;
#endif
		while(nFragm)                             // => "Delete Decayed Hadrons" case
		{
          G4QHadron* thePrev = theQHadrons[nOfOUT-2];
          nFragm = thePrev->GetNFragments();
#ifdef debug
          G4int          prevPDG = thePrev->GetPDGCode();
		  G4cout<<"G4QEnv::EvaRes:DelTheLast, nFr="<<nFragm<<", pPDG="<<prevPDG<<G4endl;
#endif
          theQHadrons.pop_back();            // the prev QHadron is excluded from OUTPUT
          delete theLast; //!!When killing, DON'T forget to delete the last QHadron as an instance!!
          theLast = thePrev;                   // Update the Last pointer (Prev instead of Last)
		  nOfOUT--;
		}
        if(nOfOUT)
		{
          if(lastBN<1&&nOfOUT>1)               // => "Skip Meson/Gams & Antibaryons" case @@ A few ?
		  {
            G4QHadron* thePrev = theQHadrons[nOfOUT-2];
            theQHadrons.pop_back();          // the last QHadron is excluded from OUTPUT
            theQHadrons.pop_back();          // the prev QHadron is excluded from OUTPUT
            theQHadrons.push_back(theLast);       // the Last becomes the Prev (first part of exchange)
            theQHadrons.push_back(thePrev);       // the Prev becomes the Last (second part of exch.)
            theLast = thePrev;                 // Update the Last pointer (Prev instead of Last)
		  }
          G4LorentzVector last4M = theLast->Get4Momentum();
          G4QContent  lastQC = theLast->GetQC();
          G4double lastM  = last4M.m();        // Mass of the BackFused Fragment
          totQC+=lastQC;                       // Update (increase) the total QC
          q4M+=last4M;                         // Update (increase) the total 4-momentum
          totMass=q4M.m();                     // Calculate new real total mass
          G4int totPDG=totQC.GetSPDGCode();    // The updated PDG for the Total Residual Nucleus
          G4int totBN=totQC.GetBaryonNumber(); // Baryon number of the Total Residual Nucleus
          G4double dM=totMass-GSMass -lastM;
#ifdef debug
		  G4cout<<"G4QE::EvR:TM="<<totMass<<"-LM="<<lastM<<lastQC<<"-GSM="<<GSMass<<"="<<dM<<G4endl;
#endif
          if(dM>-0.001)
		  {
            G4QHadron* evH = new G4QHadron(totPDG,q4M);// Create QHadron for the TotalResidNucleus
            if(dM<=0.)
            {
              theQHadrons.pop_back();        // lastQHadron is excluded from QHadrV as it's in TRN
              if(totBN==2)DecayDibaryon(evH);  // Fill dibaryon (with decay products)
              else theQHadrons.push_back(evH);    // Just Fill TRN to HVect as it's (delete equivalent)
		    }
            else                               // Decay TotalResidualNucleus in GSM + Last and Break
		    {
              G4LorentzVector r4Mom(0.,0.,0.,GSMass);
              if(!G4QHadron(q4M).DecayIn2(last4M,r4Mom))
              {
                theQHadrons.pop_back();      // lastQHadron is excluded from QHadrV as it's in TRN
                delete theLast; //!When killing, DON'T forget to delete last QHadron as an instance!
                theQHadrons.push_back(evH);       // Just Fill TRN to Vect as it is (delete equivalent)
#ifdef debug
                G4cout<<"***G4QEnv::EvaRes: DecayIn L"<<lastQC<<"+TRN"<<totQC<<" failed"<<G4endl;
#endif
	          }
              else
              {
                delete evH;                    // Delete the Hadron for the Total Residual Nucleus
                theLast->Set4Momentum(last4M); // Already exists: don't create&fill, just set 4Mom
                G4QHadron* nuclH = new G4QHadron(thePDG,r4Mom);
                if(thePDG==92000000||thePDG==90002000||thePDG==90000002)DecayDibaryon(nuclH);//DelEq
                else theQHadrons.push_back(nuclH);// Fill the Residual Nucleus (delete equivalent)
              }
              break;
		    }
		  }
          thePDG=totPDG;                       // Make a Residual Nucleus out of the TotResidNucl
		  GSMass=G4QPDGCode(thePDG).GetMass(); // Update the Total Residual Nucleus mass
          theQHadrons.pop_back();            // the last QHadron is excluded from OUTPUT
          delete theLast;//!! When killing, DON'T forget to delete the last QHadron as an instance!!
          nOfOUT--;                            // Update the value of OUTPUT entries
		}
	  }
      if(!nOfOUT || !corFlag)
      {
        G4LorentzVector h4Mom(0.,0.,0.,GSMass);// GSMass should be updated in previous while-LOOP
        G4LorentzVector g4Mom(0.,0.,0.,0.);
        if(!G4QHadron(q4M).DecayIn2(h4Mom, g4Mom))
        {
          G4cerr<<"***G4QEnv::EvaRes: h="<<thePDG<<"(GSM="<<GSMass<<")+gamma>tM="<<totMass<<G4endl;
          G4Exception("G4QEnvironment::EvaporateResidual: Initial Decay in Gamma failed");
        }
#ifdef debug
	    G4cout<<"G4QEnv::EvaRes: "<<q4M<<"->totResN="<<thePDG<<h4Mom<<" + gamma="<<g4Mom<<G4endl;
#endif
        G4QHadron* curG = new G4QHadron(22,g4Mom);
        theQHadrons.push_back(curG);              // Fill the gamma (delete equivalent)
        G4QHadron* curH = new G4QHadron(thePDG,h4Mom);
        if(thePDG==92000000||thePDG==90002000||thePDG==90000002) DecayDibaryon(curH); //(del.equiv.)
        else theQHadrons.push_back(curH);         // Fill the TotalResidualNucleus (delete equivalent)
	  }
	}
    else if(bA==2) DecayDibaryon(qH);          // Decay the residual dibaryon (delete equivalent)
    else if(bA==3&&(bZ==bA||bN==bA||bS==bA)) DecayThreeBaryon(qH); // Decay 3p,3n,3l (delete equiv.)
    else if(bA==5&&abs(bZ-bN)==1) DecayAlphaBar(qH);// Decay alpha-nucleon state (delete equivalent)
    else if(bZ==4&&bN==2&&!bS) DecayAlphaDiN(qH);  // Decay alpha+2protons state (delete equivalent)
    else if(bZ==4&&bN==4&&!bS) DecayAlphaAlpha(qH); // Decay alpha+alpha state (delete equivalent)
    else if(abs(totMass-GSMass)<.003&&!bsCond&&!dbsCond) theQHadrons.push_back(qH);//FillAsItIs(del.eq)
    else if(totMass<GSMass+.003&&(bsCond||dbsCond)) //==> " M <= GSM but decay is possible" case
    {
#ifdef debug
	  G4cout<<"G4QE::EvaRes:2B="<<dbsCond<<",BDA="<<bsCond<<",M="<<totMass<<"<GSM="<<GSMass<<G4endl;
#endif
      G4double gResM  =1000000.;               // Prototype of mass of residual for a neutron
      G4int    gResPDG=0;                      // Prototype of PDGCode of residual for a neutron
      if(bN==4&&bZ==2&&!bS)                    // It's He6 nucleus
	  {
        G4QContent resQC=totQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        gResPDG= thePDG;                       // PDG of the Residual Nucleus
        gResM  = mHel6;                        // min mass of the Residual Nucleus
	  }
      G4double nResM  =1000000.;               // Prototype of mass of residual for a neutron
      G4int    nResPDG=0;                      // Prototype of PDGCode of residual for a neutron
      if(bsCond&&bN>0&&bA>1)                   // It's nucleus and there is a neutron
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
      if(bsCond&&bZ>0&&bA>1)                   // It's nucleus and there is a proton
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
      if(bsCond&&bS>0&&bA>1)                   // It's nucleus and there is a Lambda
	  {
        G4QContent resQC=totQC-lambQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        lResPDG=resN.GetPDG();                 // PDG of the Residual Nucleus
        if     (lResPDG==90000001) lResM=mNeut;
        else if(lResPDG==90001000) lResM=mProt;
        else if(lResPDG==91000000) lResM=mLamb;
        else lResM  =resN.GetMZNS();           // min mass of the Residual Nucleus
	  }
      G4double dResM  =1000000.;               // Prototype of mass of residual for a Alpha
      G4int    dResPDG=0;                      // Prototype of PDGCode of residual for a Alpha
      if(bsCond&&bN>0&&bZ>0&&bA>2)             // It's nucleus and there is a Deuteron
	  {
        G4QContent resQC=totQC-deutQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        dResPDG=resN.GetPDG();                 // PDG of the Residual Nucleus
        if     (dResPDG==90000001) dResM=mNeut;
        else if(dResPDG==90001000) dResM=mProt;
        else if(dResPDG==91000000) dResM=mLamb;
        else dResM  =resN.GetMZNS();           // min mass of the Residual Nucleus
	  }
      G4double aResM  =1000000.;               // Prototype of mass of residual for a Alpha
      G4int    aResPDG=0;                      // Prototype of PDGCode of residual for a Alpha
      if(bsCond&&bN>1&&bZ>1&&bA>4)             // It's nucleus and there is an Alpha
	  {
        G4QContent resQC=totQC-alphQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        aResPDG=resN.GetPDG();                 // PDG of the Residual Nucleus
        if     (aResPDG==90000001) aResM=mNeut;
        else if(aResPDG==90001000) aResM=mProt;
        else if(aResPDG==91000000) aResM=mLamb;
        else aResM  =resN.GetMZNS();           // min mass of the Residual Nucleus
	  }
      G4double nnResM  =1000000.;              // Prototype of mass of residual for a dineutron
      G4int    nnResPDG=0;                     // Prototype of PDGCode of residual for a dineutron
      if(dbsCond&&bN>1&&bA>2)                  // It's nucleus and there is a dineutron
	  {
        G4QContent resQC=totQC-neutQC-neutQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        nnResPDG=resN.GetPDG();                // PDG of the Residual Nucleus
        if     (nnResPDG==90000001) nnResM=mNeut;
        else if(nnResPDG==90001000) nnResM=mProt;
        else if(nnResPDG==91000000) nnResM=mLamb;
        else nnResM  =resN.GetMZNS();          // min mass of the Residual Nucleus
	  }
      G4double ppResM  =1000000.;              // Prototype of mass of residual for a diproton
      G4int    ppResPDG=0;                     // Prototype of PDGCode of residual for a diproton
      if(dbsCond&&bZ>1&&bA>2)                  // It's nucleus and there is a diproton
	  {
        G4QContent resQC=totQC-protQC-protQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        ppResPDG=resN.GetPDG();                // PDG of the Residual Nucleus
        if     (ppResPDG==90000001) ppResM=mNeut;
        else if(ppResPDG==90001000) ppResM=mProt;
        else if(ppResPDG==91000000) ppResM=mLamb;
        else ppResM  =resN.GetMZNS();          // min mass of the Residual Nucleus
	  }
      G4double npResM  =1000000.;              // Prototype of mass of residual for a proton+neutron
      G4int    npResPDG=0;                     // Prototype of PDGCode of residual for a prot+neut
      if(dbsCond&&bN>0&&bZ>0&&bA>2)            // It's nucleus and there is a proton and a neutron
	  {
        G4QContent resQC=totQC-neutQC-protQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        npResPDG=resN.GetPDG();                // PDG of the Residual Nucleus
        if     (npResPDG==90000001) npResM=mNeut;
        else if(npResPDG==90001000) npResM=mProt;
        else if(npResPDG==91000000) npResM=mLamb;
        else npResM  =resN.GetMZNS();          // min mass of the Residual Nucleus
	  }
      G4double lnResM  =1000000.;              // Prototype of mass of residual for a lambda+neutron
      G4int    lnResPDG=0;                     // Prototype of PDGCode of residual for a lamb+neut
      if(dbsCond&&bN>0&&bS>0&&bA>2)            // It's nucleus and there is a lambda and a neutron
	  {
        G4QContent resQC=totQC-lambQC-protQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        lnResPDG=resN.GetPDG();                // PDG of the Residual Nucleus
        if     (lnResPDG==90000001) lnResM=mNeut;
        else if(lnResPDG==90001000) lnResM=mProt;
        else if(lnResPDG==91000000) lnResM=mLamb;
        else lnResM  =resN.GetMZNS();          // min mass of the Residual Nucleus
	  }
      G4double lpResM  =1000000.;              // Prototype of mass of residual for a proton+lambda
      G4int    lpResPDG=0;                     // Prototype of PDGCode of residual for a prot+lamb
      if(dbsCond&&bS>0&&bZ>0&&bA>2)            // It's nucleus and there is a proton and a lambda
	  {
        G4QContent resQC=totQC-neutQC-protQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        lpResPDG=resN.GetPDG();                // PDG of the Residual Nucleus
        if     (lpResPDG==90000001) lpResM=mNeut;
        else if(lpResPDG==90001000) lpResM=mProt;
        else if(lpResPDG==91000000) lpResM=mLamb;
        else lpResM  =resN.GetMZNS();          // min mass of the Residual Nucleus
	  }
      G4double llResM  =1000000.;              // Prototype of mass of residual for a di-lambda
      G4int    llResPDG=0;                     // Prototype of PDGCode of residual for a di-lambda
      if(dbsCond&&bS>1&&bA>2)                  // It's nucleus and there is a di-lambda
	  {
        G4QContent resQC=totQC-neutQC-protQC;
        G4QNucleus resN(resQC);                // Pseudo nucleus for the Residual Nucleus
        llResPDG=resN.GetPDG();                // PDG of the Residual Nucleus
        if     (llResPDG==90000001) llResM=mNeut;
        else if(llResPDG==90001000) llResM=mProt;
        else if(llResPDG==91000000) llResM=mLamb;
        else llResM  =resN.GetMZNS();          // min mass of the Residual Nucleus
	  }
#ifdef debug
      G4cout<<"G4QEnv::EvaRes: rP="<<pResPDG<<",rN="<<nResPDG<<",rL="<<lResPDG<<",N="
            <<bN<<",Z="<<bZ<<",nL="<<bS<<",totM="<<totMass<<",n="<<totMass-nResM-mNeut
            <<",p="<<totMass-pResM-mProt<<",l="<<totMass-lResM-mLamb<<G4endl;
#endif
      if(   thePDG==90004004 || thePDG==90002004 && totMass>mHel6+.003
         || bA>4 && bsCond && bN>1 && bZ>1 && totMass>aResM+mAlph
         || bA>1 && bsCond && (   bN>0&&totMass>nResM+mNeut
                               || bZ>0&&totMass>pResM+mProt
							   || bS>0&&totMass>lResM+mLamb)
         || bA>2 && (bN>0&&bZ>0 && (bsCond && totMass>dResM+mDeut || dbsCond && totMass>dResM+mDeut)
				     || dbsCond && (   bN>1&&totMass>nnResM+mNeut+mNeut
                                    || bZ>1&&totMass>ppResM+mProt+mProt
                                    || bS>1&&totMass>llResM+mLamb+mLamb
                                    || bN&&bS&&totMass>lnResM+mLamb+mNeut
                                    || bZ&&bS&&totMass>lpResM+mLamb+mProt)))
	  {
        G4int barPDG = 90002002;               // Just for the default case of Be8->alpha+alpha
        G4int resPDG = 90002002;
        G4int thdPDG = 0;
        G4double barM= mAlph;
        G4double resM= mAlph;
        G4double thdM= mNeut;                  // This default value is used in the IF
		if(gResPDG&&totMass>mHel6+.003)        // Can make radiative decay of He6 (priority 0)
		{
          barPDG=90002004;
          resPDG=22;
          barM  =mHel6;
          resM  =0.;
		}
		else if(nResPDG&&totMass>nResM+mNeut)  // Can radiate a neutron (priority 1)
		{
          barPDG=90000001;
          resPDG=nResPDG;
          barM  =mNeut;
          resM  =nResM;
		}
		else if(pResPDG&&totMass>pResM+mProt)  // Can radiate a proton (priority 2)
		{
          barPDG=90001000;
          resPDG=pResPDG;
          barM  =mProt;
          resM  =pResM;
		}
		else if(lResPDG&&totMass>lResM+mLamb)  // Can radiate a Lambda (priority 3)
		{
          barPDG=91000000;
          resPDG=lResPDG;
          barM  =mLamb;
          resM  =lResM;
		}
		else if(ppResPDG&&totMass>ppResM+mProt+mProt)// Can radiate a DiProton (priority 4)
		{
          barPDG=90001000;
          resPDG=ppResPDG;
          thdPDG=90001000;
          barM  =mProt;
          resM  =ppResM;
          thdM  =mProt;
		}
		else if(nnResPDG&&totMass>nnResM+mNeut+mNeut)// Can radiate a DiNeutron (priority 5)
		{
          barPDG=90000001;
          resPDG=nnResPDG;
          thdPDG=90000001;
          barM  =mNeut;
          resM  =nnResM;
		}
		else if(npResPDG&&totMass>npResM+mProt+mNeut)// Can radiate a neutron+proton (priority 6)
		{
          barPDG=90001000;
          resPDG=npResPDG;
          thdPDG=90000001;
          barM  =mProt;
          resM  =npResM;
		}
		else if(lnResPDG&&totMass>lnResM+mLamb+mNeut)// Can radiate a Lambda+neutron (priority 7)
		{
          barPDG=91000000;
          resPDG=lnResPDG;
          thdPDG=90000001;
          barM  =mLamb;
          resM  =lnResM;
		}
		else if(lpResPDG&&totMass>lpResM+mLamb+mProt)// Can radiate a Lambda+proton (priority 8)
		{
          barPDG=91000000;
          resPDG=lpResPDG;
          thdPDG=90001000;
          barM  =mLamb;
          resM  =lpResM;
          thdM  =mProt;
		}
		else if(llResPDG&&totMass>llResM+mLamb+mLamb)// Can radiate a DiLambda (priority 9)
		{
          barPDG=91000000;
          resPDG=llResPDG;
          thdPDG=91000000;
          barM  =mLamb;
          resM  =llResM;
          thdM  =mLamb;
		}
        else if(thePDG!=90004004&&bN>1&&bZ>1&&bA>4&&totMass>aResM+mAlph) // Try to decay in alpha
		{
          barPDG=90002002;
          resPDG=aResPDG;
          barM  =mAlph;
          resM  =aResM;
		}
		else if(dResPDG&&totMass>dResM+mDeut)  // Can radiate a Deuteron (lowest priority)
		{
          barPDG=90001001;
          resPDG=dResPDG;
          barM  =mDeut;
          resM  =dResM;
		}
        else if(thePDG!=90004004&&totMass>GSMass)// If it's not Be8 decay in gamma & GSM
		{
          barPDG=22;
          resPDG=thePDG;
          barM  =0.;
          resM  =GSMass;
		}
        else if(thePDG!=90004004)
		{
          G4cerr<<"***G4QEnv::EvaRes:PDG="<<thePDG<<",M="<<totMass<<"< GSM="<<GSMass<<G4endl;
          G4Exception("***G4QEnvironment::EvaporateResidual: M<GSM & can't decay in p,n,l,d,alpha");
		}
        G4LorentzVector a4Mom(0.,0.,0.,barM);
        G4LorentzVector b4Mom(0.,0.,0.,resM);
        if(!thdPDG)
        {
          if(!qH->DecayIn2(a4Mom,b4Mom))
          {
            theQHadrons.push_back(qH);              // Fill as it is (delete equivalent)
            G4cout<<"G4QEnv::EvaRes: rP="<<pResPDG<<",rN="<<nResPDG<<",rL="<<lResPDG<<",N="
                  <<bN<<",Z="<<bZ<<",L="<<bS<<",totM="<<totMass<<",n="<<totMass-nResM-mNeut
                  <<",p="<<totMass-pResM-mProt<<",l="<<totMass-lResM-mLamb<<G4endl;
            G4cerr<<"***G4QE::EvaporResid:DecayIn2 failed bPDG="<<barPDG<<",rPDG="<<resPDG<<G4endl;
	      }
          else
          {
            delete qH;
            G4QHadron* HadrB = new G4QHadron(barPDG,a4Mom);
            theQHadrons.push_back(HadrB);           // Fill the baryon (delete equivalent)
            G4QHadron* HadrR = new G4QHadron(resPDG,b4Mom);
            // @@ Self-call !!
            if(HadrR->GetBaryonNumber()>1) EvaporateResidual(HadrR);//Continue decay (delete equiv.)
            else theQHadrons.push_back(HadrR);      // Fill ResidNucleus=Baryon to Output HadronVector
          }
        }
        else
        {
          G4LorentzVector c4Mom(0.,0.,0.,thdM);
          if(!qH->DecayIn3(a4Mom,b4Mom,c4Mom))
          {
            theQHadrons.push_back(qH);              // Fill as it is (delete equivalent)
            G4cout<<"G4QEnv::EvaRes:rNN="<<nnResPDG<<",rNP="<<npResPDG<<",rPP="<<ppResPDG<<",N="
                  <<bN<<",Z="<<bZ<<",L="<<bS<<",totM="<<totMass<<",nn="<<totMass-nnResM-mNeut-mNeut
                  <<",np="<<totMass-npResM-mProt-mNeut<<",pp="<<totMass-ppResM-mProt-mProt<<G4endl;
            G4cerr<<"***G4QE::EvaporResid:DecayIn2 failed bPDG="<<barPDG<<",rPDG="<<resPDG<<G4endl;
	      }
          else
          {
            delete qH;
            G4QHadron* HadrB = new G4QHadron(barPDG,a4Mom);
            theQHadrons.push_back(HadrB);           // Fill the first baryon (delete equivalent)
            G4QHadron* HadrC = new G4QHadron(thdPDG,c4Mom);
            theQHadrons.push_back(HadrC);           // Fill the second baryon (delete equivalent)
            G4QHadron* HadrR = new G4QHadron(resPDG,b4Mom);
            // @@ Self-call !!
            if(HadrR->GetBaryonNumber()>1) EvaporateResidual(HadrR);//Continue decay (delete equiv.)
            else theQHadrons.push_back(HadrR);      // Fill ResidNucleus=Baryon to Output HadronVector
          }
        }
	  }
      else if (abs(totMass-GSMass)<.003) theQHadrons.push_back(qH);// Fill as it is (delete equivalent)
      else                                     // "System is below mass shell and can't decay" case
	  {
        G4cerr<<"***G4QEnv::EvaRes: M="<<totMass<<" < GSM="<<GSMass<<", d="<<totMass-GSMass<<G4endl;
        G4Exception("***G4QEnv::EvaporateResidual: Nuclear Mass is below the Ground State value");
	  }
    }
    else                                       // ===> Evaporation of excited system
	{
#ifdef debug
      G4cout<<"G4QEnv::EvaRes: ***EVA*** tPDG="<<thePDG<<", tM="<<totMass<<" > GS="<<GSMass
            <<qNuc.Get4Momentum()<<G4endl;
#endif
      G4LorentzVector b4M;
      G4LorentzVector r4M;
      G4bool evC=true;
      G4bool dcC=false;
      G4int bPDG=0;
      G4int rPDG=0;
      G4double bM   = 0.;                      // Prototype of Real Mass of the EvaporatedDibaryon
      G4double rM   = 0.;                      // Prototype of Real Mass of the residual nucleus
      G4int bB=0;                              // Proto of Baryon Number of the evaporated baryon
      G4int rB=0;                              // Proto of Baryon Number of the residual nucleus
      G4QHadron* bHadron = new G4QHadron;
      G4QHadron* rHadron = new G4QHadron;
      G4int evcn=0;
      //G4int evcm=27;
      G4int evcm=9;
      while(evC&&evcn<evcm)
      {
        evC=true;
        dcC=false;
        evcn++;
        if(!qNuc.EvaporateBaryon(bHadron,rHadron))
	    {
          G4cerr<<"***G4QEnv::EvaRes: ***EVA*** PDG="<<thePDG<<",tM="<<totMass<<G4endl;
          delete bHadron;
          delete rHadron;
          theQHadrons.push_back(qH);              // Fill hadron in the HadronVector as it is
          evC=false;
          dcC=true;
          //@@Temporary
          G4Exception("G4QEnvironment::EvaporateResidual: Failed to evaporate baryon/alpha");
	    }
        if(!dcC)
        {
          evC=false;
          b4M=bHadron->Get4Momentum();
          r4M=rHadron->Get4Momentum();
          bM   = b4M.m();                      // Real mass of the evaporated dibaryon
          rM   = r4M.m();                      // Real mass of the residual nucleus
          bB=bHadron->GetBaryonNumber();       // Baryon number of the evaporated baryon
          rB=rHadron->GetBaryonNumber();       // Baryon number of the residual nucleus
          G4int bC=bHadron->GetCharge();       // Baryon number of the evaporated baryon
          G4int rC=rHadron->GetCharge();       // Baryon number of the residual nucleus
          G4double bCB=qNuc.CoulombBarrier(bC,bB);
          G4double rCB=qNuc.CoulombBarrier(rC,rB);
          bPDG=bHadron->GetPDGCode();
          rPDG=rHadron->GetPDGCode();
#ifdef debug
          G4cout<<"G4QEnv::EvaRes: Attempt #"<<evcn<<" > "<<evcm<<", rPDG="<<rPDG<<", bPDG="
                <<bPDG<<", bE="<<b4M.e()-b4M.m()<<" > bCB="<<bCB<<G4endl;
#endif
          //if(b4M.e()-b4M.m()<bCB&&evcn<evcm) evC=true;
		}
	  }
      if(!dcC)
      {
#ifdef debug
        G4cout<<"G4QEnv::EvaRes: *** EVA DONE *** EvaporatedFragment="<<bPDG<<b4M<<",bB="<<bB
              <<", ResNuclPDG="<<rPDG<<r4M<<",rB="<<rB<<G4endl;
#endif
        delete qH;
        if(bB<2)theQHadrons.push_back(bHadron);   // Fill Evaporated Baryon (delete equivalent)
        else if(bB==2) DecayDibaryon(bHadron); // => "Dibaryon" case needs decay
        else if(bB==4) theQHadrons.push_back(bHadron); // "Alpha radiation" case (delete equival.)
        else if(bB==5) DecayAlphaBar(bHadron);      // "Alpha+Baryon Decay" case (delete equival.)
        else if(bPDG==90004004) DecayAlphaAlpha(bHadron);// "Alpha+Alpha Decay" case (delete equiv.)
        else
	    {
          G4cerr<<"***G4QEnv::EvaRes: bB="<<bB<<" > 2 - unexpected evaporated fragment"<<G4endl;
          G4Exception("***G4QEnvironment::EvaporateResidual: Unexpected evaporation act");
	    }
        if(rB>2) EvaporateResidual(rHadron);   // Continue evaporation (@@ Self-call)
        else if(rB==2)                         // => "Dibaryon" case needs decay
	    {
          G4double rGSM = rHadron->GetQPDG().GetMass();// Ground State mass of the dibaryon
#ifdef debug
		  G4cout<<"G4QEnv::EvaRes: ResidDibarionM="<<rM<<",GSM="<<rGSM<<", M-GSM="<<rM-rGSM<<G4endl;
#endif

          if(rM<=rGSM-0.001)
		  {
            G4cerr<<"***G4QEnv::EvaRes: <residual> M="<<rM<<" < GSM="<<rGSM<<G4endl;
            G4Exception("***G4QEnvironment::EvaporateResidual:Evaporation is below the mass shell");
		  }
          else if(abs(rM-rGSM)<0.001&&rPDG==90001001)theQHadrons.push_back(rHadron);//(delete equival.)
          else DecayDibaryon(rHadron);         // => "Dibaryon Decay" case (delete equivalent)
	    }
        else if(rB==5) DecayAlphaBar(rHadron);      // "Alpha+Baryon Decay" case (delete equival.)
        else if(rPDG==90004004) DecayAlphaAlpha(rHadron);// "Alpha+Alpha Decay" case (delete equiv.)
        else theQHadrons.push_back(rHadron);      // Fill ResidNucleus=Baryon to Output HadronVector
	  } // End of fail in decay check
	} // End of Evaporation of excited system
  }
  else                                         // => "Decay if it is impossible to evaporate" case
  {
#ifdef debug
    G4cout<<"G4QEnv::EvaRes: InputHadron4M="<<q4M<<", PDG="<<thePDG<<G4endl;
#endif
    if(thePDG)
    {
      if(thePDG==10)                           // "Chipolino decay" case 
	  {
        G4QContent totQC = qH->GetQC();        // Quark content of the hadron
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
            G4cerr<<"***G4QEnv::EvaRes: tM="<<totMass<<"-> h1M="<<m1<<" + h2M="<<m2<<G4endl;
		    G4Exception("G4QEnvironment::EvaporateResidual: Chip->h1+h2 DecayIn2 did not succeed");
	      }
          delete qH;
          G4QHadron* H2 = new G4QHadron(h2.GetPDGCode(),qe4M);
          theQHadrons.push_back(H2);              // (delete equivalent)
          G4QHadron* H1 = new G4QHadron(h1.GetPDGCode(),fq4M);
          theQHadrons.push_back(H1);              // (delete equivalent)
		}
        else
	    {
          G4cerr<<"***G4QEnv::EvaRes: M="<<totMass<<"<"<<m1<<"+"<<m2<<", d="<<m1+m2-totMass<<G4endl;
          G4Exception("***G4QEnvironment::EvaporateResidual: Chipolino is under a Mass Shell");
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
          delete qH;
          G4QHadron* h1H = new G4QHadron(nucPDG,n14M);
          theQHadrons.push_back(h1H);           // (delete equivalent)
          G4QHadron* h2H = new G4QHadron(nucPDG,n24M);
          theQHadrons.push_back(h2H);           // (delete equivalent)
          G4QHadron* piH = new G4QHadron(piPDG,pi4M);
          theQHadrons.push_back(piH);           // (delete equivalent)
		}
		else
	    {
          G4cerr<<"***G4QEnv::EvaRes: IdPDG="<<thePDG<<", q4M="<<q4M<<", M="<<totMass
                <<" < M_2N+Pi, d="<<totMass-2*nucM-mPi<<G4endl;
          G4Exception("***G4QEnvironment::EvaporateResidual: ISO-dibaryon is under the Mass Shell");
	    }
	  }
      else                                     // "Hadron" case
	  {
        G4double totM=G4QPDGCode(thePDG).GetMass();
        if(abs(totMass-totM)<0.001||abs(thePDG)-10*(abs(thePDG)/10)>2)theQHadrons.push_back(qH);
        else if (totMass>totM)                 // "Radiative Hadron decay" case
	    {
          G4LorentzVector fq4M(0.,0.,0.,0.);
          G4LorentzVector qe4M(0.,0.,0.,totM);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
		  {
            G4cerr<<"***G4QEnv::HadQEnv: tM="<<totMass<<"-> h1M="<<totM<<" + gamma"<<G4endl;
		    G4Exception("G4QEnv::HadrQEnv: H*->H+gamma DecayIn2 did not succeed");
	      }
          delete qH;
          G4QHadron* H1 = new G4QHadron(22,fq4M);
          theQHadrons.push_back(H1);              // (delete equivalent)
          G4QHadron* H2 = new G4QHadron(thePDG,qe4M);
          theQHadrons.push_back(H2);              // (delete equivalent)
	    }
        else
	    {
          G4cerr<<"***G4QEnv::EvaRes: ResNuc="<<thePDG<<theQC<<", q4M="<<q4M<<", M="<<totMass
                <<" < GSM="<<totM<<G4endl;
          G4Exception("***G4QEnvironment::EvaporateResidual: Hadron is under the Mass Shell");
	    }
	  }
	}
    else
    {
      G4cerr<<"***G4QEnv::EvaRes: ResNuc="<<thePDG<<theQC<<", q4M="<<q4M<<", qM="<<totMass<<G4endl;
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
  static const G4QPDGCode gQPDG(22);
  static const G4QPDGCode nQPDG(2112);
  static const G4QPDGCode pQPDG(2212);
  static const G4QPDGCode lQPDG(3122);
  static const G4QPDGCode aQPDG(90002002);
  static const G4QPDGCode a6QPDG(90002004);
  static const G4QPDGCode be6QPDG(90004002);
  static const G4QPDGCode b7QPDG(90005002);
  static const G4QPDGCode a8QPDG(90002006);
  static const G4QPDGCode c10QPDG(90006004);
  static const G4QPDGCode o14QPDG(90008006);
  static const G4QPDGCode o15QPDG(90008007);
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mHe6 = G4QPDGCode(2112).GetNuclMass(2,4,0);
  static const G4double mBe6 = G4QPDGCode(2112).GetNuclMass(4,2,0);
  static const G4double mB7  = G4QPDGCode(2112).GetNuclMass(5,2,0);
  static const G4double mHe8 = G4QPDGCode(2112).GetNuclMass(2,6,0);
  static const G4double mC10 = G4QPDGCode(2112).GetNuclMass(6,4,0);
  static const G4double mO14 = G4QPDGCode(2112).GetNuclMass(8,6,0);
  static const G4double mO15 = G4QPDGCode(2112).GetNuclMass(8,7,0);
  static const G4double third=1./3.;
  static const G4double nPDG=90000001;
  G4int envA=theEnvironment.GetBaryonNumber();
  G4int envC=theEnvironment.GetCharge();
#ifdef pdebug
  G4cout<<"G4QEnvironment(G4QE)::Fragment(Fr): ***called*** envA="<<envA<<G4endl;
#endif
  HadronizeQEnvironment();
  G4int nHadr=theQHadrons.size();
#ifdef pdebug
  G4cout<<"G4QEnvironment::Fragment: after HadronizeQEnvironment #ofHadrons="<<nHadr<<G4endl;
#endif
  G4int lHadr=theQHadrons[nHadr-1]->GetBaryonNumber();
  if(lHadr>1)                             // The last hadron is a nucleus: try to decay/evaporate it
  {
    G4QHadron* theLast = theQHadrons[nHadr-1];
    G4QHadron* curHadr = new G4QHadron(theLast);
#ifdef pdebug
    G4cout<<"G4QE::Fr:Before nH="<<nHadr<<", PDG="<<curHadr->GetPDGCode()<<G4endl;
#endif
    theQHadrons.pop_back();                    // the last QHadron-Nucleus is excluded from OUTPUT
    delete theLast; //**!! When killing, DON'T forget to delete the last QHadron as an instance !!**
    EvaporateResidual(curHadr);                  // Try to evaporate Hadron-Nucleus (delete equiv.)
    nHadr=theQHadrons.size();
#ifdef pdebug
    G4cout<<"G4QE::Fr:After nH="<<nHadr<<G4endl;
#endif
  }
  if(nHadr)for(G4int ipo=0; ipo<theQHadrons.size(); ipo++)// Find NuclFragm and try to decay/evap
  {
    G4int hBN  = theQHadrons[ipo]->GetBaryonNumber();
    G4int hPDG = theQHadrons[ipo]->GetPDGCode();
#ifdef pdebug
    G4int hNF  = theQHadrons[ipo]->GetNFragments();
    if(hBN==1)G4cout<<"G4QE::Fr:h#"<<ipo<<": hPDG="<<hPDG<<", hNFrag="<<hNF<<G4endl;
#endif
    if(hBN>lHadr) lHadr=hBN;                     // Now lHadr is the biggest nucleus
    G4LorentzVector h4Mom = theQHadrons[ipo]->Get4Momentum();
    G4double hE=h4Mom.e();
    if(hBN>1 || hPDG==22 && hE<.001)             // Nucleus is found: Swap with the last & evap
	{
      G4QHadron* theLast = theQHadrons[nHadr-1];
      G4QHadron* theCurr = theQHadrons[ipo];
      G4QHadron* curHadr = new G4QHadron(theCurr);// Remember the current hadron
      theCurr->SetQPDG(theLast->GetQPDG());      // the formal current is substituted by the last
      theCurr->Set4Momentum(theLast->Get4Momentum());
      theQHadrons.pop_back();                  // the last QHadron is excluded from OUTPUT
      delete theLast; //*!! When killing, DON'T forget to delete the last QHadron as an instance !!*
      if(hPDG!=22) EvaporateResidual(curHadr);   // Try to evaporate the Nucleus (delete equiv)
      nHadr=theQHadrons.size();
	}
#ifdef pdebug
    G4int           hNFrag= theQHadrons[ipo]->GetNFragments();
    G4QContent      hQC   = theQHadrons[ipo]->GetQC();
    G4cout<<"G4QE::Fr:h#"<<ipo<<": hPDG="<<hPDG<<hQC<<", h4M="<<h4Mom<<", hNFrag="<<hNFrag<<G4endl;
#endif
  }
  G4double apa=1.;
  if(envA>5) apa=pow(envA-4.,third); 
  G4double apa2=apa*apa;
  //G4double cut=8.;
  //if(envA>3) cut=27./envA;
  //G4double cut=40./envA;
  //////G4double cut=20./sqrt(envA);
  //G4double cut=2.7;
  //G4double cut=4.;
  //G4double cut=1.;
  //G4double p2cut=32400.; //(200/1.1)
  //G4double p2cut=40000.; //(200 fm*MeV)
  //G4double p2cut=24800.; //(200/1.27)
  //G4double p2cut=20000.; //(200/1.41)
  G4double p2cut=17543.; //(200/1.51) // G4 estimate
  //G4double p2cut=10000.;
  //if(envA>1) p2cut/=envA;
  //if(envA>5) p2cut*=pow(envA-4.,-.66666666);
  if(envA>5) p2cut/=apa2;
  //if(envA>1) p2cut*=pow(envA,-.3333333);
  //if(envA>1) p2cut/=sqrt(envA);
  //if(envA>5) p2cut*=pow(envA-4,-1.3333333);
  //if(envA>5&&nHadr>2) p2cut*=pow((envA-4)/sqrt(nHadr-2),-1.3333333);
  //G4double p2cut=30000.*pow(envA,-1.3333333);
  //G4double p2cut=40000.*pow(envA,-1.3333333);
  //G4double p2cut=60000.*pow(envA,-1.3333333);
  //
  //G4double p2cut2=0.; // No cut for only two baryons
  //G4double p2cut2=0.; // No cut for only two baryons
  //if(envA>4) p2cut2=2700.*apa2;
  G4double p2cut2=2700.; //cut for only two baryons
  if(envA<16) p2cut2=0.;
  //
  G4int bfCountM=1;
  if(envA>8) bfCountM=(envA-1)/4;
  G4bool bfAct = true;
  G4int bfCount= 0;
  G4LorentzVector tmp4Mom=tot4Mom;
  while(bfAct&&bfCount<bfCountM)
  {
    tot4Mom=tmp4Mom;
    bfAct=false;
    bfCount++;
    if(nHadr) for(G4int hadron=0; hadron<theQHadrons.size(); hadron++) // Back Fusion LOOP
    {
#ifdef pdebug
      G4cout<<"G4QE::Fr:LOOP START,h="<<hadron<<", PDG="<<theQHadrons[hadron]->GetPDGCode()<<G4endl;
#endif
      G4QHadron* curHadr = theQHadrons[hadron];    // Get a pointer to the current Hadron
      G4int         hPDG = curHadr->GetPDGCode();
      G4int           hS = curHadr->GetStrangeness();
      G4int           hF = curHadr->GetNFragments();
      G4LorentzVector h4m= curHadr->Get4Momentum();
      G4double hM        = h4m.m();                // Mass of the first fragment
      G4int hB           = curHadr->GetBaryonNumber();
      G4int hC           = curHadr->GetCharge();
      if(!hF&&(hPDG>80000000&&hPDG<90000000||hPDG==90000000||
               hPDG>90000000&&(hPDG%1000000>200000||hPDG%1000>300)))
        G4cerr<<"***G4QEnv::Fragment: PDG("<<hadron<<")="<<hPDG<<", M="<<hM<<G4endl;
#ifdef ppdebug
      G4cout<<"G4QE::Fr:>>>h="<<hPDG<<",hS="<<hS<<",hB="<<hB<<",h#"<<hadron<<"<nH="<<nHadr<<G4endl;
#endif
	  if(hadron&&!hF&&hB>0&&!hS&&(nHadr>3||hB<2)) // ThermonuclBackFusion condition (VIMP for gamA TotCS)
	  //if(2>3)                           // Close the ThermonuclearBackFusion (VIMP for gamA TotCS)
      {
#ifdef ppdebug
        //if(nHadr==3)
          G4cout<<"G4QE::Fr:h="<<hPDG<<",B="<<hB<<",h#"<<hadron<<" < nH="<<nHadr<<G4endl;
#endif
        G4QContent hQC = curHadr->GetQC();
        if(hadron&&!hF&&hB>0) for(G4int pt=0; pt<hadron; pt++)
	    {
          G4QHadron* backH = theQHadrons[pt]; // just get a pointer to one of the previous hadrons
          G4int bPDG = backH->GetPDGCode();
          G4int   bF = backH->GetNFragments();
          G4LorentzVector b4m= backH->Get4Momentum();
          G4double bM= b4m.m();               // Mass of the second fragment
          G4QContent bQC = backH->GetQC();
          G4int   bB = backH->GetBaryonNumber();
          G4int   bS = backH->GetStrangeness();
          G4int   bC = backH->GetCharge();
          G4QContent sQC=bQC+hQC;
          G4int sPDG=sQC.GetZNSPDGCode();
          G4QPDGCode sQPDG(sPDG);
          G4double tM=sQPDG.GetMass();
          G4LorentzVector s4M=h4m+b4m;
          G4double sM2=s4M.m2();
          G4double sM=sqrt(sM2);
          G4double dsM2=sM2+sM2;
          G4double rm=bM-hM;
          G4double sm=bM+hM;
          G4double pCM2=(sM2-rm*rm)*(sM2-sm*sm)/(dsM2+dsM2);
#ifdef ppdebug
          //if(nHadr==3)
		  G4cout<<"G4QE::F:"<<pt<<",B="<<bB<<",S="<<bS<<",p="<<pCM2<<"<"<<p2cut<<",hB="<<hB
                <<",bM+hM="<<bM+hM<<">tM="<<tM<<",tQC="<<sQC<<G4endl;
#endif
          //if(!bF&&(bB==1||hB==1)&&bM+hM>tM+.001&&pCM2<p2cut)           // Only baryons == pcut
		  if(!bF&&!bS&&(bB==1&&hB>0||hB==1&&bB>0)&&bM+hM>tM+.001
             && (pCM2<p2cut&&nHadr>3||pCM2<p2cut2&&nHadr==3))
		  //if(!bF&&(bB==1||hB==1)&&bM+hM>tM+.001&&(pCM2<p2cut&&nHadr>3 ||
		  //   pCM2<p2cut2&&nHadr==3&&bPDG>90000000))
          //if(!bF&&(bB<3||hB<3)&&bM+hM>tM+.001&&pCM2<p2cut)           // Only baryons == pcut
		  //if(!bF&&(bB==1||hB==1)&&(nHadr>3||bPDG>90000000)&&bM+hM>tM+.001&&pCM2<p2cut)//Bar+n=3,NN
		  //if(!bF&&(bB==1&&!bC||hB==1&&!hC)&&bM+hM>tM+.001&&pCM2<p2cut) // Only neuterons == pcut
		  //if(!bF&&(bB==1||hB==1)&&bM+hM>tM+.001&&sM-bM-hM<cut)         // Only baryons == ecut
		  //if(!bF&&bB&&bB<fL&&bM+hM>tM+.001&&sM-bM-hM<cut)              // Light fragments == ecut
	      {
#ifdef ppdebug
            //if(nHadr==3)
	        G4cout<<"G4QE::Fr:*FUSION*#"<<hadron<<"["<<hPDG<<"]"<<hM<<"+#"<<pt<<"["<<bPDG<<"]"<<bM
                  <<"="<<bM+hM<<", sM="<<sM<<">["<<sPDG<<"]"<<tM<<",p2="<<pCM2<<"<"<<p2cut<<G4endl;
            //    <<"="<<bM+hM<<", sM="<<sM<<">["<<sPDG<<"]"<<tM<<",t="<<sM-bM-hM<<"<"<<cut<<G4endl;
#endif
            bfAct=true;
            //@@Temporary decay in gamma
            G4bool three=false;
            G4QPDGCode fQPDG=sQPDG;
            G4QPDGCode rQPDG=gQPDG;
            G4QPDGCode tQPDG=gQPDG;
            G4LorentzVector f4Mom(0.,0.,0.,tM);
            G4LorentzVector g4Mom(0.,0.,0.,0.);
            G4LorentzVector t4Mom(0.,0.,0.,0.);
            if     (sPDG==90000002)
		    {
              fQPDG=nQPDG;
              rQPDG=nQPDG;
              f4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              g4Mom=f4Mom;
		    }
            else if(sPDG==90002000)
		    {
              fQPDG=pQPDG;
              rQPDG=pQPDG;
              f4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=f4Mom;
		    }
            else if(sPDG==92000000)
		    {
              fQPDG=lQPDG;
              rQPDG=lQPDG;
              f4Mom=G4LorentzVector(0.,0.,0.,mLamb);
              g4Mom=f4Mom;
		    }
            else if(sPDG==90002003)
		    {
              rQPDG=nQPDG;
              fQPDG=aQPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              f4Mom=G4LorentzVector(0.,0.,0.,mAlph);
		    }
            else if(sPDG==90003002)
		    {
              rQPDG=pQPDG;
              fQPDG=aQPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mAlph);
		    }
            else if(sPDG==90004002)
		    {
              tQPDG=pQPDG;
              rQPDG=pQPDG;
              fQPDG=aQPDG;
              t4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mAlph);
              three=true;
		    }
            else if(sPDG==90002005)
		    {
              rQPDG=nQPDG;
              fQPDG=a6QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              f4Mom=G4LorentzVector(0.,0.,0.,mHe6);
		    }
            else if(sPDG==90005002)
		    {
              rQPDG=pQPDG;
              fQPDG=be6QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mBe6);
		    }
            else if(sPDG==90004004)
		    {
              fQPDG=aQPDG;
              rQPDG=aQPDG;
              f4Mom=G4LorentzVector(0.,0.,0.,mAlph);
              g4Mom=f4Mom;
		    }
            else if(sPDG==90006002)
		    {
              rQPDG=pQPDG;
              fQPDG=b7QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mB7);
		    }
            else if(sPDG==90002007)
		    {
              rQPDG=nQPDG;
              fQPDG=a8QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              f4Mom=G4LorentzVector(0.,0.,0.,mHe8);
		    }
            else if(sPDG==90005004)
		    {
              rQPDG=pQPDG;
              fQPDG=aQPDG;
              tQPDG=aQPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mAlph);
              t4Mom=G4LorentzVector(0.,0.,0.,mAlph);
              three=true;
		    }
            else if(sPDG==90008004)
		    {
              tQPDG=pQPDG;
              rQPDG=pQPDG;
              fQPDG=c10QPDG;
              t4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mC10);
              three=true;
		    }
            else if(sPDG==90009006)
		    {
              rQPDG=pQPDG;
              fQPDG=o14QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mO14);
		    }
            else if(sPDG==90009007)
		    {
              rQPDG=pQPDG;
              fQPDG=o15QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mO15);
		    }
            else if(sPDG==90010006)
		    {
              tQPDG=pQPDG;
              rQPDG=pQPDG;
              fQPDG=o14QPDG;
              t4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mO14);
              three=true;
		    }
#ifdef ppdebug
            G4cout<<"G4QE::Fr:**three="<<three<<"**, r="<<rQPDG<<",f="<<fQPDG<<",t="<<tQPDG<<G4endl;
#endif
            if(!three)
            {
              if(!G4QHadron(s4M).DecayIn2(f4Mom,g4Mom))
              {
                G4cerr<<"***G4QE::Fr:(2)***FUSION*** tM["<<sPDG<<"]="<<tM<<" > sM="<<sM<<G4endl;
                G4Exception("***G4QEnvironment::Fragment:DecayIn2 didn't succeed for the fusion");
              }
#ifdef ppdebug
              G4cout<<"G4QE::Fr:*FUSION DONE*,fPDG=="<<sPDG<<",PDG1="<<hPDG<<",PDG2="<<bPDG<<G4endl;
#endif
              curHadr->SetQPDG(fQPDG);
              curHadr->Set4Momentum(f4Mom);
              backH->SetQPDG(rQPDG);
              backH->Set4Momentum(g4Mom);
		    }
            else
            {
              if(!G4QHadron(s4M).DecayIn3(f4Mom,g4Mom,t4Mom))
              {
                G4cerr<<"***G4QE::Fr:(3)***FUSION*** tM["<<sPDG<<"]="<<tM<<" > sM="<<sM<<G4endl;
                G4Exception("***G4QEnvironment::Fragment:DecayIn3 didn't succeed for the fusion");
              }
#ifdef ppdebug
              G4cout<<"G4QE::Fr:*DONE*,nH="<<nHadr<<",PDG="<<sPDG<<",1="<<hPDG<<",2="<<bPDG<<G4endl;
#endif
              curHadr->SetQPDG(fQPDG);
              curHadr->Set4Momentum(f4Mom);
              backH->SetQPDG(rQPDG);
              backH->Set4Momentum(g4Mom);
              G4QHadron* newH = new G4QHadron(tQPDG.GetPDGCode(),t4Mom);
              theQHadrons.push_back(newH); // (delete equivalent for newH)
              nHadr=theQHadrons.size();
#ifdef ppdebug
              G4cout<<"G4QE::Fr:*Products*,nH="<<nHadr<<",f="<<fQPDG<<f4Mom<<",b="<<rQPDG<<g4Mom
                    <<",new="<<tQPDG<<t4Mom<<",i="<<nHadr<<",f="<<theQHadrons.size()<<G4endl;
#endif
		    }
            tot4Mom+=b4m;               // Instead of the fused hadron
            tot4Mom-=g4Mom;             // subtract from the total the new hadron
            ///////////////////////////////break; // Make fusion only for one (?)
            // Instead the curHadr parameters should be updated ______
            hQC=fQPDG.GetQuarkContent();
		    h4m=f4Mom;
            // End of Instead ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		  } // End of the fusion check
	    } // End of the LOOP over previous hadrons
	  } // End of the FUSION check
      G4LorentzVector cH4Mom = curHadr->Get4Momentum(); // 4-mom of the current hadron
      tot4Mom-=cH4Mom;                                  // Reduce the total 4mom by the current 4mom
      if(hadron+1==nHadr)                               // The Last Hadron in the set
	  {
        G4double re=tot4Mom.e();
        G4double rpx=tot4Mom.px();
        G4double rpy=tot4Mom.py();
        G4double rpz=tot4Mom.pz();
        G4double dem=re*re+rpx*rpx+rpx*rpx+rpx*rpx;
#ifdef pdebug
        G4cout<<"G4QE::Fr: Is Energy&Mom conserved? t4M="<<tot4Mom<<",sd2="<<dem<<">10e-6?"<<G4endl;
#endif
        if(dem>0.000001)                                // Energy or momentum is not conserved
	    {
          if(dem>1.)G4cerr<<"***G4QE::Fr:E/Mom conservation="<<tot4Mom<<dem<<".Correction!"<<G4endl;
          G4QHadron* prevHadr = theQHadrons[nHadr-2]; // GetPointer to the Hadr previous to the Last
          G4LorentzVector pH4Mom = prevHadr->Get4Momentum(); // 4-mom of the previous hadron
          G4double cHM = curHadr->GetMass();          // Mass of the current hadron
          G4double pHM = prevHadr->GetMass();         // Mass of the current hadron
          tot4Mom+=cH4Mom+pH4Mom;
          G4double totRM=tot4Mom.m();
          if(cHM+pHM<=totRM)                            // *** Make the final correction ***
		  {
            if(!G4QHadron(tot4Mom).DecayIn2(pH4Mom,cH4Mom))
            {
              G4cerr<<"***G4QE::Fr:**Correction**tot4M="<<tot4Mom<<totRM<<" > sM="<<cHM+cHM<<G4endl;
              G4Exception("***G4QEnvironment::Fragment:DecayIn2 didn't succeed for the CORRECTION");
            }
#ifdef pdebug
            G4cout<<"G4QE::Fr:***CORRECTION IS DONE*** "<<G4endl;
#endif
            curHadr->Set4Momentum(cH4Mom);
            prevHadr->Set4Momentum(pH4Mom);
          }
          else
          {
            G4cerr<<"*!*G4QE::Fr: Can't correct "<<cHM<<"+"<<pHM<<"="<<cHM+pHM<<">"<<totRM<<G4endl;
            //G4Exception("***G4QEnvironment::Fragment: TEMPORARY EXCEPTION"); //@@@@@
          }
        }
      }
    }
    // >| 2     | 2  | 2     | 2     | 2      | 2 - old            | 1. If gamma: add to sum4Mom
    //  |>0->sum| 3  | 3     | 3     | 3      | 3 - old            | 2. Compare BN with the Last
    //  | 5     |>5  | 4     | 4     | 4      | 4 - old            | 3. Swap if larger, del the Last
    //  | 0     | 0  |>0->sum| 5<-sum| 5->evap| 2 - new            | 4. If the Last: add the sum
    //  | 4     | 4  | 5     | ex    |        |(0 - possible gamma)| 5. Decay/Eporate the Last
    //  | 3     | ex |                        | 3 - new
    G4LorentzVector sum(0.,0.,0.,0.);
    G4int gamCount=0;
    nHadr=theQHadrons.size();
    //if(nHadr)for(G4int h=nHadr-1; h>=0; h--)//Compress gamma & kill DecayedHadrons
    //if(nHadr>3)for(G4int h=nHadr-1; h>=0; h--)//Compress gamma & kill DecayedHadrons
	for(G4int h=nHadr-1; h>=0; h--)//Compress gamma & kill DecayedHadrons
    {
      G4QHadron* curHadr = theQHadrons[h];         // Get a pointer to the current Hadron
      G4int   hF = curHadr->GetNFragments();
      G4int hPDG = curHadr->GetPDGCode();
      if(hPDG==22) gamCount++;
      //if(hF||hPDG==22&&gamCount>1)                 // It should be compressed
	  if(hF||hPDG==22)                 // It should be compressed
	  {
        G4QHadron* theLast = theQHadrons[theQHadrons.size()-1]; // Get pointer to the Last Hadron
        if(hPDG==22) sum+=curHadr->Get4Momentum(); // Add 4Mom of gamma to the "sum"
        if(h<theQHadrons.size()-1)              // Need swap with the Last
	    {
          curHadr->SetNFragments(0);
          curHadr->Set4Momentum(theLast->Get4Momentum());
          //curHadr->SetQC(theLast->GetQC());
          curHadr->SetQPDG(theLast->GetQPDG());
	    }
        theQHadrons.pop_back();                  // the last QHadron is excluded from theQHadrons
        delete theLast; //*!!When killing, DON'T forget to delete the last QHadron as an instance!!*
        nHadr--;
      }
    }
#ifdef pdebug
    G4cout<<"G4QE::Fr: nH="<<nHadr<<", sum="<<sum<<G4endl;
#endif
    if(nHadr>1)for(G4int hdr=0; hdr<theQHadrons.size()-1; hdr++) // Order "theBiggest is theLast"
    {
#ifdef pdebug
      G4cout<<"G4QE::Fr:ORD, h="<<hdr<<"<"<<nHadr<<",hPDG="<<theQHadrons[hdr]->GetPDGCode()<<G4endl;
#endif
      G4QHadron* curHadr = theQHadrons[hdr];       // Get a pointer to the current Hadron
      G4QHadron* theLast = theQHadrons[theQHadrons.size()-1]; //Get a pointer to the Last Hadron
      G4int hB           = curHadr->GetBaryonNumber();
      G4int lB           = theLast->GetBaryonNumber();
#ifdef pdebug
      G4cout<<"G4QE::Fr: hBN="<<hB<<" < lastBN="<<lB<<", lastPDG="<<theLast->GetPDGCode()<<G4endl;
#endif
      if(lB<hB)                                    // Must be swapped
	  {
        G4QContent     hQC = curHadr->GetQC();
        G4QPDGCode   hQPDG = curHadr->GetQPDG();
        G4LorentzVector h4m= curHadr->Get4Momentum();
        curHadr->Set4Momentum(theLast->Get4Momentum());
        //curHadr->SetQC(theLast->GetQC());
        curHadr->SetQPDG(theLast->GetQPDG());
        theLast->Set4Momentum(h4m);
        //theLast->SetQC(hQC);
        theLast->SetQPDG(hQPDG);
	  }
    }
    nHadr=theQHadrons.size();
    G4QHadron* theLast = theQHadrons[nHadr-1];   // Get a pointer to the Last Hadron
    G4QHadron* theNew  = new G4QHadron(theLast); // Make a New Hadron of the Last Hadron
#ifdef pdebug
    G4cout<<"G4QE::Fr: Before LastSubtract nH="<<nHadr<<", lastPDG="<<theNew->GetPDGCode()<<G4endl;
#endif
    theQHadrons.pop_back();                  // the last QHadron is excluded from OUTPUT
    delete theLast; //*!! When killing, DON'T forget to delete the last QHadron as an instance !!*
    theNew->Set4Momentum(theNew->Get4Momentum()+sum); // Icrease 4Mom of the Last by the "sum"
#ifdef pdebug
    G4cout<<"G4QE::Fr: Before Evaporation nH="<<nHadr<<", lastPDG="<<theNew->GetPDGCode()<<G4endl;
#endif
    EvaporateResidual(theNew);   // Try to evaporate the Nucleus (delete equivalent)
    //theQHadrons.push_back(theNew);    // (delete equivalent for theNew)
    G4int onH=nHadr;
    nHadr=theQHadrons.size();
    //if(nHadr>onH) bfAct=true;
  } // End of the While-LOOOP for the Back Fusion
  // Now just fill the output theFravment vector (User is responsible to ClearAndDestroy it)
  G4QHadronVector* theFragments = new G4QHadronVector; // Inter creation. User's responsibale to del
  nHadr=theQHadrons.size();
  if(nHadr) for(G4int hd=0; hd<theQHadrons.size(); hd++)
  {
#ifdef pdebug
    G4cout<<"G4QE::Fr: LOOP starts h="<<hd<<", PDG="<<theQHadrons[hd]->GetPDGCode()<<G4endl;
#endif
    G4QHadron* curHadr = new G4QHadron(theQHadrons[hd]);
    G4int hPDG=curHadr->GetPDGCode();
#ifdef pdebug
    G4cout<<"G4QE::Fr:Copy's made with PDG="<<hPDG<<G4endl;// Just to compare with the original one
#endif
    if(hPDG==91000000) curHadr->SetQPDG(G4QPDGCode(3122)); // Move it to the next loop
    theFragments->push_back(curHadr);                         // (delete equivalent - user's respons.)
  }
#ifdef pdebug
  G4cout<<"G4QE::Fr:=== OUT ==>nH="<<theQHadrons.size()<<",nF="<<theFragments->size()<<G4endl;
#endif
  return theFragments;
} // End of "Fragment"

//The public Quasmons duplication with delete responsibility of User (!)
G4QuasmonVector* G4QEnvironment::GetQuasmons()
{//              =============================
  G4int nQ=theQuasmons.size();
#ifdef pdebug
  G4cout<<"G4QEnvironment::GetQuasmons is called nQ="<<nQ<<G4endl;
#endif
  G4QuasmonVector* quasmons = new G4QuasmonVector; // Intermediate
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
#ifdef pdebug
    G4cout<<"G4QEnv::GetQuasm:Q#"<<iq<<",QQPDG="<<theQuasmons[iq]->GetQPDG()<<",QQC="
          <<theQuasmons[iq]->GetQC()<<",M="<<theQuasmons[iq]->Get4Momentum().m()<<G4endl;
#endif
    G4Quasmon* curQ = new G4Quasmon(theQuasmons[iq]);
    quasmons->push_back(curQ);                       // (delete equivalent - user)
  }
#ifdef pdebug
  G4cout<<"G4QEnvironment::GetQuasmons ===OUT==="<<G4endl;
#endif
  return quasmons;
} // End of GetQuasmons

//The public Quasmons duplication with delete responsibility of User (!)
G4QHadronVector* G4QEnvironment::GetQHadrons()
{//              =============================
  G4int nH=theQHadrons.size();
#ifdef pdebug
  G4cout<<"G4QEnvironment::GetQHadrons is called nH="<<nH<<G4endl;
#endif
  G4QHadronVector* hadrons = new G4QHadronVector; // Intermediate
  if(nH) for(G4int ih=0; ih<nH; ih++)
  {
#ifdef pdebug
    G4cout<<"G4QEnv::GetQHadrons:H#"<<ih<<",HQPDG="<<theQHadrons[ih]->GetQPDG()<<",HQC="
          <<theQHadrons[ih]->GetQC()<<",HM="<<theQHadrons[ih]->GetMass()<<G4endl;
#endif
    G4QHadron* curH = new G4QHadron(theQHadrons[ih]);
    hadrons->push_back(curH);                       // (delete equivalent - user's responsibility)
  }
#ifdef pdebug
  G4cout<<"G4QEnvironment::GetQHadrons ===OUT==="<<G4endl;
#endif
  return hadrons;
} // End of GetQHadrons

//Decay of the excited dibayon in two baryons
void G4QEnvironment::DecayDibaryon(G4QHadron* qH)
{//  ============================================
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  G4LorentzVector q4M = qH->Get4Momentum();  // Get 4-momentum of the Dibaryon
  G4int          qPDG = qH->GetPDGCode();    // PDG Code of the decayin dybaryon
#ifdef pdebug
  G4cout<<"G4QEnv::DecayDibaryon: *Called* PDG="<<qPDG<<",4M="<<q4M<<G4endl;
#endif
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
  else if(qPDG==90001001)                    // "exited deutron" case
  {
    G4double qM=q4M.m();
    if(abs(qM-mDeut)<0.003)
	{
      theQHadrons.push_back(qH);                // Fill as it is (delete equivalent)
      return;
	}
    else if(mProt+mNeut<qM)
    {
      fPDG = 2112;
      sPDG = 2212;
      fMass= mNeut;    
      sMass= mProt;
    }
    else
    {
      fPDG = 22;
      sPDG = 90001001;
      fMass= 0.;
      sMass= mDeut;    
    }
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
          <<"(sM="<<sMass<<")"<<" >? TotM="<<q4M.m()<<q4M<<G4endl;
    G4Exception("***G4QEnv::DecayDibaryon: DecayIn2 didn't succeed for dibaryon");
  }
#ifdef pdebug
  G4cout<<"G4QEnv::DecayDibaryon: *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
        <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
  //qH->SetNFragments(2);                      // Fill a#of fragments to decaying Dibaryon
  //theQHadrons.push_back(qH);                    // Fill hadron with nf=2 (delete equivalent)
  // Instead
  delete qH;
  //
  G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
  theQHadrons.push_back(H1);                    // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
  theQHadrons.push_back(H2);                    // Fill "H2" (delete equivalent)
} // End of DecayDibaryon

//Decay of the excited 3p or 3n systems in three baryons
void G4QEnvironment::DecayThreeBaryon(G4QHadron* qH)
{//  ===============================================
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  G4LorentzVector q4M = qH->Get4Momentum();  // Get 4-momentum of the Dibaryon
  G4int          qPDG = qH->GetPDGCode();    // PDG Code of the decayin dybaryon
#ifdef pdebug
  G4cout<<"G4QEnv::DecayDibaryon: *Called* PDG="<<qPDG<<",4M="<<q4M<<G4endl;
#endif
  G4int          fPDG = 3122;                // Prototype for 3 lambdas case
  G4double       fMass= mLamb;
  if     (qPDG==90000003)                    // "three-neutron" case
  {
    fPDG = 2112;
    fMass= mNeut;
  }
  else if(qPDG==90003000)                    // "exited deutron" case
  {
    fPDG = 2212;
    fMass= mProt;
  }
  else if(qPDG!=93000000)                    // "Bad call" case
  {
    G4cerr<<"***G4QEnv::DecayThreeBaryon: PDG="<<qPDG<<G4endl;
    G4Exception("***G4QEnv::DecayThreeBaryon: Can not decay this PDG Code");
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,fMass);
  G4LorentzVector t4Mom(0.,0.,0.,fMass);
  if(!G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
  {
    G4cerr<<"***G4QEnv::DecayThreeBaryon: fPDG="<<fPDG<<"(fM="<<fMass<<")*3 = "<<3*fMass
          <<" >? TotM="<<q4M.m()<<q4M<<G4endl;
    G4Exception("***G4QEnv::DecayThreeBaryon: DecayIn3 didn't succeed for the ThreeBaryon");
  }
#ifdef pdebug
  G4cout<<"G4QEnv::Decay3Bar: *DONE* fPDG="<<fPDG<<",f="<<f4Mom<<",s="<<s4Mom<<",t="<<t4Mom<<G4endl;
#endif
  //qH->SetNFragments(3);                      // Fill a#of fragments to decaying Dibaryon
  //theQHadrons.push_back(qH);                    // Fill hadron with nf=2 (delete equivalent)
  // Instead
  delete qH;
  //
  G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
  theQHadrons.push_back(H1);                    // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(fPDG,s4Mom); // Create a Hadron for the 2-nd baryon
  theQHadrons.push_back(H2);                    // Fill "H2" (delete equivalent)
  G4QHadron* H3 = new G4QHadron(fPDG,t4Mom); // Create a Hadron for the 3-d baryon
  theQHadrons.push_back(H3);                    // Fill "H3" (delete equivalent)
} // End of DecayThreeBaryon

//Decay of the excited alpha+2p or alpha+2n systems
void G4QEnvironment::DecayAlphaDiN(G4QHadron* qH)
{//  ============================================
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mHel6= G4QPDGCode(2112).GetNuclMass(2,4,0);
  G4LorentzVector q4M = qH->Get4Momentum();  // Get 4-momentum of the Dibaryon
  G4int          qPDG = qH->GetPDGCode();    // PDG Code of the decayin dybaryon
#ifdef pdebug
  G4cout<<"G4QEnv::DecayDibaryon: *Called* PDG="<<qPDG<<",4M="<<q4M<<G4endl;
#endif
  G4int          fPDG = 2212;                // Prototype for alpha+pp case
  G4double       fMass= mProt;
  G4int          sPDG = 90002002;
  G4double       sMass= mAlph;
  if     (qPDG==90002004)                    // "alpha+2neutrons" case
  {
    G4double qM=q4M.m();
    if(abs(qM-mHel6)<0.003)
	{
      theQHadrons.push_back(qH);                // Fill as it is (delete equivalent)
      return;
	}
    else if(mNeut+mNeut+mAlph<qM)
    {
      fPDG = 2112;
      fMass= mNeut;
	}
    else
    {
      G4cerr<<"***G4QEnv::DecayAlphaDiN: M(He6="<<mHel6<<")="<<qM<<" < "<<mNeut+mNeut+mAlph<<G4endl;
      G4Exception("***G4QEnv::DecayAlphaDiN: Can not decay excited He6 with such a mass");
    }
  }
  else if(qPDG!=90004002)                    // "Bad call" case
  {
    G4cerr<<"***G4QEnv::DecayAlphaDiN: PDG="<<qPDG<<G4endl;
    G4Exception("***G4QEnv::DecayAlphaDiN: Can not decay this PDG Code");
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,fMass);
  G4LorentzVector t4Mom(0.,0.,0.,sMass);
  if(!G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
  {
    G4cerr<<"***G4QEnv::DecayAlphaDiN: fPDG="<<fPDG<<"(fM="<<fMass<<")*2+mAlpha = "
          <<fMass+fMass+sMass<<" >? TotM="<<q4M.m()<<q4M<<G4endl;
    G4Exception("***G4QEnv::DecayAlphaDiN: DecayIn3 didn't succeed for the Alpha+N+N");
  }
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAl2N: *DONE* fPDG="<<fPDG<<",f="<<f4Mom<<",s="<<s4Mom<<",t="<<t4Mom<<G4endl;
#endif
  //qH->SetNFragments(3);                      // Fill a#of fragments to decaying Dibaryon
  //theQHadrons.push_back(qH);                    // Fill hadron with nf=2 (delete equivalent)
  // Instead
  delete qH;
  //
  G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
  theQHadrons.push_back(H1);                    // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(fPDG,s4Mom); // Create a Hadron for the 2-nd baryon
  theQHadrons.push_back(H2);                    // Fill "H2" (delete equivalent)
  G4QHadron* H3 = new G4QHadron(sPDG,t4Mom); // Create a Hadron for the alpha
  theQHadrons.push_back(H3);                    // Fill "H3" (delete equivalent)
} // End of DecayThreeBaryon

//Decay of the excited alpha+bayon state in alpha and baryons
void G4QEnvironment::DecayAlphaBar(G4QHadron* qH)
{//  ============================================
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  G4LorentzVector q4M = qH->Get4Momentum();  // Get 4-momentum of the Alpha-Baryon
  G4double         qM = q4M.m();             // Mass of Alpha-Baryon
  G4int          qPDG = qH->GetPDGCode();    // PDG Code of the decayin Alpha-Baryon
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAlphaBar: *Called* PDG="<<qPDG<<",4M="<<q4M<<G4endl;
#endif
  G4int          fPDG = 90002002;            // Prototype for "alpha+n" case
  G4int          sPDG = 2112;
  G4double       fMass= mAlph;
  G4double       sMass= mNeut;
  if     (qPDG==90003002)                    // "alpha+p" case
  {
    sPDG = 2212;
    sMass= mProt;    
  }
  else if(qPDG==9100202)                     // "alpha+l" case
  {
    sPDG = 3122;
    sMass= mLamb;    
  }
  else if(qPDG!=90002003)
  {
    theQHadrons.push_back(qH);                  // Fill hadron as it is (delete equivalent)
    return;
  }
  G4double dM=fMass+sMass-qM;
  if(dM>0.&&dM<1.)
  {
    G4cerr<<"***Corrected***G4QEnv::DecayAlphaBar:fPDG="<<fPDG<<"(fM="<<fMass<<") + sPDG="<<sPDG
          <<"(sM="<<sMass<<")="<<fMass+sMass<<" > TotM="<<qM<<q4M<<G4endl;
    G4double hdM=dM/2;
    fMass-=hdM;
    sMass-=hdM;
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,sMass);     // Mass is random since probab. time
  if(!G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
  {
    G4cerr<<"***G4QEnv::DecayAlphaBar:fPDG="<<fPDG<<"(fM="<<fMass<<") + sPDG="<<sPDG
          <<"(sM="<<sMass<<")="<<fMass+sMass<<" > TotM="<<q4M.m()<<q4M<<G4endl;
    G4Exception("***G4QEnv::DecayAlphaBar: DecayIn2 didn't succeed for alpha+baryon");
  }
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAlphaBar: *DONE* alpha4M="<<f4Mom<<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
  //qH->SetNFragments(2);                      // Fill a#of fragments to decaying Dibaryon
  //theQHadrons.push_back(qH);                    // Fill hadron with nf=2 (delete equivalent)
  // Instead
  delete qH;
  //
  G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the alpha
  theQHadrons.push_back(H1);                    // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the baryon
  theQHadrons.push_back(H2);                    // Fill "H2" (delete equivalent)
} // End of DecayAlphaBar

//Decay of the excited alpha+alpha state in 2 alphas
void G4QEnvironment::DecayAlphaAlpha(G4QHadron* qH)
{//  ==============================================
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  G4int          qPDG = qH->GetPDGCode();    // PDG Code of the decayin dialpha
  if(qPDG!=90004004)
  {
    G4cerr<<"***G4QEnv::DecayAlphaAlpha: qPDG="<<qPDG<<G4endl;
    G4Exception("***G4QEnv::DecayAlphaAlpha: Not Be8 state decais in 2 alphas");
  }
  G4LorentzVector q4M = qH->Get4Momentum();  // Get 4-momentum of the Dibaryon
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAlphaAlpha: *Called* PDG="<<qPDG<<",4M="<<q4M<<G4endl;
#endif
  G4int          fPDG = 90002002;
  G4int          sPDG = 90002002;
  G4double       fMass= mAlph;
  G4double       sMass= mAlph;
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,sMass);     // Mass is random since probab. time
  if(!G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
  {
    G4cerr<<"***G4QEnv::DecayAlphaAlpha:fPDG="<<fPDG<<"(fM="<<fMass<<") + sPDG="<<sPDG
          <<"(sM="<<sMass<<")"<<" >? TotM="<<q4M.m()<<q4M<<G4endl;
    G4Exception("***G4QEnv::DecayAlphaAlpha: DecayIn2 didn't succeed for dialpha");
  }
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAlphaAlpha: *DONE* fal4M="<<f4Mom<<", sal4M="<<s4Mom<<G4endl;
#endif
  //qH->SetNFragments(2);                      // Fill a#of fragments to decaying Dibaryon
  //theQHadrons.push_back(qH);                    // Fill hadron with nf=2 (delete equivalent)
  // Instead
  delete qH;
  //
  G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st alpha
  theQHadrons.push_back(H1);                    // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd alpha
  theQHadrons.push_back(H2);                    // Fill "H2" (delete equivalent)
} // End of DecayAlphaAlpha

// Check that it's still possible to decay the Total Residual Nucleus in Quasmon + Environ & correct
G4bool G4QEnvironment::CheckGroundState(G4Quasmon* quasm, G4bool corFlag) // Cor's forbidden by def.
{ //   ==================================================================
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  ///@@@///
  ///////////////corFlag=true;
  ///@@@///
  G4QContent valQ=quasm->GetQC();                // Quark content of the Quasmon
  G4int    resQPDG=valQ.GetSPDGCode();           // Reachable in a member function
  G4double resQMa=G4QPDGCode(resQPDG).GetMass(); // GS Mass of the Residual Quasmon
  G4double resEMa=0.;                            // GS Mass of the Empty Residual Environment
  G4bool   bsCond=false;                         // BaryonSeparetionCondition for Quasmon in vacuum
  G4LorentzVector enva4M=G4LorentzVector(0.,0.,0.,0.);
  G4LorentzVector reTLV=quasm->Get4Momentum();   // Prototyoe of the 4-Mom of the Residual Nucleus
  G4double resSMa=resQMa;                        // Prototype of MinimalSplitMass of ResidualNucleus
  G4int envPDG=theEnvironment.GetPDG();
  if(envPDG!=90000000)                           // => "Environment is not vacuum" case
  { // @@@@@@@@@@@@@@@@@@@ CALL SUBROUTINE @@@@@@@@@
    resEMa=theEnvironment.GetMZNS();             // GSMass of the Residual Environment
    enva4M=G4LorentzVector(0.,0.,0.,resEMa);     // 4-Mom of the Residual Environment
    reTLV+=enva4M;                               // 4-Mom of Residual Nucleus
    resSMa+=resEMa;                              // Minimal Split Mass of Residual Nucleus
  }
  else                                           // Calculate BaryonSeparetionCondition for vacQuasm
  {
    G4double resQM=reTLV.m();                    // CM Mass of the Residual vacQuasmon
    G4int   baryn=valQ.GetBaryonNumber();        // Baryon Number of the Residual vacQuasmon
    if(baryn>1)                                  // => "Can split baryon" case
	{
      if(valQ.GetN())                            // ===> "Can split neutron" case
	  {
        G4QContent resQC=valQ-neutQC;                 // QC of Residual for the Neutron
        G4int    resPDG=resQC.GetSPDGCode();          // PDG of Residual for the Neutron
        G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
        if(resMas+mNeut<resQM) bsCond=true;
	  }
      else if(valQ.GetP())                       // ===> "Can split proton" case
	  {
        G4QContent resQC=valQ-protQC;                 // QC of Residual for the Proton
        G4int    resPDG=resQC.GetSPDGCode();          // PDG of Residual for the Proton
        G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
        if(resMas+mProt<resQM) bsCond=true;
	  }
      else if(valQ.GetL())                       // ===> "Can split lambda" case
	  {
        G4QContent resQC=valQ-lambQC;                 // QC of Residual for the Lambda
        G4int    resPDG=resQC.GetSPDGCode();          // PDG of Residual for the Lambda
        G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
        if(resMas+mLamb<resQM) bsCond=true;
	  }
	}
  }
  G4double resTMa=reTLV.m();                     // CM Mass of the Residual Nucleus (Quasm+Environ)
  G4int nOfOUT = theQHadrons.size();          // Total #of QHadrons at this point
#ifdef pdebug
  G4cout<<"G4QEnv::CheckGS:tM="<<resTMa<<"<rQM+rEM="<<resSMa<<"&n="<<nOfOUT<<">0. Correct?"<<G4endl;
#endif
  if((resTMa<resSMa || !resEMa&&!bsCond) && nOfOUT>0 && corFlag) // *** CORRECTION ***
  {
    G4QHadron*  theLast = theQHadrons[nOfOUT-1];
    if(theLast->GetNFragments() || theLast->GetPDGCode()==22) return false;// Decayed H or gamma
    G4LorentzVector hadr4M = theLast->Get4Momentum();
    G4double  hadrMa=hadr4M.m();
    G4LorentzVector tmpTLV=reTLV+hadr4M;          // Tot (ResidNucl+LastHadron) 4-Mom
#ifdef pdebug
	G4cout<<"G4QEnv::CheckGS: YES,4M/M="<<tmpTLV<<tmpTLV.m()<<" > rSM+hM="<<resSMa+hadrMa<<G4endl;
#endif
    if(tmpTLV.m()>resSMa+hadrMa)                  // resMa contains 2 Hadrons: resQ and Environ
    {
      if(resEMa)                                  // => "The not vacuum Environment exist" case
      {
        G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resQMa); // GS Mass of Quasmon
        G4QHadron* quasH = new G4QHadron(valQ, quas4M);
        G4QHadron* envaH = new G4QHadron(theEnvironment.GetQCZNS(),enva4M);
        if(!G4QHadron(tmpTLV).DecayIn3(hadr4M,quas4M,enva4M))
        {
          delete quasH;                            // Delete "Quasmon Hadron"
          delete envaH;                            // Delete "Environ Hadron"
          G4cerr<<"***G4QEnv::CheckGS: Decay in Fragm+ResQ+ResEnv did not succeeded"<<G4endl;
          return false;
        }
        else
        {
          //@@CHECK CoulBar (only for ResQuasmon in respect to ResEnv) and may be evaporate instead
          theLast->Set4Momentum(hadr4M);
          quasH->Set4Momentum(quas4M);
          if(resQPDG==92000000||resQPDG==90002000||resQPDG==90000002)DecayDibaryon(quasH); // DelEqu
          else if(resQPDG==90002003||resQPDG==90003002) DecayAlphaBar(quasH); //DelEqu
          else if(resQPDG==90004004) DecayAlphaAlpha(quasH); //DelEqu
          else theQHadrons.push_back(quasH);          // Fill ResidQuasm Hadron (delete equivalent)
          envaH->Set4Momentum(enva4M);
          if(envPDG==92000000||envPDG==90002000||envPDG==90000002) DecayDibaryon(envaH); // DelEqu
          else if(envPDG==90002003||envPDG==90003002) DecayAlphaBar(envaH); //DelEqu
          else if(envPDG==90004004) DecayAlphaAlpha(envaH); //DelEqu
          else theQHadrons.push_back(envaH);          // Fill 2nd Hadron (delete equivalent)
		}
	  }
      else                                         // => "The Environment is vacuum" case (DecayIn2)
      {
        G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resQMa); // GS Mass of Quasmon
        G4QHadron* quasH = new G4QHadron(valQ, quas4M);
        if(!G4QHadron(tmpTLV).DecayIn2(hadr4M,quas4M))
        {
          delete quasH;                            // Delete "Quasmon Hadron"
          G4cerr<<"***G4QEnv::CheckGS: Decay in Fragm+ResQ did not succeeded"<<G4endl;
          return false;
        }
        else
        {
          //@@CHECK CoulBar (only for ResQuasmon in respect to ResEnv) and may be evaporate instead
          theLast->Set4Momentum(hadr4M);
          quasH->Set4Momentum(quas4M);
          if(resQPDG==92000000||resQPDG==90002000||resQPDG==90000002)DecayDibaryon(quasH); // DelEq
          else if(resQPDG==90002003||resQPDG==90003002) DecayAlphaBar(quasH); //DelEqu
          else if(resQPDG==90004004) DecayAlphaAlpha(quasH); //DelEqu
          else theQHadrons.push_back(quasH);               // Fill ResidQuasmHadron (delete equivalent)
		}
	  }
    }
    else                                                // "CORRECTION" !!!
    {
      if(nOfOUT>1 && corFlag)
	  {
        G4QHadron*  thePrev = theQHadrons[nOfOUT-2];    // Get pointer to previos before last hadron
        if(thePrev->GetNFragments()||thePrev->GetNFragments()) return false; // Decayed H or gamma
        G4LorentzVector prev4M = thePrev->Get4Momentum();
        G4double  prevMa=prev4M.m();                    // Mass of previous hadron
        tmpTLV+=prev4M;                                 // Increment Total 4-Mom of TotalResidNucl
        G4QContent totQC=valQ+theEnvironment.GetQCZNS();// Quark Content of the ResidNucl=ResidQ+Env
        G4int      totPDG=totQC.GetSPDGCode();          // PDG Code of Total Residual Nucleus 
        G4double   tQMa=G4QPDGCode(totPDG).GetMass();   // GS Mass of the Residual Nucleus
#ifdef pdebug
	    G4cout<<"G4QEnv::CheckGS:NO, T="<<tmpTLV<<tmpTLV.m()<<">r+p+l="<<tQMa+hadrMa+prevMa<<G4endl;
#endif
        if(tmpTLV.m()>tQMa+hadrMa+prevMa)
        {
          G4LorentzVector nuc4M = G4LorentzVector(0.,0.,0.,tQMa); // 4-Mom of ResidNucleus at rest
          G4QHadron* nucH = new G4QHadron(totQC, nuc4M);
          if(!G4QHadron(tmpTLV).DecayIn3(hadr4M,prev4M,nuc4M))
          {
            delete nucH;                                // Delete "Residual Nucleus Hadron"
            G4cerr<<"***G4QEnv::CheckGS: Decay in ResidNucl+LastH+PrevH did not succeeded"<<G4endl;
            return false;
          }
          else
          {
            theLast->Set4Momentum(hadr4M);
            thePrev->Set4Momentum(prev4M);
            nucH->Set4Momentum(nuc4M);
#ifdef pdebug
	        G4cout<<"G4QEnv::CheckGS:***SUCCEED***>CHECK, D4M="<<tmpTLV-hadr4M-prev4M-nuc4M<<G4endl;
#endif
            if(totPDG==92000000||totPDG==90002000||totPDG==90000002) DecayDibaryon(nucH); //DelEqu
            else if(totPDG==90002003||totPDG==90003002) DecayAlphaBar(nucH); //DelEqu
            else if(totPDG==90004004) DecayAlphaAlpha(nucH); //DelEqu
            else theQHadrons.push_back(nucH);               // Fill ResidNuclHadron (delete equivalent)
	      }
		}
        else return false;
	  }
      else return false;
    }
  }
  else return false;
  return true;
} // End of "CheckGroundState"

// Try to decay the Total Residual Nucleus in Environ+Quasmon
G4bool G4QEnvironment::DecayInEnvQ(G4Quasmon* quasm)
{ //   =============================================
  G4QContent valQ=quasm->GetQC();                // Quark content of the Quasmon
  G4int    resQPDG=valQ.GetSPDGCode();           // Reachable in a member function
  G4double resQMa=G4QPDGCode(resQPDG).GetMass(); // GS Mass of the Residual Quasmon
  G4LorentzVector enva4M=G4LorentzVector(0.,0.,0.,0.);
  G4LorentzVector reTLV=quasm->Get4Momentum();   // Prototyoe of the 4-Mom of the Residual Nucleus
  G4double resSMa=resQMa;                        // Prototype of MinimalSplitMass of ResidualNucleus
  G4int envPDG=theEnvironment.GetPDG();          // PDG Code of the Environment
  if(envPDG!=90000000)                           // => "Environment is not vacuum" case
  {
    G4double resEMa=theEnvironment.GetMZNS();    // GSMass of the Residual Environment
    enva4M=G4LorentzVector(0.,0.,0.,resEMa);     // 4-Mom of the Residual Environment
    reTLV+=enva4M;                               // 4-Mom of Residual Nucleus
    resSMa+=resEMa;                              // Minimal Split Mass of Residual Nucleus
    G4double resTMa=reTLV.m();                   // CM Mass of the Residual Nucleus (Quasm+Environ)
	//#ifdef pdebug
    G4cout<<"G4QEnv::DecayInEnvQ: totM="<<reTLV<<resTMa<<" > rQM+rEM="<<resSMa<<G4endl;
	//#endif
    if(resTMa>resSMa)
    {
      G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resQMa); // GS Mass of Quasmon
      G4QHadron* quasH = new G4QHadron(valQ, quas4M);
      G4QHadron* envaH = new G4QHadron(envPDG, enva4M);
      if(!G4QHadron(reTLV).DecayIn2(enva4M,quas4M))
      {
        delete quasH;                            // Delete "Quasmon Hadron"
        delete quasH;                            // Delete "Enviromnent Hadron"
        G4cerr<<"***G4Quasm::DecayInEnvQ: Decay in Envuronment+ResQuasm did not succeeded"<<G4endl;
        return false;
      }
      else
      {
        quasH->Set4Momentum(quas4M);
        if(resQPDG==92000000||resQPDG==90002000||resQPDG==90000002) DecayDibaryon(quasH); // DelEqu
        else if(resQPDG==90002003||resQPDG==90003002) DecayAlphaBar(quasH); //DelEqu
        else if(resQPDG==90004004) DecayAlphaAlpha(quasH); //DelEqu
        else theQHadrons.push_back(quasH);          // Fill ResidQuasm Hadron (delete equivalent)
        envaH->Set4Momentum(enva4M);
        if(envPDG==92000000||envPDG==90002000||envPDG==90000002) DecayDibaryon(envaH); // DelEquival
        else if(envPDG==90002003||envPDG==90003002) DecayAlphaBar(envaH); //DelEqu
        else if(envPDG==90004004) DecayAlphaAlpha(envaH); //DelEqu
        else theQHadrons.push_back(envaH);          // Fill Environment Hadron (delete equivalent)
	  }
    }
    else return false;
  }
  else return false;                             // => "Environment is vacuum" case
  return true;
} // End of "DecayInEnvQ"




