// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Quasmon.cc,v 1.19 2000-09-24 16:03:19 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4Quasmon ----------------
//             by Mikhail Kossov, July 1999.
//  class for an excited hadronic state used by the CHIPS Model
// ------------------------------------------------------------

//#define debug
//#define pdebug
//#define ppdebug
//#define sdebug

#include "G4Quasmon.hh"

G4Quasmon::G4Quasmon(G4QContent qQCont, G4LorentzVector q4M, G4LorentzVector ph4M):
  q4Mom(q4M), valQ(qQCont), phot4M(ph4M)
{
#ifdef pdebug
  G4cout<<"G4Quasmon: ***DIRECT CREATION*** QQC="<<qQCont<<", Q4M="<<q4M
        <<", (ph4M="<<ph4M<<")"<<G4endl;
#endif
  G4int nP= theWorld.GetQPEntries();      // A#of initialized particles in CHIPS World
  G4int          nMesons  = 27;
  if     (nP<24) nMesons  =  9;
  else if(nP<41) nMesons  = 18;
  G4int          nBaryons = 36;
  if     (nP<18) nBaryons = -1;
  else if(nP<35) nBaryons = 16;
  G4int          nClusters= nP-73;
  InitCandidateVector(nMesons,nBaryons,nClusters);
#ifdef pdebug
  G4cout<<"G4Quasmon: Hadronization candidates are initialized with nMesons="<<nMesons
        <<",nBaryons="<<nBaryons<<",nClusters="<<nClusters<<G4endl;
#endif
  status=3;                                   
}

G4Quasmon::G4Quasmon(const G4Quasmon& right)
{
  q4Mom                 = right.q4Mom;
  valQ                  = right.valQ;
  //theEnvironment        = right.theEnvironment;
  status                = right.status;
  //theQHadrons (Vector)
  G4int nQH             = right.theQHadrons.entries();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
    theQHadrons.insert(curQH);
  }
  theWorld              = right.theWorld;
  phot4M                = right.phot4M;
  nBarClust             = right.nBarClust;
  nOfQ                  = right.nOfQ;
  //theQCandidates (Vector)
  G4int nQC             = right.theQCandidates.entries();
  if(nQC) for(G4int iq=0; iq<nQC; iq++)
  {
    G4QCandidate* curQC = new G4QCandidate(right.theQCandidates[iq]);
    theQCandidates.insert(curQC);
  }
  f2all                 = right.f2all;
  rEP                   = right.rEP;
  rMo                   = right.rMo;
}

G4Quasmon::G4Quasmon(G4Quasmon* right)
{
#ifdef debug
  G4cout<<"G4Quasmon:Copy-Constructor is called E="<<right->theEnvironment<<G4endl;
#endif
  q4Mom                 = right->q4Mom;
  valQ                  = right->valQ;
  //theEnvironment        = right->theEnvironment;
  status                = right->status;
  //theQHadrons (Vector)
  G4int nQH             = right->theQHadrons.entries();
#ifdef debug
  G4cout<<"G4Quasmon:Copy-Constructor:nQH="<<nQH<<G4endl;
#endif
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
#ifdef debug
    G4cout<<"G4Quasmon:Copy-Constructor:H#"<<ih<<",QH="<<right->theQHadrons[ih]<<G4endl;
#endif
    G4QHadron* curQH    = new G4QHadron(right->theQHadrons[ih]);
    theQHadrons.insert(curQH);
  }
  theWorld              = right->theWorld;
  phot4M                = right->phot4M;
  nBarClust             = right->nBarClust;
  nOfQ                  = right->nOfQ;
  //theQCandidates (Vector)
  G4int nQC             = right->theQCandidates.entries();
#ifdef debug
  G4cout<<"G4Quasmon:Copy-Constructor:nCand="<<nQC<<G4endl;
#endif
  if(nQC) for(G4int iq=0; iq<nQC; iq++)
  {
#ifdef debug
    G4cout<<"G4Quasmon:Copy-Constructor:C#"<<iq<<",QC="<<right->theQCandidates[iq]<<G4endl;
#endif
    G4QCandidate* curQC = new G4QCandidate(right->theQCandidates[iq]);
    theQCandidates.insert(curQC);
  }
  f2all                 = right->f2all;
  rEP                   = right->rEP;
  rMo                   = right->rMo;
#ifdef debug
  G4cout<<"G4Quasmon:Copy-Constructor: >>>DONE<<<"<<G4endl;
#endif
}

G4Quasmon::~G4Quasmon()
{
  G4int nCan=theQCandidates.entries();
#ifdef sdebug
  G4cout<<"~G4Quasmon: before theQCandidates nC="<<nCan<<G4endl;
#endif
  if(nCan) for (G4int ic=0; ic<nCan; ic++)
  {
#ifdef sdebug
    G4cout<<"~G4Quasmon: delete the QCandidate #"<<ic<<G4endl;
#endif
    delete theQCandidates[ic];
#ifdef sdebug
    G4cout<<"~G4Quasmon: after delete"<<G4endl;
#endif
  }
  //theQCandidates.clearAndDestroy();
#ifdef sdebug
  G4cout<<"~G4Quasmon: before theQHadrons"<<G4endl;
#endif
  theQHadrons.clearAndDestroy();
#ifdef sdebug
  G4cout<<"~G4Quasmon: === DONE ==="<<G4endl;
#endif
}

G4double G4Quasmon::Temperature=180.;  
G4double G4Quasmon::SSin2Gluons=0.1;  
G4double G4Quasmon::EtaEtaprime=0.3;
// Fill the private static parameters
void G4Quasmon::SetParameters(G4double temperature, G4double ssin2g, G4double etaetap)
{//  =================================================================================
  Temperature=temperature; 
  SSin2Gluons=ssin2g; 
  EtaEtaprime=etaetap;
}
void G4Quasmon::SetTemper(G4double temperature) {Temperature=temperature;}
void G4Quasmon::SetSOverU(G4double ssin2g)      {SSin2Gluons=ssin2g;}
void G4Quasmon::SetEtaSup(G4double etaetap)     {EtaEtaprime=etaetap;}

const G4Quasmon& G4Quasmon::operator=(const G4Quasmon& right)
{ //=========================================================
  q4Mom                 = right.q4Mom;
  valQ                  = right.valQ;
  //theEnvironment        = right.theEnvironment;
  status                = right.status;
  //theQHadrons (Vector)
  G4int nQH             = right.theQHadrons.entries();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
    theQHadrons.insert(curQH);
  }
  theWorld              = right.theWorld;
  phot4M                = right.phot4M;
  nBarClust             = right.nBarClust;
  nOfQ                  = right.nOfQ;
  //theQCandidates (Vector)
  G4int nQC             = right.theQCandidates.entries();
  if(nQC) for(G4int iq=0; iq<nQC; iq++)
  {
    G4QCandidate* curQC = new G4QCandidate(right.theQCandidates[iq]);
    theQCandidates.insert(curQC);
  }
  f2all                 = right.f2all;
  rEP                   = right.rEP;
  rMo                   = right.rMo;
}

//Initialize a Candidate vector for the instance of a Quasmon
void G4Quasmon::InitCandidateVector(G4int maxMes, G4int maxBar, G4int maxClust)
//   ==========================================================================
{
  static const G4int nOfMesons =45;//a#of S=0,1,2,3,4 Mesons, possible candidates to out hadrons
  static const G4int nOfBaryons=72;//a#of 1/2,3/2,5/2,7/2 Baryons,possible candidates to hadrons
  // Scalar resonances   (0):           Eta,Pi0,Pi+,APi-,Ka0,Ka+,AKa0,AKa-,Eta*
  static G4int mesonPDG[nOfMesons]  =  {221,111,211,-211,311,321,-311,-321,331
  // Vector resonances   (1):           omega,Rh0,Rh+,Rho-,K0*,K+*,AK0*,AK-*,Phi
								       ,223,113,213,-213,313,323,-313,-323,333
  // Tensor D-resonances (2):            f2,a20,a2+, a2-,K20,K2+,AK20,AK2-,f2'
                                       ,225,115,215,-215,315,325,-315,-325,335
  // Tensor F-resonances (3):           om3,ro3,r3+,rh3-,K30,K3+,AK30,AK3-,Phi3
                                       ,227,117,217,-217,317,327,-317,-327,337
  // Tensor G-resonances (4):           f4,a40,a4+, a4-,K40,K4+,AK40,AK4-,f4'
								       ,229,119,219,-219,319,329,-319,-329,339};

  // Baryon octet      (1/2):            n  , an  , p  , ap  ,lamb,alamb,sig-,asig-
  static G4int baryonPDG[nOfBaryons] = {2112,-2112,2212,-2212,3122,-3122,3112,-3112
							          //sig0,asig0,sig+,asig+,ksi-,aksi-,ksi0,aksi0
                                       ,3212,-3212,3222,-3222,3312,-3312,3322,-3322
  // Baryon decuplet   (3/2):           del-,adel-,del0,adel0,del+,adel+,dl++,adl++,sis-,asis-
                                       ,1114,-1114,2114,-2114,2214,-2214,2224,-2224,3114,-3114
							          //sis0,asis0,sis+,asis+,kss-,akss-,kss0,akss0,omeg,aomeg
                                       ,3214,-3214,3224,-3224,3314,-3314,3324,-3324,3334,-3334
  // Baryon octet      (5/2):           n5/2,an5/2,p5/2,ap5/2,l5/2,al5/2,si5-,asi5-
                                       ,2116,-2116,2216,-2216,3126,-3126,3116,-3116
							          //si50,asi50,si5+,asi5+,ks5-,aks5-,ks50,aks50
                                       ,3216,-3216,3226,-3226,3316,-3316,3326,-3326
  // Baryon decuplet   (7/2):           dl5-,adl5-,dl50,adl50,dl5+,adl5+,d5++,ad5++,si5-,asi5-
                                       ,1118,-1118,2118,-2118,2218,-2218,2228,-2228,3118,-3118
							          //si50,asi50,si5+,asi5+,ks5-,aks5-,ks50,aks50,ome5,aome5
								       ,3218,-3218,3228,-3228,3318,-3318,3328,-3328,3338,-3338};
  G4int i=0;
  G4int ind=0;
  if(maxMes>nOfMesons) maxMes=nOfMesons;
  if(maxMes>=0) for (i=0; i<maxMes; i++) 
  {
    theQCandidates.insert(new G4QCandidate(mesonPDG[i]));
#ifdef sdebug
    cout<<"G4Quasmon::InitCandidateVector: "<<ind++<<", Meson # "<<i<<" with code = "
        <<mesonPDG[i]<<", QC="<<theQCandidates[i]->GetQC()<<" is initialized"<<endl;
#endif
  }
  if(maxBar>nOfBaryons) maxBar=nOfBaryons;
  if(maxBar>=0) for (i=0; i<maxBar; i++) 
  {
    theQCandidates.insert(new G4QCandidate(baryonPDG[i]));
#ifdef sdebug
    cout<<"G4Quasmon::InitCandidateVector: "<<ind++<<", Baryon # "<<i<<" with code = "
        <<baryonPDG[i]<< ", QC="<<theQCandidates[i]->GetQC()<<" is initialized"<<endl;
#endif
  }
  if(maxClust>=0) for (i=0; i<maxClust; i++) 
  {
    G4int clustQCode = i+72;          // Q-code of the cluster in the CHIPS world
    G4QPDGCode clustQPDG;
    clustQPDG.InitByQCode(clustQCode);
    G4int clusterPDG=clustQPDG.GetPDGCode();
    theQCandidates.insert(new G4QCandidate(clusterPDG));
#ifdef sdebug
    cout<<"G4Quasmon::InitCandidateVector:"<<ind++<<", Cluster # "<<i<<" with code = "
        <<clusterPDG<<", QC="<<clustQPDG.GetQuarkContent()<<" is initialized"<<endl;
#endif
  }
}

// Fragmentation of the Quasmon (the main member function)
G4QHadronVector G4Quasmon::HadronizeQuasmon(G4QNucleus& qEnv)
//              =============================================
{
  static const G4int NUCPDG  = 90000000;
  static const G4double np   = 1877.9; //@@ temporary = a mass of np pair > m_n + m_p
  static const G4int maxRet  = 7;
  static const G4int mEv     = 0;
  static const G4double BIG  = 100000.;
  static const G4double BIG2 = BIG*BIG;
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QContent PiQC(0,1,0,1,0,0);
  static const G4QContent K0QC(1,0,0,0,0,1);
  static const G4QContent KpQC(0,1,0,0,0,1);
  static const G4QContent zeroQC(0,0,0,0,0,0);
  static const G4LorentzVector nothing(0.,0.,0.,0.);
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mCarb= G4QPDGCode(2112).GetNuclMass(6,6,0);
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();
  static const G4double mEta = G4QPDGCode(221).GetMass();
  static const G4double mEtaP= G4QPDGCode(331).GetMass();
  static const G4double mPi2 = mPi * mPi;
  static const G4double mPi02= mPi0* mPi0;
  static const G4double mK2  = mK  * mK;
  static const G4double mK02 = mK0 * mK0;
  static const G4double mEta2= mEta*mEta;
  static const G4double mEP2 = mEtaP*mEtaP;
  static const G4double CBKinMin = 0.1;          // .1 MeV of KinEnergy is needed to try to pen.
  static const G4double alpha = 1./137.0359895;  // Fine-structure constant
  static const G4double third =1./3.;
  static const G4double conCon=197.327;          // Conversion constant (hc)
  static const G4double rCB   = 1.09;            // R=r_CB*(a^1/3+A^1/3) - CoulBarrierRadius(fm)
#ifdef pdebug
  G4cout<<"G4Quasm::HadrQuasm: ***===>>> START QUASMON HADRONIZATION <<<===***"<<G4endl;
#endif
  theEnvironment=qEnv;
  G4double momPhoton=phot4M.rho();
  G4double addPhoton=phot4M.e();
  if(!status||q4Mom==nothing)                    // This Quasmon is done (Sould not be here)
  {
#ifdef ppdebug
    G4cout<<"G4Quasm::HadrQuasm:NOTHING-TO-DO: Q4M="<<q4Mom<<", QEnv="<<theEnvironment<<G4endl;
#endif
    if(addPhoton)
	{
	  G4cerr<<"***G4Quasmon::HQ: Q4M="<<q4Mom<<",status="<<status<<", phE="<<addPhoton<<G4endl;
      G4Exception("G4Quasmon::HadronizeQuasmon: Overhead photon for the Zero Quasmon");
	}
    KillQuasmon();                               // This Quasmon is done
    qEnv=theEnvironment;                         // Update QEnvironment
    return theQHadrons;
  }
  status=2;                                      // Nothing is done yet
  G4int j=0;                                     // counter of Quasmon Hadronization acts
  G4int first=true;                              // First act flag
  G4int retC=0;                                  // Counter of repetitions for the first act
  G4int sPDG=0;                                  // Prototype of PDG of a Selected Candidate
  G4int pPDG=0;                                  // ProtTemporary PDG Code of the Parent Cluster
  G4double fmh=false;                            // Flag of hadronization in nuclear matter
  G4double rmM=0.;                               // Prototype of coalescence mass
  G4double npqp2=0;                              // A#of quark-partons -2 in a selected fragment
  G4double sMass=0.;                             // Mass of selected candidate
  G4double sM2=0.;                               // Squared mass of selected candidate
  G4int    pBaryn=0;                             // Parent cluster'c Baryon Num for sel.fragment
  G4double pMass=0.;                             // Bounded parent cluster mass for sel.fragment
  G4double pM2=0.;                               // Sq. Bounded par.clust. mass for sel.fragment
  G4double delta=0.;                             // Binding energy = (sM2-pM2)/(2*pMass)
  G4double minSqT=0.;                            // Minimal mass of Residual Quasmon
  G4double hili=0.;                              // High limit of quark exchange randomization
  G4double loli=0.;                              // Low limit of quark exchange randomization
  G4QContent curQ(0,0,0,0,0,0);                  // ProtTemporary copy of valQ to estimate MinM2
  G4QContent memQ(0,0,0,0,0,0);                  // ProtTemporary copy of valQ to remember state
  G4QContent pQC(0,0,0,0,0,0);                   // ProtTemporary Quark Content of ParentCluster
  G4QContent sQC(0,0,0,0,0,0);                   // ProtTemporary Quark Content of the fragment
  G4QContent transQC(0,0,0,0,0,0);               // ProtTemporary Quark Content of ExchangeMeson
  G4LorentzVector m4Mom(0.,0.,0.,0.);            // 4Momentum to memorize a Quasmon's 4-momentum
  G4LorentzVector kp4Mom(0.,0.,0.,0.);           // 4-mom prototype for kappa (recoil q)
  G4LorentzVector check=-theEnvironment.Get4Momentum()-q4Mom;//4Mom sum to check conservation
  G4int ccheck=-theEnvironment.GetZ()-valQ.GetCharge();      //For checking charge conservation
  G4bool start=true;                             //@@ start=true in case of j=0 decision @@
  while (q4Mom.m2()>-.000001&&(theEnvironment.GetPDG()==NUCPDG||start))//=***=The MAIN LOOP=***=
  {
    start         =false;                        //@@ In matter - only one hadronization attempt
    G4bool   quexf=false;                        // Flag of successful quark exchange
    G4double qM2  = q4Mom.m2();                  // Current squared mass of Quasmon 
    G4double quasM= sqrt(qM2);                   // Current mass of Quasmon 
    G4ThreeVector qltb = q4Mom.boostVector();    // Boost vector for backward Lor.Trans. in LS
    G4double b2=qltb.mag2();                     // beta^2 of Quasmon
    CalculateNumberOfQPartons(quasM);            // Fills private parameter nOfQ (#of QPartons)
    //=======================
    G4int envPDG=theEnvironment.GetPDG();        // PDGCode of the current Nuclear Environment
    G4int envN  =theEnvironment.GetN();          // N of current Nuclear Environment
    G4int envZ  =theEnvironment.GetZ();          // Z of current Nuclear Environment
    G4double envM=theEnvironment.GetMass();      // mass of the current Nuclear Environment
    G4QContent envQC=theEnvironment.GetQCZNS();  // QuarkContent of the current Nuclear Environ
#ifdef pdebug
	G4cout<<"G4Quasmon::HadrQ:"<<j<<",ePDG="<<envPDG<<",eM="<<envM<<",eQC="<<envQC<<G4endl;
#endif
    G4int envBarN=envQC.GetBaryonNumber();       // Baryon Number of the current Nuclear Environ
    //theEnvironment.UpdateClusters(nBarClust);    // @@ new Clusters for reduced nucleus @@
    G4QContent totQC=valQ+envQC;                 // Total Quark Content
    G4int      totBN=totQC.GetBaryonNumber();    // Total Baryon Number of the Total System
    G4int      totS=totQC.GetStrangeness();      // Total Strangeness of the Total System
    G4QNucleus totN(totQC);                      // Pseudo nucleus for the Total System
    G4int      totNeut=totN.GetN();              // Total number of neutrons in the system
    G4int      totProt=totN.GetZ();              // Total number of protons in the system
    G4double totM  =totN.GetMZNS();              // min mass of the Total System
#ifdef pdebug
	G4cout<<"G4Quasmon::HadrQ:"<<j<<",tN="<<totN<<",tGSM="<<totM<<",tQC="<<totQC<<G4endl;
#endif
    G4double protCB=CoulombBarrier(totProt,totBN,1.,1.);
    G4int    resNPDG=0;
    G4double resNM =10000000.;                   // Prototype of residual mass after n separated
    if(totNeut>0)
	{
      G4QContent resNQC=totQC-G4QContent(2,1,0,0,0,0);
      G4QNucleus resNN(resNQC);
      resNM  = resNN.GetMZNS();
      resNPDG= resNN.GetPDG();
	}
    G4int    totPDG=totN.GetPDG();               // Total PDG Code for the Current compound
    G4LorentzVector tot4M =q4Mom+G4LorentzVector(0.,0.,0.,envM);
    G4double totMass=tot4M.m();
#ifdef pdebug
	G4cout<<"G4Quasmon::HadrQ:"<<j<<",tPDG="<<totPDG<<",tM="<<totMass<<",rPDG="<<resNPDG<<G4endl;
#endif
    G4ThreeVector totBoost = tot4M.boostVector();// Boost vector for Total System (backward)
    G4ThreeVector totRBoost= -totBoost;          // Boost vector for Total System (forward)
    G4int    iniPDG =valQ.GetSPDGCode();
    G4int    iniBarN=valQ.GetBaryonNumber();
    G4int    iniQChg=valQ.GetCharge();
    G4int    iniN   =valQ.GetN();
    G4int    iniP   =valQ.GetP();
    G4int    iniS   =valQ.GetL();
#ifdef pdebug
	G4cout<<"G4Quasmon::HadrQ:"<<j<<",iniPDG="<<iniPDG<<",Z="<<iniP<<",N="<<iniN<<",S="<<iniS<<G4endl;
#endif
    G4QNucleus iniRN(iniP,iniN-1,iniS);
#ifdef pdebug
	G4cout<<"G4Quasmon::HadrQ:"<<j<<",iniRN="<<iniRN<<G4endl;
#endif
    G4double iniRM = iniRN.GetMZNS();            // Mass of Residual Quasmon when neutron is rad
    if(iniBarN<2||envBarN>0) iniRM=0.;
#ifdef pdebug
	G4cout<<"G4Quasmon::HadrQ:"<<j<<",iniRM="<<iniRM<<G4endl;
#endif
    G4double iniQM =G4QPDGCode(iniPDG).GetMass();// Minimum mass of Quasmon
#ifdef pdebug
	G4cout<<"G4Quasmon::HadrQ:"<<j<<",iniQM="<<iniQM<<G4endl;
#endif
    G4double iniQM2= iniQM*iniQM;
    G4double bndQM = totM-envM;
    if(envPDG==NUCPDG) bndQM=iniQM;
    G4double bndQM2= bndQM*bndQM;
    G4double quen  =iniQM+envM;
    G4double freeE =(totMass-totM)*iniQM;
#ifdef pdebug
    cout<<"G4Quasm::HadrQuasm:j="<<j<<",mQ="<<quasM<<valQ<<bndQM2<<",nQ="<<nOfQ<<", Env="<<envM
        <<theEnvironment<<envQC<<",Q+E="<<quen<<",tM="<<totM<<totQC<<totPDG<<"<"<<totMass<<endl;
#endif
    if(envPDG>80000000&&envPDG!=90000000&&totM>totMass&&quen<totMass)// Decay in Quasm and Environment    
	{ // @@ Can not evaporate as it is under the ground state of the nucleus (?)
#ifdef ppdebug
      G4cout<<"G4Quasm::HadrQuasm:UNDER-MASS-SHELL: Q4M="<<q4Mom<<", QEnv="<<theEnvironment<<G4endl;
#endif
      // @@ Temporary put the panic flag
      status=-1;
      // @@ End of temporary solution
      if(addPhoton)
      {
        q4Mom+=phot4M;
        addPhoton=0.;                              // the Photon is adopted: clean it up
        momPhoton=0.;
        phot4M=G4LorentzVector(0.,0.,0.,0.);
	  }
      qEnv=theEnvironment;
      return theQHadrons;
	}
    else if(totBN>1&&iniQM>quasM&&totMass>totM&&iniPDG!=1114&&iniPDG!=2224&&totS>=0&&envPDG>80000000&&envPDG!=90000000)
	{
#ifdef ppdebug
      G4cout<<"G4Quasm::HadrQuasm:NEEDS-EVAPORATION-1: Q4M="<<q4Mom<<", QEnv="<<theEnvironment<<G4endl;
#endif
      if(addPhoton)
      {
        q4Mom+=phot4M;
        addPhoton=0.;                              // the Photon is adopted: clean it up
        momPhoton=0.;
        phot4M=G4LorentzVector(0.,0.,0.,0.);
	  }
      qEnv=theEnvironment;
      return theQHadrons;
	}
    G4int tQ    = valQ.GetTot();                 // Total number of quarks for current Quasmon
    G4int bQ    = abs(valQ.GetBaryonNumber());   // Baryon number of the current Quasmon
    G4QContent cQ = valQ;                        // Temporary copy of Quasmon QC
    G4int   s   = 4;                             // Mesonic
    if (bQ) s   = 3*bQ + 2;                      // Barionic
    if (tQ> s) cQ.DecQAQ((tQ-s)/2);              // Reduce QC to minimum QC
    G4int rsPDG = cQ.GetSPDGCode();              // PDG for the lowest residual Quasmon state
#ifdef pdebug
    cout<<"G4Quasmon::HadrQuasmon:j="<<j<<",eN="<<envN<<",eZ="<<envZ<<",QPDG="<<rsPDG<<cQ<<endl;
#endif
    PrepareCandidates(j);                        // Calculate Preprobabilities
    ModifyInMatterCandidates();                  // Calculate InMediaMasses of Cand. & Possibil.
    G4double kMom = 0.;                          // Energy of primary qParton in Q-CMS
    G4double kLS  = 0;                           // Energy of primary qParton in LS
    G4double cost = 0.;                          // Cos(theta) of k in QS InRespecTo Q direction
    G4bool   cond = true;                        // Adoptable k condition
    G4bool   fskip=false;                        // Flag to skip when sucked
    G4bool   fred =false;                        // Flag of Environment reduction
    G4LorentzVector k4Mom(0.,0.,0.,0.);          // 4-momentum prototype for k
    G4LorentzVector cr4Mom(0.,0.,0.,0.);         // 4-momentum prototype for the ColResQuasmon
    G4int counter =0;                            // Counter of attempts of k for hadronization
    //G4int maxCount=27;
    G4int maxCount=7;
    //G4int maxCount=3;
    //G4int maxCount=1+static_cast<int>(pow(2,j));
    G4double kTmp=0.;
    //if(addPhoton>0&&j&&fmh)
	if(2>3)                                      // @@@@@@@@@@@@@@@@@@@@@@@@@@@
	{
      cond=false;
      k4Mom=kp4Mom;
      kLS=k4Mom.e();
      if(b2>0.000001)                            // Not zero speed of Quasmon
	  {
        G4double beta=sqrt(b2);
        G4ThreeVector qv = k4Mom.vect();         // 3D vector for q
        G4double ct=qv.dot(qltb)/beta/qv.mag();  // cos(theta) between Q & q
        kMom = kLS*(1.-beta*ct)/sqrt(1.-b2);     // k* in Quasmon CMS
	  }
#ifdef pdebug
	  G4cout<<"G4Quasmon::HadrQuasmon:j="<<j<<",k=q="<<kLS<<",kQCM="<<kMom<<G4endl;
#endif
	}
    while(counter<maxCount&&cond)
	{
#ifdef pdebug
	  G4cout<<"G4Quasmon::HadrQuasmon: START LOOP j="<<j<<G4endl;
#endif
      cond=true;
      //if(addPhoton>0.) nOfQ=valQ.GetQ()-valQ.GetAQ();
      //if(addPhoton>0.&&quasM<1500.&&G4UniformRand()<f2all)
      if(addPhoton>0.&&quasM<1500.&&G4UniformRand()<0.5)
      {
        kMom=(1.-pow(G4UniformRand(),1./static_cast<double>(nOfQ-2)))*quasM/2.;
	    cost=addPhoton/momPhoton;
	  }
      else
      {
        kMom = GetQPartonMomentum(0.,0.);          // Calculate value of primary qParton
        cost = 1.-2.*G4UniformRand();
	  }
      kTmp=kMom;
      kLS = kMom;
      if(b2>0.000001)                            // Not zero speed of Quasmon
	  {
        G4double beta=sqrt(b2);
        kLS = kMom*(1.+beta*cost)/sqrt(1.-b2);   // k in Laboratory System
	  }
      if(addPhoton>0.) kLS=sqrt(momPhoton*(2*kLS*cost+momPhoton)+kLS*kLS);
      G4double cQM2=qM2-(kTmp+kTmp)*quasM;
      G4double cQM=sqrt(cQM2);                   // Mass of the coloured residual Quasmom
      k4Mom=G4LorentzVector(0.,0.,0.,0.);
      cr4Mom=G4LorentzVector(0.,0.,0.,cQM);
      if(!G4QHadron(q4Mom).RelDecayIn2(k4Mom, cr4Mom, q4Mom, cost, cost))// Q->ColResidQ+k
      {
#ifdef pdebug
         G4cerr<<"***G4Quasmon::HadrQ: QM="<<quasM<<",cQM="<<cQM<<",cost="<<cost<<G4endl;
#endif
         cond=true;
      }
      else
      {
#ifdef pdebug
	    cout<<"G4Quasmon::HadrQuasm: j="<<j<<",k="<<kMom<<",tk="<<kTmp<<",k4M="<<k4Mom<<endl;
#endif
        G4double quasMt=0.;
        if(addPhoton>0.)
        {
          G4LorentzVector tmp4M=q4Mom+phot4M;    // temporary Quasmon
          quasMt= sqrt(tmp4M.m2());              // Temporary Current mass of Quasmon 
        }
        G4double rEn=cr4Mom.e();
        rMo=cr4Mom.rho();                        // p for the ResidualColouredQuasmon in LS
        G4double rEmP=rEn-rMo;                   // E-p for the ResidualColouredQuasmon in LS
        rEP=rEn+rMo;                             // E+p for the ResidualColouredQuasmon in LS
        G4int totCand = theQCandidates.entries();// Total number of candidates
#ifdef sdebug
	    cout<<"G4Quasmon::HadrQuasm: j="<<j<<":it#"<<counter<<",k="<<kMom<<k4Mom<<endl;
#endif
        for (G4int index=0; index<totCand; index++)
        {
          G4QCandidate* curCand=theQCandidates[index];
          G4int cPDG = curCand->GetPDGCode();
          if(cPDG==90000001||cPDG==90001000||cPDG==91000000||cPDG<80000000)//@@ k-PreAttempts(Acc)
          {
            G4bool poss= curCand->GetPossibility();
#ifdef pdebug
		    if(cPDG==90000001 || cPDG==90001000 || cPDG==90000002 || cPDG==90001001)
	          cout<<"G4Quasmon::HadrQuasm:pos="<<poss<<",cPDG="<<cPDG<<",iQC="<<iniQChg<<endl;
#endif
		    if(poss)
            {
              G4double   cMs=curCand->GetMass();  // Bound mass of ParentCluster(=Fragment)
              G4QContent cQC=curCand->GetQC();
              G4double   cfM=curCand->GetQPDG().GetMass();
              G4QContent rtQC=curQ+envQC-cQC;     // TotalResidualQuarkContent - OutFragm
              G4QNucleus rtN(rtQC);               // Create a pseudo-nucleus for residual
              G4double totr = rtN.GetMZNS();      // Mass of total Residual
              G4double bnM  = totr-envM+cMs;      // Bound mass of residual Quasmon
              G4double pmk  = rMo*cMs/kLS;
#ifdef pdebug
		      if(cPDG==90000001 || cPDG==90001000 || cPDG==90000002 || cPDG==90001001)
			    cout<<"G4Quasmon::HadrQ:cfM="<<cfM<<",cMs="<<cMs<<",ind="<<index<<endl;
#endif
              G4double bnM2=bnM*bnM;
              G4double k = kMom;
              if(cPDG>80000000&&cPDG!=90000000) k=kLS; // ===> Nuclear case (Lab System)
              G4double kMin=(cfM*cfM-cMs*cMs)/(cMs+cMs);
              if(kMin<0.) kMin=0.;
              //G4double dR=(bnM*bnM-cQM2)/(rEP+rEP);
              //if(b2<.000001) dR=((bnM*bnM-cQM2)/2.+pmk*(kLS-kMin))/(rEP+pmk);
              //if (dR<0.) dR=0.;
              G4double dR=0.;
#ifdef pdebug
     		  if(cPDG==90000001 || cPDG==90001000 || cPDG==90000002 || cPDG==90001001)
	            cout<<"G4Quasmon::HadrQ:i="<<index<<",cPDG="<<cPDG<<",k="<<kMom<<",kLS="<<kLS
                    <<",bM="<<bnM<<",rEP="<<rEP<<",d="<<kMin<<",dR="<<dR<<",c="<<cond<<endl;
#endif
              ///***Cap***???
			  //if(cPDG<80000000&&k>kMin&&quasM>cMs+iniQM||cPDG>80000000&&dR<kLS-kMin)cond=false;
			  if(cPDG<80000000&&k>kMin&&quasMt>cMs+iniQM||cPDG>80000000&&cPDG!=90000000&&dR<kLS-kMin
                 && kLS-freeE/(iniQM+cMs)<cMs/2.)cond=false;
		    }
		  }
	    }
	  }
      counter++;
      //if(counter>=maxCount)cerr<<"***G4Quasmon::HadronizeQuasmon:j="<<j<<",c="<<counter<<endl;
	}
    if(addPhoton>0.&&first)
    {
      first=false;
#ifdef pdebug
	  cout<<"G4Quasmon::HadrQuasmon:j="<<j<<",k="<<k4Mom<<",Q="<<q4Mom<<",Ph="<<phot4M<<endl;
#endif
      k4Mom+= phot4M;                            // add photon energy to radiated quark
      q4Mom+= phot4M;                            // tmp Quasmon with added photon
      qM2   = q4Mom.m2();                        // NEW Current squared mass of Quasmon 
      quasM = sqrt(qM2);                         // NEW Current mass of Quasmon 
      tot4M+= phot4M;
      totMass = tot4M.m();                       // NEW Total mass of compound
      totBoost  = tot4M.boostVector();           // NEW Boost vector for Total System (backward)
      totRBoost = -totBoost;                     // NEW Boost vector for Total System (forward)
      addPhoton=0.;                              // the Photon is adopted: clean it up
      momPhoton=0.;
      phot4M=G4LorentzVector(0.,0.,0.,0.);
    }
    if(envPDG>80000000&&envPDG!=90000000&&counter>=maxCount&&j&&totS>=0)
	{
#ifdef ppdebug
      G4cout<<"G4Quasm::HadrQuasm:NEEDS-EVAPORATION-2: Q4M="<<q4Mom<<", QEnv="<<theEnvironment<<G4endl;
#endif
      if(addPhoton)
      {
        q4Mom+=phot4M;
        addPhoton=0.;                              // the Photon is adopted: clean it up
        momPhoton=0.;
        phot4M=G4LorentzVector(0.,0.,0.,0.);
	  }
      qEnv=theEnvironment;
      return theQHadrons;
	}
    j++;                                         // Increment Quasmon Decay Counter
    G4double dk = kMom+kMom;                     // Double QCM k-value
#ifdef pdebug
    cout<<"G4Quasm::HadrQ:===>"<<j<<": kMom="<<kMom<<k4Mom<<",t="<<theQCandidates[64]->GetMass()
        <<endl;
#endif
    CalculateHadronizationProbabilities(kMom,kLS,j); //Mass of resonance is randomized here
    G4int nCandid = theQCandidates.entries();
    G4bool fprob = true;                         // Flag of existing decay probability
    ////G4bool fdul  = true;                        // Prototype of flag of resonance decay
    G4bool fdul  = false;                        // Prototype of flag of resonance decay
    int i=0;                                     // "i" will point to the selected candidate
    G4double maxP = theQCandidates[nCandid-1]->GetIntegProbability();
#ifdef pdebug
    cout << "G4Quasmon::HadronizeQuasmon: NofCandidates="<<nCandid<<", maxP="<<maxP<<endl;
#endif
    if (maxP<=0.)                                // No possible channels for this k value
    {
      if(j<2)                                    // First act
	  {
        j=0;
        fprob = false;                           // Probab selection did not succeed
      }
      else
      {
        G4double smP = theQCandidates[nCandid-1]->GetSecondIntProb();
        if (smP<=0.)                               // Impossible to hadronize in anything
        {
#ifdef pdebug
          cout<<"***G4Quasmon::HadronizeQuasmon:fMaxP="<<maxP<<",sMaxP="<<smP<<",QM="<<quasM
              <<",QC="<<valQ<<endl;
#endif
          fprob = false;                           // Probab selection did not succeed
	    }
        G4double totP = smP * G4UniformRand();
        while(theQCandidates[i]->GetSecondIntProb() < totP) i++;
	  }
    }
    else
	{
      G4double totP = maxP * G4UniformRand();
      while(theQCandidates[i]->GetIntegProbability() < totP) i++;
	}
#ifdef pdebug
    cout<<"G4Quasmon::HQ:fp="<<fprob<<",i="<<i<<",QM="<<quasM<<",QQC="<<valQ<<",k="<<kMom<<endl;
#endif
    if (i>=nCandid) cerr << "***G4Quasmon::HadronizeQuasmon: Cand#"<<i<<">=Tot#"<<nCandid<<endl;
    G4bool hsflag=false;                         // Prototype of H+S decay flag
    G4bool fchipo=false;                         // Final decay of Quasmon-Chipolino
    curQ = valQ;                                 // Temporary copy of valQ to estimate MinM2
    G4int rPDG=0;                                // Prototype of the residual Quasmon PDG
    G4double rMass=0.;                           // Prototype of the residual Quasmon mass
    G4LorentzVector rQ4Mom(0.,0.,0.,0.);         // 4-momentum of residual Quasmon
    G4LorentzVector fr4Mom(0.,0.,0.,0.);         // 4-momentum prototype for the fragment
    G4double totCM=0.;
    G4double totCGSM=BIG;                        // Prototype of Min Mass Of the Total Residual
    G4double totCPDG=0;                          // Prototype of PDG Code Of the Total Residual
    if (!fprob)                                  // Prepare final decay (nothing was selected)
	{
      if(j)                                      // normally j>0 except for "nothing is found"
      {
        if(totBN>1&&totMass>totM&&totS>=0&&envPDG>80000000&&envPDG!=90000000)
	    {
#ifdef ppdebug
          G4cout<<"G4Quasm::HadrQuasm:NEEDS-EVAPORATION-3: Q4M="<<q4Mom<<", QEnv="<<theEnvironment<<G4endl;
#endif
          if(addPhoton)
          {
            q4Mom+=phot4M;
            addPhoton=0.;                              // the Photon is adopted: clean it up
            momPhoton=0.;
            phot4M=G4LorentzVector(0.,0.,0.,0.);
	      }
          qEnv=theEnvironment;
          return theQHadrons;
	    }
        else
        {
          memQ=curQ;                             // Remembe QC for the case of 2H decay
          sPDG=curQ.GetSPDGCode();               // PDG of current Quasmon as S-hadron      
          if(!sPDG)
          {
		    G4cerr<<"***G4Quasmon::HQ: QuasmQC="<<curQ<<",MQ="<<quasM<<",sPDG=0,j="<<j<<G4endl;
            G4Exception("G4Quasmon::HadronizeQuasmon: Can not decay the Quasmon");
          }
          else if(sPDG==10)                      // ---> Chipolino case
	      {
            fchipo=true;
            G4QChipolino chipQ(valQ);
            G4QPDGCode QPDG1=chipQ.GetQPDG1();
            sPDG = QPDG1.GetPDGCode();
            sMass= QPDG1.GetMass();
            G4QPDGCode QPDG2=chipQ.GetQPDG2();
            rPDG = QPDG2.GetPDGCode();
            rMass= QPDG2.GetMass();
	      }
          else                                   // ---> @@ Final decay in MinHadr + PI0
          {
            sMass=G4QPDGCode(sPDG).GetMass();
            rPDG=envPDG;
            if (rPDG>80000000&&rPDG!=90000000)
		    {
              rMass=theEnvironment.GetMZNS();
              q4Mom+=G4LorentzVector(0.,0.,0.,rMass);
              valQ +=theEnvironment.GetQC();
              quasM=q4Mom.m();
              if(envPDG>80000000&&envPDG!=90000000)
              {
                theEnvironment.Reduce(rPDG);     //Update NuclEnvironment
                fred=true;                       // Flag of reduction
              }
		    }
            else
            {
#ifdef pdebug
              G4cout<<"G4Quasmon::HadronizeQuasmon: DECAY IN PI0 : rPDG="<<rPDG<<G4endl;
#endif
              rMass=mPi0;                        //@@ Safety decay in pi0 ??
              rPDG=111;
            }
          }
          hsflag=true;                           // Two particle decay is forced ???
#ifdef pdebug
          G4cout<<"G4Quas::HadrQ:sPDG="<<sPDG<<",sM="<<sMass<<",rPDG="<<rPDG<<",rM="<<rMass<<G4endl;
#endif
		}
	  }
	}
    else
	{
      G4QCandidate* curCand = theQCandidates[i]; // Pointer to selected candidate (hadr./fragm.)
      sPDG  = curCand->GetPDGCode();             // PDG of the selected candidate 
      G4double prpr=curCand->GetPreProbability();
      //@@ For clusters should be another randomization & a loop over possible parent clusters
      if(sPDG>80000000&&sPDG!=90000000)          // ===> Fragment
	  {
        int ip=0;
        G4int nParCandid = curCand->GetPClustEntries();
        G4double sppm  = curCand->TakeParClust(nParCandid-1)->GetProbability();
        if (sppm<=0)                             // Impossible to find a parent cluster
        {
          cerr<<"***G4Quasmon::HadronizeQuasmon: cP="<<theQCandidates[i]->GetIntegProbability()
              <<",nC="<<nParCandid<<",pP="<<sppm<<",QM="<<quasM<<",QC="<<valQ;
          for(int ipp=0; ipp<nParCandid; ipp++)
            cerr<<", "<<ipp<<": "<<curCand->TakeParClust(ip)->GetProbability();
          cerr<<endl;
          G4Exception("G4Quasmon::HadronizeQuasmon: No parent cluster for the fragment");
	    }
        else
	    {
          G4double totPP = sppm * G4UniformRand();
          while(curCand->TakeParClust(ip)->GetProbability() < totPP) ip++;
#ifdef pdebug
          cout<<"G4Quasmon::HadrQuasm:p#ip="<<ip<<",f#i="<<i<<",tP="<<totPP<<",sP="<<sppm<<endl;
#endif
	    }
        G4QParentCluster* parCluster=curCand->TakeParClust(ip);
        pPDG  = parCluster->GetPDGCode();
	    G4QPDGCode pQPDG(pPDG);
        pQC   = pQPDG.GetQuarkContent();
        pBaryn= pQC.GetBaryonNumber();
	    pMass = parCluster->GetMass();
        pM2   = pMass*pMass;
	    transQC = parCluster->GetTransQC();
	    delta = parCluster->GetBind();
	    loli  = parCluster->GetLow();
	    hili  = parCluster->GetHigh();
        //G4double dhil=.0001*(hili-loli);         // Safety factor
        //loli += dhil;
        //hili -= dhil;
	    npqp2 = parCluster->GetNQPart2();
        // @@ Here one can get other useful parameters of the parent cluster for hadronization
	    G4QPDGCode sQPDG(curCand->GetPDGCode());
        sQC   = sQPDG.GetQuarkContent();
        //if(sPDG==90001001 && G4UniformRand()>0.75) sMass=np; //@@ n-p pair
        if(sPDG==90001001 && G4UniformRand()>1.0) sMass=np; //@@ no n-p pair
        else                                      sMass = sQPDG.GetMass();
        sM2   = sMass*sMass;
        curQ += transQC;                         // Subtract ExchangeMesonQC from QC of Quasmon
#ifdef pdebug
        cout<<"G4Quasmon::HadrQuasm:valQ="<<valQ<<" + tQ="<<transQC<<"("<<pPDG<<" to "<<sPDG
            <<") = curQ="<<curQ<<", bind="<<delta<<",pM="<<pMass<<pQC<<endl;
#endif
	  }
      else                                       // ===> Hadron
      {
        pBaryn=1;
        sMass = theQCandidates[i]->GetMass();    // Mass is randomized on probability level
        sM2=theQCandidates[i]->GetMass2();       // Sq. Mass is randomized on probability level
        curQ-= theQCandidates[i]->GetQC();       // Subtract outHadron QC from QC of Quasmon
#ifdef pdebug
        cout<<"G4Quasmon::HadronizeQuasmon: valQ="<<valQ<<" - sQ="<<theQCandidates[i]->GetQC()
            <<" (PDG="<<theQCandidates[i]->GetPDGCode()<<") = curQ="<<curQ<<",sM2="<<sM2<<endl;
#endif
	  }
      memQ=curQ;                                 // Remembe QC for the case of 2H decay
      G4double rtM=0.;
      G4double reM=0.;
      G4int rqPDG = curQ.GetSPDGCode();
      G4int rQQ=G4QPDGCode(curQ).GetQCode();
      if(rqPDG==10) minSqT=G4QChipolino(curQ).GetMass2();// Minimal mass of DoubleHadron
      else if(!rqPDG||rQQ<-1)
	  {
        cerr<<"***G4Quasmon::HadronizeQuasmon:ResQuasm PDG="<<rqPDG<<curQ<<rQQ<<endl;
        minSqT=1000000.;
	  }
      else
      {
        G4int baryn=curQ.GetBaryonNumber();
        G4double minT=G4QPDGCode(rqPDG).GetMass();
        //if(envPDG>pPDG&&rqPDG>80000000&&sPDG>80000000)
        //if(((sPDG<80000000&&envPDG>80000000)||(sPDG>80000000&&envPDG>pPDG))&&rqPDG>80000000)
        if((sPDG<80000000&&envPDG>80000000&&envPDG!=90000000||sPDG>80000000&&sPDG!=90000000&&envPDG>pPDG)
            && (rqPDG>80000000&&rqPDG!=90000000||rqPDG==2112||rqPDG==2212||rqPDG==3122) )
        {
          G4double newT=0.;
          if(sPDG<80000000)
		  {
            G4QContent rtQC=curQ+envQC;           // Total Residual Quark Content
            G4QNucleus rtN(rtQC);                 // Create a pseudo-nucleus for residual
            rtM=rtN.GetMZNS();                    // Min Mass of total residual Nucleus
            newT=rtM-envM;
#ifdef pdebug
		    G4cout<<"G4Quasmon::HadrQuasm:vacFrag:M="<<newT<<",rM="<<rtM<<rtQC
                  <<",eM="<<envM<<",mM="<<minT<<G4endl;
#endif
		  }
          else
		  {
            G4QContent reQC=envQC-pQC;            // Total Residual Quark Content
            G4QNucleus reN(reQC);                 // Create a pseudo-nucleus for residual
            reM=reN.GetMZNS();                    // Min Mass of residual Environment
            G4QContent rtQC=curQ+reQC;            // Total Residual Quark Content
            G4QNucleus rtN(rtQC);                 // Create a pseudo-nucleus for residual
            rtM=rtN.GetMZNS();                    // Min Mass of total residual Nucleus
            newT=rtM-envM+pMass;
#ifdef pdebug
		    G4cout<<"G4Quasmon::HadrQuasm:nucFrag:M="<<newT<<",rM="<<rtM<<rtQC
                  <<",eM="<<envM<<",pM="<<pMass<<G4endl;
#endif
		  }
          if(newT<minT) minT=newT;
	    }
        minSqT=minT*minT;
      }
#ifdef pdebug
	  G4cout<<"G4Quasmon::HadrQuasm:rqPDG="<<rqPDG<<",mST="<<minSqT<<",rtM="<<rtM<<G4endl;
#endif
      if(minSqT==0.)
      {
        cerr<<"***G4Quasmon::HadronizeQuasmon: minSqT=0(!), curQC="<<curQ<<endl;
        G4Exception("G4Quasmon::HadronizeQuasmon: Minimal Residual Mass can not be calculated");
      }
      G4double m2 = BIG2; //@@ just a BIG number // Prototype/Squared Mass of Residual Quasmon
      G4double kt=0.;                            // Prototype of (Max Residual Quasmon Mass)^2
      if(sPDG>80000000&&sPDG!=90000000)          // ===> Nuclear fragment candidate hadronization
	  {
#ifdef pdebug
	    cout<<"G4Quasmon::HadrQuasm: BoundM="<<pMass<<",FragM="<<sMass<<",QM="<<quasM<<endl;
#endif
        // = = = =  P u r e   k i n e m a t i c a l   c a l c u l a t i o n s:  = = = = =
        // Fusion of k + parentCluster => colouredCluster (cc)
        G4LorentzVector cl4Mom(0.,0.,0.,pMass);       // 4-momentum prototype for parent cluster
        G4LorentzVector tot4Mom=q4Mom+cl4Mom;         // @@ Just for checking
#ifdef pdebug
	    cout<<"G4Quasmon::HadrQuasm:Q("<<quasM<<")->k("<<k4Mom<<")+CRQ("<<cr4Mom.m()<<")"<<endl;
#endif
        G4LorentzVector cc4Mom=k4Mom+cl4Mom;          // 4-mom of ColoredFragment (before kappa)
        G4double ccM2=cc4Mom.m2();
        G4double frM2=sMass*sMass;
        if(ccM2<=frM2)
        {
#ifdef pdebug
	      cout<<"***G4Q::HQ:Failed with Fragment Mass:"<<ccM2<<"< "<<frM2<<", pM="<<pMass
              <<" + k="<<k4Mom.m()<<k4Mom<<" = "<<sqrt(ccM2)<<cc4Mom<<" < fm="<<sMass<<endl;
#endif
          hsflag=true;
		}
        else
        {
          G4double crMass2 = cr4Mom.m2();             // SquredMass of ColouredResidualQuasmon
          G4double newh=0.;
          if(hili<loli) hili=loli;
          G4double pw= 1./static_cast<double>(npqp2);
          // ----------->>>>>Decay of the ColouredCluster in a Fragment + kappa (fixed ctc)<<<
          G4QContent totCQC=curQ+envQC-pQC;           // Total Residual Quark Content
          G4QNucleus totCN(totCQC);                   // Pseudo Nucleus for the Total Residual
          totCGSM=totCN.GetMZNS();                    // Min Mass Of the Total Residual
          totCPDG=totCN.GetPDG();                     // PDGCode Of the Total Residual
#ifdef pdebug
	      cout<<"G4Quasmon::HadronizeQuasmon:kt="<<kt<<" < minSqT="<<minSqT<<"? totCM="<<totCM
          <<" < totCGSM="<<totCGSM<<"?"<<endl;
#endif
          G4double z = pow(loli+(hili-loli)*G4UniformRand(),pw);
	      G4double ex=kLS-delta;                      // Excess of parton energy in LS
          G4double ctkk=1.-((ex+ex)/(1.-z)-pMass)/kLS;// cos(theta_k,kappa) in LS
#ifdef pdebug
	      cout<<"G4Q::HQ:ct="<<ctkk<<",pM="<<pMass<<",z="<<z<<",zl="<<pow(loli,pw)<<",zh="
              <<pow(hili,pw)<<",ex="<<ex<<",li="<<loli<<",hi="<<hili<<endl;
#endif
          if(abs(ctkk)>1.)
          {
#ifdef pdebug
            G4cout<<"***G4Q:HadQ:ctkk="<<ctkk<<",ex="<<ex<<",z="<<z<<",pM="<<pMass
                  <<",new="<<newh<<",hi="<<hili<<",lo="<<loli<<",n="<<npqp2<<G4endl;
#endif
            if(ctkk> 1.)ctkk= 1.;
            if(ctkk<-1.)ctkk=-1.;
          }
          G4double kappa=ex/(1.-kLS*ctkk/pMass);    // Energy of the recoil quark
          G4double cen=kLS+pMass;                   // LS Energy of k+parentCluster CompSystem
          G4double ctc=(cen*ctkk-kLS)/(cen-kLS*ctkk);//cos(theta_k,kap) in k+pClast CompSystem
          if(abs(ctc)>1.)cerr<<"**G4Quasmon:HadrQuasm:e="<<cen<<",k="<<kLS<<",c="<<ctc<<endl;
          kp4Mom=G4LorentzVector(0.,0.,0.,0.);      // 4-mom prototype update for kappa (q)
          fr4Mom=G4LorentzVector(0.,0.,0.,sMass);   // 4-momentum prototype for the fragment
          if(!G4QHadron(cc4Mom).RelDecayIn2(kp4Mom, fr4Mom, k4Mom, ctc, ctc))
          {
            cerr<<"***G4Quasmon::HadrQuasm: c4M="<<cc4Mom<<",sM="<<sMass<<",ct="<<ctc<<endl;
            G4Exception("G4Quasmon::HadronizeQuasm: Can't decay FakeColClust in fr+ kappa");
          }
          G4double ccM=sqrt(ccM2);
#ifdef pdebug
	      cout<<"G4Q::HQ: CF("<<ccM<<")->F("<<sMass<<fr4Mom<<")+kp("
              <<kp4Mom<<", "<<kp4Mom.m()<<")"<<endl;
#endif
          fmh=true;
          // Fusion of the ColouredResidualQuasmon + kappa (LS) -> get residual Quasmon mass
          rQ4Mom=cr4Mom+kp4Mom;                      // 4-momentum of residual Quasmon
          quexf=true;                                // Quark Exchange is successfully done
#ifdef pdebug
	      cout<<"G4Q::HQ:RQ("<<rQ4Mom.m()<<")=CRQ("<<cr4Mom.m()<<")+kp("<<kp4Mom<<")"<<endl;
#endif
          kt=rQ4Mom.m2();
          G4LorentzVector totC4Mom=rQ4Mom+G4LorentzVector(0.,0.,0.,envM-pMass); // TotatResidual
          totCM=totC4Mom.m();
          m2=kt;
          tot4Mom-=rQ4Mom+fr4Mom;
#ifdef pdebug
	      cout<<"G4Quasm::HQ: t4M="<<tot4Mom<<".Is "<<kt<<">"<<minSqT<<"?"
              <<",m2="<<m2<<endl;
#endif
          if(kt<minSqT) hsflag=true;
		}
	  }
      else kt = (quasM-dk)*(quasM-sM2/dk);             // ===> "Hadronic candidate" case
      G4LorentzVector tL=rQ4Mom+G4LorentzVector(0.,0.,0.,reM);
      G4double tM=tL.m();
      // Residual S+S limit (absolute low limit for corresponding P-resonance) for R->S+S Decay
#ifdef pdebug
      cout<<"G4Quasmon::HadronizeQuasmon:k="<<kMom<<".Is r2="<<kt<<"> T2="<<minSqT<<" & tM="<<tM
          <<" > rtM="<<rtM<<" to avoid H+S(sPDG="<<sPDG<<"),QM="<<quasM<<",dk="<<dk<<endl;
#endif
      
      if(kt>minSqT+.01&&(sPDG<80000000||tM>rtM)&&!hsflag)// Mass(RQ)>M_min && Mass(RN)>M_min
	  {
        if(sPDG<80000000)                              // Hadronic candidate:finish calculations
	    {
          G4double np = nOfQ - 3;                      // Power for an residual quasmon mass
          G4double cM = pow(minSqT/kt,np);             // Cut for PossibleQuasmon residual mass
          G4double uR = G4UniformRand();
          G4double rn = pow(cM +(1.-cM)*uR, 1./np);
#ifdef pdebug
          cout<<"G4Quasmon::HadronizeQuasmon: YES it is big enough: t="<<kt<<" > T="<<minSqT
              <<", np="<<np<<", cM="<<cM<<", uR="<<uR<<", rn="<<rn<<endl;
#endif
          m2 = kt*rn;                                  // SquaredMass of ResidualQuasmon (hadr)
	    }
        else m2 = kt;                                  // SquaredMass of ResidualQuasmon (fragm)
      }
      else hsflag=true;                                // Decay in ThisHadron/Fragment+(SHadron)
#ifdef pdebug
	  cout<<"G4Quasmon::HadronizeQuasmon:****** rMass="<<rMass<<",sqm2="<<sqrt(m2)<<endl;
#endif
      rMass=sqrt(m2);                                  // if(hsflag)(TwoPartDecay) it's fake 0.
      G4int cB    = abs(curQ.GetBaryonNumber());       // Baryon Number of residual
      G4int cS    = abs(curQ.GetStrangeness());        // Strangenes of residual
      rPDG        = curQ.GetSPDGCode();                // PDG of lowest residual hadronic state
      G4int aPDG  = abs(rPDG);
      G4int rb    = abs(curQ.GetBaryonNumber());       // BaryNum of residual hadronic state
      G4double rcMass=-BIG;     //@@ just a BIG number // Prototype of minimal mass for residual
      if (!rPDG)
      {
        cerr<<"***G4Quasmon::HadronizeQuasmon: rQ ="<<curQ<<", rPDG="<<rPDG<<"(b="<<rb
            <<") + sPDG="<<sPDG<<"(sM="<<sMass<<")"<<endl;
        G4Exception("G4Quasmon::HadronizeQuasmon: unidentifiable residual Hadron");
      }
      G4double sB = 1473.;                             // @@ Mean of DELTA(3/2) & N(1680)(5/2)
      if(rPDG!=10) rcMass=G4QPDGCode(rPDG).GetMass();  // residual mass for not Chipolino case
      else         sB=0.;                              // Chipolino never decays in hadrons
      if(rPDG==221 || rPDG==331) rcMass=mPi0;
      G4double bs = rcMass+mPi0;
      G4bool rFl=false;                                // true: ResidualResonance,false: Quasmon
      if(sPDG<80000000&&envPDG==NUCPDG)rFl=G4UniformRand()<bs*bs/m2;// ProbabFunct: m_min^2/m^2
#ifdef pdebug
	  cout<<"G4Quasmon::HadronizeQuasmon: sPDG="<<sPDG<<", hsflag="<<hsflag<<", rPDG="<<rPDG
          <<curQ<<",rM="<<rMass<<",rb="<<rb<<",rFl="<<rFl<<endl;
#endif
      if(!hsflag && rFl && rPDG && rPDG!=10 && rb<2 && aPDG!=1114 && aPDG!=2224 && aPDG!=3334)
	  { // ------------------------------->>>>>>>>>>>>>>> Hadron-Parton Duality decay
        G4int regPDG = 0;                              // PDG prototype for G-meson
        G4int refPDG = 0;                              // PDG prototype for F-meson
        G4int redPDG = 0;                              // PDG prototype for D-meson
        G4int repPDG = 0;                              // PDG prototype for P-meson
        if(rPDG && rPDG!=10)                           // Can the residual be a Hadron ?
	    {
          if     (rPDG== 3122) rPDG= 3212;             // LAMBDA* converted to SIGMA*
          else if(rPDG==-3122) rPDG=-3212;
          if(rPDG>0)repPDG=rPDG+2;                     // YES: Make P-state out of S-state
          else      repPDG=rPDG-2;                     // Subtract 2 for the negative PDG
          if(repPDG>0)redPDG=repPDG+2;                 // YES: Make D-state out of P-state
          else        redPDG=repPDG-2;                 // Subtract 2 for the negative PDG
          if(redPDG>0)refPDG=redPDG+2;                 // YES: Make F-state out of D-state
          else        refPDG=redPDG-2;                 // Subtract 2 for the negative PDG
          if(refPDG>0)regPDG=refPDG+2;                 // YES: Make G-state out of F-state
          else        regPDG=refPDG-2;                 // Subtract 2 for the negative PDG
          if(rPDG==221)      rPDG=111;                 // eta=>Pi0
#ifdef pdebug
	      cout<<"G4Quasmon::HadronizeQuasmon: H("<<sPDG<<")+RGS("<<rPDG<<") decay of M="<<quasM
              <<"(QC="<<valQ<<")"<<endl;
#endif
          G4double resM2  = G4QPDGCode(  rPDG).GetMass2(); // Mass^2 of the S-resonance
          G4double repM2  = G4QPDGCode(repPDG).GetMass2(); // Mass^2 of the P-resonance
          G4double redM2  = G4QPDGCode(redPDG).GetMass2(); // Mass^2 of the D-resonance
          G4double refM2  = G4QPDGCode(refPDG).GetMass2(); // Mass^2 of the F-resonance
          sB              = sqrt((resM2+repM2)/2.);    // Boundary between S&P resonances
          G4double     pB = sqrt((repM2+redM2)/2.);    // Boundary between P&D resonances
          G4double     dB = sqrt((redM2+refM2)/2.);    // Boundary between D&F resonances
          G4double     fB = sqrt(refM2)+150.;    // Fake boundary between F&G(no BR) resonances
          if(!cB)      fB=  sqrt((refM2+G4QPDGCode(regPDG).GetMass2())/2.);
          G4double dif=quasM-sMass;
          G4double rM = GetRandomMass(repPDG,dif);     // Randomize Mass of P-resonance
          G4double dM = GetRandomMass(redPDG,dif);     // Randomize Mass of D-resonance
          G4double fM = GetRandomMass(refPDG,dif);     // Randomize Mass of F-resonance
#ifdef pdebug
          cout<<"G4Quasmon::HadronizeQuasmon: rM="<<rMass<<", sB="<<sB<<" (bQ="<<bQ<<"), pB="
              <<pB<<", dB="<<dB<<", fB="<<fB<<endl;
#endif
          if(((rM>0 && rMass<pB && rMass>sB) || (dM>0 && rMass>pB && rMass<dB) ||
             (fM>0 && rMass>dB && rMass<fB)) && theEnvironment.GetPDG()==NUCPDG)
		  { // Final H+R decay of QUASMON in vacuum case (should not exist if Environ exists)
            if     (rMass>pB && rMass<dB && dM>0)      // D-resonance case
		    {
              repPDG=redPDG;
              rM=dM;
		    }
            else if(rMass>dB && rMass<fB && dM>0)      // F-resonance case
		    {
              repPDG=refPDG;
              rM=fM;
		    }
#ifdef pdebug
            G4cout<<"G4Quasmon::HadronizeQuasmon: sPDG="<<sPDG<<",Q => rM="<<rMass<<"(rPDGmin="
                  <<rPDG<<curQ<<") + sB="<<sB<<G4endl;
#endif
            if(quasM<rM+sMass &&(sPDG==221||sPDG==331))// Change eta-Candidate to pi-Candidate
	        {
              sPDG=111;
              sMass=mPi0;
  	        }
            G4LorentzVector r4Mom(0.,0.,0.,rM);        // P/D-resonance with a random Mass
            G4LorentzVector s4Mom(0.,0.,0.,sMass);     // Mass's random since probability time
            if(!G4QHadron(q4Mom).DecayIn2(r4Mom, s4Mom))
            {
              G4cerr<<"**G4Quasmon::HadronizeQuasmon: rPD="<<repPDG<<"(rM="<<rMass<<") + sPD="
                  <<sPDG<<"(sM="<<sMass<<")"<<G4endl;
    	      G4Exception("G4Quasmon::HadronizeQuasmon: Hadron+Resonance decay didn't succeed");
            }
#ifdef ppdebug
	        G4cout<<"G4Quasmon::HadronizeQuasmon:=== 1 ===> HVec, Q="<<q4Mom<<" -> s4M="
                  <<s4Mom<<"("<<sPDG<<"), r4M="<<r4Mom<<"("<<repPDG<<")"<<G4endl;
#endif
            G4QHadron* curHadr1 = new G4QHadron(repPDG,r4Mom,curQ);// Create Hadron for Residual
            FillHadronVector(curHadr1);                // Fill "new curHadr1" (delete equivalent)
            G4QHadron* curHadr2 = new G4QHadron(sPDG,s4Mom); // Creation Hadron for a Candidate
            FillHadronVector(curHadr2);                // Fill "new curHadr2" (delete equivalent)
            KillQuasmon();                             // This Quasmon is done
            if(addPhoton)
            {
              q4Mom+=phot4M;
              addPhoton=0.;                              // the Photon is adopted: clean it up
              momPhoton=0.;
              phot4M=G4LorentzVector(0.,0.,0.,0.);
	        }
            qEnv=theEnvironment;                       // Return the QEnvironment
            return theQHadrons;                        // The last decay of the quasmon...
	      }
	    }
      }
      curQ = memQ;                                     // Recover original curQ=valQ-candQ
      fdul = (rMass<sB&&rFl&&rPDG!=10);
    }
    if(j)                                              // j>0 except for "nothing is found"
    {
      if(fprob)rPDG=curQ.GetSPDGCode();                // PDG of residual Quasmon as S-hadron
      G4double reMass=sqrt(minSqT);
      if (!rPDG)
      {
        cerr<<"***G4Q::HQ: rQ="<<curQ<<",rPDG="<<rPDG<<"+sPDG="<<sPDG<<"(sM="<<sMass<<")"<<endl;
        G4Exception("G4Quasmon::HadronizeQuasmon: unidentifiable residual Hadron");
      }
      if(rPDG==221)reMass=mPi0;
      G4double aMass=mPi0;
      //if(envPDG>80000000&&(sPDG<80000000||envPDG!=pPDG))aMass=0.;
      if((sPDG<80000000&&envPDG>80000000&&envPDG!=90000000||sPDG>80000000&&sPDG!=90000000&&envPDG>pPDG)&&iniBarN>0
        ||iniBarN>1||rPDG==10) aMass=0.;               // No Pi0 condition if NucEnviron exists
#ifdef pdebug
      G4cout<<"G4Quasmon::HadrQuasmon: Is [hsf="<<hsflag<<"] or [rdf="<<fdul<<"] or [rM="<<rMass
            <<" < "<<reMass<<" + "<<aMass<<"] to decay in sH + minResQ?"<<G4endl;
#endif
      if(hsflag||rMass<reMass+aMass||fdul) //=====> Decay Q->S+H or Q/C->H1+H2 or suck or evapor
	  {
#ifdef pdebug
        cout<<"G4Quasmon::HadronizeQuasmon: Yes, hsf="<<hsflag<<",sPDG="<<sPDG<<",pM="<<pMass
            <<",Env="<<envPDG<<",QM="<<quasM<<valQ<<", fpr="<<fprob<<endl;
#endif
        G4LorentzVector q4tmp=q4Mom+G4LorentzVector(0.,0.,0.,pMass); // LorVectorQ+parentMass
        G4QContent tmpQC =valQ+pQC;                  // QuasmonQC+parentQC
        G4int      tmpPDG=tmpQC.GetSPDGCode();       // PDGCode of Compound(Quasmon+Parent)
        G4QPDGCode tmpQPDG(tmpPDG);                  // QPDG of CompoundQuasmon
        G4QContent tmpQ=envQC-pQC;                   // Residual for resid Environ !!
        G4QNucleus tmpN(tmpQ);                       // Pseudo nucleus for resid Environ
        G4int      tmpNPDG=tmpN.GetPDG();            // PDG Code of residual Environment
        G4double   tmpNM=tmpN.GetMZNS();             // Mass of residual environment !!
        G4QContent tmpRQ=valQ+transQC;               // QContent of Residual Quasmon !!
        G4QNucleus tmpR(tmpRQ);                      // Nucleus for Residual Quasmon
        G4int      tmpRPDG=tmpR.GetPDG();            // PDG for Residual Quasmon
        G4double   tmpRM=tmpR.GetMZNS();             // Mass of Residual Quasmon   !!
        G4QContent tmpTQ=tmpQ+tmpRQ;                 // Residual Nucleus for a Fragment
        G4QNucleus tmpT(tmpTQ);                      // Nucleus for Residual for a Fragment
        G4double   tmpTM=tmpT.GetMZNS();             // Mass of Residual Nucleus for Fragment !!
        G4int      tmpTPDG=tmpT.GetPDG();            // PDGC of Residual Nucleus for Fragment !!
        G4double   tmpQM=totM-tmpNM;                 // Bound Mass of newQuasmon
        G4double   tmpM =tmpQPDG.GetMass();          // Minimum free mass for newQ
        if(tmpQM<tmpM) tmpM=tmpQM;
        if(totBN>1&&totMass>totM&&totS>=0&&envPDG>80000000&&envPDG!=90000000)
		{
#ifdef ppdebug
          G4cout<<"G4Quasm::HadrQuasm:NEEDS-EVAPORATION-4: Q4M="<<q4Mom<<", QEnv="<<theEnvironment<<G4endl;
#endif
          if(addPhoton)
          {
            q4Mom+=phot4M;
            addPhoton=0.;                              // the Photon is adopted: clean it up
            momPhoton=0.;
            phot4M=G4LorentzVector(0.,0.,0.,0.);
	      }
          qEnv=theEnvironment;
          return theQHadrons;
		}
        G4double dm=quasM-sMass;
#ifdef pdebug
        G4cout<<"G4Quasm::HQ: f="<<fprob<<",d="<<dm<<",r="<<rPDG<<G4endl;
#endif
        // First decay suppression for the final decay in 2 particles
        if(j==1&&((rPDG==111||rPDG==221||rPDG==331) &&
                  (sPDG==221||sPDG==331||sPDG==223||sPDG==333||sPDG==225||sPDG==335) ||
                  (rPDG==221||rPDG==331) && (sPDG==111||sPDG==113||sPDG==115)
                 )&& G4UniformRand()>EtaEtaprime && dm>782.) rPDG=223;
        // Final state eta/eta' suppression
        if(rPDG==331&&dm<958.)rPDG=221;
        if(rPDG==221&&dm<548.)rPDG=111;
        if(rPDG<80000000)                            // => ResidualQuasmon is not NuclearCluster
        {
          reMass=GetRandomMass(rPDG,dm);             // Randomize mass of the Hadron
          if(reMass==0.)
	      {
            if(sPDG==221 || sPDG==331)               // Change eta-Candidate to pi-Candidate
	        {
              if     (sPDG==221) dm+=mEta-mPi0;
              else if(sPDG==331) dm+=mEtaP-mPi0;
              sPDG=111;
              sMass=mPi0;
              reMass=GetRandomMass(rPDG,dm);       // Rerandomize mass of resonance
              if(reMass==0.) cerr<<"***G4Quasmon::HQ: (2) rPDG="<<rPDG<<", dm="<<dm<<endl;
		    }
            else if(rPDG=111)
			{
              rPDG=22;
              reMass=0.;
		    }
          }
		}
        //else if(rPDG!=90000000) reMass=G4QPDGCode(rPDG).GetMass(); // Get Nuclear Cluster mass
        else reMass=G4QPDGCode(rPDG).GetMass(); // Get Nuclear Cluster mass
        //if(reMass+sMass>quasM)                   // Cann't decay Quasmon (mass is too small)
		//if(totBN>1&&totMass>totM&&totS>=0)
		if((totBN>1&&totMass>totM&&totS>=0||reMass+sMass>quasM)&&envPDG>80000000&&envPDG!=90000000)
        {
#ifdef ppdebug
          G4cout<<"G4Quasm::HadrQuasm:NEEDS-EVAPORATION-5: Q4M="<<q4Mom<<", QEnv="<<theEnvironment<<G4endl;
#endif
          if(addPhoton)
          {
            q4Mom+=phot4M;
            addPhoton=0.;                          // the Photon is adopted: clean it up
            momPhoton=0.;
            phot4M=G4LorentzVector(0.,0.,0.,0.);
	      }
          qEnv=theEnvironment;
          return theQHadrons;
		}
        G4LorentzVector r4Mom(0.,0.,0.,reMass);
        G4LorentzVector s4Mom(0.,0.,0.,sMass);     // Mass is random since probab. time
        if(!G4QHadron(q4Mom).DecayIn2(r4Mom, s4Mom))
        {
          G4cerr<<"***G4Quasmon::HadrQuasmon: M="<<q4Mom.m()<<" => rPDG="<<rPDG<<"(rM="
                <<reMass<<") + sPDG="<<sPDG<<"(sM="<<sMass<<")"<<G4endl;
    	  G4Exception("***G4Quasmon::HadrQuasmon:Hadron+SHadron DecayIn2 didn't succeed");
        }
#ifdef ppdebug
	    G4cout<<"G4Quasmon::HadrQuasmonon:=== 2.3 ===>HdVect s4M="<<s4Mom<<",sPDG="<<sPDG
              <<", r4M="<<r4Mom<<",rPDG="<<rPDG<<G4endl;
#endif
        if(rPDG==90000000)
        {
          G4cerr<<"***G4Quasmon::HadrQuasmon: MV="<<reMass<<G4endl;
          G4Exception("***G4Quasmon::HadrQuasm: Particle-Vacuum");
        }
        G4QHadron* curHadr1 = new G4QHadron(rPDG,r4Mom);// Create RealHadron for ResidEnv
        FillHadronVector(curHadr1);                // Fill "new curHadr1" (delete equivalent)
        G4QHadron* curHadr2 = new G4QHadron(sPDG,s4Mom);// Creation Hadron for theCandidate
        FillHadronVector(curHadr2);                // Fill "new curHadr2" (delete equivalent)
        KillQuasmon();                             // This Quasmon is done
        if(addPhoton)
        {
          q4Mom+=phot4M;
          addPhoton=0.;                            // the Photon is adopted: clean it up
          momPhoton=0.;
          phot4M=G4LorentzVector(0.,0.,0.,0.);
	    }
        qEnv=theEnvironment;                       // Return the QEnvironment
        return theQHadrons;                        // This is the last decay of the quasmon
	  }
#ifdef pdebug
	  G4cout<<"G4Quasmon::HQ:MQ="<<quasM<<" -> resQM="<<reMass<<" + sPDG="<<sPDG
            <<",M="<<sMass<<",frM="<<fr4Mom.m()<<",fskip="<<fskip<<G4endl;
#endif
	}
	else if(retC<maxRet) 
	{
      fskip=true;
      retC++;
	}
	if(!fskip)                                       // Continue search for fragmentation
	{
	  G4int ePDG=theEnvironment.GetPDG();
#ifdef pdebug
	  G4cout<<"G4Quasmon::HadrQ: sPDG="<<sPDG<<", fr4m="<<fr4Mom<<", ePDG="<<ePDG<<G4endl;
#endif
      //if(sPDG<80000000&&ePDG==NUCPDG)// ====>>>> Vacuum case @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(sPDG<80000000)                                // ====>>>> "Hadron candidate" case
      {
        G4int SQ=totQC.GetStrangeness();
#ifdef pdebug
	    G4cout<<"G4Quasmon::HadrQ: sPDG="<<sPDG<<", sM="<<sMass<<", SQ="<<SQ<<G4endl;
#endif
        if(!sPDG&&SQ<0)           // decay in K+/aK0 & residual
		{
          sPDG=321;
          sMass=mK;
          G4QContent resKPQC=totQC-G4QContent(0,1,0,0,0,1); // Residual Quark Content for K+
          G4QNucleus rKPN(resKPQC);                  // Pseudo nucleus for the Resid System
          G4double rKPM = rKPN.GetMZNS();            // min mass of the Residual System
          G4int  rKPPDG = rKPN.GetPDG();             // PDG of Residual
          G4QContent resK0QC=totQC-G4QContent(1,0,0,0,0,1); // Residual Quark Content for K0
          G4QNucleus rK0N(resK0QC);                  // Pseudo nucleus for the Resid System
          G4int  rK0PDG = rK0N.GetPDG();             // PDG of Residual
          G4double rK0M = rK0N.GetMZNS();            // min mass of the Residual System
          if(rKPM+mK > totMass && rK0M+mK0 > totMass || rKPPDG==NUCPDG || rK0PDG==NUCPDG)
		  {
#ifdef ppdebug
	        G4cout<<"G4Quasmon::HadrQuasm: ***PANIC*2* tM="<<totMass<<" < KM="<<mK<<","<<mK0
            <<", rM="<<rKPM<<","<<rK0M<<", d="<<mK+rKPM-totMass<<","<<mK0+rK0M-totMass<<G4endl;
#endif
            status =-1;                              // Panic exit
            if(addPhoton)
            {
              q4Mom+=phot4M;
              addPhoton=0.;                          // the Photon is adopted: clean it up
              momPhoton=0.;
              phot4M=G4LorentzVector(0.,0.,0.,0.);
    	    }
            qEnv=theEnvironment;                     // Return the QEnvironment
            return theQHadrons;
		  }
          if(rKPM + mK > rK0M + mK0)
		  {
            rPDG  = rK0PDG;                          // PDG of the Residual System to K0
            rMass = rK0M;
            sPDG  = 311;
            sMass = mK0;
		  }
          else
		  {
            rPDG  = rKPPDG;                          // PDG of the Residual System to K+
            rMass = rKPM;
            sPDG  = 321;
            sMass = mK;
		  }
          G4LorentzVector r4Mom(0.,0.,0.,rMass);
          G4LorentzVector s4Mom(0.,0.,0.,sMass);     // Mass is random since probab. time
          if(!G4QHadron(tot4M).DecayIn2(r4Mom, s4Mom))
          {
            G4cerr<<"***G4Quasmon::HadrQuasmon: tM="<<tot4M.m()<<totQC<<" => rPDG="<<rPDG
                  <<"(rM="<<rMass<<") + sPDG="<<sPDG<<"(sM="<<sMass<<")"<<G4endl;
    	    G4Exception("***G4Quasmon::HadrQuasmon:K+ResidNucl DecayIn2 didn't succeed");
          }
#ifdef ppdebug
	      G4cout<<"G4Quasmon::HadrQuasmonon:=== 2.4 ===>HdVect s4M="<<s4Mom<<",sPDG="<<sPDG
              <<", r4M="<<r4Mom<<",rPDG="<<rPDG<<G4endl;
#endif
          G4QHadron* curHadr1 = new G4QHadron(rPDG,r4Mom);// Create RealHadron for ResidEnv
          FillHadronVector(curHadr1);                // Fill "new curHadr1" (delete equivalent)
          G4QHadron* curHadr2 = new G4QHadron(sPDG,s4Mom);// Creation Hadron for theCandidate
          FillHadronVector(curHadr2);                // Fill "new curHadr2" (delete equivalent)
          KillQuasmon();                             // This Quasmon is done
          if(addPhoton)
          {
            q4Mom+=phot4M;
            addPhoton=0.;                              // the Photon is adopted: clean it up
            momPhoton=0.;
            phot4M=G4LorentzVector(0.,0.,0.,0.);
	      }
          qEnv=theEnvironment;                       // Return the QEnvironment
          return theQHadrons;                        // This is the last decay of the quasmon
		}
        G4bool ffin=false;                           // Flag of final decay in gamma-residual
        if(quasM<rMass+sMass&&(sPDG==221||sPDG==331))// Change eta-Candidate or ? to pi
	    {
          sPDG = 111;
          sMass=mPi0;
	    }
        else if(!sPDG)
		{
          if     (iniS<0&&iniQChg+iniQChg>=iniBarN)  // Try to decay in K+
		  {
            sPDG = 321;
            sMass= mK;
            G4QNucleus totN(valQ+KpQC);              // Nucleus Residual after Kp sub.
            rPDG = totN.GetPDG();
		    rMass= totN.GetMZNS();
		  }
          else if(iniS<0)                            // Try to decay in K0
		  {
            sPDG = 311;
            sMass= mK0;
            G4QNucleus totN(valQ+K0QC);              // Nucleus Residual after K0 sub.
            rPDG = totN.GetPDG();
		    rMass= totN.GetMZNS();
		  }
          else if(iniQChg>iniBarN-iniS)              // Try to decay in Pi+
		  {
            sPDG = 211;
            sMass= mPi;
            G4QNucleus totN(valQ-PiQC);              // Nucleus Residual after Pi+ sub.
            rPDG = totN.GetPDG();
		    rMass= totN.GetMZNS();
		  }
          else if(iniQChg<0)                         // Try to decay in Pi-
		  {
            sPDG = -211;
            sMass= mPi;
            G4QNucleus totN(valQ+PiQC);              // Nucleus Residual after Pi- sub.
            rPDG = totN.GetPDG();
		    rMass= totN.GetMZNS();
		  }
		  else
		  {
            sPDG = 22;
            sMass= 0.;
            rPDG = iniPDG;
		    rMass= iniQM;
		  }
          ffin = true;
		}
#ifdef pdebug
	    G4cout<<"G4Quasmon::HadrQ: MQ="<<q4Mom.m()<<"-> sPDG="<<sPDG<<", sM="<<sMass
              <<" + rPDG="<<rPDG<<",rM="<<rMass<<G4endl;
#endif
        if(q4Mom.m()<rMass+sMass)
		{
#ifdef ppdebug
	      G4cout<<"G4Quasmon::HadrQuasmonon: ***PANIC*1* tM="<<q4Mom.m()<<" < rM="<<rMass
              <<", sM="<<sMass<<", d="<<rMass+sMass-q4Mom.m()<<G4endl;
#endif
          status =-1;                                // Panic exit
          if(addPhoton)
          {
            q4Mom+=phot4M;
            addPhoton=0.;                              // the Photon is adopted: clean it up
            momPhoton=0.;
            phot4M=G4LorentzVector(0.,0.,0.,0.);
    	  }
          qEnv=theEnvironment;                       // Return the QEnvironment
          return theQHadrons;
		}
        G4LorentzVector resQ4Mom(0.,0.,0.,rMass);    // 4-mom of residual Quasmon in CMS
        G4LorentzVector s4Mom(0.,0.,0.,sMass);       // Mass has been rand. at probab. level
        if(!G4QHadron(q4Mom).DecayIn2(resQ4Mom, s4Mom))
        {
          G4cerr<<"***G4Quasmon::HadronizeQuasmon: MQ="<<q4Mom.m()<<", rPDG="<<rPDG<<", (M="
                <<rMass<<") + sPDG="<<sPDG<<"(M="<<sMass<<")"<<G4endl;
	      G4Exception("G4Quasmon::HadronizeQuasmon: Quasmon+Hadron DecayIn2 did not succeed");
        }
#ifdef pdebug
	    G4cout<<"G4Quasm::HadrQ:Decay Q="<<q4Mom<<"->s="<<sPDG<<s4Mom<<"+RQ="<<resQ4Mom<<G4endl;
#endif
        G4QHadron* candHadr = new G4QHadron(sPDG,s4Mom);// Creation a Hadron for the Candidate
        if(ffin)
		{
          theQHadrons.insert(candHadr);              // Fill the emergency photon (delete equivalent)
          G4QHadron* candHRes = new G4QHadron(rPDG,resQ4Mom);// Creation Hadron for the QResid
          FillHadronVector(candHRes);                // Fill "new candHRes" (delete equivalent)
#ifdef ppdebug
          G4cout<<"G4Quasmon::HadrQ: QEnv="<<theEnvironment<<G4endl;
#endif
          KillQuasmon();                             // This Quasmon is done
          if(addPhoton)
          {
            q4Mom+=phot4M;
            addPhoton=0.;                              // the Photon is adopted: clean it up
            momPhoton=0.;
            phot4M=G4LorentzVector(0.,0.,0.,0.);
	      }
          qEnv=theEnvironment;                       // @@ can be directly vacuum
          return theQHadrons;
		}
        else
		{
          FillHadronVector(candHadr);                // Fill "new candHadr"
          check+= s4Mom;                             // @@ Just for checking
          ccheck+=G4QPDGCode(sPDG).GetCharge();      // @@ Just for checking
          q4Mom = resQ4Mom;                          // Update ResidualQuasmon LorentzVector
          valQ    = curQ;                            // Update the Quark Content of Quasmon
		}
      }
      else                                           // ==> "Nuclear media || NuclCand" case
	  {
        G4QHadron* candHadr = new G4QHadron(sPDG,fr4Mom);// Createa a Hadron for the Candidate
        FillHadronVector(candHadr);                  // Fill the RadiatedHadron (delete equivalent)
        G4LorentzVector sumL=theEnvironment.Get4Momentum()+q4Mom; //@@ For Check Print Only
        check += fr4Mom;                             //@@ Just for checking
        ccheck+=G4QPDGCode(sPDG).GetCharge();        //@@ Just for checking
        q4Mom  = rQ4Mom;                             // Update the 4Mom of Quasmon
        valQ  += transQC;                            // Update the Quark Content of Quasmon
#ifdef pdebug
        cout<<"G4Quasmon::HadronizeQuasm:QuarkExchange QQC="<<valQ<<",tranQC="<<transQC<<endl;
#endif
        theEnvironment.Reduce(pPDG);                 // Update NuclearEnviron after Q->RQ+F
        sumL-=theEnvironment.Get4Momentum()+q4Mom+fr4Mom;
#ifdef pdebug
        cout<<"G4Quasmon::HadronizeQuasmon: >>>>>>>>>NM SUBCHECK>>>>>>>>>>>>>"<<sumL<<endl;
#endif
      }
	}
    G4LorentzVector sumLor=theEnvironment.Get4Momentum()+q4Mom+check;
    G4int eZ   = theEnvironment.GetZ();
    G4int sumC = eZ+valQ.GetCharge()+ccheck;
    G4int curPDG=valQ.GetSPDGCode();
#ifdef debug
    cout<<"G4Quasm::HadrQuasm:j="<<j<<",eZ="<<eZ<<valQ<<">>CHEK>>"<<sumLor<<",C="<<sumC<<endl;
    if(!curPDG)cerr<<"***G4Q::HQ: Quasmon-Tripolino j="<<j<<",QC="<<valQ<<endl;
    cout<<"G4Quasmon::HadronizeQuasmon: ResidualQuasmon LorV="<<q4Mom<<", QCont="<<valQ<<endl;
#endif
  }
#ifdef pdebug
  G4cout<<"G4Quasmon::HadrQuasm:***OUT***MQ2="<<q4Mom.m2()<<",EPDG="<<theEnvironment.GetPDG()<<G4endl;
#endif
  if(addPhoton)
  {
    q4Mom+=phot4M;
    addPhoton=0.;                                    // the Photon is adopted: clean it up
    momPhoton=0.;
    phot4M=G4LorentzVector(0.,0.,0.,0.);
  }
  qEnv=theEnvironment;                               // Return the QEnvironment
  return theQHadrons;
}

// Test for Decay, Decay, Filling of hadron vector
void G4Quasmon::FillHadronVector(G4QHadron* qH)
{//  ==========================================
  // SHORT-RANGE CORRELATIONS
  // ========================
  // dN/kdk=C/(A+k^2)^2, where A=-1/2ar. {NN=pp,nn,np}
  // Randomization: p_NN^2=A_NN/(1/rndm-1), E_NN=sqrt(m_NN^2+p_NN^2)
  // A_pp=888.3
  // A_nn=410.1
  // A_pn=297.5
  // Breite-Wigner representation:
  // m_pp=-0.9450 MeV; G_pp=0
  // m_nn=-0.4336 MeV; G_nn=0
  // m_pn=-0.3139 MeV; G_pn=0
  // ----------------------------------------------------------------
  static const G4double mAlpha = G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QContent PiQC(0,1,0,1,0,0);
  static const G4QContent K0QC(1,0,0,0,0,1);
  static const G4QContent KpQC(0,1,0,0,0,1);
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();

  status=1;                                    // Something is going to be filled by Quasmon
  G4int thePDG      = qH->GetPDGCode();        // Get PDG code of the Hadron to switch
  G4LorentzVector t = qH->Get4Momentum();      // 4-Mom of Chipolino
#ifdef ppdebug
  G4cout<<"G4Quasmon::FillHadronVector:Hadron's PDG="<<thePDG<<",4Mom="<<qH->Get4Momentum()<<G4endl;
#endif
  // In the decay scheme only one resonance can be present & it should be the last
  if (thePDG==90000000||thePDG==90999999||thePDG==90999000||thePDG==90000999||thePDG==89999001)
  {
    G4cerr<<"***G4Quasmon::FillHadronVector:PDG="<<thePDG<<",M="<<qH->Get4Momentum().m()<<G4endl;
  }
  if (thePDG==10) // Chipolino decays (@@always - Chipolino is not kept in HadV (*Example*))
  {
    G4double rM = t.m();                       // Mass of Chipolino
    G4QContent chipQC = qH->GetQC();           // QC of Chipolino
    G4QContent h1QC = chipQC.SplitChipo(rM);   // Extract QC of one of the Hadrons of Chipolino
    G4QContent h2QC = chipQC - h1QC;           // Define QC of the second Hadron
    G4int h1PDG = h1QC.GetSPDGCode();          // PDGCode of the First Hadron
    G4int h2PDG = h2QC.GetSPDGCode();          // PDGCode of the First Hadron
    if(!h1PDG || !h2PDG)
	{
      cerr<<"***FillHV:h1QC="<<h1QC<<"(PDG="<<h1PDG<<"),h2QC="<<h2QC<<"(PDG="<<h2PDG<<")"<<endl;
	  G4Exception("G4Quasmon::FillHadronVector: Cipolino cann't be defragmented");
	}
    G4QHadron* fHadr = new G4QHadron(h1PDG);   // the First Hadron is created
    G4QHadron* sHadr = new G4QHadron(h2PDG);   // the Second Hadron is created
    G4LorentzVector f4Mom = fHadr->Get4Momentum();
    G4LorentzVector s4Mom = sHadr->Get4Momentum();
    if(!qH->DecayIn2(f4Mom,s4Mom))
    {
      delete fHadr;                            // Delete "new fHadr"
      delete sHadr;                            // Delete "new sHadr"
      cerr<<"***G4Q::FillHadronVector:ChipQC"<<chipQC<<":PDG1="<<h1PDG<<",PDG2="<<h2PDG<<endl;
	}
    else
    {
      qH->SetNFragments(2);
      theQHadrons.insert(qH);                  // (delete equivalent)
      fHadr->Set4Momentum(f4Mom);              // Put the randomized 4Mom to 1-st Hadron
      FillHadronVector(fHadr);                 // Fill 1st Hadron (delete equivalent)
      sHadr->Set4Momentum(s4Mom);              // Put the randomized 4Mom to 2-nd Hadron
      FillHadronVector(sHadr);                 // Fill 2nd Hadron (delete equivalent)
    }
  }
  else if(thePDG>80000000&&thePDG!=90000000)   // === Decay-Evaporation ===
  {
    G4double totMass=qH->GetMass();            // Real Mass of the nuclear fragment
    G4QNucleus qNuc(t,thePDG);                 // Make a Nucleus out of the Hadron
    G4double GSMass =qNuc.GetGSMass();         // GrState Mass of the nuclear fragment
    G4QContent totQC=qNuc.GetQCZNS();          // Total Quark Content of Residual Nucleus
    G4int    nN     =qNuc.GetN();              // A#of neutrons in the Nucleus
    G4int    nZ     =qNuc.GetZ();              // A#of protons in the Nucleus
    G4int    nS     =qNuc.GetS();              // A#of protons in the Nucleus
    G4int    bA     =qNuc.GetA();              // A#of baryons in the Nucleus
#ifdef pdebug
	G4cout<<"G4Quasm::FillHadrVect: Nucl="<<qNuc<<",nPDG="<<thePDG<<",GSM="<<GSMass<<G4endl;
#endif
    if((nN<0||nZ<0||nS<0)&&bA>0)               // => "Anti-strangeness or ISOBAR" case
	{
      G4double m1=mPi;         
      G4int  PDG1=-211;
      G4QNucleus  newNpm(totQC+PiQC);
      G4int  PDG2=newNpm.GetPDG();
      G4double m2=newNpm.GetMZNS();
      if(nN<0)
	  {
        PDG1       =211;
        G4QNucleus  newNpp(totQC-PiQC);
        PDG2       =newNpp.GetPDG();
        m2         =newNpp.GetMZNS();
	  }
      if(nS<0)
      {
        m1         =mK;         
        PDG1       =321;
        G4QNucleus  newNp(totQC-KpQC);
        PDG2       =newNp.GetPDG();
        m2         =newNp.GetMZNS();
        G4QNucleus  newN0(totQC-K0QC);
        G4double m3=newN0.GetMZNS();
        if (m3+mK0<m2+mK)                      // => "aK0+ResA is better" case
        {
          m1  =mK0;
          PDG1=311;
          m2  =m3;
          PDG2=newN0.GetPDG();
        }
	  }
      if(totMass>m1+m2)                        // => "can decay" case
      {
        G4LorentzVector fq4M(0.,0.,0.,m1);
        G4LorentzVector qe4M(0.,0.,0.,m2);
        if(!qH->DecayIn2(fq4M,qe4M))
	    {
          G4cerr<<"***G4Quasm::FillHadrV: QM="<<t.m()<<"-> Mes="<<PDG1<<"(M="
	  	        <<m1<<") + ResA="<<PDG2<<"(M="<<m2<<")"<<G4endl;
	      G4Exception("G4Quasm::FillHadrV: Mes+ResA DecayIn2 did not succeed");
	    }
        qH->SetNFragments(2);                  // Put a#of Fragments=2
        theQHadrons.insert(qH);
        G4QHadron* H1 = new G4QHadron(PDG1,fq4M);
        theQHadrons.insert(H1);                // (delete equivalent)
        G4QHadron* H2 = new G4QHadron(PDG2,qe4M);
        FillHadronVector(H2);                  // (delete equivalent)
	  }
      else if(m1+m2-totMass<0.01)              // Split the 4-momentum
	  {
        G4double r1=m1/totMass;
        G4double r2=1.-r1;
        qH->SetNFragments(2);                  // Put a#of Fragments=2
        theQHadrons.insert(qH);
        G4QHadron* H1 = new G4QHadron(PDG1,r1*t);
        theQHadrons.insert(H1);                // (delete equivalent)
        G4QHadron* H2 = new G4QHadron(PDG2,r2*t);
        FillHadronVector(H2);                  // (delete equivalent)
	  }
      else
	  {
        G4cerr<<"***G4Quasm::FillHadrVec: PDG="<<thePDG<<"("<<t.m()<<") < Mes="<<PDG1
              <<"("<<m1<<") + ResA="<<PDG2<<"("<<m2<<"), d="<<totMass-m1-m2<<G4endl;
	    G4Exception("G4Quasm::FillHadrVec: mass of decaying hadron is too small");
	  }
	}
    else if(abs(totMass-GSMass)<.1)            // the Nucleus is too close the Ground State
    {
#ifdef pdebug
	  G4cout<<"G4Quasm::FillHadrVect: Ground state decay"<<G4endl;
#endif
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
#ifdef pdebug
      G4cout<<"G4Quasm::FillHadrVect: rP="<<pResPDG<<",rN="<<nResPDG<<",rL="<<lResPDG<<",nN="
            <<nN<<",nZ="<<nZ<<",nL="<<nS<<",totM="<<totMass<<",n="<<totMass-nResM-mNeut<<",p="
            <<totMass-pResM-mProt<<",l="<<totMass-lResM-mLamb<<G4endl;
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
          barPDG = 90000001;
          resPDG = nResPDG;
          barM= mNeut;
          resM= nResM;
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
          G4cerr<<"***G4Q::FillHadrVect:PDG="<<thePDG<<",M="<<totMass<<"< GSM="<<GSMass<<G4endl;
          G4Exception("***G4Quasmon::FillHadronVector: Below GSM but can not decay in p or n");
		}
        G4LorentzVector a4Mom(0.,0.,0.,barM);
        G4LorentzVector b4Mom(0.,0.,0.,resM);
        if(!qH->DecayIn2(a4Mom,b4Mom))
        {
          theQHadrons.insert(qH);              // No decay (delete equivalent)
          G4cerr<<"***G4Q::FillHadronVector: Be8 decay did not succeeded"<<G4endl;
	    }
        else
        {
          qH->SetNFragments(2);                // Fill a#of fragments to decaying Hadron
          theQHadrons.insert(qH);              // Fill hadron with nf=2 (delete equivalent)
          G4QHadron* HadrB = new G4QHadron(barPDG,a4Mom);
          FillHadronVector(HadrB);             // Fill 1st Hadron (delete equivalent)
          G4QHadron* HadrR = new G4QHadron(resPDG,b4Mom);
          FillHadronVector(HadrR);             // Fill 2nd Hadron (delete equivalent)
        }
	  }
      else
      {
#ifdef pdebug
	    G4cout<<"G4Quasm::FillHadrVect: Decay as it is"<<G4endl;
#endif
        theQHadrons.insert(qH);             // No decay  (delete equivalent)
	  }
    }
    else if (totMass<GSMass)
	{
      cerr<<"***G4Quasm::FillHV: tM="<<totMass<<" > GSM="<<GSMass<<", d="<<totMass-GSMass<<endl;
      G4Exception("***G4Quasmon::FillHadronVector: Out Mass is below the Ground State value");
	}
    else                                       // ===> Evaporation of excited system
	{
#ifdef pdebug
      G4cout<<"G4Quasm::FillHadrVect:Evaporate "<<thePDG<<",tM="<<totMass<<" > GS="<<GSMass
            <<qNuc.Get4Momentum()<<", m="<<qNuc.Get4Momentum().m()<<G4endl;
#endif
      G4QHadron* bHadron = new G4QHadron;
      G4QHadron* rHadron = new G4QHadron;
      if(!qNuc.EvaporateBaryon(bHadron,rHadron))
	  {
        cerr<<"***G4Quasmon::FillHadronVector:Evaporation, PDG="<<thePDG<<",tM="<<totMass<<endl;
        delete bHadron;
        delete rHadron;
        theQHadrons.insert(qH);                // Fill hadron in the HadronVector as it is
	  }
#ifdef pdebug
      cout<<"G4Quasm::FillHadrVec:Done b="<<bHadron->GetQPDG()<<",r="<<rHadron->GetQPDG()<<endl;
#endif
      qH->SetNFragments(2);                    // Fill a#of fragments to decaying Hadron
      theQHadrons.insert(qH);                  // Fill hadron with nf=2 (delete equivalent)
      FillHadronVector(bHadron);               // Fill Evapor. Baryon (delete equivalent)
      FillHadronVector(rHadron);               // Fill Residual Nucl. (delete equivalent)
	}
  }
  else if(!DecayOutHadron(qH))                 // (delete equivalent)
  {
#ifdef pdebug
    G4cout<<"***G4Quasmon::FillHadronVector: Emergency OUTPUT, PDGcode= "<<thePDG<<G4endl;
#endif
  }
}

// Calculate a momentum of quark-parton greater then minimum value kMin
G4double G4Quasmon::GetQPartonMomentum(G4double mR2, G4double mC2)
//       =========================================================
{
  //gives k>kMin QParton Momentum for the current Quasmon
#ifdef pdebug
  cout<<"G4Quasmon::GetQPartonMomentum: is called with mR="<<sqrt(mR2)<<",mC="<<sqrt(mC2)<<endl;
#endif
  G4double qMass = q4Mom.m();
  G4double kLim  = qMass/2.;               //Kinematikal limit for "k"
  G4double twM   = qMass+qMass;
  G4double kMax  = kLim;
  G4double kMin  = mC2/twM;
  if(mR2!=0.)
  {
    G4double qM2   = qMass*qMass;          // Squared Mass of Quasmon
    G4double frM   = twM+twM;              // 4*M_Quasmon
    G4double fM2m2 = frM*qMass*mC2;        // 4*MQ**2*mC2=0.
    G4double Mmum  = qM2+mC2-mR2;          // QM**2-mR2
    G4double Mmum2 = Mmum*Mmum;            // (QM**2-mR2)**2
    if(Mmum2<fM2m2)G4Exception("G4Quasmon::QPartMom: mR & mC are too big in comp. with mQ");
    G4double sqM   = sqrt(Mmum2-fM2m2);    // QM**2-mR2
    kMin  = (Mmum-sqM)/frM;                // kMin=0.
    kMax  = (Mmum+sqM)/frM;                // kMax=2*(QM**2-mR2)/4QM
  }
  if (kMin<0 || kMin>kLim || kMax<0 || kMax>kLim || qMass<=0. || nOfQ<3)
  {
    G4cerr<<"***G4Q::GetQPartMom: kMax="<<kMax<<", kMin="<<kMin<<", MQ="<<qMass
          <<", n="<<nOfQ<<G4endl;
	G4Exception("G4Quasmon::GetQPartonMomentum: Can not generate quark-parton");   
  }
#ifdef pdebug
  G4cout<<"G4Q::GetQPartonMom: kLim="<<kLim<<",kMin="<<kMin<<",kMax="<<kMax
        <<", nQ="<<nOfQ<<G4endl;
#endif
  G4int n=nOfQ-2;
  G4double fn=n;
  G4double vRndm = G4UniformRand();
  if (kMin>0.) // If there is a minimum cut for the QpMomentum 
  {
    G4double xMin=kMin/kLim;
    if (kMax>=kLim) vRndm = vRndm*pow((1.-xMin),n)*(1.+n*xMin);
    else
	{
      G4double xMax=kMax/kLim;
      G4double vRmin = pow((1.-xMin),n)*(1.+n*xMin);
      G4double vRmax = pow((1.-xMax),n)*(1.+n*xMax);
      vRndm = vRmax + vRndm*(vRmin-vRmax);
    }
  }
  else if (kMax<kLim)
  {
    G4double xMax=kMax/kLim;
    G4double vRmax = pow((1.-xMax),n)*(1.+n*xMax);
    vRndm = vRmax + vRndm*(1.-vRmax);
  }
  if (n==1) return kLim*sqrt(1.-vRndm); // Direct solution for nOfQ==3
  else                                  // Needs iterations
  {
    G4double x  = 1. - pow(vRndm*(1+n*vRndm)/(fn+1.),1/fn);
    G4double ox = x;
    G4int    it = 0;
    G4double d  = 1.;                   // Difference prototype
    G4double df = 1./static_cast<double>(nOfQ);
    G4double f  = df*(static_cast<int>(nOfQ*nOfQ*n*x/5.)+(nOfQ/2));
    if(f>99.)
    {
#ifdef pdebug
      cout<<"G4Quasmon::GetQPartonMomentum: f="<<f<<" is changer to 99"<<endl;
#endif
      f  = 99.;
    }
    G4double r  = 1.-x;
    G4double p  = r;
    if (n>2) p  = pow(r,n-1);
    G4double nx = n*x;
    G4double c  = p*r*(1.+nx);
    G4double od = c - vRndm;
#ifdef pdebug
	cout<<"G4Q::QPMom:>>>>>>First x="<<x<<", n="<<n<<", f="<<f<<", d/R(first)="<<od/vRndm<<endl;
#endif
    G4int nitMax=n+n+n;
    if(nitMax>100)nitMax=100;
    while ( abs(d/vRndm) > 0.001 && it <= nitMax)
	{
      x  = x + f*od/(r*nx*(fn+1.));
      r  = 1.-x;
      if (n>2) p  = pow(r,n-1);
      else     p  = r;
      nx = n*x;
      c  = p*r*(1.+nx);
      d  = c - vRndm;
      if ((od>0&&d<0)||(od<0&&d>0))
      {
        if (f>1.0001) f=1.+(f-1.)/2.;
        if (f<0.9999) f=1.+(1.-f)*2;
        x = (x + ox)/2.;
        r  = 1.-x;
        if (n>2) p  = pow(r,n-1);
        else     p  = r;
        nx = n*x;
        c  = p*r*(1.+nx);
        d  = c - vRndm;
      }
      else
	  {
        if (f>1.0001&&f<99.) f=1.+(f-1.)*2;
        if (f<0.99999999999) f=1.+(1.-f)/2.;
	  }
#ifdef pdebug
	  cout << "G4Q::QPMom: Iteration#" << it << ": (c=" << c << " - R=" << vRndm 
           << ")/R =" << d/vRndm << ", x=" << x << ", f=" << f << endl;
#endif
      od = d;
      ox = x;
      it++;
	};
    if(it>nitMax)
	{
#ifdef pdebug
      cout<<"...G4Quasmon::GetQPartonMomentum: a#of iterations > nitMax="<<nitMax<<endl;
#endif
    }
    if(x>1.) x=0.9;
    if(x<0.) x=.01;
    G4double kCand=kLim*x;
    if(kCand>=kMax)kCand=kMax-.001;
    if(kCand<=kMin)kCand=kMin+.001;
    return kCand;
  }
  return 0.;
} 

// For the given quasmon mass calculate a number of quark-partons in the system
void G4Quasmon::CalculateNumberOfQPartons(G4double qMass)
//   =====================================================
{
  static const G4double mK0  = G4QPDGCode(311).GetMass();
  // @@ Temporary here. To have 3 quarks in Nucleon Temperature should be < M_N/4 (234 MeV)
  // M^2=4*n*(n-1)*T^2 => n = M/2T + 1/2 + T/4M + o(T^3/16M^3)
  // @@ Genius (better than 10**(-3) even for n=2!) but useless
  //G4double qMOver2T = qMass/(Temperature+Temperature);
  //G4double est = qMOver2T+1.+0.125/qMOver2T;
  // @@ Longer but exact
  G4double qMOverT = qMass/Temperature;
  G4int valc = valQ.GetTot();
  // .................................
  // Exponent, Double Split, Poisson 1 ============
  //G4int b = valQ.GetBaryonNumber();
  //G4int mq= 3*b;
  //if (!b) mq=2;
  //G4double mean = ((1.+sqrt(1.+qMOverT*qMOverT))/2. - mq)/2.;
  //if(mean<0.) nOfQ=mq;
  // Exponent ------
  //else nOfQ=mq-2*mean*log(G4UniformRand());
  // Poisson 1 ------
  //else nOfQ=mq+2*RandomPoisson(mean);
  // Double Split ------
  //else
  //{
  //  G4int imean = static_cast<int>(mean);
  //  G4double dm = mean - imean;
  //  if(G4UniformRand()>dm) nOfQ=mq+imean+imean;
  //  else nOfQ=mq+imean+imean+2;
  //}
  // .........
  // Poisson 2 =============
  if(valc%2==0)nOfQ = 2*RandomPoisson((1.+sqrt(1.+qMOverT*qMOverT))/4.); // abs(b) is even
  else   nOfQ = 1+2*RandomPoisson((1.+sqrt(1.+qMOverT*qMOverT))/4.-0.5); // abs(b) is odd
  //
#ifdef pdebug
  cout<<"G4Quasm::Calc#ofQP:QM="<<qMass<<",T="<<Temperature<<",QC="<<valQ<<",nOfQ="<<nOfQ<<endl;
#endif
  G4int absb = abs(valQ.GetBaryonNumber());
  G4int tabn = 0;
  if(absb>0)tabn=3*absb;    // Minimal QC for fragmentation
  else if(tabn<4) tabn=4;   // Mesonic system
  if (nOfQ<tabn) nOfQ=tabn;
  G4int nSeaPairs = (nOfQ-valc)/2;
  G4int stran = abs(valQ.GetS());
  G4int astra = abs(valQ.GetAS());
  if(astra>stran) stran=astra;
  G4int nMaxStrangeSea=static_cast<int>((qMass-stran*mK0)/(mK0+mK0));//KK is a minimum for s-sea
  if (absb) nMaxStrangeSea=static_cast<int>((qMass-absb)/672.); //LambdaK is a minimum for s-sea
#ifdef pdebug
  G4cout<<"G4Q::Calc#ofQPart:"<<valQ<<",INtot="<<valc<<",FinTot="<<nOfQ<<",Sea="<<nSeaPairs<<G4endl;
#endif
  if (nSeaPairs)            // Add sea to initial quark content
  {
    G4int morDec=0;
    if(nSeaPairs>0)valQ.IncQAQ(nSeaPairs,SSin2Gluons);
    else    morDec=valQ.DecQAQ(nSeaPairs);
    if(morDec)
	{
#ifdef pdebug
      G4cout<<"***G4Q::Calc#ofQPart: not totally reduced #pairs="<<morDec<<G4endl;
#endif
      nOfQ+=morDec+morDec;
	}
    G4int sSea=valQ.GetS(); // Content of strange quarks
    G4int asSea=valQ.GetAS();
    if(asSea<sSea) sSea=asSea;
    if(sSea>nMaxStrangeSea) // @@@@@@@ Too many strange sea ??????
	{
      sSea-=nMaxStrangeSea; // Strange sea excess
      valQ.DecS(sSea);      // Reduce strange sea to adoptable limit
      valQ.DecAS(sSea);
      valQ.IncQAQ(sSea,0.); // Add notstrange sea ????????
	}
  }
#ifdef pdebug
  G4cout<<"G4Quasmon::Calc#ofQPart: FinalQC="<<valQ<<G4endl;
#endif
} 

// Calculate a#of clusters in the nucleus
void G4Quasmon::PrepareClusters()
//   ============================
{
  G4double ze             = theEnvironment.GetZ();        // For bn<3 (in all nucleus)
  G4double ne             = theEnvironment.GetN();
  G4double se             = theEnvironment.GetS();
  G4double ae             = ze + ne + se;
  G4double ze1            = theEnvironment.GetDZ() + 1;   // For bn>2 (nly in a dense region)
  G4double ne1            = theEnvironment.GetDN() + 1;
  G4double se1            = theEnvironment.GetDS() + 1;
  G4double pos            = theEnvironment.GetProbability(); // Vacuum value
  for (G4int index=0; index<theQCandidates.entries(); index++)
  {
    G4QCandidate* curCand=theQCandidates[index];
    G4int cPDG  = curCand->GetPDGCode();
    if(cPDG>80000000&&cPDG!=90000000)              // ===> Cluster case
	{
      G4QNucleus cN(cPDG);
      G4int zc = cN.GetZ();                        // "Z" of the cluster
      G4int nc = cN.GetN();                        // "N" of the cluster
      G4int sc = cN.GetS();                        // "S" of the cluster
      G4int ac = cN.GetA();                        // "A" of the cluster
      pos      = theEnvironment.GetProbability(ac);// Get a cluster prob. normalization factor
      G4double dense=1.;
      if(ac==1)dense=theEnvironment.GetProbability(254)/pos;
      if(ac==2)dense=theEnvironment.GetProbability(255)/pos;
#ifdef sdebug
	  cout<<"G4Quasmon::PrepareClusters: cPDG="<<cPDG<<",norm="<<pos<<",zc="<<zc<<",nc="<<nc
          <<",sc="<<sc<<",ze1="<<ze1<<",ne1="<<ne1<<",se1="<<se1<<endl;
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
          else cerr<<"***G4Quasmon::PrepareClusters:z="<<zc<<",n="<<nc<<",s="<<sc<<endl;
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
	  cout<<"G4Quasmon::PrepareClusters:ClastPDG="<<cPDG<<",prePro="<<pos<<",d="<<dense<<endl;
#endif
	}
	else
    {
      curCand->SetPreProbability(pos);            // ===> Hadronic case in Vacuum     
      curCand->SetDenseProbability(0.);           // ===> Hadronic case in Vacuum
    }
	curCand->SetPossibility(true);                // All candidates are possible at this point
  }
}// End of "PrepareClusters"

// Calculate the hadronization pre-probabiliies for all candidates,
// by default all candidates are possible, calculate k_min(=BindEnergy)
void G4Quasmon::PrepareCandidates(G4int j)
//   =====================================
{
  G4double hadronicProbab = theEnvironment.GetProbability();
  G4double ze             = theEnvironment.GetZ();       // For bn<3
  G4double ne             = theEnvironment.GetN();
  G4double se             = theEnvironment.GetS();
  G4double ae             = ze + ne + se;
  G4double ze1            = theEnvironment.GetDZ() + 1;  // For bn>2
  G4double ne1            = theEnvironment.GetDN() + 1;
  G4double se1            = theEnvironment.GetDS() + 1;
  if(hadronicProbab==0.) hadronicProbab=1.;       // If ther is no environment (vacuum)
  G4double quasmM         = q4Mom.m();            // Mass of quasmon
  G4double dQuasmM        = quasmM+quasmM;
  G4int    maxB           = theEnvironment.GetMaxClust();// A of maxCluster of the Environment
  G4double totZ           = theEnvironment.GetZ() + valQ.GetCharge();       // Z of the Nucleus
  G4double totA           = theEnvironment.GetA() + valQ.GetBaryonNumber(); // A of the Nucleus
  G4QContent totQC        = theEnvironment.GetQCZNS()+valQ;
  G4QNucleus totN(totQC);                         // pseudo nucleus for the Total System
  G4double totM           = totN.GetMZNS();       // min mass Of the Total System
  G4int nCand=theQCandidates.entries();
#ifdef pdebug
  G4cout<<"G4Quasmon::PrepareCandidates:hP="<<hadronicProbab<<",A="<<ae<<",mB="<<maxB<<",nC="
	    <<nCand<<",p(1)="<<theEnvironment.GetProbability(1)<<",zE="<<ze<<",nE="<<ne<<G4endl;
#endif
  if(nCand) for (G4int index=0; index<nCand; index++)
  {
    G4QCandidate*  curCand=theQCandidates[index];
    //G4double realM=curCand->GetMass();          // Real mass (BoundMass for a Cluster)
    // Calculate bound mass taking into account Quasmon+Environ
    G4QContent canQC=curCand->GetQC();
    G4QContent resQC=totQC-canQC;
#ifdef sdebug
    G4cout<<"G4Quasmon::PrepareCandidates: #"<<index<<",QC="<<canQC<<", rQC="<<resQC<<G4endl;
#endif
    G4QNucleus resN(resQC);                       // pseudo nucleus for Resid System
#ifdef sdebug
    G4cout<<"G4Quasmon::PrepareCandidates: resN="<<resN<<G4endl;
#endif
    G4double resM   = resN.GetMZNS();             // min mass Of the Residual System
    G4double realM  = totM-resM;
#ifdef sdebug
    G4cout<<"G4Quasmon::PrepareCandidates: realM="<<realM<<G4endl;
#endif
    G4int cPDG      = curCand->GetPDGCode();
    G4double PDGM   = G4QPDGCode(cPDG).GetMass(); // Vacuum (PDG) mass for the cluster
#ifdef sdebug
    G4cout<<"G4Quasmon::PrepareCandidates: PDG="<<cPDG<<",PDGM="<<PDGM<<G4endl;
#endif
    if (realM>PDGM) realM=PDGM;
    if(cPDG>80000000&&cPDG!=90000000)             // Nuclear candidate case
	{
      G4QNucleus cN(cPDG);
      G4int zc = cN.GetZ();                       // "Z" of the cluster
      G4int nc = cN.GetN();                       // "N" of the cluster
      G4int sc = cN.GetS();                       // "S" of the cluster
      G4int ac = cN.GetA();                       // "A" of the cluster
      //if(ac>maxB||ze<zc||ne<nc||se<sc) curCand->SetPossibility(false);
      if(ac>maxB) curCand->SetPossibility(false);
      else
	  {
        G4double pos = theEnvironment.GetProbability(ac);
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
            else G4cerr<<"***G4Quasmon::PrepareCandidates:z="<<zc<<",n="<<nc<<",s="<<sc<<G4endl;
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
        curCand->SetMass(realM);                 // Keep real (bounded) mass for a Candidate
        G4double curKMin=(PDGM*PDGM-realM*realM)/(realM+realM);
        if(totZ>zc) curKMin+=CoulombBarrier(totZ,totA,zc,ac);
        curCand->SetKMin(curKMin);
	    curCand->SetPossibility(true);
	  }
	}
	else                                         // Hadronic case     
	{
      curCand->SetPreProbability(hadronicProbab);
      G4double cz=curCand->GetCharge();
#ifdef sdebug
      G4cout<<"G4Quasmon::PrepareCandidates: Chg"<<cz<<", dQM="<<dQuasmM<<G4endl;
#endif
      G4double curKMin=realM*realM/dQuasmM;
      if(totZ>cz)curKMin+=CoulombBarrier(totZ,totA,cz,1.);
      curCand->SetKMin(curKMin);
	  curCand->SetPossibility(true);
	}
  }
#ifdef pdebug
  G4cout<<"G4Quasmon::PrepareCandidates: >>>OUT<<<"<<G4endl;
#endif
}

// Modify Candidate masses in nuclear matter and set possibilities
void G4Quasmon::ModifyInMatterCandidates()
//   ======================================
{
  G4double envM = theEnvironment.GetMass();      // Mass of the Current Environment
  G4int envPDGC = theEnvironment.GetPDGCode();   // PDG Code of Current Environment
  G4QContent envQC=theEnvironment.GetQCZNS();    // QuarkContent of the current Nuclear Environ.
  G4int eP = theEnvironment.GetZ();              // A#of protons in the Current Environment
  G4int eN = theEnvironment.GetN();              // A#of neutrons in the Current Environment
  G4int eL = theEnvironment.GetS();              // A#of lambdas in the Current Environment
  for (G4int ind=0; ind<theQCandidates.entries(); ind++)
  {
    G4QCandidate* curCand=theQCandidates[ind];   // Pointer to the Candidate
    G4int  cPDG = curCand->GetPDGCode();         // PDGC of the Candidate
    if(cPDG>80000000&&cPDG!=90000000)            //Fragments are modified to take into account CurNuclEnv
	{
      G4QNucleus cNuc(cPDG);                     // Fake nucleus for the Candidate
      G4int cP = cNuc.GetZ();                    // A#of protons in the Current Environment
      G4int cN = cNuc.GetN();                    // A#of neutrons in the Current Environment
      G4int cL = cNuc.GetS();                    // A#of lambdas in the Current Environment
      G4QPDGCode cQPDG(cPDG);                    // QPDG of the Current Cluster
      if(eP>=cP&&eN>=cN&&eL>=cL)                 // Cluster exists
      {
        G4double                  clM = 0.;      // Prototype for the BoundCluster mass
        G4double renvM = 0;                      // Prototype of the residual Environ mass
        if(cP==eP&&cN==eN&&cL==eL)clM = cQPDG.GetMass();//The only(not bound) Cluster of Environ
        else                                     // Bound Cluster
		{
          renvM = cQPDG.GetNuclMass(eP-cP,eN-cN,eL-cL);
          clM   = envM-renvM;
        }
        curCand->SetParPossibility(true);
        curCand->SetMass(clM);
#ifdef sdebug
        cout<<"G4Quasmon:CalcHP:C="<<cPDG<<cNuc<<clM<<",Env="<<envPDGC<<",eM="<<renvM<<endl;
#endif
	  }
      else curCand->SetParPossibility(false);
	} // @@ Modification of hadron masses in nuclear matter are not implemented yet
  }
}

// Randomize the Resonance masses and calculate probabilities of hadronization for them
void G4Quasmon::CalculateHadronizationProbabilities(G4double kVal, G4double kLS, G4int j)
//   ====================================================================================
{
  G4int    vap = nOfQ-3;                         // Vacuum power
  G4double mQ2 = q4Mom.m2();                     // Squared Mass of the Quasmon
  G4double eQ  = q4Mom.e();                      // LS Energy of the Quasmon
  G4double mQ  = sqrt(mQ2);                      // Mass of the decaying Quasmon
  G4double dk  = kVal + kVal;                    // Double momentu of quark-parton in QCM
  G4double rQ2 = mQ2-dk*mQ;                      // Min Residual Colored Quasmon Squared Mass
  G4double rQ  = sqrt(rQ2);
  G4double mQk = mQ-dk;                          // For acceleration
  //G4double vaf = theEnvironment.GetProbability()*mQk/kVal/vap; //@@ !! Vacuum factor
  G4double vaf = mQk/kVal/vap;                   //@@ Vacuum factor (instead of in G4QNucleus)
  G4double accumulatedProbability = 0.;
  G4double secondAccumProbability = 0.;
  G4int    iniPDG =valQ.GetSPDGCode();
  G4double iniQM = G4QPDGCode(iniPDG).GetMass(); // Not bounded quasmon GS mass
  G4double freeE = (mQ-iniQM)*iniQM;
  G4int qBar =valQ.GetBaryonNumber();
  G4int absb = abs(qBar);                        // Abs BaryNum of Quasmon
  G4int nofU = valQ.GetU()- valQ.GetAU();        // A#of u-quarks
  G4int dofU = nofU+nofU;
  G4int nofD = valQ.GetD()- valQ.GetAD();        // A#of d-quarks
  G4int dofD = nofD+nofD;
  G4int qChg = valQ.GetCharge();
  G4int qIso = qBar-qChg-qChg;
  G4int maxC = theQCandidates.entries();         // A#of candidates
  G4int maxB = theEnvironment.GetMaxClust();     // Maximum BaryNum for clusters
  G4double totZ = theEnvironment.GetZ() + valQ.GetCharge();       // Z of the Nucleus
  G4double totA = theEnvironment.GetA() + valQ.GetBaryonNumber(); // A of the Nucleus
  G4double envM = theEnvironment.GetMass();      // Mass of the Current Environment
  G4int envPDGC = theEnvironment.GetPDGCode();   // PDG Code of Current Environment
  G4int envN    = theEnvironment.GetN();         // N of current Nuclear Environment
  G4QContent envQC=theEnvironment.GetQCZNS();    // QuarkContent of the current Nuclear Environ.
  G4double addPhoton=phot4M.e();                 // Photon energy for capture by quark-parton
#ifdef sdebug
  cout<<"G4Quasmon::CalculateHadrProbab:QM="<<mQ<<valQ<<",vaf="<<vaf<<",mC="<<maxB<<endl;
#endif
  // ================= Calculate probabilities for candidates
  for (G4int index=0; index<theQCandidates.entries(); index++)
  {
    G4QCandidate* curCand=theQCandidates[index];
    curCand->ClearParClustVector();             // Clear Parent Cluster Vector for the Fragment
    G4double probability = 0.;
    G4double secondProbab = 0.;
    G4int    resPDG=0;
  	G4double comb=0.;
    G4int cPDG = curCand->GetPDGCode();
    G4QContent candQC = curCand->GetQC();
    G4int cU=candQC.GetU()-candQC.GetAU();
    G4int dU=cU+cU;
    G4int cD=candQC.GetD()-candQC.GetAD();
    G4int dD=cD+cD;
    G4int dUD=abs(cU-cD);
    G4int cS=candQC.GetS()-candQC.GetAS();
    G4bool pos=curCand->GetPossibility();
	if(pos&&(cPDG<80000000||(cPDG>80000000&&cPDG!=90000000&&dUD<2))) //@@=>  2 || 3
	{
      G4QContent curQ = valQ;                   // Make a current copy of the Quasmon QuarkCont
      G4int aPDG = abs(cPDG);
      G4int baryn= candQC.GetBaryonNumber();    // Baryon number of the Candidate
      G4int cC   = candQC.GetCharge();          // Charge of the Candidate
      G4int cI   = baryn-cC-cC;
      G4int cNQ  = candQC.GetTot()-2;           // A#of quark-partons in a candidate -2
      G4double resM=0.;                         // Prototype for minMass of residual hadron
      if(cPDG>80000000&&cPDG!=90000000&&baryn<=maxB)// ===>Nuclear Case (QUark EXchange mechanism)
      {
        G4int      pc=0;                        // Parents counter
        G4double   pcomb=0.;                    // Summed probability of parent clusters
        G4QPDGCode cQPDG(cPDG);                 // QPDG for the nuclear candidate to a fragment
        G4double   frM=cQPDG.GetMass();         // Vacuum mass of the candidate to nuc.fragment
        G4double   cZ=candQC.GetCharge();
        G4double   CB=CoulombBarrier(totZ,totA,cZ,baryn);
        G4double   qMax=frM+CB-kLS;
        G4double   qM2=qMax+qMax;
        G4int iQmin=0;
        G4int iQmax=3;
        iQmax=2;                                //@@@@@@@@@@ TEMPORARY for PICAPTURE
        G4int oQmin=0;
        G4int oQmax=3;
        oQmax=2;                                //@@@@@@@@@@ TEMPORARY for PICAPTURE
		if     (dofU<=nofD) iQmax=1;            // d-quark (u-quark is forbiden)
        else if(dofD<=nofU) iQmin=1;            // u-quark (d-quark is forbiden)
#ifdef pdebug
        if(cPDG==90001000||cPDG==90000001||cPDG==90000002||cPDG==90001001)
        cout<<"G4Quasmon::CHP:***>>>c="<<cPDG<<",m="<<frM<<",imi="<<iQmin<<",ima="<<iQmax<<endl;
#endif
        if(oQmax>oQmin) for(int iq=iQmin; iq<iQmax; iq++) // Entering (from Quasmon) quark
		{
          G4double qFact=1.;
          if(j==1&&iq==1&&addPhoton>0.) qFact=4.;
          G4double nqInQ=0.;                    // A#of quarks of this sort in a Quasmon
          if     (!iq)   nqInQ=valQ.GetD();
          else if(iq==1) nqInQ=valQ.GetU();
          else if(iq==2) nqInQ=valQ.GetS();
          comb=0.;                              // Local summ for the i-quark of the Quasmon
#ifdef sdebug
          cout<<"G4Q::CHP:i="<<iq<<",cU="<<cU<<",cD="<<cD<<",omi="<<oQmin<<",oma="<<oQmax<<endl;
#endif
          if(oQmax>oQmin) for(int oq=oQmin; oq<oQmax; oq++) // Exiting (to Quasmon) quark
		  {
            G4int shift= cQPDG.GetRelCrossIndex(iq, oq);
            G4QContent ioQC=cQPDG.GetExQContent(iq, oq);
            G4QContent resQC=valQ+ioQC;         // Quark Content of the residual Quasmon
#ifdef sdebug
			cout<<"G4Quasmon::CHadProb:iq="<<iq<<",oq="<<oq<<",QC="<<ioQC<<",rQC="<<resQC<<endl;
#endif
            G4QPDGCode resQPDG(resQC);          // QPDG of the residual Quasmon
            resPDG=resQPDG.GetPDGCode();        // just for the final print (?)
#ifdef pdebug
            if(cPDG==90001000||cPDG==90000001||cPDG==90000002||cPDG==90001001)
			cout<<"G4Quasm::CalculateHadronizationProbability:i="<<iq<<",o="<<oq<<ioQC
                <<",shift="<<shift<<", cQPDG="<<cQPDG<<", residQC="<<resQC<<resQPDG<<endl;
#endif
            if(resQPDG.GetQCode()>-2)           // Such residual Quasmon is possible
			{
              G4int is=index+shift;
              if(shift!=7&&is<maxC)                         // This quark exchange is possible
			  {
                G4QCandidate* parCand=theQCandidates[is];   // Pointer to ParentCluster of Cand.
                G4int possib=parCand->GetParPossibility();
#ifdef pdebug
                if(cPDG==90001000||cPDG==90000001||cPDG==90000002||cPDG==90001001)
			    cout<<"G4Quasm::CalculateHadronizationProbability: possibility="<<possib<<endl;
#endif
                if(possib)
                {
                  G4QContent rQQC = valQ+ioQC;              // Quark Content of Residual Quasmon
                  G4int rQU=rQQC.GetU()-rQQC.GetAU();
                  G4int rQD=rQQC.GetD()-rQQC.GetAD();
                  G4QContent parQC = parCand->GetQC();      // Quark Content of Parent Cluster
                  G4int barot = parQC.GetBaryonNumber();    // Baryon Number of Parent Cluster
                  G4int charge= parQC.GetCharge();          // Charge of Parent Cluster
                  G4int isos  = barot-charge-charge;        // Isospin of Parent Cluster
                  G4double pUD= abs(isos)+1.;
                  if(barot!=baryn) cerr<<"*G4Q::CHPr:cand="<<candQC<<",par="<<parQC<<",shift="
                                       <<shift<<",ind="<<index<<",is="<<is<<endl;
                  G4int    zZ=qChg*charge;                  // ChargeProduct for ParClust & ResQ
                  G4int    dI=qIso-isos;                    // IsotopicShiftSum for ParCl & ResQ
                  //********* ISOTOPIC  FOCUSING MATRIX *******************
				  if(abs(dI)< 1 || abs(dI)<2 && abs(cC-charge)<2  ||
					 ((barot==1 || barot>=4) && (dI>=2&&cC<charge || dI<=-2&&cC>charge)) || 
					 ((barot==2 || barot==3) && (dI> 2&&cC<charge || dI<=2&&cC<=charge || dI<=-barot))
					)
                  {
                    G4double pPP=parCand->GetPreProbability();// Probability of Parent Cluster
                    G4int    parPDG=parCand->GetPDGCode();    // PDGCode of the Parent Clucter
                    G4double boundM=parCand->GetMass();       // Bound mass of Parent Cluster
#ifdef pdebug
                    if(cPDG==90001000||cPDG==90000001||cPDG==90000002||cPDG==90001001)
					  cout<<"G4Q::CHP:c="<<cPDG<<",p="<<parPDG<<",pM="<<boundM<<",i="<<is
                          <<",a="<<parCand<<endl;
#endif
                    // Kinematical analysis of decay possibility
                    G4double   minM  =0.;                     // Prototype of minM of ResidQuasm
                    if (resPDG==10)minM=G4QChipolino(resQC).GetMass();// ResidualQuasmonChipo
                    else if(resPDG)minM=G4QPDGCode(resPDG).GetMass(); // ResidualQuasmonHadron
                    //G4LorentzVector comp=q4Mom+G4LorentzVector(0.,0.,0.,boundM);// Compound LV
                    G4double kCut=boundM/2.+freeE/(iniQM+boundM);
#ifdef pdebug
                    if(cPDG==90001000||cPDG==90000001||cPDG==90000002||cPDG==90001001)
					  cout<<"G4Quasm::CHP:cPDG="<<cPDG<<",miM="<<minM<<",kLS="<<kLS<<",c="<<kCut
                          <<",fE="<<freeE<<",iM="<<iniQM<<",PDG="<<parPDG<<",bM="<<boundM<<endl;
#endif
                    //if(resPDG&&minM>0.&&comp.m()>minM) //Kinematical analysis of hadronization
                    ////***Cap***???/////
                    if(resPDG&&minM>0.&&kLS<kCut) // Kinematical analysis
					//if(resPDG&&minM>0.) // Kinematical analysis of hadronization
                    {
#ifdef sdebug
                      cout<<"G4Quasmon::HadrProbab: fM="<<frM<<", bM="<<boundM<<endl;
#endif
                      G4double pmk=rMo*boundM/kLS;
                      G4double delta=(frM*frM-boundM*boundM)/(boundM+boundM);
                      G4double kd =kLS-delta;
                      G4double dkd=kd+kd;
                      G4double hz=1.-dkd/(boundM+kLS+kLS);
                      G4double kf=0.;
                      if(hz>0.)kf=pow(hz,cNQ);
                      if(envM>boundM)
                      {
                        G4QContent rtQC=valQ+envQC-parQC; // Total Residual Quark Content
                        G4QNucleus rtN(rtQC);             // Create PseudoNucleus for totResid.
                        G4double rtM=rtN.GetGSMass();     // Mass of residual Q+E system
                        minM=rtM-envM+boundM;
                      }
                      G4double minM2=minM*minM;
                      G4double ph=kf;
					  G4double newl=0.;
					  G4double newh=1.;
                      if(minM2>rQ2)                       // ==> Check of residual Quasm. Mass
                      {
                        G4double nz=1.-(minM2-rQ2)/(boundM*rEP);
                        G4bool atrest=(eQ-mQ)/mQ<.001;
                        if(atrest) nz=1.-(minM2-rQ2+pmk*dkd)/(boundM*(rEP+pmk));
                        if(nz>0.&&nz<1.)
                        {
                          newh=pow(nz,cNQ);
                          //if(atrest && newh<kf) kf=newh;
                          if(newh<kf)kf=newh;
                        }
                      }
                      if(kf<0.)kf=0.;
                      if(kf>1.)kf=1.;
                      G4double high = kf;
#ifdef pdebug
                      if(cPDG==90001000||cPDG==90000001||cPDG==90000002||cPDG==90001001)
						cout<<"G4Quasmon::CHProb:kf="<<kf<<",mM2="<<minM2<<",rQ2="<<rQ2<<endl;
#endif
                      G4double lz=1.-dkd/boundM;
                      G4double low=0.;
                      if(lz>0.)low=pow(lz,cNQ);// Low limit for randomization
                      else
                      {
#ifdef sdebug
                        cerr<<"***G4Q::HP:l="<<low<<",d="<<dkd<<",m="<<boundM<<",n="<<cNQ<<endl;
#endif
                        low=0.;
                      }
                      G4double pl=low;
                      if(totZ>cZ)                         // ==> Check of Coulomb Barrier
                      {
                        G4double nz=(qMax+qMax)/boundM-1.;
                        if(nz>0.&&nz<1.)
                        {
                          newl=pow(nz,cNQ);
                          if(newl>low) low=newl;
                        }
                        else if(nz>1.) low=10.;
                      }
                      kf-=low;
#ifdef pdebug
                      if(cPDG==90001000||cPDG==90000001||cPDG==90000002||cPDG==90001001)
                      cout<<"G4Q::HP:>>cPDG="<<cPDG<<",l="<<low<<",h="<<high<<",ol="<<pl<<",oh="
                          <<ph<<",nl="<<newl<<",nh="<<newh<<",kf="<<kf<<",delta="<<delta<<endl;
#endif
                      G4double probab=0.;
                      if(kf>0)
                      {
                        kf*=boundM/kLS/cNQ;      // Final value of kinematical (i,o) factor
                        G4int noc=cQPDG.GetNumOfComb(iq, oq);
                        probab=qFact*kf*nqInQ*pPP*noc/pUD;
#ifdef sdebug
                        cout<<"G4Q::HP:p="<<probab<<",qF="<<qFact<<",iq="<<iq<<",oq="<<oq<<",j="
                            <<j<<",Ph="<<phot4M<<endl;
#endif
                        if(probab<0.) probab=0.;
                      }
                      pcomb += probab;           // Update integrated probability for ParntClust
                      G4QParentCluster* curParC = new G4QParentCluster(parPDG,pcomb);
                      curParC->SetTransQC(ioQC); // Keep Quark Content of the Exchange Meson
                      curParC->SetLow(low);      // Keep the Low limit of randomization
                      curParC->SetHigh(high);    // Keep the High limit of randomization
                      curParC->SetMass(boundM);  // Keep bounded mass for future calculations
                      curParC->SetBind(delta);   // Keep BindingEnerergy for future calculations
                      curParC->SetNQPart2(cNQ);  // Keep a#of quark-partons in the fragment
#ifdef sdebug
		 	          cout<<"G4Quasm::CalculateHadrProbab: FillParentClaster="<<*curParC<<endl;
#endif
                      curCand->FillPClustVec(curParC);//Fill possible ParentClust to FragmVector
                      delete curParC;
                      comb += probab;
#ifdef sdebug
		 	          cout<<"G4Quasmon::CalcHP:i="<<index<<",cC="<<cPDG<<",pc"<<pc<<parQC<<",E="
						 <<theEnvironment<<",p="<<parCand->GetParPossibility()<<","<<comb<<endl;
#endif
                      pc++;
					}
				  }
#ifdef sdebug
				  else cout<<"*G4Q::HPr:z="<<zZ<<",dI="<<dI<<",cC="<<cC<<",rQ="<<rQ<<endl;
#endif
				}
			  }                                // >>>> End of if of Quark Exchange possibility
			}                                  // +++> End of if of existinr residual Q Code
            probability+=comb;                 // Collect the probability for the fragment
#ifdef sdebug
	 	    cout<<"G4Quasmon::CHPr:pr="<<probability<<"(+"<<comb<<"),iq="<<iq<<",oq="<<oq<<endl;
#endif
		  }                                    // ...> End of Quark Exchange "oq" Test LOOP
		}                                      // ...> End of Quark Exchange "iq" Test LOOP
	  }                                        // ---> End of Nuclear Case of fragmentation
      else if(cPDG<80000000)                     // ===> Hadronic case (QUark FUsion mechanism)
      {
        comb = valQ.NOfCombinations(candQC);
        if(cPDG==111||aPDG==211||cPDG==221||cPDG==331||cPDG==113||aPDG==213||cPDG==223||
           cPDG==333||cPDG==115||aPDG==215||cPDG==225||cPDG==335||cPDG==110||cPDG==220||
           cPDG==330)
        {
          G4QContent tQCd(1,0,0,1,0,0);
          G4QContent tQCu(0,1,0,0,1,0);
          G4QContent tQCs(0,0,1,0,0,1);
          G4double cmd=valQ.NOfCombinations(tQCd);
          G4double cmu=valQ.NOfCombinations(tQCu);
          G4double cms=valQ.NOfCombinations(tQCs);
          if(cPDG!=333&&cPDG!=335&&aPDG!=211&&aPDG!=213&&aPDG!=215) comb=(cmd+cmu)/2.;
          if(cPDG==331||cPDG==221)comb=(comb + cms)/2.; //eta,eta'
          // Hadronization in eta/eta' suppression (conversion to rho/omega)
          if(cPDG==111&&j==1||cPDG==221||cPDG==331) comb*=EtaEtaprime;
          if(cPDG==113||aPDG==213) comb*=2.-EtaEtaprime;
          if(j==1&&cPDG==223) comb*=5.-4*EtaEtaprime;
          else if (cPDG==223) comb*=4.-3*EtaEtaprime;
          if(j>1&&cPDG==111) comb*=2-EtaEtaprime;
        }
        curQ -= candQC;                           // This is a quark content of residual quasmon
        resPDG = curQ.GetSPDGCode();              // PDG for the lowest residual hadronic state
#ifdef sdebug
        //if(cPDG==2112||cPDG==2212||cPDG==3122)
	    cout<<"G4Q:CHProb:cQC="<<candQC<<",iQC="<<valQ<<",rPDG=0,rQC="<<curQ<<endl;
#endif
        if (resPDG==221 && (mQ<683.|| absb)) resPDG=111;// pi0 minimum residual instead of eta
        if (resPDG==111 && (mQ>683.&&!absb)) resPDG=221;// eta minimum residual instead of pi0
#ifdef sdebug
        //if(cPDG==2112||cPDG==2212||cPDG==3122)
        cout<<"G4Quasmon::CalcHadrProbability:comb="<<comb<<",rPDG="<<resPDG<<curQ<<endl;
#endif
        if(comb && resPDG && (resPDG>80000000&&resPDG!=90000000||resPDG<10000))
	    {
#ifdef sdebug
          //if(cPDG==2112||cPDG==2212||cPDG==3122)
          cout<<"G4Q:HPr:NEW"<<index<<",Q="<<valQ<<",cPDG="<<cPDG<<",rPDG="<<resPDG<<curQ<<endl;
#endif
          if(resPDG!=10)resM=G4QPDGCode(resPDG).GetMass();// PDG mean mass for the resid. hadron
          else resM=G4QChipolino(curQ).GetMass();  // Chipolino mass for the residual hadron
          G4int resQCode=G4QPDGCode(curQ).GetQCode();
#ifdef sdebug
          //if(cPDG==2112||cPDG==2212||cPDG==3122)
          cout<<"G4Q:HPr:ResidQMass="<<resM<<curQ<<",envPDG="<<envPDGC<<",rQC="<<resQCode<<endl;
#endif
          if(envPDGC>80000000&&envPDGC!=90000000&&resM>0.&&cPDG>1000&&resPDG!=10) // => Environment
          { 
            G4QContent rtQC=curQ+envQC;            // Total Residual Quark Content
            G4QNucleus rtN(rtQC);                  // Create a pseudo-nucleus for residual
            G4double rtM =rtN.GetMZNS();           // Min Mass of total residual Nucleus
            G4double bnRQ=rtM-envM;                // Bound mass of residual Quasmon
#ifdef sdebug
            //if(cPDG==2112||cPDG==2212||cPDG==3122)
            cout<<"G4Q:HPr: Recalculated RQMass="<<bnRQ<<",envM="<<envM<<",rtM="<<rtM<<endl;
#endif
            if(bnRQ<resM) resM=bnRQ;
          }
          if(resM>0. && resQCode>-2)
	      {
            G4double rndM=GetRandomMass(cPDG,mQ-resM); // Candidate's Mass randomization
#ifdef sdebug
            //if(cPDG==2112||cPDG==2212||cPDG==3122)
            cout<<"G4Q:HPr: CandMassRNDM="<<rndM<<",R+R="<<resM+rndM<<" < mQ="<<mQ<<endl;
#endif
            // --- Kinematical Factors ---
            if(rndM>0. && resM+rndM<mQ)
	        {
              curCand->SetMass(rndM);              // Randomized a mass value of the Candidate
              G4double mH2 = rndM*rndM;            // Squared mass of the candidate (Mu2)
              G4double rHk = mH2/dk;
              G4double zMax = 1.-rHk/mQ;           // z_max
              G4double mR2 = resM*resM;
              //G4double dSM = mQ2+mH2-mR2;
              //G4double det = sqrt(dSM*dSM-4.*mQ2*mH2);
              ////////G4double zMin= mR2/mQ/(mQ-dk);       // z_min //@@ ??
              ////////zMin=0.;
			  G4double possibility=zMax;
			  ////////G4double possibility=zMax-zMin;
#ifdef sdebug
              //if(cPDG==2112||cPDG==2212||cPDG==3122)
              cout<<"G4Q::QHadProb: M(PDG="<<cPDG<<", QC="<<curCand->GetQC()<<")="<<rndM
                  <<",pos="<<possibility<<",rHk="<<rHk<<",dk="<<dk<<",resM="<<resM<<endl;
#endif
              if (resPDG==10)                      // Chipolino case - check minimum
		      {
                G4double rM2 = mQk*(mQ-rHk);
                if(rM2<resM*resM) possibility = 0.;
              }
              if (possibility > 0. && vap > 0)
	          {
                probability = vaf*pow(zMax, vap);
                //probability = vaf*(pow(zMax, vap)-pow(zMin, vap));
#ifdef sdebug
                //if(cPDG==2112||cPDG==2212||cPDG==3122)
                cout<<"G4Q::HProb:#"<<index<<",m2="<<mH2<<",n="<<nOfQ<<",p="<<probability
                    <<",vaf="<<vaf<<",vap="<<vap<<",zMax="<<zMax<<endl;
#endif
                
                if(qBar>1&&baryn>0)                //---> High baryon number ("nuclear") case
                {
                  G4QContent rtQC=curQ+envQC;            // Total Residual Quark Content
                  G4QNucleus rtN(rtQC);                  // Create a pseudo-nucleus for residual
                  G4double rtM =rtN.GetMZNS();           // Min Mass of total residual Nucleus
                  G4double bnRQ=rtM-envM;                // Bound mass of residual Quasmon
                }
                else                               //---> Low baryon number case (tuned on p-ap)
                {
                  if(cPDG==110||cPDG==220||cPDG==330) 
                       probability*=comb;                // f_0's (sigma) have spin 0
			      else probability*=comb*(abs(cPDG)%10); // Spin of resonance
				}
	          }
            }
            else
		    {
#ifdef sdebug
              //if(cPDG==2112||cPDG==2212||cPDG==3122)
              cout<<"G4Q::QHadProb:crM=0["<<cPDG<<"],mQ="<<mQ<<valQ<<",resM="<<resM<<curQ<<endl;
#endif
            }
	      }
          else
	      {
#ifdef sdebug
            //if(cPDG==2112||cPDG==2212||cPDG==3122)
            cout<<"***G4Q:HadProb:M=0,#"<<index<<valQ<<",cP="<<cPDG<<"+rP="<<resPDG<<curQ<<endl;
#endif
	      }
        }
        else
	    {
#ifdef sdebug
          //if(cPDG==2112||cPDG==2212||cPDG==3122)
          cout<<"G4Q:HadProb:comb=0 @ #"<<index<<valQ<<",cP="<<cPDG<<"+rP="<<resPDG<<curQ<<endl;
#endif
	    }
        //if(cPDG==111) secondProbab = 1.;
      }// ---> End of Hadronic Case of fragmentation
      else probability=0.;
#ifdef pdebug
      if(cPDG==90001000||cPDG==90000001||cPDG==90000002||cPDG==90001001)
      cout<<"QHadProb:^^cPDG="<<cPDG<<",pos="<<pos<<",rPDG="<<resPDG<<curQ<<resM<<", p="
          <<probability<<", s="<<accumulatedProbability<<",sp="<<secondProbab<<endl;
#endif
    }                                                 // ===> End of possibility check
    curCand->SetRelProbability(probability);
    accumulatedProbability += probability;
    curCand->SetIntegProbability(accumulatedProbability);
    curCand->SetSecondRelProb(secondProbab);
    secondAccumProbability += secondProbab;
    curCand->SetSecondIntProb(secondAccumProbability);
  }                                                   // ***> End of LOOP over candidates}
} // End of "CalculateHadronizationProbability"

// Decay of outgoing hadron
G4bool G4Quasmon::DecayOutHadron(G4QHadron* qHadron)
{ //   =============================================
  G4QPDGCode theQPDG  = qHadron->GetQPDG();
  G4int      thePDG   = theQPDG.GetPDGCode();         // Get the PDG code of decaying hadron
  G4int        pap = 0;                               // --- particle
  if(thePDG<0) pap = 1;                               // --- anti-particle
  G4LorentzVector t = qHadron->Get4Momentum();        // Get 4-momentum of decaying hadron
  G4double m = t.m();                                 // Get the mass value of decaying Hadron
  // --- Randomize a channel of decay
  G4QDecayChanVector decV = theWorld.GetQParticle(theQPDG)->GetDecayVector();
  G4int nChan = decV.entries();
#ifdef pdebug
  G4cout<<"G4Quasmon::DecayOutHadron: PDG="<<thePDG<<", m="<<m<<",("<<nChan<<" channels)"<<G4endl;
#endif
  if(nChan)
  {
    G4int i=0;
    if(nChan>1)
	{
      G4double rnd = G4UniformRand();                 // Random value to select a Decay Channel
      for(i=0; i<nChan; i++)
	  {
        G4QDecayChan* dC = decV[i];                   // The i-th Decay Channel
#ifdef pdebug
  cout<<"G4Quasmon::DecayOutHadron: i="<<i<<",r="<<rnd<<" < dl="<<dC->GetDecayChanLimit()
      <<", mm="<<dC->GetMinMass()<<endl;
#endif
        if(rnd<dC->GetDecayChanLimit() && m>dC->GetMinMass()) break;
	  }
	}
    G4QPDGCodeVector cV=decV[i]->GetVecOfSecHadrons();// PDGVector of the selected Decay Channel
    G4int nPart=cV.entries();                         // A#of particles to decay in
#ifdef pdebug
	cout<<"G4Quasmon::DecayOutHadron: resi="<<i<<",nP="<<nPart<<":"<<cV[0]->GetPDGCode()
        <<","<<cV[1]->GetPDGCode();
    if(nPart>2) cout<<","<<cV[2]->GetPDGCode();
    cout<<endl;
#endif
    if(nPart<2||nPart>3)
    {
      cerr<<"***G4Q::DecayOutH: n="<<nPart<<", ch#"<<i<<",PDG="<<thePDG<<endl;
      return false;
	}
#ifdef debug
    G4cout<<"G4Quasmon::DecayOutHadron: Fill PDG= "<<thePDG<<t<<m<<"(nP="<<nPart<<")>>>>>"<<G4endl;
#endif
    if(nPart==2)
	{
      G4QHadron* fHadr;
      G4QHadron* sHadr;
      G4int fPDG=cV[0]->GetPDGCode();
      G4int sPDG=cV[1]->GetPDGCode();
      if(cV[0]->GetWidth()==0.)
	  { // Randomize only the second Hardon or noone
        fHadr = new G4QHadron(cV[0]->GetPDGCode());   // the First Hadron is created *1*
        if(cV[1]->GetWidth()==0.)sHadr = new G4QHadron(sPDG);// the Second Hadron is created *2*
        else
		{
          G4QParticle* sPart=theWorld.GetQParticle(cV[1]);// Particle for the Second Hadron
          G4double sdm = m - fHadr->GetMass();          // MaxMassLimit for the 2-nd Hadron
          sHadr = new G4QHadron(sPart,sdm);             // the Second Hadron is created *2*
          if(sPDG<0) sHadr->MakeAntiHadron();
        }
	  }
      else
	  {  // Randomize masses of both Hadrons
        G4QParticle* sPart=theWorld.GetQParticle(cV[1]);// Particle for the Second Hadron
        G4double mim = sPart->MinMassOfFragm();         // MinMassLimit for the Second Hadron
        G4double fdm = m - mim;                         // MaxMassLimit for the First Hadron
        G4QParticle* fPart=theWorld.GetQParticle(cV[0]);// Particle for the First Hadron
        fHadr = new G4QHadron(fPart,fdm);               // the 1-st Hadron is initialized *1*
        if(fPDG<0) fHadr->MakeAntiHadron();
        G4double sdm = m - fHadr->GetMass();            // MaxMassLimit for the Second Hadron
        sHadr = new G4QHadron(sPart,sdm);               // the 2-nd Hadron is initialized *2*
        if(sPDG<0) sHadr->MakeAntiHadron();
      }
#ifdef pdebug
  cout<<"G4Quasmon::DecayOutHadron:2Dec. m1="<<fHadr->GetMass()<<",m2="<<sHadr->GetMass()<<endl;
#endif
      if(pap)
      {
        fHadr->MakeAntiHadron();
        sHadr->MakeAntiHadron();
      }
      G4LorentzVector f4Mom = fHadr->Get4Momentum();   // Get 4-mom of the First Hadron (mass) 
      G4LorentzVector s4Mom = sHadr->Get4Momentum();   // Get 4-mom of the Second Hadron (mass) 
      if(!qHadron->DecayIn2(f4Mom,s4Mom))
      {
        delete fHadr;                                  // Delete "new fHadr"
        delete sHadr;                                  // Delete "new sHadr"
        G4cerr<<"***G4Q::DecayOutHadron: PDGC="<<thePDG<<", ch#"<<i<<G4endl;
        theQHadrons.insert(qHadron);                   // Fill as it is (delete equivalent)
        return false;
	  }
      else
      {
        qHadron->SetNFragments(2);
        theQHadrons.insert(qHadron);                   // Fill with NFr=2 (delete equivalent)
        fHadr->Set4Momentum(f4Mom);                    // Put the randomized 4Mom to 1-st Hadron
        FillHadronVector(fHadr);                       // Fill 1st Hadron (delete equivalent)
        sHadr->Set4Momentum(s4Mom);                    // Put the randomized 4Mom to 2-nd Hadron
        FillHadronVector(sHadr);                       // Fill 2nd Hadron (delete equivalent)
      }
    }
    else if(nPart==3)
	{
      G4QHadron* fHadr;
      G4QHadron* sHadr;
      G4QHadron* tHadr;
      G4int fPDG=cV[0]->GetPDGCode();
      G4int sPDG=cV[1]->GetPDGCode();
      G4int tPDG=cV[2]->GetPDGCode();
      if(cV[0]->GetWidth()==0.)                        // Don't randomize the First Hardon
	  {
        fHadr = new G4QHadron(fPDG);                   // the First Hadron is created  *1*
        if(cV[1]->GetWidth()==0.)
        {
          sHadr = new G4QHadron(sPDG);                         // the Second Hadron is created *2*
          if(cV[2]->GetWidth()==0.)tHadr = new G4QHadron(tPDG);//the Third Hadron is created *3*
          else
		  {
            G4QParticle* tPart=theWorld.GetQParticle(cV[2]);   // Particle for the 3-rd Hadron
            G4double tdm = m-fHadr->GetMass()-sHadr->GetMass();// MaxMassLimit for the 2d Hadron
            tHadr = new G4QHadron(tPart,tdm);                  // the 3-rd Hadron is created *3*
            if(tPDG<0) tHadr->MakeAntiHadron();
		  }
        }
        else                                              // Randomize 2nd & 3rd Hadrons
		{
          m-=fHadr->GetMass();                            // Reduce the residual MaxMassLimit
          G4QParticle* tPart=theWorld.GetQParticle(cV[2]);// Particle for the 3-rd Hadron
          G4double mim = tPart->MinMassOfFragm();         // MinMassLimit for the 3rd Hadron
          G4double sdm = m - mim;                         // MaxMassLimit for the 2nd Hadron
          G4QParticle* sPart=theWorld.GetQParticle(cV[1]);// Particle for the 2-nd Hadron
          sHadr = new G4QHadron(sPart,sdm);               // the Second Hadron is created *2*
          if(sPDG<0) sHadr->MakeAntiHadron();
          G4double tdm = m - sHadr->GetMass();            // MaxMassLimit for the 3-rd Hadron
          tHadr = new G4QHadron(tPart,tdm);               // the Third Hadron is created *3*
          if(tPDG<0) tHadr->MakeAntiHadron();
        }
	  }
      else  // Randomize masses of all three Hadrons
	  {
        G4QParticle* sPart=theWorld.GetQParticle(cV[1]);   // Particle for the Second Hadron
        G4double smim = sPart->MinMassOfFragm();           // MinMassLimit for the Second Hadron
        G4QParticle* tPart=theWorld.GetQParticle(cV[2]);   // Particle for the Third Hadron
        G4double tmim = tPart->MinMassOfFragm();           // MinMassLimit for the Third Hadron
        G4double fdm = m - smim - tmim;                    // MaxMassLimit for the First Hadron
        G4QParticle* fPart=theWorld.GetQParticle(cV[0]);   // Particle for the First Hadron
        fHadr = new G4QHadron(fPart,fdm);                  // the First Hadron is created *1*
        if(fPDG<0) fHadr->MakeAntiHadron();
        m-=fHadr->GetMass();                               // Reduce the residual MaxMassLimit !
        G4double  sdm = m - tmim;                          // MaxMassLimit for the Second Hadron
        sHadr = new G4QHadron(sPart,sdm);                  // the Second Hadron is created *2*
        if(sPDG<0) sHadr->MakeAntiHadron();
        G4double  tdm = m - sHadr->GetMass();              // MaxMassLimit for the Third Hadron
        tHadr = new G4QHadron(tPart,sdm);                  // the Third Hadron is created *3*
        if(tPDG<0) tHadr->MakeAntiHadron();
      }     
#ifdef pdebug
      G4cout<<"G4Quasmon::DecayOutHadron:3Dec. m1="<<fHadr->GetMass()<<",m2="<<sHadr->GetMass()
            <<",m3="<<tHadr->GetMass()<<G4endl;
#endif
      if(pap)
      {
        fHadr->MakeAntiHadron();
        sHadr->MakeAntiHadron();
        tHadr->MakeAntiHadron();
      }
      G4LorentzVector f4Mom = fHadr->Get4Momentum();   // Get 4-mom of the First Hadron (mass) 
      G4LorentzVector s4Mom = sHadr->Get4Momentum();   // Get 4-mom of the Second Hadron (mass) 
      G4LorentzVector t4Mom = tHadr->Get4Momentum();   // Get 4-mom of the Third Hadron (mass) 
      if(!qHadron->DecayIn3(f4Mom,s4Mom,t4Mom))
      {
        delete fHadr;                                  // Delete "new fHadr"
        delete sHadr;                                  // Delete "new sHadr"
        delete tHadr;                                  // Delete "new tHadr"
        G4cerr<<"***G4Q::DecayOutHadron: PDGC="<<thePDG<<", ch#"<<i<<G4endl;
        theQHadrons.insert(qHadron);                   // Fill as it is (delete equivalent)
        return false;
	  }
      else
      {
        qHadron->SetNFragments(3);
        theQHadrons.insert(qHadron);                   // Fill with NFr=3 (delete equivalent)
        fHadr->Set4Momentum(f4Mom);                    // Put the randomized 4Mom to 1-st Hadron
        FillHadronVector(fHadr);                       // Fill 1st Hadron (delete equivalent)
        sHadr->Set4Momentum(s4Mom);                    // Put the randomized 4Mom to 2-nd Hadron
        FillHadronVector(sHadr);                       // Fill 2nd Hadron (delete equivalent)
        tHadr->Set4Momentum(t4Mom);                    // Put the randomized 4Mom to 3-rd Hadron
        FillHadronVector(tHadr);                       // Fill 3rd Hadron (delete equivalent)
      }
	}
  }
  else
  {
#ifdef pdebug
    G4cout<<"G4Quasmon::DecayOutHadron: Fill PDG= "<<thePDG<<t<<m<<"***0***>>>>>"<<G4endl;
#endif
    theQHadrons.insert(qHadron);                       // Fill as it is (delete equivalent)
  }
  return true;                                         // FAKE (impossible to reach)
}

// Random integer value for the Poiasson Distribution with meanValue
G4int G4Quasmon::RandomPoisson(G4double meanValue)
//               =================================
{
  if (meanValue<=0.)
  {
    cerr<<"***G4Quasmon::RandomPoisson: negative or zero Mean Value="<<meanValue<<endl;
    //G4Exception("***G4Quasmon::RandomPoisson: negative 0r zero Mean Value");
    return -1;
  }
  G4double r=G4UniformRand();
  G4double t=exp(-meanValue);
  G4double s=t;
  if (r<s) return 0;
  t*=meanValue; // To avoid /1
  s+=t;
  if (r<s) return 1;
  G4int i=1;
  while ( s<r && i<100 )
  {
    i++;
    t*=meanValue/i;
    s+=t;
  }
  return i;
}

//Coulomb Barrier
G4double G4Quasmon::CoulombBarrier(const G4double& tZ, const G4double& tA,
                                   const G4double& cZ, const G4double& cA)
{//                 ======================================================
  static const G4double third=1./3.;
  G4double rZ=tZ-cZ;
  G4double rA=tA-cA;
  G4double zz=rZ*cZ;
  //G4double r=(pow(rA,third)+pow(cA,third))*(1.51+.00921*zz)/(1.+.009443*zz);
  G4double r=1.*(pow(rA,third)+pow(cA,third)); // {1.44=200(MeV)/137}*z*Z/R0*(a**1/3+A**1/3)
  return   zz/r;
  //return   exp(-cA*rA/2000.)*zz/r;
}

//The public Hadronisation routine with delete responsibility of User (!)
G4QHadronVector* G4Quasmon::Fragment(G4QNucleus& nucEnviron)
{//              ===========================================
#ifdef pdebug
  G4cout<<"G4Quasmon::Fragment is called theEnviron="<<nucEnviron<<G4endl;
#endif
  HadronizeQuasmon(nucEnviron);
  G4int nHadrs=theQHadrons.entries();
#ifdef ppdebug
  G4cout<<"G4Quasmon::Fragment after HadronizeQuasmon nH="<<nHadrs<<G4endl;
#endif
  G4QHadronVector* theFragments = new G4QHadronVector;
  for (int hadron=0; hadron<nHadrs; hadron++)
  {
    G4QHadron* curHadr = new G4QHadron(theQHadrons[hadron]);
    theFragments->insert(curHadr);
  }
  return theFragments;
}
