// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Quasmon.cc,v 1.33 2001-10-31 17:34:04 mkossov Exp $
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
//#define tdebug
//#define sdebug

#include "G4Quasmon.hh"

G4Quasmon::G4Quasmon(G4QContent qQCont, G4LorentzVector q4M, G4LorentzVector ph4M):
  q4Mom(q4M), valQ(qQCont), phot4M(ph4M)
{
#ifdef debug
  G4cout<<"G4Quasmon: **DIRECT CREATION** QC="<<qQCont<<",Q4M="<<q4M<<",photonE="<<ph4M.e()<<G4endl;
#endif
  if(phot4M.e()>0.) q4Mom+=phot4M; // In case of capture by a quark it will be subtracted back
  G4int nP= theWorld.GetQPEntries();      // A#of initialized particles in CHIPS World
  //@@ Make special parametyer to cut high resonances !!
  G4int          nMesons  = 45;
  if     (nP<24) nMesons  =  9;
  else if(nP<41) nMesons  = 18;
  else if(nP<55) nMesons  = 27;
  else if(nP<72) nMesons  = 36;
  G4int          nBaryons = 72;
  if     (nP<35) nBaryons = 16;
  else if(nP<49) nBaryons = 36;
  else if(nP<66) nBaryons = 52;
  //G4int          nClusters= nP-72;
  G4int          nClusters= nP-80;       // "IsoNuclei"
  InitCandidateVector(nMesons,nBaryons,nClusters);
#ifdef debug
  G4cout<<"G4Quasmon: Candidates are initialized. nMesons="<<nMesons<<",nBaryons="<<nBaryons
        <<",nClusters="<<nClusters<<G4endl;
#endif
  valQ.DecQAQ(-1);
  status=4;                                   
}

G4Quasmon::G4Quasmon(const G4Quasmon& right)
{
  q4Mom                 = right.q4Mom;
  valQ                  = right.valQ;
  //theEnvironment        = right.theEnvironment;
  status                = right.status;
  //theQHadrons (Vector)
  G4int nQH             = right.theQHadrons.size();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
    theQHadrons.push_back(curQH);
  }
  theWorld              = right.theWorld;
  phot4M                = right.phot4M;
  nBarClust             = right.nBarClust;
  nOfQ                  = right.nOfQ;
  //theQCandidates (Vector)
  G4int nQC             = right.theQCandidates.size();
  if(nQC) for(G4int iq=0; iq<nQC; iq++)
  {
    G4QCandidate* curQC = new G4QCandidate(right.theQCandidates[iq]);
    theQCandidates.push_back(curQC);
  }
  f2all                 = right.f2all;
  rEP                   = right.rEP;
  rMo                   = right.rMo;
}

G4Quasmon::G4Quasmon(G4Quasmon* right)
{
#ifdef sdebug
  G4cout<<"G4Quasmon::Copy-Constructor: ***CALLED*** E="<<right->theEnvironment<<G4endl;
#endif
  q4Mom                 = right->q4Mom;
  valQ                  = right->valQ;
  //theEnvironment        = right->theEnvironment;
  status                = right->status;
  //theQHadrons (Vector)
  G4int nQH             = right->theQHadrons.size();
#ifdef sdebug
  G4cout<<"G4Quasmon::Copy-Constructor:nQH="<<nQH<<G4endl;
#endif
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
#ifdef debug
    G4cout<<"G4Quasmon:Copy-Constructor:H#"<<ih<<",QH="<<right->theQHadrons[ih]<<G4endl;
#endif
    G4QHadron* curQH    = new G4QHadron(right->theQHadrons[ih]);
    theQHadrons.push_back(curQH);
  }
  theWorld              = right->theWorld;
  phot4M                = right->phot4M;
  nBarClust             = right->nBarClust;
  nOfQ                  = right->nOfQ;
  //theQCandidates (Vector)
  G4int nQC             = right->theQCandidates.size();
#ifdef sdebug
  G4cout<<"G4Quasmon:Copy-Constructor: nCand="<<nQC<<G4endl;
#endif
  if(nQC) for(G4int iq=0; iq<nQC; iq++)
  {
#ifdef sdebug
    G4cout<<"G4Quasmon:Copy-Constructor:C#"<<iq<<",QC="<<right->theQCandidates[iq]<<G4endl;
#endif
    G4QCandidate* curQC = new G4QCandidate(right->theQCandidates[iq]);
    theQCandidates.push_back(curQC);
  }
  f2all                 = right->f2all;
  rEP                   = right->rEP;
  rMo                   = right->rMo;
#ifdef sdebug
  G4cout<<"G4Quasmon:Copy-Constructor: >>>DONE<<<"<<G4endl;
#endif
}

G4Quasmon::~G4Quasmon()
{
  G4int nCan=theQCandidates.size();
#ifdef sdebug
  G4cout<<"G4Quasmon::Destructor before theQCandidates nC="<<nCan<<G4endl;
#endif
  if(nCan) for (G4int ic=0; ic<nCan; ic++)
  {
#ifdef sdebug
    G4cout<<"G4Quasmon::Destructor delete the QCandidate #"<<ic<<G4endl;
#endif
    delete theQCandidates[ic];
#ifdef sdebug
    G4cout<<"G4Quasmon::Destructor after delete"<<G4endl;
#endif
  }
  //theQCandidates.clearAndDestroy();
#ifdef sdebug
  G4cout<<"G4Quasmon::Destructor before theQHadrons"<<G4endl;
#endif
  G4std::for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron());
#ifdef sdebug
  G4cout<<"G4Quasmon::Destructor === DONE ==="<<G4endl;
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
  G4int nQH             = right.theQHadrons.size();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
    theQHadrons.push_back(curQH);
  }
  theWorld              = right.theWorld;
  phot4M                = right.phot4M;
  nBarClust             = right.nBarClust;
  nOfQ                  = right.nOfQ;
  //theQCandidates (Vector)
  G4int nQC             = right.theQCandidates.size();
  if(nQC) for(G4int iq=0; iq<nQC; iq++)
  {
    G4QCandidate* curQC = new G4QCandidate(right.theQCandidates[iq]);
    theQCandidates.push_back(curQC);
  }
  f2all                 = right.f2all;
  rEP                   = right.rEP;
  rMo                   = right.rMo;

  return *this;
} // End of "="

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
#ifdef sdebug
  G4int ind=0;
#endif
  if(maxMes>nOfMesons) maxMes=nOfMesons;
  if(maxMes>=0) for (i=0; i<maxMes; i++) 
  {
    theQCandidates.push_back(new G4QCandidate(mesonPDG[i]));
#ifdef sdebug
    G4cout<<"G4Quasmon::InitCandidateVector: "<<ind++<<", Meson # "<<i<<" with code = "
          <<mesonPDG[i]<<", QC="<<theQCandidates[i]->GetQC()<<" is initialized"<<G4endl;
#endif
  }
  if(maxBar>nOfBaryons) maxBar=nOfBaryons;
  if(maxBar>=0) for (i=0; i<maxBar; i++) 
  {
    theQCandidates.push_back(new G4QCandidate(baryonPDG[i]));
#ifdef sdebug
    G4cout<<"G4Quasmon::InitCandidateVector: "<<ind++<<", Baryon # "<<i<<" with code = "
          <<baryonPDG[i]<< ", QC="<<theQCandidates[i]->GetQC()<<" is initialized"<<G4endl;
#endif
  }
  if(maxClust>=0) for (i=0; i<maxClust; i++) 
  {
    //G4int clustQCode = i+72;          // Q-code of the cluster in the CHIPS world
    G4int clustQCode = i+80;          // Q-code of the cluster in the CHIPS world "IsoNuclei"
    G4QPDGCode clustQPDG;
    clustQPDG.InitByQCode(clustQCode);
    G4int clusterPDG=clustQPDG.GetPDGCode();
    theQCandidates.push_back(new G4QCandidate(clusterPDG));
#ifdef sdebug
    G4cout<<"G4Quasmon::InitCandidateVector:"<<ind++<<", Cluster # "<<i<<" with code = "
          <<clusterPDG<<", QC="<<clustQPDG.GetQuarkContent()<<" is initialized"<<G4endl;
#endif
  }
} // End of "InitCandidateVector"

// Fragmentation of the Quasmon (the main member function)
G4QHadronVector G4Quasmon::HadronizeQuasmon(G4QNucleus& qEnv, G4int nQuasms)
//              ============================================================
{
  //////////static const G4double kpMax= 10.;      // Maximum Momentum/M_RQ for the Residual Quasmon
  static const G4int NUCPDG  = 90000000;
  static const G4double np   = 1877.9; //@@ temporary = a mass of np pair > m_n + m_p
  static const G4double BIG  = 1000000.;
  static const G4double BIG2 = BIG*BIG;
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QContent PiQC(0,1,0,1,0,0);
  static const G4QContent K0QC(1,0,0,0,0,1);
  static const G4QContent KpQC(0,1,0,0,0,1);
  static const G4QContent zeroQC(0,0,0,0,0,0);
  static const G4LorentzVector zeroLV(0.,0.,0.,0.);
  static const G4QNucleus vacuum= G4QNucleus(zeroLV,NUCPDG);
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mTrit= G4QPDGCode(2112).GetNuclMass(1,2,0);
  static const G4double mHel3= G4QPDGCode(2112).GetNuclMass(2,1,0);
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  ///////////static const G4double mLit6= G4QPDGCode(2112).GetNuclMass(3,3,0);
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  ///////////static const G4double m2N  = mNeut*mNeut;
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();
  static const G4double mEta = G4QPDGCode(221).GetMass();
  static const G4double mEtaP= G4QPDGCode(331).GetMass();
  static const G4double diPiM= mPi0 + mPi0;
  static const G4double PiNM = mPi + mNeut;
  ////////////static const G4double mPi2 = mPi * mPi;
  static const G4double mPi02= mPi0* mPi0;
  ////////////static const G4double mK2  = mK  * mK;
  ////////////static const G4double mK02 = mK0 * mK0;
  ////////////static const G4double mEta2= mEta*mEta;
  ////////////static const G4double mEP2 = mEtaP*mEtaP;
  static const G4double mP2  = mProt*mProt;
  ////////////static const G4double CBKinMin = 0.1; // .1 MeV of KinEnergy is needed to try to pen.
  ////////////static const G4double alpha = 1./137.0359895;    // Fine-structure constant
  ////////////static const G4double third =1./3.;
  ////////////static const G4double conCon=197.327;            // Conversion constant (hc)
  ////////////static const G4double rCB   = 1.09;    // R=r_CB*(a^1/3+A^1/3) - CoulBarrierRadius(fm)
#ifdef debug
  G4cout<<"G4Quasmon(G4Q)::HadronizeQuasmon(HQ): *==>>> START QUASMON HADRONIZATION <<<==*"<<G4endl;
#endif
  G4double momPhoton=phot4M.rho();
  G4double addPhoton=phot4M.e();
  G4bool piF = false;                              // FirstFragmentationFlag for PiCaptureAtRest
  G4bool gaF = false;                              // FirstFragmentationFlag for GammaNucleus
  if(addPhoton<0.)                                 // "PionCapture at rest" case
  {
#ifdef debug
	G4cout<<"G4Q::HQ: PionAtRest, addP="<<addPhoton<<", momP="<<momPhoton<<G4endl;
#endif
    addPhoton=0.;
    momPhoton=0.;
    piF=true;
  }
  else if(addPhoton>0.) gaF=true;
  theEnvironment=qEnv;
  if(!status||q4Mom==zeroLV)                       // This Quasmon is done (Sould not be here)
  {
#ifdef debug
    G4cout<<"G4Q::HQ:NOTHING-TO-DO: Q4M="<<q4Mom<<", QEnv="<<theEnvironment<<G4endl;
#endif
    if(addPhoton)
	{
	  G4cerr<<"***G4Quasmon::HQ: Q4M="<<q4Mom<<",status="<<status<<", phE="<<addPhoton<<G4endl;
      G4Exception("G4Quasmon::HadronizeQuasmon: Overhead photon for the Zero Quasmon");
	}
    KillQuasmon();                                 // This Quasmon is done
    qEnv=theEnvironment;                           // Update QEnvironment
    return theQHadrons;
  }
  status=2;                                        // Nothing is done (yet)
  //if(status!=4) status=2;                          // Nothing is done (yet)
  G4int sPDG=0;                                    // Prototype of PDG of a Selected Candidate
  G4int pPDG=NUCPDG;                               // ProtTemporary PDG Code of the Parent Cluster
  G4bool fmh=false;                                // Flag of hadronization in nuclear matter
  /////////////G4double rmM=0.;                                 // Prototype of coalescence mass
  G4double npqp2=0;                                // A#of quark-partons -2 in a selected fragment
  G4double sMass=0.;                               // Mass of selected candidate
  G4double sM2=0.;                                 // Squared mass of selected candidate
  G4int    pBaryn=0;                               // Parent cluster'c Baryon Num for sel.fragment
  G4double pMass=0.;                               // EnvironBoundedParentCluster Mass for sel.fragm
  G4double pNMass=0.;                              // NucleusBoundedParentCluster Mass for sel.fragm
  G4double pM2=0.;                                 // Sq. Bounded par.clust. mass for sel.fragment
  G4double delta=0.;                               // Binding energy
  G4double deltaN=0.;                              // Binding energy in Total Nucleus
  G4double minSqT=0.;                              // MinimalSqMass of ResidualQuasmon
  G4double minSqB=0.;                              // MinimalSqMass of BoundedResidualQuasmon
  G4double minSqN=0.;                              // MinimalSqMass of ResidQuasmon+ResidEnvironment
  G4double hili=0.;                                // High limit of quark exchange randomization
  G4double loli=0.;                                // Low limit of quark exchange randomization
  G4double tmpTM2=BIG2;                            // GSMass of TotalResidualNucleus for Fragment
  G4double reTNM2=0.;                              // Real mass of TotalResidNucleus for Fragment
  G4QContent curQ=zeroQC;                          // ProtTemporary copy of valQ to estimate MinM2
  G4QContent memQ=zeroQC;                          // ProtTemporary copy of valQ to remember state
  G4QContent pQC=zeroQC;                           // ProtTemporary Quark Content of ParentCluster
  G4QContent sQC=zeroQC;                           // ProtTemporary Quark Content of the fragment
  G4QContent transQC=zeroQC;                       // ProtTemporary Quark Content of ExchangeMeson
  G4LorentzVector m4Mom=zeroLV;                    // 4Momentum to memorize a Quasmon's 4-momentum
  G4LorentzVector kp4Mom=zeroLV;                   // 4-mom prototype for kappa (recoil q)
  G4LorentzVector check=-theEnvironment.Get4Momentum()-q4Mom;//4Mom sum to check conservation
  G4int ccheck=-theEnvironment.GetZ()-valQ.GetCharge();      //For checking charge conservation
  G4bool start=true;
  while(theEnvironment.GetPDG()==NUCPDG || start)  // ***=The MAIN LOOP (a LOOP only for vacuum)=***
  {
    start=false;
    G4bool   quexf=false;                          // Flag of successful quark exchange
    G4double qM2  = q4Mom.m2();                    // Current squared mass of Quasmon
    if(qM2<0. && qM2>-.000001)
	{
      qM2=0.;
      G4double p=q4Mom.rho();
      if(p<.01) 
	  {
#ifdef debug
        G4cout<<"G4Q::HQ:NothingToDo after correction: Q4M="<<q4Mom<<", E="<<theEnvironment<<G4endl;
#endif
        KillQuasmon();                                 // This Quasmon is done
        qEnv=theEnvironment;                           // Update QEnvironment
        return theQHadrons;
	  }
      else q4Mom.setE(p);
    }
    G4double quasM= sqrt(qM2);                     // Current mass of Quasmon
    G4double qurF=quasM/(q4Mom.e()-q4Mom.rho());   // Factor for k Lor.Trans. to LS
    G4ThreeVector qltb = q4Mom.boostVector();      // Boost vector for backward Lor.Trans. in LS
    //////////G4double b2=qltb.mag2();                       // beta^2 of Quasmon
#ifdef debug
	G4cout<<"G4Q::HQ: Quasm="<<q4Mom<<",qM="<<quasM<<",qQC="<<valQ<<G4endl;
#endif
    CalculateNumberOfQPartons(quasM);              // Fills PrivateParameter nOfQ (a#OfQPartonsInQ)
    //===============================
    G4int envPDG=theEnvironment.GetPDG();          // PDGCode of the current Nuclear Environment
    G4int envN  =theEnvironment.GetN();            // N of current Nuclear Environment
    G4int envZ  =theEnvironment.GetZ();            // Z of current Nuclear Environment
    //////////////G4int envS  =theEnvironment.GetS();            // S of current Nuclear Environment
    G4int envA  =theEnvironment.GetA();            // A of current Nuclear Environment
    if(envA>16&&envA<256)
	{
	  G4int oEn=static_cast<int>(27./sqrt(static_cast<double>(envA)));
	  //G4int oEn=static_cast<int>(21./sqrt(static_cast<double>(envA)));
	  //G4int oEn=static_cast<int>(16./sqrt(static_cast<double>(envA)));
      G4int zEn=oEn+oEn;                           // Charge of the Limit Active Nuclear Environment
      G4int nEn     = zEn;                         // Charge of the Limit Active Nuclear Environment
      bEn  = zEn+nEn;                              // BaryoN of the Limit Active Nuclear Environment
      mbEn = G4QPDGCode(2112).GetNuclMass(zEn,nEn,0); // Mass of the LimActNucEnv
      bEn4M=G4LorentzVector(0.,0.,0.,mbEn);        // 4-momentum of the LimitActiveNuclearEnviron
      bEnQC=G4QContent(bEn+nEn,bEn+zEn,0,0,0,0);   // QuarkContent of the LimitActiveNuclEnv
	}
    else                                           // @@ Can be made STATIC CONSTANTS
	{
      bEn   = 256;                                 // BaryoN of the Limit Active Nuclear Environment
      mbEn  = G4QPDGCode(2112).GetNuclMass(128,128,0); // Mass of the LimActNucEnv
      bEn4M = G4LorentzVector(0.,0.,0.,mbEn);      // 4-momentum of the LimitActiveNuclearEnviron
      bEnQC = G4QContent(384,384,0,0,0,0);         // QuarkContent of the LimitActiveNuclEnv
    }
    G4double envM=theEnvironment.GetMass();        // mass of the current Nuclear Environment
    G4QContent envQC=theEnvironment.GetQCZNS();    // QuarkContent of the current Nuclear Environ
#ifdef debug
	G4cout<<"G4Q::HQ: ePDG="<<envPDG<<",eM="<<envM<<",eQC="<<envQC<<G4endl;
#endif
    G4QContent totQC=valQ+envQC;                   // Total Quark Content
    G4int      totBN=totQC.GetBaryonNumber();      // Total Baryon Number of the Total System
    G4int      totS=totQC.GetStrangeness();        // Total Strangeness of the Total System
    G4int      totZ=totQC.GetCharge();             // Total Charge of the Total System
    G4QNucleus totN(totQC);                        // Pseudo nucleus for the Total System
    G4int      totNeut=totN.GetN();                // Total number of neutrons in the system
    /////////G4int      totProt=totN.GetZ();                // Total number of protons in the system
    G4double totM  =totN.GetMZNS();                // Minimum Mass of the Total System
#ifdef debug
	G4cout<<"G4Q::HQ: tN="<<totN<<",tGSM="<<totM<<",tQC="<<totQC<<G4endl;
#endif
    /////////////////////G4double protCB=theEnvironment.CoulombBarrier(1.,1.);
    G4int    resNPDG=0;
    G4double resNM =10000000.;                     // Prototype of residual mass after n separated
    if(totNeut>0)
	{
      G4QContent resNQC=totQC-G4QContent(2,1,0,0,0,0);
      G4QNucleus resNN(resNQC);
      resNM  = resNN.GetMZNS();
      resNPDG= resNN.GetPDG();
	}
    G4LorentzVector env4M =G4LorentzVector(0.,0.,0.,envM);
    G4LorentzVector tot4M =q4Mom+env4M;
    totMass=tot4M.m();
#ifdef debug
    G4int    totPDG=totN.GetPDG();                 // Total PDG Code for the Current compound
	G4cout<<"G4Q::HQ: tPDG="<<totPDG<<",tM="<<totMass<<",rPDG="<<resNPDG<<G4endl;
#endif
    G4ThreeVector totBoost = tot4M.boostVector();  // Boost vector for Total System (backward)
    G4ThreeVector totRBoost= -totBoost;            // Boost vector for Total System (forward)
    G4int    iniPDG =valQ.GetSPDGCode();
    G4int    iniBN  =valQ.GetBaryonNumber();
    G4int    iniQChg=valQ.GetCharge();
    G4int    iniN   =valQ.GetN();
    G4int    iniP   =valQ.GetP();
    G4int    iniS   =valQ.GetL();
#ifdef debug
	G4cout<<"G4Q::HQ: iniPDG="<<iniPDG<<", Z="<<iniP<<", N="<<iniN<<", S="<<iniS<<G4endl;
#endif
    G4QNucleus iniRN(iniP,iniN-1,iniS);
    G4double iniRM = iniRN.GetMZNS();              // Mass of Residual Quasmon when neutron is rad
    if(iniBN<2||envA>0) iniRM=0.;
    G4double iniQM =G4QPDGCode(iniPDG).GetMass();  // Minimum mass of Quasmon
#ifdef debug
	G4cout<<"G4Q::HQ: iniRN="<<iniRN<<", iniRM="<<iniRM<<", iniQM="<<iniQM<<G4endl;
#endif
    /////////////////G4double iniQM2= iniQM*iniQM;
    G4double bndQM = totM-envM;
    if(envPDG==NUCPDG) bndQM=iniQM;
    G4double excE  = totMass-totM;
    /////////////////G4double freeE = excE*iniQM;
#ifdef debug
    G4double bndQM2= bndQM*bndQM;
    G4double quen  = iniQM+envM;
    G4cout<<"G4Q::HQ:mQ="<<quasM<<valQ<<bndQM2<<",nQ="<<nOfQ<<",Env="<<envM
          <<envQC<<",Q+E="<<quen<<",tM="<<totPDG<<totQC<<totM<<"<"<<totMass<<G4endl;
#endif
    G4int tQ    = valQ.GetTot();                   // Total number of quarks for current Quasmon
    G4int bQ    = abs(valQ.GetBaryonNumber());     // Baryon number of the current Quasmon
    G4QContent cQ = valQ;                          // Temporary copy of Quasmon QC
    G4int   s   = 4;                               // Mesonic
    if (bQ) s   = 3*bQ + 2;                        // Barionic
    if (tQ> s) cQ.DecQAQ((tQ-s)/2);                // Reduce QC to minimum QC
#ifdef debug
    G4int rsPDG = cQ.GetSPDGCode();                // PDG for the lowest residual Quasmon state
    G4cout<<"G4Q::HQ:eN="<<envN<<",eZ="<<envZ<<",Q="<<rsPDG<<cQ<<",piF="<<piF<<",gaF="<<gaF<<G4endl;
#endif
    ////////////////G4int nCL=theEnvironment.UpdateClusters(false); // New A-clusters are calculated
    //theEnvironment.PrepareCandidates(theQCandidates,false,false);// Calc.PrePreprob's of Candid's
    //G4bool fF=piF||gaF;
	// piF,piF or fF,fF is correct. It's made to avoid the specific interaction time cluster choice.
    //theEnvironment.PrepareCandidates(theQCandidates,fF,fF);// Calculate PrePreprobab of Candid's
    theEnvironment.PrepareCandidates(theQCandidates,piF,piF);// Calculate PrePreprobab of Candid's
    ModifyInMatterCandidates();                    // Calculate InMediaMasses of Cand. & Possibil.
    G4double kMom = 0.;                            // Energy of primary qParton in Q-CMS
    G4double minK = 0.;                            // k_min for randomization
    G4double maxK = quasM/2.;                      // k_max for randomization
    G4double kLS  = 0;                             // Energy of primary qParton in LS
    G4double cost = 0.;                            // Cos(theta) of k in QS InRespecTo Q direction
    G4bool   kCond = true;                         // Not adoptable k condition
    G4bool   qCond = true;                         // Not adoptable q condition
    G4bool   fskip=false;                          // Flag to skip when sucked
    G4bool   fred =false;                          // Flag of Environment reduction
    G4bool   ffdc =true;                           // Flag of not successful decay in fragment
    G4LorentzVector k4Mom=zeroLV;                  // 4-momentum prototype for k
    G4LorentzVector cr4Mom=zeroLV;                 // 4-momentum prototype for the ColResQuasmon
    G4int kCount =0;                               // Counter of attempts of k for hadronization
    //G4int kCountMax=27;
    //G4int kCountMax=9;
    G4int kCountMax=3;
    //G4int kCountMax=0;                              //@@ *** Close search for the minimum k ***
    //G4int kCountMax=1;                              //@@ No reson to increase it (@@ change text)
    //G4int qCountMax=27;
    //G4int qCountMax=7;
    //G4int qCountMax=3;
    G4int qCountMax=1; // No reson to increase it (it can changes probability in phi-distribution)
	G4bool gintFlag=false;                          // Flag of gamma interaction with one quark
	while(kCount<kCountMax&&kCond)
	{
      kCond=true;
      G4double miM2=0.;
      if(envPDG==NUCPDG)
	  {
        if(excE>mPi0) miM2=mPi02;
        else          miM2=mP2;
	  }
      else                                          // "Env. exists" case - find k_min & k_max
	  {
        minK=100000.;
        if(piF&&quasM>iniQM) maxK=mNeut/2.+(quasM-iniQM)*iniQM/(iniQM+mNeut);
        else if(quasM>iniQM)
		{
          G4double limK=envM/2.+(quasM-iniQM)*iniQM/(iniQM+envM);
          if(limK<maxK) maxK=limK;
		}
        //@@ for acceleration check it only for P,N,d (& alpfa ?), can be done for all clusters(?)
        //if(iniN>0&&iniBN>1)                    // "Neutron-hadron" estimate
        //{
        //  G4double iK=1000000.;
        //  G4double dqm=quasM+quasM;
        //  if(envA>0)
        //  {
        //    G4QContent rtQC=valQ-neutQC+envQC;      // Total Residual Quark Content
        //    G4QNucleus rtN(rtQC);                   // Create a pseudo-nucleus for the TotalResidual
        //    G4double rtM =rtN.GetMZNS();            // Min Mass of total residual Nucleus
        //    G4double bnRQ=rtM-envM;                 // Bound mass of residual Quasmon
        //    G4double sm2=qM2+m2N-bnRQ*bnRQ;
        //    G4double fqm=dqm+dqm;
        //    G4double aK=sm2/fqm;
#ifdef debug
        //    G4double kts=.135;                      // Test value of k
        //    G4double dkts=kts+kts;
        //    G4double fu=dkts*(dkts*quasM-sm2)+quasM*m2N;
	    //    G4cout<<"G4Q::HQ:MQ="<<quasM<<",rM="<<bnRQ<<",f("<<kts<<")="<<fu<<"<0 ?"<<G4endl;
#endif
        //    if(quasM>=mNeut+bnRQ)
        //    {
        //      G4double srm=sqrt(sm2*sm2-4.*qM2*m2N)/fqm;
        //      iK=aK-srm;
        //      aK+=srm;
        //    }
        //    else iK=aK;
        //    if(aK<maxK) maxK=aK;
        //  }
        //  else iK=mNeut*mNeut/(quasM+quasM);
        //  if(iK<minK) minK=iK;
        //}
        G4int envA=envZ+envN;
        G4int totN=totBN-totZ;
        if(totN>0&&totBN>1)                   // "Neutron-cluster" estimate
        {
		  G4QNucleus tmpNN(totZ,totN-1,0);
          G4double delN=tmpNN.GetMZNS()+mNeut-totM;
          if(envN>0&&envA>1)
          {
			G4QNucleus envNN(envZ,envN-1,0);
            G4double delEN=envNN.GetMZNS()+mNeut-envM;
            if(delEN>delN) delN=delEN;
		  }
          delN*=qurF;
          if(delN<minK) minK=delN;
        }
        if(totZ>0&&totBN>1)                   // "Proton-cluster" estimate
        {
          G4double proCB=theEnvironment.CoulombBarrier(1,1);
		  G4QNucleus tmpPN(totZ-1,totN,0);
          G4double delP=tmpPN.GetMZNS()+mProt-totM+proCB;
          if(envZ>0&&envA>1)
          {
			G4QNucleus envPN(envZ-1,envN,0);
            G4double delEP=envPN.GetMZNS()+mProt-envM+proCB;
            if(delEP>delP) delP=delEP;
		  }
          delP*=qurF;
          if(delP<minK) minK=delP;
        }
        if(totN>0&&totZ>0&&totBN>2)           // "Deuteron-cluster" estimate
        {
          G4double proCB=theEnvironment.CoulombBarrier(1,2);
		  G4QNucleus tmpDN(totZ-1,totN-1,0);
          G4double delD=tmpDN.GetMZNS()+mDeut-totM+proCB;
          if(envN>0&&envZ>0&&envA>2)
          {
		    G4QNucleus envDN(envZ-1,envN-1,0);
            G4double delED=envDN.GetMZNS()+mDeut-envM+proCB;
            if(delED>delD) delD=delED;
		  }
          delD*=qurF;
          if(delD<minK) minK=delD;
        }
        if(totN>1&&totZ>0&&totBN>3)           // "Triton-cluster" estimate
        {
          G4double proCB=theEnvironment.CoulombBarrier(1,3);
		  G4QNucleus tmpTN(totZ-1,totN-2,0);
          G4double delT=tmpTN.GetMZNS()+mTrit-totM+proCB;
          if(envN>1&&envZ>0&&envA>3)
          {
		    G4QNucleus envTN(envZ-1,envN-2,0);
            G4double delET=envTN.GetMZNS()+mTrit-envM+proCB;
            if(delET>delT) delT=delET;
		  }
          delT*=qurF;
          if(delT<minK) minK=delT;
        }
        if(totN>0&&totZ>1&&totBN>3)           // "He3-cluster" estimate
        {
          G4double proCB=theEnvironment.CoulombBarrier(2,3);
		  G4QNucleus tmpRN(totZ-2,totN-1,0);
          G4double delR=tmpRN.GetMZNS()+mHel3-totM+proCB;
          if(envN>0&&envZ>1&&envA>3)
          {
		    G4QNucleus envRN(envZ-2,envN-1,0);
            G4double delER=envRN.GetMZNS()+mHel3-envM+proCB;
            if(delER>delR) delR=delER;
		  }
          delR*=qurF;
          if(delR<minK) minK=delR;
        }
        if(totN>1&&totZ>1&&totBN>4)            // "Alpha-cluster" estimate
        {
          G4double proCB=theEnvironment.CoulombBarrier(2,4);
		  G4QNucleus tmpAN(totZ-2,totN-2,0);
          G4double delA=tmpAN.GetMZNS()+mAlph-totM+proCB;
          if(envN>1&&envZ>1&&envA>4)
          {
		    G4QNucleus envAN(envZ-2,envN-2,0);
            G4double delEA=envAN.GetMZNS()+mAlph-envM+proCB;
            if(delEA>delA) delA=delEA;
		  }
          delA*=qurF;
          if(delA*qurF<minK) minK=delA*qurF;
        }
        //if(totN>2&&totZ>2&&totBN>6)                 // "Li6-cluster" estimate ??????
        //{
        //  G4double proCB=theEnvironment.CoulombBarrier(3,6);
		//  G4QNucleus tmpLN(totZ-3,totN-3,0);
        //  G4double delL=tmpLN.GetMZNS()+mLit6-totM+proCB;
        //  if(envN>2&&envZ>2&&envA>6)
        //  {
		//    G4QNucleus envLN(envZ-3,envN-3,0);
        //    G4double delEL=envLN.GetMZNS()+mLit6-envM+proCB;
        //    if(delEL>delL) delL=delEL;
		//  }
        //  delL*=qurF;
        //  if(delL<minK) minK=delL;
        //}
        if(minK<0.||minK+minK>quasM) minK=0.;
	  }
      //if(addPhoton>0.) nOfQ=valQ.GetQ()-valQ.GetAQ();
      //if(addPhoton>0.&&quasM<1500.&&G4UniformRand()<f2all)
	  //if(addPhoton>0.&&iniBN<2)// PhotonAbsorbDiagramContrib(@@HiddenPar)
	  if(addPhoton>0.&&(iniBN<2||G4UniformRand()<0.5))// PhotonAbsorbDiagramContrib(@@HiddenPar)
	  //if(addPhoton>0.&&(G4UniformRand()<.5||iniBN<2))// PhotonAbsorbDiagramContrib(@@HiddenPar)
	  //if(addPhoton>0.&&G4UniformRand()<0.5)// PhotonAbsorbDiagramContrib(@@HiddenPar)
	  //if(addPhoton>0.&&G4UniformRand()<0.75)// PhotonAbsorbDiagramContrib(@@HiddenPar)
	  //if(addPhoton>0.)// PhotonAbsorbDiagramContrib(@@ A photon is always captured by quark)
	  //if(2>3)                                        // Photon capture by quark is closed
      {
        // ===OR=== 1/(p_g,p_q) interaction
	    gintFlag=true;
        q4Mom-= phot4M;                              // recover Quasmon-Cluster without the Photon
        qM2   = q4Mom.m2();                          // Update the Current squared mass of Quasmon 
        quasM = sqrt(qM2);                           // Update the Current Mass of Quasmon 
        G4double kpow=static_cast<double>(nOfQ-2);
        if(kpow<1.) kpow=1.;
        G4bool frat=true;                            // k/(k+p) factor of QuarkExchange
		while (frat)
        {
          G4double rn=pow((1.-(minK+minK)/quasM),kpow)*G4UniformRand();
          kMom=(1.-pow(rn,1./kpow))*quasM/2.;
		  frat=G4UniformRand()>addPhoton/(addPhoton+kMom);
		}
	    cost=addPhoton/momPhoton;                    // =1 for the real photon
        // Takes into account that the mass of the quark-parton is small but not zero
        //G4double hms=12.5;                            // m_q^2/2 (12.5 means m_q=5 MeV)
        //G4double hms=2.;                            // m_q^2/2 (2. means m_q=2 MeV)
        G4double hms=50.;                            // m_q^2/2 (50. means m_q=10 MeV)
        G4double x=kMom*kMom/hms;
        G4double dc=(pow(x,G4UniformRand())-1.)/x;   // smearing of the "delta-function" for m=0
        cost-=dc;                                    // @@ smear in both sides for virtual photons
#ifdef debug
	    G4cout<<"G4Q::HQ:**PHOTON is SUB'ED from Q** k="<<kMom<<",ct="<<cost<<",QM="<<quasM<<G4endl;
#endif
	  }
      else
	  {
        //gaF=false;                                   // GammaFirstAct flag is only for the gamma-q
	    gintFlag=false;
        // ==== Probabiliti proportional to k (without 1/k factor of hadronization)
        if(!miM2) miM2=(minK+minK)*q4Mom.m();        // Make minimum mass for randomization
        if(qM2<.0001) kMom=0.;
        else kMom = GetQPartonMomentum(maxK,miM2);    // Calculate value of primary qParton
        // ==== Direct calculation of the quark spectrum
        //G4double kpow=static_cast<double>(nOfQ-2);
        //G4double kst=0.;
		//if(maxK+maxK<quasM+.1) kst=pow((1.-(maxK+maxK)/quasM),kpow);
        //if(maxK>minK)
        //{
        //  G4double rn=(pow((1.-(minK+minK)/quasM),kpow)-kst)*G4UniformRand();
        //  kMom=(1.-pow(kst+rn,1./kpow))*quasM/2.;
        //}
        //else kMom=(minK+maxK);
        //^^^ Direct calculation
        G4double rnc=G4UniformRand();
        cost = 1.-rnc-rnc;
	  }
      G4double cQM2=qM2-(kMom+kMom)*quasM;
      G4double cQM=sqrt(cQM2);                       // Mass of the coloured residual Quasmom
      k4Mom=zeroLV;
      cr4Mom=G4LorentzVector(0.,0.,0.,cQM);
      if(!G4QHadron(q4Mom).RelDecayIn2(k4Mom, cr4Mom, q4Mom, cost, cost))// Q->ColResidQ+k
      {
#ifdef debug
        G4cerr<<"*G4Q::HQ:(PhBack)QM="<<quasM<<",cQM="<<cQM<<",ct="<<cost<<",gF="<<gintFlag<<G4endl;
#endif
        kCond=true;
        if(addPhoton&&gintFlag)
        {
          q4Mom+= phot4M;                            // Recover Full Quasmon with Photon, if filed
          qM2   = q4Mom.m2();                        // Update the Current squared mass of Quasmon 
          quasM = sqrt(qM2);                         // Update the Current Mass of Quasmon 
		}
      }
      else
      {
        if(addPhoton&&gintFlag)
        {
          q4Mom+= phot4M;                            // Recover the Full Quasmon with the Photon
          k4Mom+= phot4M;                            // add photon energy to a quark of the Cluster
          qM2   = q4Mom.m2();                        // Update the Current squared mass of Quasmon 
          quasM = sqrt(qM2);                         // Update the Current Mass of Quasmon 
		}
#ifdef debug
	    G4cout<<"G4Q::HQ:(PhBack)k="<<kMom<<",k4M="<<k4Mom<<",ct="<<cost<<",gF="<<gintFlag<<G4endl;
#endif
        kLS=k4Mom.e();
        G4double rEn=cr4Mom.e();
        rMo=cr4Mom.rho();                            // p for the ResidualColouredQuasmon in LS
        rEP=rEn+rMo;                                 // E+p for the ResidualColouredQuasmon in LS
        G4int totCand = theQCandidates.size();    // Total number of candidates
#ifdef sdebug
	    G4cout<<"G4Q::HQ: ***>>>K-ITERATION #"<<kCount<<", k="<<kMom<<k4Mom<<G4endl;
#endif
        for (G4int index=0; index<totCand; index++)
        {
          G4QCandidate* curCand=theQCandidates[index];
          G4int cPDG = curCand->GetPDGCode();
          if(cPDG==90000001||cPDG==90001000||cPDG==91000000||cPDG<80000000)//@@ k-PreAttempts(Acc)
          {
            G4bool poss= curCand->GetPossibility();
#ifdef debug
		    if(cPDG==90000001 || cPDG==90001000 || cPDG==90000002 || cPDG==90001001)
	          G4cout<<"G4Q::HQ:pos="<<poss<<",cPDG="<<cPDG<<",iQC="<<iniQChg<<G4endl;
#endif
		    if(poss)
            {
              G4double   cMs=curCand->GetMass();     // Bound mass of ParentCluster(=Fragment)
              G4QContent cQC=curCand->GetQC();
              G4double   cfM=curCand->GetQPDG().GetMass();
              G4QContent rtQC=curQ+envQC-cQC;        // TotalResidualQuarkContent - OutFragm
              G4QNucleus rtN(rtQC);                  // Create a pseudo-nucleus for residual
              ///in print///G4double totr = rtN.GetMZNS();         // Mass of total Residual
              ///in print///G4double bnM  = totr-envM+cMs;         // Bound mass of residual Quasmon
              //////////////////G4double pmk  = rMo*cMs/kLS;
#ifdef debug
		      if(cPDG==90000001 || cPDG==90001000 || cPDG==90000002 || cPDG==90001001)
			    G4cout<<"G4Q::HQ:cfM="<<cfM<<",cMs="<<cMs<<",ind="<<index<<G4endl;
#endif
              G4double k = kMom;
              if(cPDG>80000000&&cPDG!=NUCPDG) k=kLS; // ===> Nuclear case (Lab System)
              G4double kMin=(cfM*cfM-cMs*cMs)/(cMs+cMs);
              if(kMin<0.) kMin=0.;
#ifdef debug
              G4double totr = rtN.GetMZNS();         // Mass of total Residual
              G4double bnM  = totr-envM+cMs;         // Bound mass of residual Quasmon
              //G4double dR=(bnM*bnM-cQM2)/(rEP+rEP);
              //if(b2<.000001) dR=((bnM*bnM-cQM2)/2.+pmk*(kLS-kMin))/(rEP+pmk);
              //if (dR<0.) dR=0.;
              G4double dR=0.;
     		  if(cPDG==90000001 || cPDG==90001000 || cPDG==90000002 || cPDG==90001001)
	            G4cout<<"G4Q::HQ:i="<<index<<",cPDG="<<cPDG<<",k="<<kMom<<",kLS="<<kLS
                      <<",bM="<<bnM<<",rEP="<<rEP<<",d="<<kMin<<",dR="<<dR<<",kCon="<<kCond<<G4endl;
#endif
		    }
		  }
	    }
	  }
      kCount++;
	}
    // End of search for "k": *** Here we forget about the initial photon forever ***???***
#ifdef debug
    if(addPhoton)
	  G4cout<<"G4Q::HQ:***PHOTON IS ADOPTED*** k="<<k4Mom<<",Q="<<q4Mom<<",PhE="<<addPhoton<<G4endl;
#endif
    addPhoton=0.;                                    // the Photon is adopted: prepare to kill
    momPhoton=0.;
    //phot4M=zeroLV;
    G4double dk = kMom+kMom;                         // Double QCM k-value
#ifdef debug
    G4cout<<"G4Q::HQ:=K's FOUND=>k="<<kMom<<",ki="<<minK<<",ka="<<maxK<<",kLS="<<k4Mom<<kLS<<G4endl;
#endif
    CalculateHadronizationProbabilities(excE,kMom,kLS,piF,gaF); //ResonanceMass is randomized here
    G4int nCandid = theQCandidates.size();
    G4bool fprob = true;                             // Flag of existing decay probability
    ////G4bool fdul  = true;                            // Prototype of flag of resonance decay
    G4bool fdul  = false;                            // Prototype of flag of resonance decay
    int i=0;                                         // "i" will point to the selected candidate
    G4double maxP = theQCandidates[nCandid-1]->GetIntegProbability();
#ifdef debug
    G4cout<<"G4Q::HQ: ***RANDOMIZE A CANDIDATE*** #OfCandidates="<<nCandid<<", maxP="<<maxP<<G4endl;
#endif
    if (maxP<=0.)                                    // No possible channels for this k value
    {
#ifdef debug
	  if(status==4) G4Exception("***G4Q::HadrQ: TEMPORARY EXCEPTION - 0 chanals of fragmentation");
#endif
      ///////G4double qCB=theEnvironment.CoulombBarrier(iniP,iniBN);
      //////G4double pCB=theEnvironment.CoulBarPenProb(qCB,excE,iniP,iniBN);//@@ excE is wrong (tmp)
      //////G4bool qenv=envPDG!=90000000&&envPDG!=90000002&&envPDG!=90002000&&iniPDG!=90002000
      //////              &&iniPDG!=90000002; // @@@@ Should be improved for nn & pp FinalStateInter
      //if(excE>qCB&&envM+iniQM<totMass&&nQuasms==1&&qenv&&G4UniformRand()<pCB) // Ddecay in GSQ+GSE
	  if(2>3)                                        // Decay in GSQuasmon+GSEnviron **Is closed**
	  {
        G4QHadron* resNuc = new G4QHadron(valQ,q4Mom);// Createa Hadron for Quasmon
        resNuc->CorMDecayIn2(iniQM,env4M);           // Recalculate the fr4Mom for ResidualGSMass
        //@@CHECK CoulBar and may be evaporate instead
        FillHadronVector(resNuc);                    // Fill the OutgoingQuasmon (delete equivalent)
#ifdef debug
        G4LorentzVector oQ4Mom=resNuc->Get4Momentum();
        G4cout<<"G4Q::HQ: outQM="<<oQ4Mom.m()<<oQ4Mom<<",GSQM="<<iniQM<<G4endl;
#endif
        if(nQuasms==1) qEnv = G4QNucleus(env4M,envPDG); // Createa Nucleus for the MovingEnvironment
        else
        {
          G4QHadron* envH = new G4QHadron(envPDG,env4M); //@@ Moving Environment !
          FillHadronVector(envH);                    // Fill Moving Environment (delete equivalent)
          qEnv = vacuum;
		}
#ifdef debug
        G4double eM=env4M.m();
        G4LorentzVector dif=tot4M-oQ4Mom-env4M;
        G4cout<<"G4Q::HQ: envM="<<envM<<"=="<<eM<<", envT="<<env4M.e()-eM<<dif<<G4endl;
#endif
        KillQuasmon();                               // This Quasmon is done
        return theQHadrons;                          // The last decay of the quasmon...
	  }
      else                                           // May be evaporation is possible
	  {
        if(CheckGroundState()) KillQuasmon();    // This Quasmon is done
        //if(CheckGroundState(true)) KillQuasmon();    // This Quasmon is done
		//else if(envPDG==NUCPDG && quasM>iniQM && iniBN>1) // Convert Quasmon to Environment
		//{
#ifdef debug
        //  G4cerr<<"***G4Q::HQ:Q2E:E="<<theEnvironment<<", valQ="<<valQ<<", tot4M="<<tot4M<<G4endl;
#endif
        //  theEnvironment=G4QNucleus(valQ,tot4M);
        //  KillQuasmon();
		//}
		else if(envPDG==NUCPDG && quasM>iniQM)       // Emergency decay in Gamma & GSQuasmon
		{
#ifdef debug
          G4cout<<"***G4Q::HQ: Emergency Decay in Gamma/Pi0 + Residual GSQuasmon"<<G4endl;
		  //G4Exception("***G4Q::HadrQ: Emergency Decay in Gamma + Residual GSQuasmon");
#endif
          G4int gamPDG=22;
          G4double gamM=0.;
          if(quasM>mPi0+iniQM)
		  {
            gamPDG=111;
            gamM=mPi0;
		  }
          G4LorentzVector r4Mom(0.,0.,0.,gamM);      // mass of the photon/Pi0
          G4LorentzVector s4Mom(0.,0.,0.,iniQM);     // mass of the GSQuasmon
          if(!G4QHadron(q4Mom).DecayIn2(r4Mom, s4Mom))
          {
            G4cerr<<"**G4Q::HQ:Q="<<q4Mom<<quasM<<"->g/pi0+GSQ="<<iniPDG<<"(M="<<iniQM<<")"<<G4endl;
    	    G4Exception("G4Quasmon::HadronizeQuasmon:Gam/Pi0+GSQuasmon decay didn't succeed,Env=0");
          }
#ifdef debug
	      G4cout<<"G4Q::HQ:=== 0 ===> HadrVec, Q="<<q4Mom<<quasM<<"->g/pi0("<<gamPDG<<")="<<r4Mom
                <<gamM<<"+GSQ="<<iniPDG<<r4Mom<<iniQM<<G4endl;
#endif
          G4QHadron* curHadr2 = new G4QHadron(gamPDG,r4Mom); // Creation Hadron for Gamma
          FillHadronVector(curHadr2);                // Fill "new curHadr2" (delete equivalent)
          G4QHadron* curHadr1 = new G4QHadron(iniPDG,s4Mom);// Create Hadron for Residual GSQ
          FillHadronVector(curHadr1);                // Fill "new curHadr1" (delete equivalent)
          KillQuasmon();                             // This Quasmon is done
		}
#ifdef pdebug
	    G4cout<<"***G4Q::HQ:dE="<<excE<<">minK="<<minK<<",Env="<<theEnvironment<<",k="<<kLS<<",Q="
              <<valQ<<quasM<<G4endl;
        //if(excE>8.)G4Exception("G4Quas::HQ: Why a candidate is not found?");
        //if(totMass>envM+iniQM&&excE>minK&&excE>1.&&minK)G4Exception("G4Q::HQ:Why it'sn't found?");
        //G4Exception("G4Quasmon::HadronizeQuasmon: Why a candidate is not found?");
#endif
        qEnv=theEnvironment;
        return theQHadrons;
	  }
    }
    else
	{
      G4double totP = maxP * G4UniformRand();
      while(theQCandidates[i]->GetIntegProbability() < totP) i++;
	}
#ifdef debug
    G4cout<<"G4Q::HQ: fp="<<fprob<<",QM="<<quasM<<",QQC="<<valQ<<",k="<<kMom<<G4endl;
#endif
    if (i>=nCandid)
    {
      G4cerr<<"***G4Q::HQ: Cand#"<<i<<" >= Tot#"<<nCandid<<G4endl;
      G4Exception("G4Quasmon::HadronizeQuasmon: Too big number of the candidate");
	}
    G4bool nucflag=false;                            // NuclearBinding Flag (can't reduce Environ)
    G4bool hsflag=false;                             // Prototype of H+S decay flag
    G4bool fchipo=false;                             // Final decay of Quasmon-Chipolino
    curQ = valQ;                                     // Temporary copy of valQ to estimate MinM2
    G4int rPDG=0;                                    // Prototype of the residual Quasmon PDG
    G4double rMass=0.;                               // Prototype of the residual Quasmon mass
    G4LorentzVector rQ4Mom=zeroLV;                   // 4-momentum of residual Quasmon
    G4LorentzVector fr4Mom=zeroLV;                   // 4-momentum prototype for the fragment
    if (!fprob)                                      // Prepare final decay, "nothing was selected"
	{
      memQ=curQ;                                     // Remembe QC for the case of 2H decay
      sPDG=curQ.GetSPDGCode();                       // PDG of current Quasmon as S-hadron      
      if(!sPDG&&theEnvironment.GetPDG()!=NUCPDG&&totBN>1&&totMass>totM&&totS>=0)
      {
#ifdef ppdebug
        G4cout<<"G4Quasm::HadrQuasm:NEEDS-EVAP-1:Q="<<q4Mom<<valQ<<",QEnv="<<theEnvironment<<G4endl;
        G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (1)"); //@@ TMP
#endif
        qEnv=theEnvironment;
        return theQHadrons;
      }
      else if(!sPDG)
	  {
        G4cerr<<"***G4Q::HQ:sPDG=0,QE="<<theEnvironment<<",tB="<<totBN<<",tS="<<totS<<G4endl;
        G4Exception("G4Quasmon::HadronizeQuasmon: Decay particle isn't selected & can't evaporate");
	  }
      else if(sPDG==10)                              // ---> Chipolino case
	  {
        fchipo=true;
        G4QChipolino chipQ(valQ);
        G4QPDGCode QPDG1=chipQ.GetQPDG1();
        sPDG = QPDG1.GetPDGCode();
        sMass= QPDG1.GetMass();
        G4QPDGCode QPDG2=chipQ.GetQPDG2();
        rPDG = QPDG2.GetPDGCode();
        rMass= QPDG2.GetMass();
        if(sMass+rMass>quasM)                        //@@ Evaporate but try 3Ptcl decay with Environ
		{
		  if(totBN>1&&totMass>totM&&totS>=0)
		  {
#ifdef ppdebug
            G4cout<<"G4Quasm::HadrQ:NEEDS-EVAP-2:Q4M="<<q4Mom<<valQ<<",QE="<<theEnvironment<<G4endl;
            G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (2)"); //@@ TMP
#endif
            qEnv=theEnvironment;
            return theQHadrons;
		  }
          else
		  {
            G4cerr<<"***G4Q::HQ:QM="<<quasM<<" < M1="<<sMass<<" + M2="<<rMass<<" = "<<sMass+rMass
                  <<",tB="<<totBN<<",tS="<<totS<<", tM="<<totMass<<" > miM="<<totM<<G4endl;
            G4Exception("G4Quasmon::HadronizeQuasmon: Virtual Chipolino can't evaporate Nucleus");
		  }
		}
      }
      else                                           // --->@@ Final decay in MinHadr+(PI0 or GAMMA)
      {
        sMass=G4QPDGCode(sPDG).GetMass();
        rPDG=envPDG;
        if (rPDG>80000000&&rPDG!=NUCPDG)
	    {
          rMass=theEnvironment.GetMZNS();
          q4Mom+=G4LorentzVector(0.,0.,0.,rMass);
          valQ +=theEnvironment.GetQC();
          quasM=q4Mom.m();
          KillEnvironment();
          fred=true;                                 // Flag of environment reduction
	    }
        else if(rPDG!=NUCPDG&&totBN>1&&totMass>totM&&totS>=0)  // Evaporation
	    {
#ifdef ppdebug
          G4QContent nTotQC=totQC-neutQC;
          G4QNucleus nTotN(nTotQC);                  // PseudoNucleus for TotalResidual to neutron
          G4double nTotM =nTotN.GetMZNS();           // MinimumMass of the Total Residual to neutron
          G4double dMnT=totMass-nTotM-mNeut;         // Energy excess for the neutron
          G4QContent pTotQC=totQC-protQC;
          G4QNucleus pTotN(pTotQC);                  // PseudoNucleus for TotalResidual to proton
          G4double pTotM =pTotN.GetMZNS();           // Minimum Mass of the Total Residual to proton
          G4double dMpT=totMass-pTotM-mProt;         // Energy excess for neutron
          if(dMpT>dMnT)dMnT=dMpT;
          G4cout<<"G4Q::HQ:NEEDS-EVAP-3: sP="<<sPDG<<",Q4M="<<q4Mom<<valQ<<",rPDG="<<rPDG<<G4endl;
          G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (3)"); //@@ TMP
#endif
          qEnv=theEnvironment;
          return theQHadrons;
	    }
        else if(rPDG==NUCPDG)
        {
#ifdef debug
          G4cout<<"G4Quasmon::HadronizeQuasmon: SAFETY DECAY IN PI0/GAMMA, rPDG="<<rPDG<<G4endl;
#endif
          if(totMass-totM>mPi0)
          {
            rMass=mPi0;                              // Safety decay in pi0
            rPDG=111;
		  }
          else
          {
            rMass=0.;                                // Safety decay in gamma
            rPDG=22;
		  }
        }
      }
      hsflag=true;                                   // Two particle decay is forced ???
#ifdef debug
      G4cout<<"G4Q::HQ:*hsflag!* sPDG="<<sPDG<<",sM="<<sMass<<",rPDG="<<rPDG<<",rM="<<rMass<<G4endl;
#endif
	}
    else                                             // ---> "Something was selected" case
	{
      G4QCandidate* curCand = theQCandidates[i];     // Pointer to selected candidate (hadr./fragm.)
      sPDG  = curCand->GetPDGCode();                 // PDG of the selected candidate 
      //////////////////G4double prpr=curCand->GetPreProbability();
      //@@ For clusters should be another randomization & a loop over possible parent clusters
      if(sPDG>80000000&&sPDG!=NUCPDG)                // ===> "Fragment" case
	  {
        int ip=0;
        G4int nParCandid = curCand->GetPClustEntries();
        G4double sppm  = curCand->TakeParClust(nParCandid-1)->GetProbability();
        if (sppm<=0)                                 // Impossible to find a parent cluster
        {
          G4cerr<<"***G4Quasmon::HadronizeQuasmon: cP="<<theQCandidates[i]->GetIntegProbability()
                <<",nC="<<nParCandid<<",pP="<<sppm<<",QM="<<quasM<<",QC="<<valQ;
          for(int ipp=0; ipp<nParCandid; ipp++)
            G4cerr<<", "<<ipp<<": "<<curCand->TakeParClust(ip)->GetProbability();
          G4cerr<<G4endl;
          G4Exception("G4Quasmon::HadronizeQuasmon: No parent cluster for the fragment");
	    }
        else                                         // ---> "Parent cluster was found" case
	    {
          G4double totPP = sppm * G4UniformRand();
          while(curCand->TakeParClust(ip)->GetProbability() < totPP) ip++;
#ifdef debug
          G4cout<<"G4Q::HQ:p#ip="<<ip<<",f#i="<<i<<",tP="<<totPP<<",sP="<<sppm<<G4endl;
#endif
	    }
        G4QParentCluster* parCluster=curCand->TakeParClust(ip);
        pPDG  = parCluster->GetPDGCode();
	    G4QPDGCode pQPDG(pPDG);
        pQC   = pQPDG.GetQuarkContent();
        pBaryn= pQC.GetBaryonNumber();
	    pMass = parCluster->GetEBMass();             //Environment Bounded Mass
	    pNMass = parCluster->GetNBMass();            // NucleusBoundedMass @@ Env.Reduce dosn't work
        pM2   = pMass*pMass;
	    transQC = parCluster->GetTransQC();
	    delta = parCluster->GetEBind();              // Environmental Binding
	    deltaN = parCluster->GetNBind();             // Nuclear Binding
	    loli  = parCluster->GetLow();
	    hili  = parCluster->GetHigh();
        //G4double dhil=.0001*(hili-loli);           // Safety factor
        //loli += dhil;
        //hili -= dhil;
	    npqp2 = parCluster->GetNQPart2();
        // @@ Here one can get other useful parameters of the parent cluster for hadronization
	    G4QPDGCode sQPDG(curCand->GetPDGCode());
        sQC   = sQPDG.GetQuarkContent();
        if(sPDG==90001001 && G4UniformRand()>0.75) sMass=np; //@@ n-p pair
        //if(sPDG==90001001 && G4UniformRand()>1.0) sMass=np; //@@ no n-p pair
        else                                      sMass = sQPDG.GetMass();
        sM2   = sMass*sMass;                         // Squared mass of the fragment
        curQ += transQC;                             // Subtract ExchangeMesonQC from QC of Quasmon
#ifdef debug
        G4cout<<"G4Q::HQ:valQ="<<valQ<<"+transQ="<<transQC<<"("<<pPDG<<" to "<<sPDG
              <<") = curQ="<<curQ<<",InvBinding="<<delta<<",pM="<<pMass<<pQC<<G4endl;
#endif
	  }
      else                                           // ===> "Hadron" case
      {
        pBaryn=1;
        sMass = theQCandidates[i]->GetMass();        // Mass is randomized on probability level
        sM2=theQCandidates[i]->GetMass2();           // Sq. Mass is randomized on probability level
        curQ-= theQCandidates[i]->GetQC();           // Subtract outHadron QC from QC of Quasmon
#ifdef debug
        G4cout<<"G4Q::HQ:valQ="<<valQ<<" - sQ="<<theQCandidates[i]->GetQC()<<",sM="<<sMass<<"(PDG="
              <<theQCandidates[i]->GetPDGCode()<<") = curQ="<<curQ<<", sM2="<<sM2<<G4endl;
#endif
	  }
      memQ=curQ;                                     // Remembe QC for the case of DecayIn2 Hadrons
      G4double rtM=0.;
      G4double reM=0.;
      // ---> Calculate the Squared Ground State Mass of the ResidualQuasmon+ResidualEnvironment
      if(envPDG>pPDG)
      {
        // *** LIM ***
		G4QContent RNQC=curQ+envQC;
        if(sPDG>80000000&&sPDG!=NUCPDG)              // ==> "Nuclear Fragment Radiation" case
        {
          if(envA-pBaryn<=bEn) RNQC-=pQC;
          else RNQC-=envQC-bEnQC;
		}
        G4int RNPDG = RNQC.GetSPDGCode();
        if(RNPDG==10) minSqN=G4QChipolino(RNQC).GetMass2();// MinSqMass of DoubleHadron of Chipolino
        else if(!RNPDG)
	    {
#ifdef debug
          G4cout<<"***G4Q::HQ: *** ResidualNucleus *** PDG="<<RNPDG<<curQ<<",QC="<<RNQC<<G4endl;
#endif
          minSqN=1000000.;
	    }
        else
		{
          G4double minT=G4QPDGCode(RNPDG).GetMass();
          minSqN=minT*minT;                          // SquaredGSMass of ResidQuasmon+ResidEnviron
        }
	  }
      // ---> Calculate the Minimum Squared Mass of the Residual Quasmon
      G4int rqPDG = curQ.GetSPDGCode();
      G4int rQQ=G4QPDGCode(curQ).GetQCode();
      if(rqPDG==10)
      {
        minSqT=G4QChipolino(curQ).GetMass2();        // MinSqMass of DoubleHadron of Chipolino
        minSqB=minSqT;
	  }
      else if(!rqPDG||rQQ<-1)
	  {
#ifdef debug
        G4cerr<<"***G4Q::HQ: *** ResidualQuasmon *** PDG="<<rqPDG<<curQ<<",Q="<<rQQ<<G4endl;
#endif
        minSqT=1000000.;
        minSqB=1000000.;
	  }
      else
      {
        //////////G4int baryn=curQ.GetBaryonNumber();
        G4double minT=G4QPDGCode(rqPDG).GetMass();
        minSqT=minT*minT;                            // Squared minimum mass of residual quasmon
        if((sPDG<80000000&&envPDG>80000000&&envPDG!=NUCPDG||sPDG>80000000&&sPDG!=NUCPDG&&
          envPDG>pPDG) && (rqPDG>80000000&&rqPDG!=NUCPDG||rqPDG==2112||rqPDG==2212||rqPDG==3122))
        {
          G4double newT=0.;                          // Prototype of minimal Quasmon mass
          if(sPDG<80000000)
		  {
            // *** LIM ***
            G4QContent rtQC=curQ;                    // Total Residual Quark Content
            if(envA<=bEn) rtQC+=envQC;
            else rtQC+=bEnQC;
            G4QNucleus rtN(rtQC);                    // Create a pseudo-nucleus for residual
            rtM=rtN.GetMZNS();                       // MinMass of Total Residual Nucleus (hadron)
            if(envA<=bEn) newT=rtM-envM;
            else newT=rtM-mbEn;
#ifdef debug
		    G4cout<<"G4Q::HQ: ***VacuumFragmentation** M="<<newT<<",rM="<<rtM<<rtQC<<",eM="<<envM
                  <<",mM="<<minT<<G4endl;
#endif
		  }
          else
		  {
            G4QContent reQC=envQC-pQC;               // Total Residual Quark Content
            G4QNucleus reN(reQC);                    // Create a pseudo-nucleus for residual
            reM=reN.GetMZNS();                       // Min Mass of residual Environment
            // *** LIM ***
            G4QContent rtQC=curQ;                    // Total Residual Quark Content
            if(envA-pBaryn<=bEn) rtQC+=reQC;
            else rtQC+=bEnQC;
            G4QNucleus rtN(rtQC);                    // Create a pseudo-nucleus for residual
            rtM=rtN.GetMZNS();                       // MinMass of Total Residual Nucleus (fragment)
            // *** LIM ***
            if(envA-pBaryn<=bEn) newT=rtM-envM+pMass;
            else newT=rtM-mbEn;
            
#ifdef debug
		    G4cout<<"G4Q::HQ: *** NuclearFragmentation*** M="<<newT<<",rM="<<rtM<<rtQC<<",eM="<<envM
                  <<",pM="<<pMass<<G4endl;
#endif
		  }
          if(newT<minT) minT=newT;
	    }
        minSqB=minT*minT;                            // MinSquaredMass of BoundedResidualQuasmon
      }
#ifdef debug
	  G4cout<<"G4Q::HQ:rqPDG="<<rqPDG<<",mST="<<minSqT<<",mSB="<<minSqB<<",rtM="<<rtM<<G4endl;
#endif
      if(!minSqT)
      {
        G4cerr<<"***G4Quasmon::HadronizeQuasmon: minSqT=0(!), curQC="<<curQ<<G4endl;
        G4Exception("G4Quasmon::HadronizeQuasmon: Minimal Residual Mass can not be calculated");
      }
      G4double m2 = BIG2; //@@ just a BIG number     // Prototype/Squared Mass of Residual Quasmon
      G4double kt=0.;                                // Prototype of SquaredRealResidualQuasmon Mass
      G4double kp=0.;                                // 3-Mpm/Mass for the residual Quasmon
      G4double kn=0.;                                // Prototype of SquaredRealResQuasm+ResEnv Mass
      if(sPDG>80000000&&sPDG!=NUCPDG)                // ==> Nuclear fragment candidate hadronization
	  {
#ifdef debug 
	    G4cout<<"G4Q::HQ: BoundM="<<pMass<<",FreeM="<<sMass<<",QM="<<quasM<<G4endl;
#endif
        // = = = =  P u r e   k i n e m a t i c a l   c a l c u l a t i o n s:  = = = = =
        // Fusion of k + parentCluster => colouredCluster (cc)
        G4LorentzVector cl4Mom(0.,0.,0.,pMass);      // 4-momentum prototype for parent cluster
        G4LorentzVector tot4Mom=q4Mom+cl4Mom;        // @@ Just for checking
#ifdef debug
	    G4cout<<"G4Q::HQ:Q("<<quasM<<")->k("<<k4Mom<<")+CRQ("<<cr4Mom.m()<<")"<<G4endl;
#endif
        G4LorentzVector cc4Mom=k4Mom+cl4Mom;         // 4-mom of ColoredFragment (before kappa)
        G4double ccM2=cc4Mom.m2();                   // SquaredMass of the ColoredFragment
        G4double frM2=sMass*sMass;                   // MinSuaredMass of the OutgoingFragment
        if(ccM2<=frM2)
        {
#ifdef debug
	      G4cout<<"***G4Q::HQ: Failed tofind FragmentMass:"<<ccM2<<"<"<<frM2<<", pM="<<pMass
                <<" + k="<<k4Mom<<" = "<<sqrt(ccM2)<<cc4Mom<<" < fm="<<sMass<<",miK="<<minK<<G4endl;
#endif
          pMass=pNMass;
          delta=deltaN;
          cl4Mom=G4LorentzVector(0.,0.,0.,pMass);    // 4-momentum prototype for parent cluster
          tot4Mom=q4Mom+cl4Mom;                      // @@ Just for checking
          cc4Mom=k4Mom+cl4Mom;                       // 4-mom of ColoredFragment (before kappa)
          ccM2=cc4Mom.m2();
          if(ccM2<=frM2)
          {
#ifdef debug
	        G4cout<<"G4Q::HQ: ***NUCLEAR BINDING*** ccM2="<<ccM2<<" < frM2="<<frM2<<G4endl;
#endif
            nucflag=true;  // Nuclear binding was used
		  }
          else
          {
#ifdef debug
	        G4cout<<"G4Q::HQ:****hsflag!**** NUCLEAR BINDING ccM2="<<ccM2<<" > frM2="<<frM2<<G4endl;
#endif
            hsflag=true;
          }
		}
        else
        {
#ifdef debug
          G4double crMass2 = cr4Mom.m2();            // SquredMass of ColouredResidualQuasmon
	      G4cout<<"G4Q::HQ:cM2="<<crMass2<<"="<<rEP*(rEP-rMo-rMo)<<",h="<<hili<<",l="<<loli<<G4endl;
#endif
          //////////G4double newh=0.;
          if(hili<loli) hili=loli;
          G4double pw= 1./static_cast<double>(npqp2);
          // ----------->>>>>Decay of the ColouredCluster in a Fragment + kappa (fixed ctc)<<<
          qCond=true;
          G4int qCount=0;
          G4double dM=0.;
          if(sPDG==90001001) dM=2.25;                // Binding energy of the deuteron
          if(excE<sMass-pMass-dM)
          //if(excE<sMass-pNMass-dM)
		  {
#ifdef ppdebug
            G4cout<<"***G4Quasmon::HadronizeQuasmon: ***PANIC#1*** TotalDE="<<excE<<"<bindE="
                  <<sMass-pMass-dM<<", dM="<<dM<<", sM="<<sMass<<", bM="<<pMass<<G4endl;
            //G4Exception("G4Quasmon::HadronizeQuasmon: Why PANIC? (1)"); //@@ TMP
#endif
            status =-1;                              // Panic exit
            qEnv=theEnvironment;                     // Update the QEnvironment
            return theQHadrons;
          }
	      G4double   ex=kLS-delta;                   // EnvironmentExcess of parton energy in LS
          G4QContent tmpEQ=envQC-pQC;                // Quark content for Residual Environment
          G4QNucleus tmpN(tmpEQ);                    // Pseudo nucleus for Residual Environment
          G4double   tmpNM=tmpN.GetMZNS();           // Mass of Residual Environment
#ifdef debug
		  G4cout<<"G4Q::HQ:envQC="<<envQC<<",pQC="<<pQC<<",tmpNM="<<tmpNM<<G4endl;
#endif
          G4QContent tmpRQ=valQ+transQC;             // QContent of Residual Quasmon
          G4QContent tmpTQ=tmpEQ+tmpRQ;              // QC for TotalResidualNucleus for a Fragment
          G4QNucleus tmpT(tmpTQ);                    // Nucleus for TotalResidNucleus for Fragment
          G4double   tmpTM=tmpT.GetMZNS();           // GSMass of TotalResidualNucleus for Fragment
          tmpTM2=tmpTM*tmpTM;                        // SquaredGSMass of TotResNuc for Fragment (g)
          G4LorentzVector ResEnv4Mom(0.,0.,0.,tmpNM);
          G4LorentzVector tCRN=cr4Mom+ResEnv4Mom;
#ifdef debug
          G4double tcEP=tCRN.e()+tCRN.rho();
          // Test that the edge of the Phase space is reachable
          G4double cta=1.-((ex+ex)/(1.-pow(loli,pw))-pMass)/kLS;// cos(theta_k,kappa_max) in LS
          //G4double kap_a=ex/(1.-kLS*cta/pMass);      // Energy of the recoil quark (q_max)
          G4double kap_a=ex/(1.+kLS*(1.-cta)/pMass); // Energy of the recoil quark (q_max)
          G4double cti=1.-((ex+ex)/(1.-pow(hili,pw))-pMass)/kLS;// cos(theta_k,kappa_max) in LS
          //G4double kap_i=ex/(1.-kLS*cti/pMass);      // Energy of the recoil quark (q_min)
          G4double kap_i=ex/(1.+kLS*(1.-cti)/pMass); // Energy of the recoil quark (q_min)
          G4double q_lim=(tmpTM2-tCRN.m2())/(tcEP+tcEP);
          G4cout<<"G4Q::HQ:qi="<<kap_i<<",ci="<<cti<<",qa="<<kap_a<<",ca="<<cta<<",q="<<q_lim
                <<",M2="<<tmpTM*tmpTM<<",R2="<<tCRN.m2()<<",EP="<<tcEP<<G4endl;
#endif
          G4LorentzVector MEMkp4M(0.,0.,0.,0.);
          G4LorentzVector MEMfr4M(0.,0.,0.,0.);
          G4LorentzVector MEMrQ4M(0.,0.,0.,0.);
          G4double MEMrQM2=0.;
          G4double MEMreM2=0.;
          while(qCount<qCountMax&&qCond)
          {
            G4double z = pow(loli+(hili-loli)*G4UniformRand(),pw); // *** q-RANDOMIZATION ***
            //**PSLtest**//G4double z = pow(hili-(hili-loli)*.01,pw); // *** q-RANDOMIZATION ***@@@@
            G4double ctkk=1.-((ex+ex)/(1.-z)-pMass)/kLS;// cos(theta_k,kappa) in LS
#ifdef ppdebug
	        if(qCount)
              G4cout<<"G4Q::HQ:"<<qCount<<",ct="<<ctkk<<",M="<<pMass<<",z="<<z<<",zl="<<pow(loli,pw)
                    <<",zh="<<pow(hili,pw)<<",dE="<<totMass-totM<<",bE="<<sMass-pMass<<G4endl;
#endif
#ifdef debug
            G4cout<<"G4Q::HQ:ct="<<ctkk<<",pM="<<pMass<<",z="<<z<<",zl="<<pow(loli,pw)
                  <<",zh="<<pow(hili,pw)<<",ex="<<ex<<",li="<<loli<<",hi="<<hili<<G4endl;
#endif
            if(abs(ctkk)>1.)
            {
#ifdef debug
              G4cout<<"***G4Q:HQ:ctkk="<<ctkk<<",ex="<<ex<<",z="<<z<<",pM="<<pMass
                    <<",kLS="<<kLS<<",hi="<<hili<<",lo="<<loli<<",n="<<npqp2<<G4endl;
#endif
              if(ctkk> 1.)ctkk= 1.;
              if(ctkk<-1.)ctkk=-1.;
            }
            G4double cen=kLS+pMass;                  // LS Energy of k+parentCluster CompSystem
            G4double ctc=(cen*ctkk-kLS)/(cen-kLS*ctkk);//cos(theta_k,kap) in k+pClast CompSystem
            if(abs(ctc)>1.)
            {
              G4cout<<"***G4Quasm:HadronizeQuasmon: e="<<cen<<", k="<<kLS<<", cost="<<ctc<<G4endl;
              if(ctc>1.) ctc=1.;
              else if(ctc<-1.) ctc=-1.;
			}
            kp4Mom=zeroLV;                           // 4-mom prototype update for RecoilQuark (q)
            fr4Mom=G4LorentzVector(0.,0.,0.,sMass);  // 4-momentum prototype for the fragment
            if(!G4QHadron(cc4Mom).RelDecayIn2(kp4Mom, fr4Mom, k4Mom, ctc, ctc))
            {
              G4cerr<<"*G4Quasmon::HadronizeQuasm:c4M="<<cc4Mom<<",sM="<<sMass<<",ct="<<ctc<<G4endl;
              G4Exception("G4Quasmon::HadronizeQuasmon: Can't decay FakeColClust in fr+ kappa");
            }
#ifdef debug
            G4double ccM=sqrt(ccM2);
            G4double kappa=ex/(1.+kLS*(1.-ctkk)/pMass);// Energy of the RecoilQuark in LS
	        G4cout<<"G4Q::HQ:>ColDec>>>CF("<<ccM<<")->F("<<sMass<<")+q"<<kp4Mom<<"="<<kappa<<G4endl;
#endif
            fmh=true;
            // Fusion of the ColouredResidualQuasmon + kappa (LS) -> get residual Quasmon mass
            rQ4Mom=cr4Mom+kp4Mom;                    // 4-momentum of residual Quasmon
            /////G4double rQM2=rQ4Mom.m2();               // @@ Temporary before the "kt" is defined
            G4LorentzVector retN4Mom=rQ4Mom+ResEnv4Mom;
            reTNM2=retN4Mom.m2();                    // Real mass of TotalResidNucleus for Fragment
            //@@ === One can add here anothe cuts for q randomization ===
            //if(rQM2>=minSqT)                         // tmpTM2 = SqrdGSMass of TotResidualNucleus
            // ***VTN***
            if(reTNM2>=tmpTM2)                       // tmpTM2 = SqrdGSMass of TotResidualNucleus
            // ***VRQ***
			//if(rQM2>=minSqB) // tmpTM2 = SqrdGSMass of TotResidualNucleus
            {
               qCond=false;
               ffdc=false;
#ifdef debug
			   //G4cout<<"G4Q::HQ:Attempt#"<<qCount<<" succeeded.rM2="<<rQM2<<">"<<minSqT<<G4endl;
               // ***VTN***
			   G4cout<<"G4Q::HQ:Attempt#"<<qCount<<" succeeded.rM2="<<reTNM2<<">"<<tmpTM2<<G4endl;
               // ***VRQ***
			   //G4cout<<"G4Q::HQ:Attempt#"<<qCount<<" succeeded.rM2="<<rQM2<<">"<<minSqB<<G4endl;
#endif
			}
			else
            {
#ifdef debug
			  //G4cout<<"G4Q::HQ:Attempt#"<<qCount<<" failed. rM2="<<rQM2<<"<"<<minSqT<<G4endl;
              // ***VTN***
			  G4cout<<"G4Q::HQ:Attempt#"<<qCount<<" failed. rM2="<<reTNM2<<"<"<<tmpTM2<<G4endl;
              // ***VRQ***
              //G4cout<<"G4Q::HQ:Attempt#"<<qCount<<" failed. rM2="<<rQM2<<"<"<<minSqB<<G4endl;
#endif
              //if(rQM2>MEMrQM2)
              // ***VTN***
              if(reTNM2>MEMrQM2)
			  {
                //MEMrQM2=rQM2;
                // ***VTN***
                MEMrQM2=reTNM2;
                //------------
                MEMkp4M=kp4Mom;
                MEMfr4M=fr4Mom;
                MEMrQ4M=rQ4Mom;
                MEMreM2=reTNM2;
              }
              else
			  {
                //rQM2=MEMrQM2;
                // ***VTN***
                reTNM2=MEMrQM2;
                //-----------
                kp4Mom=MEMkp4M;
                fr4Mom=MEMfr4M;
                rQ4Mom=MEMrQ4M;
                reTNM2=MEMreM2;
			  }
			}
            qCount++;
		  }
          quexf=true;                                // Quark Exchange is successfully done
#ifdef debug
	      G4cout<<"G4Q::HQ:RQ("<<rQ4Mom.m()<<")=CRQ("<<cr4Mom.m()<<")+q"<<kp4Mom<<G4endl;
#endif
          kt=rQ4Mom.m2();
          kp=rQ4Mom.rho()/sqrt(kt);
          // *** LIM ***
          G4LorentzVector totC4Mom=rQ4Mom; // TotatResidualNucleus prototype
          if(envA-pBaryn<=bEn) totC4Mom+=G4LorentzVector(0.,0.,0.,envM-pMass);
          else totC4Mom+=bEn4M;
          kn=totC4Mom.m2();
          m2=kt;
          tot4Mom-=rQ4Mom+fr4Mom;
#ifdef debug
	      G4cout<<"G4Q::HQ: t4M="<<tot4Mom<<".Is "<<kt<<">"<<minSqT<<"?"<<",m2="<<m2<<G4endl;
#endif
          //if(kt<minSqT&&ffdc)                        // Max enreachment kinematic limit (VirtualQ)
		  //if(kt<minSqT)
		  //if(kt<minSqB)
		  if(kn<minSqN)
		  //if(kn<minSqN&&kt<minSqT)
          // ***VTN***
		  //if(kt<minSqT&&sPDG!=90000001||kn<minSqN&&sPDG==90000001) // NuclEnreach for OnlyNeutrons
		  //if(kt<minSqT&&sMass>1200.||kn<minSqN&&sMass<1200.) // NuclearEnreach for OnlyNucleons
          // ***VRQ*** Close to get just GS value
		  //if(kt<minSqT&&sMass>1200.||kt<minSqB&&sMass<1200.) // BindingEnreachment for OnlyNucleon
          // @@@@ Priveleged Nuteron !! ??? @@@@@
		  //if(kn<minSqN&&sPDG==90000001||kt<minSqB&&sPDG!=90000001) // BindEnreach for OnlyNuteron
		  //if(kt<minSqB) // BindingEnreachment for OnlyNucleon
          {
            hsflag=true;
#ifdef debug
	        G4cout<<"G4Q::HQ:***hsflag***No t="<<kt<<"<"<<minSqT<<" or n="<<kn<<"<"<<minSqN<<G4endl;
#endif
		  }
		}
	  }
      else                                             // ==> "HadronicCandidate hadronization" case
      {
        kt = (quasM-dk)*(quasM-sM2/dk);
        //@@@@@@@@ rQ4Mom (4-momentum of residual Quasmon) could be defined here?! How?
        if(envPDG>80000000&&envPDG!=NUCPDG)
        {
          // *** LIM ***
          G4LorentzVector TCRN=q4Mom-k4Mom;
          if(envA<=bEn) TCRN+=env4M;
          else TCRN+=bEn4M;
          kn=TCRN.m2();
		}
        else kn=kt;
      }
      // *** LIM ***
      G4LorentzVector tL=rQ4Mom;                     // @@ Is rQ4Mom calculated for hadrons??
      if(envA-pBaryn<=bEn) tL+=G4LorentzVector(0.,0.,0.,reM);
      else tL+=bEn4M;
      G4double tM=tL.m();                            // Real Residual Total Nucleus Mass (hadr/frag)
      // Residual S+S limit (absolute low limit for corresponding P-resonance) for R->S+S Decay
#ifdef debug
      G4cout<<"G4Q::HQ:k="<<kMom<<".Is r2="<<kt<<">T2="<<minSqT<<" or "<<minSqB<<" & tM="<<tM
            <<">rtM="<<rtM<<" to avoid H+S("<<sPDG<<"),QM="<<quasM<<",dk="<<dk<<G4endl;
#endif
      //@@@???@@@ tM is not refined fpr hadrons as rQ4Mom is not defined (try to refine!) @@@???@@@
      // In following minSqT is for not bounded recoil Quasmon, minSqB is for the bounded Quasmon
      if(kn>minSqN+.01&&(sPDG<80000000||tM>rtM)&&!hsflag)// Mass(RQ)>M_min && Mass(RNucleus)>M_min
	  //if(kt>minSqB+.01&&(sPDG<80000000||tM>rtM)&&!hsflag)// Mass(RQ)>M_min && Mass(RNucleus)>M_min
      //if(kt>minSqT+.01&&(sPDG<80000000||tM>rtM)&&!hsflag)// Mass(RQ)>M_min && Mass(RNucleus)>M_min
      //if((kt>minSqT+.01||sMass<1200.)&&tM>rtM&&!hsflag)// Mass(RQ)>M_min && Mass(RNucleus)>M_min
      //if((kt>minSqT+.01||sPDG==90000001)&&tM>rtM&&!hsflag)// Mass(RQ)>M_min && Mass(RNucl)>M_min
      //if(tM>rtM&&!hsflag)                               // Mass(RNucleus)>M_min
      //if(tM>rtM&&kp<kpMax&&!hsflag)                       // Mass(RNucleus)>M_min
	  //if(2>3)                                      // *** Close k-correction ***
	  {
        // ***VTN***
		if(sPDG<80000000)                            // Hadronic candidate: finish calculations
		//if(2>3)                                      // Close crrection for hadrons @@??p-pbar??@@
	    {
          G4double np = nOfQ - 3;                    // Power for an residual quasmon mass
          G4double cM = pow(minSqT/kt,np);           // Cut for PossibleQuasmon residual mass
          G4double uR = G4UniformRand();
          G4double rn = pow(cM +(1.-cM)*uR, 1./np);
#ifdef debug
          G4cout<<"G4Q::HQ: YES it is big enough: t="<<kt<<" > T="<<minSqT<<", np="<<np<<", cM="
                <<cM<<", uR="<<uR<<", rn="<<rn<<G4endl;
#endif
          m2 = kt*rn;                                // SquaredMass of ResidualQuasmon (hadr)
	    }
        else m2 = kt;                                // SquaredMass of ResidualQuasmon (fragm)
      }
	  else
      {
        hsflag=true;                                 // Decay in ThisHadron/Fragment+(SHadron)
#ifdef debug
		G4cout<<"G4Q::HQ:***hsflag!*** tM="<<tM<<"<"<<rtM<<G4endl;
#endif
	  }
#ifdef debug
	  G4cout<<"G4Q::HQ:******>>>>> rMass="<<rMass<<", sqm2="<<sqrt(m2)<<", hsflag="<<hsflag<<G4endl;
#endif
      rMass=sqrt(m2);                                // if(hsflag)(TwoPartDecay) it's fake 0.
      G4int cB    = abs(curQ.GetBaryonNumber());     // Baryon Number of residual
      ///////////G4int cS    = abs(curQ.GetStrangeness());      // Strangenes of residual
      rPDG        = curQ.GetSPDGCode();              // PDG of lowest ResiduaiQuasm hadronic state
      G4int aPDG  = abs(rPDG);
      G4int rb    = abs(curQ.GetBaryonNumber());     // BaryNum of residual hadronic state
      G4double rcMass=-BIG;   //@@ just a BIG number // Prototype of minimal mass for residual
      if (!rPDG)
      {
        G4cerr<<"***G4Quasmon::HadronizeQuasmon: rQ ="<<curQ<<", rPDG="<<rPDG<<"(b="<<rb
              <<") + sPDG="<<sPDG<<"(sM="<<sMass<<")"<<G4endl;
        G4Exception("G4Quasmon::HadronizeQuasmon: unidentifiable residual Hadron");
      }
      G4double sB = 1473.;                           // @@ Mean of DELTA(3/2) & N(1680)(5/2)
      if(rPDG!=10) rcMass=G4QPDGCode(rPDG).GetMass();// residual mass for not Chipolino case
      else         sB=0.;                            // Chipolino never decays in hadrons
      if(rPDG==221 || rPDG==331) rcMass=mPi0;
      G4double bs = rcMass+mPi0;
      G4bool rFl=false;                              // true: ResidualResonance,false: Quasmon
      if(sPDG<80000000&&envPDG==NUCPDG)rFl=G4UniformRand()<bs*bs/m2;// ProbabFunct: m_min^2/m^2
#ifdef debug
	  G4cout<<"G4Q::HQ: sPDG="<<sPDG<<", hsflag="<<hsflag<<", rPDG="<<rPDG<<curQ<<",rM="<<rMass
            <<",rb="<<rb<<",rFl="<<rFl<<G4endl;
#endif
      if(!hsflag && rFl && rPDG && rPDG!=10 && rb<2 && aPDG!=1114 && aPDG!=2224 && aPDG!=3334)
	  { // ------------------------------->>>>>>>>>>>>>>> Hadron-Parton Duality decay
        G4int regPDG = 0;                            // PDG prototype for G-meson
        G4int refPDG = 0;                            // PDG prototype for F-meson
        G4int redPDG = 0;                            // PDG prototype for D-meson
        G4int repPDG = 0;                            // PDG prototype for P-meson
        if(rPDG && rPDG!=10)                         // Can the residual be a Hadron ?
	    {
          if     (rPDG== 3122) rPDG= 3212;           // LAMBDA* converted to SIGMA*
          else if(rPDG==-3122) rPDG=-3212;
          if(rPDG>0)repPDG=rPDG+2;                   // Make P-state out of S-state
          else      repPDG=rPDG-2;                   // Subtract 2 for the negative PDG
          if(repPDG>0)redPDG=repPDG+2;               // Make D-state out of P-state
          else        redPDG=repPDG-2;               // Subtract 2 for the negative PDG
          if(redPDG>0)refPDG=redPDG+2;               // Make F-state out of D-state
          else        refPDG=redPDG-2;               // Subtract 2 for the negative PDG
          if(refPDG>0)regPDG=refPDG+2;               // Make G-state out of F-state
          else        regPDG=refPDG-2;               // Subtract 2 for the negative PDG
          if(rPDG==221)      rPDG=111;               // eta=>Pi0
#ifdef debug
	      G4cout<<"G4Q::HQ:MQ="<<quasM<<"(QQC="<<valQ<<")->H("<<sPDG<<")+RGS("<<rPDG<<")"<<G4endl;
#endif
          G4double resM2  = G4QPDGCode(  rPDG).GetMass2(); // Mass^2 of the S-resonance
          G4double repM2  = G4QPDGCode(repPDG).GetMass2(); // Mass^2 of the P-resonance
          G4double redM2  = G4QPDGCode(redPDG).GetMass2(); // Mass^2 of the D-resonance
          G4double refM2  = G4QPDGCode(refPDG).GetMass2(); // Mass^2 of the F-resonance
          sB              = sqrt((resM2+repM2)/2.);  // Boundary between S&P resonances
          G4double     pB = sqrt((repM2+redM2)/2.);  // Boundary between P&D resonances
          G4double     dB = sqrt((redM2+refM2)/2.);  // Boundary between D&F resonances
          G4double     fB = sqrt(refM2)+150.;        // Fake boundary between F&G resonances
          if(!cB)      fB=  sqrt((refM2+G4QPDGCode(regPDG).GetMass2())/2.);
          G4double dif=quasM-sMass;
          G4double rM = GetRandomMass(repPDG,dif);   // Randomize Mass of P-resonance
          G4double dM = GetRandomMass(redPDG,dif);   // Randomize Mass of D-resonance
          G4double fM = GetRandomMass(refPDG,dif);   // Randomize Mass of F-resonance
#ifdef debug
          G4cout<<"G4Q::HQ: rM="<<rMass<<", sB="<<sB<<" (bQ="<<bQ<<"), pB="<<pB<<", dB="<<dB
                <<", fB="<<fB<<G4endl;
#endif
          if(((rM>0 && rMass<pB && rMass>sB) || (dM>0 && rMass>pB && rMass<dB) ||
             (fM>0 && rMass>dB && rMass<fB)) && theEnvironment.GetPDG()==NUCPDG)
		  { // Final H+R decay of QUASMON in vacuum case (should not exist if Environ exists)
            if     (rMass>pB && rMass<dB && dM>0)    // D-resonance case
		    {
              repPDG=redPDG;
              rM=dM;
		    }
            else if(rMass>dB && rMass<fB && dM>0)    // F-resonance case
		    {
              repPDG=refPDG;
              rM=fM;
		    }
#ifdef debug
            G4cout<<"G4Q::HQ:s="<<sPDG<<",Q=>rM="<<rMass<<"(minQ="<<rPDG<<curQ<<")+sB="<<sB<<G4endl;
#endif
            if(quasM<rM+sMass &&(sPDG==221||sPDG==331))// Change eta-Candidate to pi-Candidate
	        {
              sPDG=111;
              sMass=mPi0;
  	        }
            G4LorentzVector r4Mom(0.,0.,0.,rM);      // P/D-resonance with a random Mass
            G4LorentzVector s4Mom(0.,0.,0.,sMass);   // Mass's random since probability time
            if(!G4QHadron(q4Mom).DecayIn2(r4Mom, s4Mom))
            {
              G4cerr<<"**G4Quasmon::HadronizeQuasmon: rPD="<<repPDG<<"(rM="<<rMass<<") + sPD="
					<<sPDG<<"(sM="<<sMass<<"), Env="<<theEnvironment<<G4endl;
    	      G4Exception("G4Quasmon::HadronizeQuasmon: Hadron+Resonance decay didn't succeed");
            }
#ifdef debug
	        G4cout<<"G4Q::HQ:=== 1 ===> HadronVec, Q="<<q4Mom<<" -> s4M="<<s4Mom<<"("<<sPDG
                  <<"), r4M="<<r4Mom<<"("<<repPDG<<")"<<G4endl;
#endif
            //@@CHECK CoulBar and may be evaporate instead
            G4QHadron* curHadr1 = new G4QHadron(repPDG,r4Mom,curQ);// Create Hadron for Residual
            FillHadronVector(curHadr1);              // Fill "new curHadr1" (delete equivalent)
            // @@@ Renaming correction to distinguish from the evaporation (for BF/FSI purposes)
            if     (sPDG==2112) sPDG=90000001;       // rename PDG of the neutron
            else if(sPDG==2212) sPDG=90001000;       // rename PDG of the proton
            else if(sPDG==3122) sPDG=91000000;       // rename PDG of the lambda
            // @@@ ^^^^^^^^^^^^^^^^^^^^
            G4QHadron* curHadr2 = new G4QHadron(sPDG,s4Mom); // Creation Hadron for a Candidate
            FillHadronVector(curHadr2);              // Fill "new curHadr2" (delete equivalent)
            KillQuasmon();                           // This Quasmon is done
            qEnv=theEnvironment;                     // Update the QEnvironment
            return theQHadrons;                      // The last decay of the quasmon...
	      }
	    }
      }
      curQ = memQ;                                   // Recover original curQ=valQ-candQ
      fdul = (rMass<sB&&rFl&&rPDG!=10);
    }
    if(fprob)rPDG=curQ.GetSPDGCode();                // PDG of residual Quasmon as S-hadron
    G4double reMass=sqrt(minSqT);                    // Min Residual Quasmon Mass after decay
    if (!rPDG)
    {
      G4cerr<<"***G4Q::HQ: rQ="<<curQ<<",rPDG="<<rPDG<<"+sPDG="<<sPDG<<"(sM="<<sMass<<")"<<G4endl;
      G4Exception("G4Quasmon::HadronizeQuasmon: unidentifiable residual Hadron");
    }
    if(rPDG==221) reMass=mPi0;
    G4double aMass=mPi0;
    //if(envPDG>80000000&&(sPDG<80000000||envPDG!=pPDG))aMass=0.;
    if((sPDG<80000000&&envPDG>80000000&&envPDG!=NUCPDG||sPDG>80000000&&sPDG!=NUCPDG
        &&envPDG>pPDG) && iniBN>0||iniBN>1||rPDG==10) aMass=0.; // No Pi0 condition (eg in NucE)
#ifdef debug
    G4cout<<"G4Q::HQ:rPDG="<<rPDG<<",QB="<<iniBN<<". Is [hsf="<<hsflag<<"] or [rdf="<<fdul
          <<"] or [rM="<<rMass<<"<"<<reMass<<"+"<<aMass<<" and rM2="<<reTNM2<<"<miM2="<<tmpTM2
          <<"] to fail?"<<G4endl;
#endif
    //        *** This is a condition, when fragmentation can be continued ***
    if(hsflag||!(rMass>reMass+aMass||envPDG>NUCPDG&&envPDG>pPDG&&reTNM2>tmpTM2)||fdul)
      // >>>>>>>>>>> Decay Q->S+H or Q/C->H1+H2 or suck in or evaporate or slow down or decay etc.
	{
      // === At this point decide what to do, if fragmentation in this Candidate is impossible ===
#ifdef debug
      G4cout<<"G4Q::HQ: Yes, hsf="<<hsflag<<",sPDG="<<sPDG<<",pM="<<pMass<<",Env="<<envPDG<<",QM="
            <<quasM<<valQ<<", fpr="<<fprob<<G4endl;
#endif
      G4QPDGCode rQPDG=G4QPDGCode(rPDG);
      if(hsflag)rMass=rQPDG.GetMass();               // Get Nuclear Cluster mass
      if(sPDG>80000000&&sPDG!=NUCPDG)                // "Quark exchange hadronization" case
	  {
        G4QContent tmpEQ=envQC-pQC;                  // Quark content for Residual Environment
        G4QNucleus tmpN(tmpEQ);                      // Pseudo nucleus for Residual Environment
        G4double   tmpNM=tmpN.GetMZNS();             // GS Mass of Residual Environment
        G4QContent tmpRQ=valQ+transQC;               // QContent of Residual Quasmon
        G4QNucleus tmpR(tmpRQ);                      // Nucleus for Residual Quasmon
        G4double   tmpRM=tmpR.GetMZNS();             // GS Mass of Residual Quasmon
        G4LorentzVector ResEnv4Mom(0.,0.,0.,tmpNM);  // 4-Mom for the Residual Environment
        if(rQ4Mom==zeroLV)
		{
#ifdef ppdebug
          G4cout<<"G4Q::HQ: ***NEEDS-EVAP-5*** Q="<<q4Mom<<valQ<<", QEnv="<<theEnvironment<<G4endl;
          G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (5)"); //@@ TMP
#endif
          qEnv=theEnvironment;
          return theQHadrons;
		}
        G4LorentzVector retN4Mom=rQ4Mom+ResEnv4Mom;// 4-Mom for the Total Residual Nucleus
        G4double retNM2=retN4Mom.m2();             // SqMass of the Total Residual Nucleus @@ double
        G4double retNM=sqrt(retNM2);               // Mass of the Total Residual Nucleus @@ double ?
        G4QContent tmpTQ=tmpEQ+tmpRQ;              // QC for TotalResidualNucleus to the Fragment
        G4QNucleus tmpT(tmpTQ);                    // Nucleus for the Total Residual Nucleus
        //G4QNucleus tmpT(tmpTQ,retN4Mom);           // Nucleus for the Total Residual Nucleus
        G4double   tmpTM=tmpT.GetMZNS();           // MinMass of the Total Residual Nucleus @@double
        if(tmpTM>retNM) tmpT=G4QNucleus(tmpTQ,retN4Mom);
        /////////G4double fM2=fr4Mom.m2();         // Squared mass of the outgoing fragment
        G4QPDGCode sQPDG(sPDG);
        // tmpNM - mass of residual environment (M), retNM - mass of TotalNucl [MN=sqrt((E+M)^2-p^2]
        // tmpRM - mass GS mass of ResidQuasmon (m_GS), rQ4Mom - 4-momentum of ResidQuasmon (E,p,m)
        ////G4double m2=rQ4Mom.m2();                     // Real Squared Mass of the ResidualQuasmon
        // potE=[sqrt(E^2*M^2-m^2*M^2+m_GS^2*MN*2)-m^2-E*M]/MN=-U (binding energy should be cut off)
        //G4double pEc=2*(tmpRM+tmpNM-retNM);        // DoubledBindingEnergy (virial theorem)
        /////////G4double pEc=tmpRM+tmpNM-retNM;              // BindingEnergy (relativistic effect)
#ifdef debug
        G4int sB=sQPDG.GetBaryNum();
        G4double sCB=theEnvironment.CoulombBarrier(sQPDG.GetCharge(),sB);
        G4bool bSplit=tmpT.SplitBaryon(); // Possibility to split a baryon from the TotResN
        G4int rB=rQPDG.GetBaryNum();
        G4double rCB=theEnvironment.CoulombBarrier(rQPDG.GetCharge(),rB);
        G4double EQ=rQ4Mom.e();                      // EnergyOfResidualQuasmon (E)
        G4double em=tmpNM*EQ;                        // ResEnvM * EnergyOfResidualQuasmon (M*E)
        G4double mM=retNM*tmpRM;                     // TotResNuclM*GSMassResidualQuasmon (m_GS*MN)
        G4double pEn=(sqrt(em*em-m2*tmpNM*tmpNM+mM*mM)-m2-em)/retNM; //Real BindingEnergy
        G4double pEt=tmpNM*(EQ+tmpNM-retNM)/retNM;   // Energy Transfer to nucleus
        G4double PQ=rQ4Mom.rho();                    // 3MomentumOfResidualQuasmon (p)
        G4double pPt=tmpNM*PQ/retNM;                 // 3Momentum Transfer to nucleus
        G4cout<<"G4Q::HQ:rNM="<<retNM<<" > mMRE="<<tmpNM<<" + mMRQ="<<tmpRM<<" = "<<tmpNM+tmpRM
              <<", CB="<<rCB<<", bS="<<bSplit<<", pEn="<<pEn<<", pEt="<<pEt<<G4endl;
#endif
        /////////G4bool epC=3<2; // ** Decay in Fragment & ResidNucleus is ByDefault closed **
        //G4bool epC=pEt<pEc;
        //G4bool epC=pEn<pEc+pEc;                    // RealBindE<normalBindE
        //G4bool epC=pEn<40.;                        // RealBindE<normalBindE
        //G4bool epC=pPt<400.;
        //G4bool epC=pPt<720./pow(envA,.1236);       //@@ Temporary solution for splitting nuclei
        //G4bool epC=2<3;              // ** Decay in Fragment & ResidNucleus is ByDefault opened **
        //if(sPDG==90000001&&pEt>pEc) epC=3<2;
        //if(sPDG==90000001&&(pPt>270.||envA-sB<4)) epC=3<2;
        //if(pPt<420.&&envA-sB>4) epC=2<3;
        //if(sB==1&&pPt<450.&&envA>5) epC=2<3;
        //if(sPDG==90000001&&pPt<480.&&envA>5) epC=2<3;
        //if(sB==1&&pPt<470.&&envA>5) epC=2<3; // !!! Was important for E-dep of tot gammaA c.s.
        //////////////////////////////////////////////////////if(sB==1&&pPt<460.&&envA>5) epC=2<3;
        //if(sB<3&&envA>5) epC=2<3;
        //if(pPt>270.||envA-sB<4) epC=3<2;
        //if(envA-sB<4) epC=3<2;
        //G4bool ppC=pPt<250.;                       // RealBindE<normalBindE
        ////////G4bool ppC=3<2;                              // *** Momentum reduction is closed ***
        //G4bool ppC=2<3;                            //  *** Momentum reduction is opened ***
        //if(sB>2&&sB!=4||rB>2&&rB!=4) ppC=3<2;      // Bad idea of stupid symmetry
        //if(sB>2&&sB!=4||rB>2) ppC=3<2;
        //if(rB>2) ppC=3<2;
        //if(sB>1&&rB>1) ppC=3<2;
        //if(rB>1&&rB!=4) ppC=3<2;
        //if(rB>1) ppC=3<2;
        //G4bool ppC=pEn<pEc+pEc;                    // RealBindE<normalBindE
		//retNM=realMOfTotResN>tmpNM=minMaOfResidEnv+minMOfResidQ == Comment for condition below ==
        //if(nQuasms==1&&(retNM>tmpNM+tmpRM+rCB||bSplit)&&ppC)//=> "Go out & drop ExcitedNucleus"
        //if(nQuasms==1&&(retNM>tmpNM+tmpRM||bSplit)&&ppC) //=> "Go out & drop ExcitedNucleus"
        // The condition below is just that the Residual Nucleus is heavier than the GS of the ResN
		if(2>3)                                      // *** Attempt GoOutDropExcNucl is closed ***
		{
          G4LorentzVector retQ4Mom=retN4Mom-G4LorentzVector(0.,0.,0.,tmpNM);
#ifdef debug
          G4double retQ3M=retQ4Mom.rho();
          G4cout<<"G4Q::HQ:*CheckFermiMotion* rQ3M="<<retQ3M<<",rQB="<<iniBN<<",E="<<excE<<G4endl;
#endif
          //if(retNM>tmpNM+tmpRM&&retQ3M>100000.)    // FermiMotionCut for DecayIn3 (@@dir) *Closed*
          //if(retNM>tmpNM+tmpRM)                    // Decay in Q+E+F (@@don't conserve dir of Fr.)
		  if(2>3)                                    // *** Closed ***
		  {
            G4LorentzVector re4M = G4LorentzVector(0.,0.,0.,tmpNM); // GS Mass of Residual Environ
            G4LorentzVector rq4M = G4LorentzVector(0.,0.,0.,tmpRM); // GS Mass of Residual Quasmon
#ifdef debug
            G4cout<<"G4Q::HQ: >***>>>YES<<<***< rQM="<<tmpRM<<",rEM="<<tmpNM<<",fM="<<sMass<<G4endl;
#endif
            if(!G4QHadron(tot4M).DecayIn3(rq4M,re4M,fr4Mom))
            {
              G4cerr<<"***G4Quasm::HadrQuasm: Decay in Fragm+ResQ+ResEnv did not succeeded"<<G4endl;
              G4Exception("***G4Quasmon::HadrQuasm: Failed to decay in Frag+ResQ+ResEnv");
            }
            else
            {
              //@@CHECK CoulBar (only for ResQuasmon in respect to ResEnv) and may be evaporate
              G4QHadron* resQH = new G4QHadron(tmpRQ,rq4M); // Create a Hadron for the ResidQuasmon
              if(nQuasms==1) qEnv=G4QNucleus(tmpEQ,re4M);// Create Nucleus for MovingResidualEnviron
              else
              {
                G4QHadron* envH = new G4QHadron(tmpEQ,re4M); //@@ Moving Environment !
                FillHadronVector(envH);              // Fill MovingEnvironment (delete equiv)
                qEnv = vacuum;
			  }
              G4QHadron* candH = new G4QHadron(sPDG,fr4Mom);// Create a Hadron for the Candidate
              FillHadronVector(candH);               // Fill CandiFragm Hadron (delete equivalent)
              FillHadronVector(resQH);               // Fill ResidQuasm Hadron (delete equivalent)
              KillQuasmon();                         // This Quasmon is done
              return theQHadrons;                    // The last decay of the quasmon...
		    }
		  }
          else                                       // ==> "Normal decay" !!!!!
          {
            G4QHadron* candHadr = new G4QHadron(sPDG,fr4Mom);// Createa a Hadron for the Candidate
            FillHadronVector(candHadr);              // Fill the RadiatedHadron (delete equivalent)
            theEnvironment.Reduce(pPDG);             // Update NuclearEnviron after Q->RQ+F
            q4Mom  = retQ4Mom;                       // Update Quasmon's 4-mom
            valQ  += transQC;                        // Update the Quark Content of Quasmon
#ifdef debug
            G4cout<<"G4Q::HQ:***LetAnyway>>s="<<sPDG<<fr4Mom<<",E="<<theEnvironment.GetPDG()<<",MQ="
                  <<q4Mom.m()<<valQ<<G4endl;
#endif
		  }
          status=1;                                  // Something was successfuly done           
          phot4M=zeroLV;
          // === Check of boundary to escape not existing state of residual nucleus ===
          if(CheckGroundState()) KillQuasmon();  // This Quasmon is done
          //if(CheckGroundState(true)) KillQuasmon();  // This Quasmon is done
          qEnv=theEnvironment;
          return theQHadrons;                        // The last decay of the quasmon...
		}
        //else if(nQuasms==1&&tmpTM+sMass<totMass&&fM2>0.) //Lose momentum,go out,leave GSResNuc
		//else if(nQuasms==1&&tmpTM+sMass+sCB<totMass&&fM2>0.) //Lose momentum,go out,leave GSResNuc
		//else if(nQuasms==1&&tmpTM+sMass<totMass&&fM2>0.&&epC)//"Lose 3-mom. to GSTotNM & go out"
		else if(2>3) // *************** Decay in Fragm + GSResidNucleus is Closed ****************
		{
          G4LorentzVector dif1=tot4M-retN4Mom-fr4Mom;
          qEnv = G4QNucleus(tmpTQ,retN4Mom);         // Createa Nucleus for the TotalResidualNucleus
#ifdef debug
          G4cout<<"***G4Q::HQ: ***KeepEnvironmentMoving***, nQ="<<nQuasms<<",pPt="<<pPt<<G4endl;
#endif
          //G4QHadron* envH = new G4QHadron(tmpTQ,retN4Mom); //@@ Moving Environment !
          //FillHadronVector(envH);                  // Fill Moving Environment (delete equivalent)
          //qEnv = vacuum;
#ifdef debug
          G4cout<<"G4Q::HQ: rnM="<<retN4Mom.m()<<",> GSM="<<tmpTM<<",dif4M="<<dif1<<G4endl;
#endif
          G4QHadron* candHadr = new G4QHadron(sPDG,fr4Mom);// Createa a Hadron for the Candidate
          FillHadronVector(candHadr);                // Fill the RadiatedHadron (delete equivalent)
#ifdef debug
          G4double frM=fr4Mom.m();
          G4LorentzVector dif2=tot4M-retN4Mom-fr4Mom;
          G4cout<<"G4Q::HQ:sM="<<sMass<<"="<<frM<<", fT="<<fr4Mom.e()-frM<<",dif24M="<<dif2<<G4endl;
#endif
          KillQuasmon();                             // This Quasmon is done
          return theQHadrons;                        // The last decay of the quasmon...
		}
        else if(totBN>1&&totMass>totM&&totS>=0&&envPDG>80000000&&envPDG!=NUCPDG) // Evaporation
		{
#ifdef ppdebug
          //@@ May be recalculate hadronization ??
          G4double fraM=fr4Mom.m();
          G4double kinE=fr4Mom.e()-fraM;
          G4double sumM=tmpTM+fraM;
          G4cout<<"G4Q::HQ:EV-6:s="<<sPDG<<",T="<<kinE<<",RM="<<retN4Mom.m()<<"<"<<tmpTM<<",tQC="
                <<transQC<<",E="<<excE<<",sM="<<sumM<<">tM="<<totMass<<",sCB="<<sCB<<",nQ="<<nQuasms
                <<G4endl;
          //if(excE>.027*fraM) G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (6)"); //@@ TMP
          //if(excE>60.) G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (6)"); //@@ TMP
          //if(totMass<sumM) G4Exception("G4Quasmon::HadrQuasm: Why Fail? (6). Prod's>TotM");//@@TMP
          G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (6) Products are more than total");
#endif
#ifdef debug
          G4cout<<"G4Quasm::HQ:Q="<<q4Mom<<quasM<<",Env="<<theEnvironment<<",Ph="<<phot4M<<G4endl;
#endif
          qEnv=theEnvironment;
          return theQHadrons;
		}
        else if(totBN==1&&nQuasms==1) // Decay of the baryonic Total Nucleus state
		{
#ifdef debug
          G4cout<<"G4Q::HQ:tB=1,nQ=1,tZ="<<totZ<<",tS="<<totS<<",T="<<totQC<<",M="<<totMass<<G4endl;
#endif
          G4double nucM= mProt;
          G4double piM = 0.;
          G4int nucPDG = 2212;
          G4int piPDG  = 22;
          if(abs(totS)==1)
          {
            if(totS==1)               // Decay of the strange hyperstate
			{
              if(!totZ&&totMass>mLamb+mPi0)
			  {
                nucM  = mLamb;
                nucPDG= 3122;
                piM   = mPi0;
                piPDG = 111;
			  }
              else if(abs(totZ)==1&&totMass>mLamb+mPi)
			  {
                nucM  = mLamb;
                nucPDG= 3122;
                piM   = mPi;
                if(totZ>0) piPDG = 211;
                else       piPDG =-211;
			  }
              else
			  {
                G4cerr<<"***G4Q::HQ:tZ="<<totZ<<",tS="<<totS<<",T="<<totQC<<",M="<<totMass<<G4endl;
    	        G4Exception("G4Quasmon::HadronizeQuasmon: Pion + Lambda decay did not succeed");
              }
			}
			else                      // Decay of the anti-strange hyperstate
			{
              if(!totZ&&totMass>mNeut+mK0)
			  {
                nucM  = mNeut;
                nucPDG= 2112;
                piM   = mK0;
                piPDG = 311;
			  }
              else if(totZ==2&&totMass>mProt+mK)
			  {
                piM   = mK;
                piPDG = 321;
			  }
              else if(totZ==1&&totMass>mProt+mK0&&G4UniformRand()>0.5)
			  {
                piM   = mK0;
                piPDG = 311;
			  }
              else if(totZ==1&&totMass>=mNeut+mK)
			  {
                nucM  = mNeut;
                nucPDG= 2112;
                piM   = mK;
                piPDG = 321;
			  }
              else
			  {
                G4cerr<<"***G4Q::HQ:tZ="<<totZ<<",tS="<<totS<<",T="<<totQC<<",M="<<totMass<<G4endl;
    	        G4Exception("G4Quasmon::HadronizeQuasmon: Kaon + Nucleon decay did not succeed");
              }
			}
		  }
          else if(totMass>PiNM&&!totS) // Decay in nucleon & pion
          {
            if(!totZ&&totMass>mProt+mPi&&G4UniformRand()<0.5)
		    {
              piM   = mPi;
              piPDG = -211;
            }
            else if(!totZ&&totMass>mNeut+mPi0)
		    {
              nucM  = mNeut;
              nucPDG= 2112;
              piM   = mPi0;
              piPDG = 111;
            }
            else if(totZ==1&&totMass>mNeut+mPi&&G4UniformRand()<0.5)
		    {
              nucM  = mNeut;
              nucPDG= 2112;
              piM   = mPi;
              piPDG = 211;
            }
            else if(totZ==1&&totMass>mProt+mPi0)
		    {
              piM   = mPi0;
              piPDG = 111;
            }
            else if(totZ==-1)
		    {
              nucM  = mNeut;
              nucPDG= 2112;
              piM   = mPi;
              piPDG = -211;
            }
            else if(totZ==2)
		    {
              piM   = mPi;
              piPDG = 211;
            }
            else
			{
              G4cerr<<"***G4Q::HQ:totZ="<<totZ<<",totB="<<totBN<<",E="<<envQC<<",Q="<<valQ<<G4endl;
    	      G4Exception("G4Quasmon::HadronizeQuasmon: Pion + Nucleon decay did not succeed");
            }
		  }
          else if(!totS)
		  {
            if(!totZ)
			{
              nucM=mNeut;
              nucPDG=2112;
            }
            else if(totZ<0||totZ>1)
			{
              G4cerr<<"***G4Q::HQ:totZ="<<totZ<<",totB="<<totBN<<",E="<<envQC<<",Q="<<valQ<<G4endl;
    	      G4Exception("G4Quasmon::HadronizeQuasmon: Photon + Nucleon decay did not succeed");
            }
          }
          G4LorentzVector pi4M(0.,0.,0.,piM);       // mass of the kaon/pion/photon
          G4LorentzVector nuc4M(0.,0.,0.,nucM);     // mass of the nucleon
          if(!G4QHadron(tot4M).DecayIn2(pi4M, nuc4M))
          {
            G4cerr<<"***G4Q::HQ:T="<<tot4M<<totMass<<"->gam/pi/K+N="<<nucPDG<<", MN="<<nucM<<G4endl;
    	    G4Exception("G4Quasmon::HadronizeQuasmon: gamma/Pi + Nucleon decay did not succeed");
          }
#ifdef debug
	      G4cout<<"G4Q::HQ:T="<<tot4M<<totMass<<"->GPK="<<piPDG<<pi4M<<"+B="<<nucPDG<<nuc4M<<G4endl;
#endif
          G4QHadron* piH = new G4QHadron(piPDG,pi4M);// Create Hadron for gamma/Pion
          FillHadronVector(piH);                     // Fill "new piH" (delete equivalent)
          G4QHadron* nucH = new G4QHadron(nucPDG,nuc4M); // Creation Hadron for the nucleon
          FillHadronVector(nucH);                    // Fill "new nucH" (delete equivalent)
          KillQuasmon();                             // This Quasmon is done          
          qEnv=vacuum;
          return theQHadrons;                        // The last decay of the total nucleus...
        }
#ifdef debug
        else G4cerr<<"***G4Q::HQ: B="<<totBN<<",tM="<<totMass<<" > M="<<totM<<",S="<<totS<<G4endl;
#endif
	  }
      G4double dm=quasM-sMass;
#ifdef debug
      G4cout<<"G4Q::HQ:f="<<fprob<<",d="<<dm<<",rPDG="<<rPDG<<",rM="<<rMass<<",M="<<reMass<<G4endl;
#endif
      G4double rWi=0.;
      if(rPDG!=10) rWi=G4QPDGCode(rPDG).GetWidth();
      if(rPDG!=10&&rMass>dm&&!rWi) // Try to use a width of the h-resonance or reduce its spin state
	  {
        G4double sWi=G4QPDGCode(sPDG).GetWidth();
        if(sWi)                                      // Hadron is a resonance
		{
          G4double ddm=quasM-rMass;
          if(sMass-ddm<1.5*sWi-.01)
          {
            sMass=GetRandomMass(sPDG,ddm);           // Randomize mass of the Reson-Hadron
#ifdef debug
            G4cout<<"G4Q::HQ: "<<sPDG<<" mass's changed to "<<sMass<<",d="<<ddm<<",W="<<sWi<<G4endl;
#endif
		  }
          //else if(G4int ds1=abs(sPDG)%10>2)        // @@ Make it C-style ?: operator ?
		  //{
          //  G4int oPDG=sPDG;
          //  if   (sPDG>0&&sPDG%2==0)   sPDG-=ds1-2;
          //  else if(sPDG>0)            sPDG-=ds1-1;
          //  else if(sPDG<0&&sPDG%2==0) sPDG+=ds1-2;
          //  else if(sPDG<0)            sPDG+=ds1-1;
          //  sMass=G4QPDGCode(sPDG).GetMass();
          //  G4cout<<"**G4Q::HQ: sPDG="<<oPDG<<" is changed to "<<sPDG<<", sM="<<sMass<<G4endl;
		  //}
		}
      }
      // First decay suppression for the final decay in 2 particles
      if(( (rPDG==111||rPDG==221||rPDG==331)
               && (sPDG==221||sPDG==331||sPDG==223||sPDG==333||sPDG==225||sPDG==335) ||
           (rPDG==221||rPDG==331)
               && (sPDG==111||sPDG==113||sPDG==115)
         ) && G4UniformRand()>EtaEtaprime && dm>782.) rPDG=223;
      // Final state eta/eta' suppression
      if(rPDG==331&&dm<958.)rPDG=221;
      if(rPDG==221&&dm<548.)rPDG=111;
      if(rPDG<80000000&&rPDG!=2212&&rPDG!=2112&&rPDG!=3122&&rPDG!=10) // => ResidQ isn't a NuclClust
      {
        reMass=GetRandomMass(rPDG,dm);               // Randomize mass of the RResidQuasmon-Hadron
#ifdef debug
        G4cout<<"G4Q::HQ:(reMass) ResidQuasmonMass="<<reMass<<" is changed for rPDG="<<rPDG<<G4endl;
#endif
        if(reMass==0.)
	    {
          if(sPDG==221 || sPDG==331)                 // Change eta-Candidate to pi/gamma
	      {
            if     (sPDG==221) dm+=mEta-mPi0;
            else if(sPDG==331) dm+=mEtaP-mPi0;
            if(dm<0)
			{
              dm+=mPi0;
              sPDG=22;
              sMass=0.;
			}
            else
            {
              sPDG=111;
              sMass=mPi0;
			}
            if(dm<mPi0-.00001&&rPDG==111)
			{
              rPDG=22;
              reMass=0.;
			}
            else reMass=GetRandomMass(rPDG,dm);      // Rerandomize mass of ResidQuasmon-Resonance
            if(reMass==0.) G4cerr<<"***G4Quasmon::HadronizeQ: (2) rPDG="<<rPDG<<", dm="<<dm<<G4endl;
		  }
          else if(rPDG==111)                         // Make a photon out of the Resid Quasmon
		  {
            rPDG=22;
            reMass=0.;
		  }
          else
		  {
            if(CheckGroundState()) KillQuasmon();// This Quasmon is done
            //if(CheckGroundState(true)) KillQuasmon();// This Quasmon is done
#ifdef ppdebug
            G4cout<<"G4Q::HQ:NEEDS-EVAP-7: sPDG="<<sPDG<<",Q4M="<<q4Mom<<valQ<<",rP="<<rPDG<<G4endl;
            G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (7)"); //@@ TMP
#endif
            qEnv=theEnvironment;
            return theQHadrons;
	      }
        }
	  }
      else if(rPDG==NUCPDG)
	  {
        if(dm>mPi0)
        {
		   rPDG=111;
           reMass=mPi0;
        }
        else
        {
  		   rPDG=22;
           reMass=0.;
        }
  	  }
      G4double freeRQM=rQPDG.GetMass();
      G4int    RQB = rQPDG.GetBaryNum();
      G4double fRQW= 3*rQPDG.GetWidth();
      if(fRQW<.001) fRQW=.001;
      G4QPDGCode sQPDG(sPDG);
      G4int sChg=sQPDG.GetCharge();
      G4int sBaryn=sQPDG.GetBaryNum();
      G4double sCB=theEnvironment.CoulombBarrier(sChg,sBaryn);
#ifdef debug
	  G4cout<<"G4Q::HQ:frCB="<<sCB<<",frC="<<sChg<<",frB="<<sBaryn<<",E="<<theEnvironment<<G4endl;
#endif
      G4int rChg=rQPDG.GetCharge();
      G4int rBaryn=rQPDG.GetBaryNum();
      G4double rCB=theEnvironment.CoulombBarrier(rChg,rBaryn);
#ifdef debug
	  G4cout<<"G4Q::HQ:rqCB="<<rCB<<",rqC="<<rChg<<",rqB="<<sBaryn<<",rQPDG="<<rQPDG<<G4endl;
#endif
	  if(totBN>1&&totS>=0 && envPDG>80000000 && envPDG!=NUCPDG
         && (reMass+sMass>quasM || sCB+rCB+reMass+sMass+envM>totMass || !RQB&&quasM<diPiM))
	  //if(2>3)                                      // This Evaporation channel is closed
      {
#ifdef ppdebug
        G4cout<<"***G4Q::HQ:E-8:RQM+SM="<<reMass+sMass<<">QM="<<quasM<<", sCB="<<sCB<<" + rCB="
              <<rCB<<" + rM="<<reMass<<" + sMass="<<sMass<<" + eM="<<envM<<" = "
              <<sCB+rCB+reMass+sMass+envM<<"> tM="<<totMass<<", sM="<<reMass+sMass+envM<<G4endl;
        G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (8)"); //@@ TMP
#endif
        qEnv=theEnvironment;
        return theQHadrons;
	  }
      if(rPDG==NUCPDG)                                
      {
        G4cerr<<"***G4Quasmon::HadronizeQuasmon: rPDG=90000000, MV="<<reMass<<G4endl;
        G4Exception("***G4Quasmon::HadronizeQuasmon: Residual Particle is Vacuum");
      }
      if(rPDG==2212&&sPDG==311&&reMass+sMass>quasM)
	  {
        if(mNeut+mK<=quasM+.001)
		{
          reMass=mNeut;
          rPDG  =2112;
          rQPDG=G4QPDGCode(rPDG);
          rChg=rQPDG.GetCharge();
          rBaryn=rQPDG.GetBaryNum();
          rCB=theEnvironment.CoulombBarrier(rChg,rBaryn);
#ifdef debug
	      G4cout<<"G4Q::HQ:NCB="<<rCB<<",NC="<<rChg<<",rqB="<<sBaryn<<",rQPDG="<<rQPDG<<G4endl;
#endif
          freeRQM=mNeut;
          RQB=1;
          fRQW=0.;
          sMass =mK;
          if(mNeut+mK<=quasM) sMass=quasM-mNeut;
          sPDG  =321;
          sQPDG=G4QPDGCode(sPDG);
          sChg=sQPDG.GetCharge();
          sBaryn=sQPDG.GetBaryNum();
          sCB=theEnvironment.CoulombBarrier(sChg,sBaryn);
#ifdef debug
	      G4cout<<"G4Q::HQ:KCB="<<sCB<<",KC="<<sChg<<",frB="<<sBaryn<<",E="<<theEnvironment<<G4endl;
#endif
          curQ=neutQC;
        }
        else
        {
          G4cerr<<"***G4Quasmon::HadronizeQuasmon:(NK) QM="<<quasM<<",d="<<quasM-mNeut-mK<<G4endl;
          G4Exception("***G4Quasmon::HadronizeQuasmon: Can't decay Q in N and K");
        }
      }
      G4LorentzVector r4Mom(0.,0.,0.,reMass);
      G4LorentzVector s4Mom(0.,0.,0.,sMass);         // Mass is random since probab. time
      if(sMass>80000000)
      {
        G4cout<<"G4Q::HQ:Q->RQ+QEXF s="<<sPDG<<",pM="<<pMass<<",E="<<theEnvironment<<G4endl;
        q4Mom+=G4LorentzVector(0.,0.,0.,pMass);
      }
      if(!G4QHadron(q4Mom).DecayIn2(r4Mom, s4Mom))
      {
        G4cerr<<"***G4Q::HQ: M="<<q4Mom.m()<<" => rPDG="<<rPDG<<"(rM="<<reMass<<") + sPDG="<<sPDG
              <<"(sM="<<sMass<<")"<<G4endl;
    	G4Exception("***G4Quasmon::HadronizeQuasmon: Hadron+SHadron DecayIn2 didn't succeed");
      }
      G4double sKE=s4Mom.e()-sMass;
#ifdef debug
	  G4cout<<"G4Q::HQ:=== 2.3 ===>HdronVect s4M="<<s4Mom<<",sPDG="<<sPDG<<", r4M/M="<<r4Mom<<reMass
            <<",fRM="<<freeRQM<<",fRW="<<fRQW<<",rPDG="<<rPDG<<",rCB="<<rCB<<",sCB="<<sCB<<G4endl;
#endif
      //////////G4double rKE=r4Mom.e()-reMass;
      //if(sKE<sCB||rKE<rCB)                           // => "KinEn is below CB, try once more" case
      if(sKE<sCB)                                    // => "KinEn is below CB, try once more" case
	  {
#ifdef ppdebug
        //G4cout<<"****G4Q::HQ:E-9: sKE="<<sKE<<"<sCB="<<sCB<<" or rKE="<<rKE<<"<rCB="<<rCB<<G4endl;
        G4cout<<"****G4Q::HQ:E-9: sKE="<<sKE<<"<sCB="<<sCB<<G4endl;
        G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (9)"); //@@ TMP
#endif
        if(sMass>80000000) q4Mom-=G4LorentzVector(0.,0.,0.,pMass);
        qEnv=theEnvironment;
        return theQHadrons;
	  }
      else if(abs(reMass-freeRQM)<fRQW||envPDG==NUCPDG) // => "ResidQ is a GSHadron/Fragment" case
      {
        G4QHadron* curHadr1 = new G4QHadron(rPDG,r4Mom);// Create RealHadron for the ResidQ
        FillHadronVector(curHadr1);                  // Fill "new curHadr1" (delete equivalent)
        G4QHadron* curHadr2 = new G4QHadron(sPDG,s4Mom);// Creation Hadron for theCandidate
        FillHadronVector(curHadr2);                  // Fill "new curHadr2" (delete equivalent)
#ifdef debug
	    G4cout<<"G4Q::HQ:Decay of Quasmon in RQ="<<rPDG<<r4Mom<<" + Fragment="<<sPDG<<s4Mom
              <<", Env="<<theEnvironment<<G4endl;
#endif
        if(sMass>80000000) theEnvironment.Reduce(pPDG);  // Update NuclearEnviron after Q->RQ+QEXF
        KillQuasmon();                               // This Quasmon is done
        qEnv=theEnvironment;                         // Update the QEnvironment
        return theQHadrons;                          // This is the last decay of the quasmon
	  }
      else                                           // => "Virtual Residual Quasmon" case
	  {
        G4LorentzVector resTotN4Mom=r4Mom+G4LorentzVector(0.,0.,0.,envM);
        G4QContent resTotNQC=envQC+curQ;             // QCof possible residual nucleus
        G4QNucleus resTotN(resTotNQC);               // Pseudo nucleus for Total Residual Nucleus
        //////G4int      resTotNPDG=resTotN.GetPDG();      // PDG Code of the Total Residual Nucleus
        G4double   resTotNM=resTotN.GetMZNS();       // Mass of the Total Residual Nucleus
        if(resTotN4Mom.m()<resTotNM)                 // Needs total evaporation
		{
          //if(totBN>1&&totMass>totM&&totS>=0)         //@@ ??
#ifdef ppdebug
          G4cout<<"G4Q::HQ:*NEEDS-EVAP-10* RTNM="<<resTotN4Mom.m()<<" < minRTNM="<<resTotNM<<G4endl;
          G4Exception("G4Quasmon::HadronizeQuasmon: Why Fail? (10)"); //@@ TMP
#endif
          if(sMass>80000000) q4Mom-=G4LorentzVector(0.,0.,0.,pMass);
          qEnv=theEnvironment;
          return theQHadrons;
        }
        else                                         // Only the ResidualTotalNucleus can evaporate
		{
          G4QHadron* curHadr2 = new G4QHadron(sPDG,s4Mom);// Creation Hadron for the Out Hadron
          FillHadronVector(curHadr2);                // Fill "new curHadr2" (delete equivalent)
          q4Mom = r4Mom;
          if(sMass>80000000)
          {
            theEnvironment.Reduce(pPDG);             // Update NuclearEnviron after Q->RQ+QEXF
            valQ  += transQC;                        // Update the Quark Content of Quasmon
          }
          else valQ  = curQ;                         // Update the Quark Content of Quasmon
#ifdef debug
          G4cout<<"***>G4Q::HQ:S="<<sPDG<<s4Mom<<",E="<<theEnvironment<<",Q="<<q4Mom<<valQ<<G4endl;
#endif
          status=1;                                  // Something was successfuly done           
          phot4M=zeroLV;
          piF=false;
          gaF=false;
          if(CheckGroundState()) KillQuasmon();  // This Quasmon is done
          //if(CheckGroundState(true)) KillQuasmon();  // This Quasmon is done
          qEnv=theEnvironment;
          return theQHadrons;
        }
      }
	}
	if(!fskip)                                       // Continue search for fragmentation
	{
      // ==== From this point for nuclear fragments the procedure is just filling =====
#ifdef debug
	  G4int ePDG=theEnvironment.GetPDG();
      G4double frKin=fr4Mom.e()-sMass;
	  G4cout<<"G4Q::HQ:s="<<sPDG<<fr4Mom<<fr4Mom.m()<<"="<<sMass<<",T="<<frKin<<",E="<<ePDG<<G4endl;
#endif
      ///if(sPDG<80000000&&ePDG==NUCPDG)             // ====>>>> "Hadron candidate in Vacuum" case
	  if(sPDG<80000000)           // Hadronic decay is always as in vacuum @@??
      {
        G4int SQ=totQC.GetStrangeness();
#ifdef debug
	    G4cout<<"G4Q::HQ: sPDG="<<sPDG<<", sM="<<sMass<<", SQ="<<SQ<<G4endl;
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
	        G4cout<<"G4Q::HQ: ***PANIC#2*** tM="<<totMass<<" < KM="<<mK<<","<<mK0<<", rM="<<rKPM
                  <<","<<rK0M<<", d="<<mK+rKPM-totMass<<","<<mK0+rK0M-totMass<<G4endl;
            G4Exception("G4Quasmon::HadronizeQuasmon: Why PANIC? (2)"); //@@ TMP
#endif
            status =-1;                              // Panic exit
            qEnv=theEnvironment;                     // Update the QEnvironment
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
            G4cerr<<"***G4Q::HQ: tM="<<tot4M.m()<<totQC<<" => rPDG="<<rPDG<<"(rM="<<rMass
                  <<") + sPDG="<<sPDG<<"(sM="<<sMass<<")"<<G4endl;
    	    G4Exception("***G4Quasmon::HadrQuasmon:K+ResidNucl DecayIn2 didn't succeed");
          }
#ifdef debug
	      G4cout<<"G4Q::HQ:=== 2.4 ===>HadronVector s="<<sPDG<<s4Mom<<", r="<<rPDG<<r4Mom<<G4endl;
#endif
          //@@CHECK CoulBar and may be evaporate instead
          G4QHadron* curHadr1 = new G4QHadron(rPDG,r4Mom);// Create RealHadron for ResidEnv
          FillHadronVector(curHadr1);                // Fill "new curHadr1" (delete equivalent)
          G4QHadron* curHadr2 = new G4QHadron(sPDG,s4Mom);// Creation Hadron for theCandidate
          FillHadronVector(curHadr2);                // Fill "new curHadr2" (delete equivalent)
          KillQuasmon();                             // This Quasmon is done
          qEnv=theEnvironment;                       // Update the QEnvironment
          return theQHadrons;                        // This is the last decay of the quasmon
		}
        G4bool ffin=false;                           // Flag of final decay in gamma-residual
        if(quasM<rMass+sMass&&(sPDG==221||sPDG==331))// Change eta-Candidate or any to pi
	    {
          sPDG = 111;
          sMass=mPi0;
	    }
        else if(!sPDG)
		{
          if     (iniS<0&&iniQChg+iniQChg>=iniBN)  // Try to decay in K+
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
          else if(iniQChg>iniBN-iniS)              // Try to decay in Pi+
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
#ifdef debug
	    G4cout<<"G4Q::HQ:MQ="<<q4Mom.m()<<"->sPDG="<<sPDG<<"(M="<<sMass<<") + rPDG="<<rPDG<<"(M="
              <<rMass<<")"<<G4endl;
#endif
        if(q4Mom.m()<rMass+sMass)
		{
#ifdef ppdebug
	      G4cout<<"G4Q::HQ: ***PANIC#3*** tM="<<q4Mom.m()<<" < rM="<<rMass<<", sM="<<sMass<<", d="
                <<rMass+sMass-q4Mom.m()<<G4endl;
          G4Exception("G4Quasmon::HadronizeQuasmon: Why PANIC? (3)"); //@@ TMP
#endif
          status =-1;                                // Panic exit
          qEnv=theEnvironment;                       // Update the QEnvironment
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
#ifdef debug
	    G4cout<<"G4Q::HQ:DecayQ="<<q4Mom<<"->s="<<sPDG<<s4Mom<<"+R="<<resQ4Mom<<",f="<<ffin<<G4endl;
#endif
        G4QHadron* candHadr = new G4QHadron(sPDG,s4Mom);// Creation a Hadron for the Candidate
        if(ffin)
		{
          //@@CHECK CoulBar and may be evaporate instead
          theQHadrons.push_back(candHadr);              // Fill the emergency PHOTON (delete equiv.)
          G4QHadron* candHRes = new G4QHadron(rPDG,resQ4Mom);// Creation Hadron for the QResid
          FillHadronVector(candHRes);                // Fill "new candHRes" (delete equivalent)
		  //#ifdef debug
          G4cout<<"***>>>G4Quasmon::HadrQ:EmergencyHadron="<<sPDG<<",QEnv="<<theEnvironment<<G4endl;
		  //#endif
          KillQuasmon();                             // This Quasmon is done
          qEnv=theEnvironment;
          return theQHadrons;
		}
        else
		{
          G4QContent outQC = G4QPDGCode(sPDG).GetQuarkContent();
          G4int outBar     = outQC.GetBaryonNumber();
          G4int outChg     = outQC.GetCharge();
          G4double outCB   = theEnvironment.CoulombBarrier(outChg,outBar); // QCarge is neglected
          // Now the CoulBar reflection should be taken into account
          G4double outT    = s4Mom.e()-s4Mom.m(); 
          G4double outProb = theEnvironment.CoulBarPenProb(outCB,outT,outChg,outBar);
          //outProb=1.;
          if(G4UniformRand()<outProb)
          {
            FillHadronVector(candHadr);              // Fill "new candHadr" (delete equivalent)
            check+= s4Mom;                           // @@ Just for checking
            ccheck+=outChg;                          // @@ Just for checking
            q4Mom = resQ4Mom;                        // Update ResidualQuasmon LorentzVector
            valQ  = curQ;                            // Update the Quark Content of Quasmon
            status= 1;                               // Something was successfuly done
            phot4M=zeroLV;
            piF=false;
            gaF=false;
          }
          else delete candHadr;
		}
      }
      else                                           // ===> "Hadron in NuclMedia or NuclCand" case
	  {
        // Now the CoulBar reflection should be taken into account
        G4QContent outQC = G4QPDGCode(sPDG).GetQuarkContent();
        G4int outBar     = outQC.GetBaryonNumber();
        G4int outChg     = outQC.GetCharge();
        G4double outCB   = theEnvironment.CoulombBarrier(outChg,outBar);
        // Now the Coulomb Potential should be taken into account  //@@ How to take into account ??
        if(nucflag) rQ4Mom+=pMass-pNMass;              // Make a correction to 4Mom of ResidQuasmon
        G4QHadron tmpRQH(valQ+transQC,rQ4Mom);         // Temporary QHadron for the Residual Quasmon
        // Now the CoulBar reflection should be taken into account
        G4double outT    = fr4Mom.e()-fr4Mom.m(); 
        G4double outProb = theEnvironment.CoulBarPenProb(outCB,outT,outChg,outBar);
        if(G4UniformRand()<outProb)
        {
          theEnvironment.Reduce(pPDG);               // Update NuclearEnviron after Q->RQ+F
          G4LorentzVector sumL=theEnvironment.Get4Momentum()+q4Mom; //@@ For Check Print Only
          check += fr4Mom;                           //@@ Just for checking
          ccheck+=G4QPDGCode(sPDG).GetCharge();      //@@ Just for checking
          q4Mom  = rQ4Mom;                           // Update the 4Mom of the Quasmon
          valQ  += transQC;                          // Update the Quark Content of Quasmon
          G4QHadron* candHadr = new G4QHadron(sPDG,fr4Mom);// Createa a Hadron for the Candidate
          FillHadronVector(candHadr);                // Fill the RadiatedHadron(delete equivalent)
#ifdef debug
          G4cout<<"G4Q::HQ:QuarkExchangeHadronizationThroughCB Q="<<valQ<<",trQC="<<transQC<<G4endl;
#endif	  
          sumL-=theEnvironment.Get4Momentum()+q4Mom+fr4Mom;
#ifdef debug
          G4cout<<"G4Q::HQ: status=1, >>>>>>>>> NuclearMatter SUBCHECK >>>>>>>>>>>>>"<<sumL<<G4endl;
#endif	  
          status=1;                                  // Something was successfuly done
          phot4M=zeroLV;
          piF=false;
          gaF=false;
	    }
        else
        {
#ifdef debug
          G4cout<<"G4Q::HQ:CB ISN'T PENETRATED,P="<<outProb<<",T="<<outT<<",M="<<fr4Mom.m()<<G4endl;
#endif	  
          if(envA<20&&G4UniformRand()>envA*envA/400)
          {
            status=3;
            if(gaF)
			{
              phot4M=zeroLV;
              gaF=false;
            }
		  }
		}
      }
	} // End of skip
    // === Check of boundary to escape not existing state of residual nucleus ===
    if(CheckGroundState())
	//if(CheckGroundState(true))
	{
      KillQuasmon();                           // This Quasmon is done
      qEnv=theEnvironment;
      return theQHadrons;      
	}
    G4LorentzVector sumLor=theEnvironment.Get4Momentum()+q4Mom+check;
#ifdef debug
    G4int eZ   = theEnvironment.GetZ();
    G4int sumC = eZ+valQ.GetCharge()+ccheck;
    G4int curPDG=valQ.GetSPDGCode();
    G4cout<<"G4Q::HQ:eZ="<<eZ<<valQ<<">>***>>FinalCHEK>>***>> 4M="<<sumLor<<",Charg="<<sumC<<G4endl;
    if(!curPDG) G4cerr<<"***G4Q::HQ: Quasmon-Tripolino QC="<<valQ<<G4endl;
    G4cout<<"G4Q::HQ:=======> ResidualQ 4M="<<q4Mom<<", QC="<<valQ<<G4endl;
#endif
  }
#ifdef debug
  G4cout<<"G4Q::HQ: Q="<<q4Mom<<valQ<<",E="<<theEnvironment<<", status="<<status<<G4endl;
#endif
  qEnv=theEnvironment;                               // Update the QEnvironment
  return theQHadrons;
} // End of "HadronizeQuasmon"

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
  static const G4LorentzVector zeroLV(0.,0.,0.,0.);
  static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0);
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
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();

  status=1;                                    // Something is going to be filled by Quasmon
  phot4M=zeroLV;
  G4int thePDG      = qH->GetPDGCode();        // Get PDG code of the Hadron to switch
  G4LorentzVector t = qH->Get4Momentum();      // 4-Mom of Chipolino
#ifdef ppdebug
  G4cout<<"G4Quasmon::FillHadronVector:Hadron's PDG="<<thePDG<<",4Mom="<<qH->Get4Momentum()<<G4endl;
#endif
  if(thePDG>80000000 && (thePDG<90000000 || thePDG%1000>500 || thePDG%1000000>500000)
	 && thePDG!=90002999 && thePDG!=89999003 && thePDG!=90003998 && thePDG!=89998004
	 && thePDG!=90003999 && thePDG!=89999004 && thePDG!=90004998 && thePDG!=89998005)
  { // Translation from CHIPS encoding to PDG encoding
    if     (thePDG==90999999) thePDG=-311;  // anti-K0 === Meson OCTET
    else if(thePDG==90999000) thePDG=-321;  // K-
    else if(thePDG==89000001) thePDG=311;   // K0
    else if(thePDG==89001000) thePDG=321;   // K+
    else if(thePDG==90000999) thePDG=211;   // pi+
    else if(thePDG==89999001) thePDG=-211;  // pi-
    else if(thePDG==89999999) thePDG=-2112; // anti-neutron === anti-OCTET(reduced to SEPTUM Lam/S0)
    else if(thePDG==89999000) thePDG=-2212; // anti-proton
    else if(thePDG==89000000) thePDG=-3122; // anti-lambda
    else if(thePDG==88999002) thePDG=-3222; // anti-SIGMA+
    else if(thePDG==89000999) thePDG=-3222; // anti-SIGMA-
    else if(thePDG==88000001) thePDG=-3322; // anti-KSI0
    else if(thePDG==88001000) thePDG=-3312; // anti-KSI-
    else if(thePDG==89999002) thePDG=1114;  // Delta-(resonance) === bary-DECUPLET/OCTET
    else if(thePDG==90001999) thePDG=2224;  // Delta++(resonance) (Delta0&Delta+ a covered by n&p)
    else if(thePDG==90999001) thePDG=3112;  // Sigma-
    else if(thePDG==91000999) thePDG=3222;  // Sigma+ (Sigma0 iz covered by Lambda)
    else if(thePDG==91999000) thePDG=3312;  // Ksi-
    else if(thePDG==91999999) thePDG=3322;  // Ksi0
    else if(thePDG==92998999) thePDG=3112;  // Omega-(resonance)
    else G4cerr<<"*G4Quasmon::FillHadronVector:PDG="<<thePDG<<",M="<<qH->Get4Momentum().m()<<G4endl;
    qH->SetQPDG(G4QPDGCode(thePDG));
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
      G4cerr<<"***FillHV:h1QC="<<h1QC<<"(PDG="<<h1PDG<<"),h2QC="<<h2QC<<"(PDG="<<h2PDG<<")"<<G4endl;
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
      G4cerr<<"***G4Q::FillHadronVector:ChipQC"<<chipQC<<":PDG1="<<h1PDG<<",PDG2="<<h2PDG<<G4endl;
	}
    else
    {
      //qH->SetNFragments(2);
      //theQHadrons.push_back(qH);                  // (delete equivalent)
      // Instead
      delete qH;
      //
      fHadr->Set4Momentum(f4Mom);              // Put the randomized 4Mom to 1-st Hadron
      FillHadronVector(fHadr);                 // Fill 1st Hadron (delete equivalent)
      sHadr->Set4Momentum(s4Mom);              // Put the randomized 4Mom to 2-nd Hadron
      FillHadronVector(sHadr);                 // Fill 2nd Hadron (delete equivalent)
    }
  }
  // pD++,nD-,ppD++,nnD-,D++D++,D-D-,pD++D++,nD-D- "nuclear states"
  else if(thePDG>80000000 && 
          thePDG==90002999&&thePDG==89999003&&thePDG==90003999&&thePDG==89999004&&thePDG==90003998&&
          thePDG==89998004&&thePDG==90004998&&thePDG==89998005 && !DecayOutHadron(qH))
  {
#ifdef pdebug
    G4cout<<"***G4Quasmon::FillHadronVector: Emergency OUTPUT, PDGcode= "<<thePDG<<G4endl;
#endif
  }
  else if(thePDG>80000000&&thePDG!=90000000)   // === Decay-Evaporation of the BarionicFragment ===
  {
    G4double fragMas=qH->GetMass();            // Real Mass of the nuclear fragment
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
      if(fragMas>m1+m2)                        // => "can decay" case
      {
        G4LorentzVector fq4M(0.,0.,0.,m1);
        G4LorentzVector qe4M(0.,0.,0.,m2);
        if(!qH->DecayIn2(fq4M,qe4M))
	    {
          G4cerr<<"***G4Quasm::FillHadrV: QM="<<t.m()<<"-> Mes="<<PDG1<<"(M="
	  	        <<m1<<") + ResA="<<PDG2<<"(M="<<m2<<")"<<G4endl;
	      G4Exception("G4Quasm::FillHadrV: Mes+ResA DecayIn2 did not succeed");
	    }
        //qH->SetNFragments(2);                  // Put a#of Fragments=2
        //theQHadrons.push_back(qH);
        // Instead
        delete qH;
        //
        G4QHadron* H1 = new G4QHadron(PDG1,fq4M);
        theQHadrons.push_back(H1);                // (delete equivalent)
        G4QHadron* H2 = new G4QHadron(PDG2,qe4M);
        FillHadronVector(H2);                  // (delete equivalent)
	  }
      else if(m1+m2-fragMas<0.01)              // Split the 4-momentum
	  {
        G4double r1=m1/fragMas;
        G4double r2=1.-r1;
        //qH->SetNFragments(2);                  // Put a#of Fragments=2
        //theQHadrons.push_back(qH);
        // Instead
        delete qH;
        //
        G4QHadron* H1 = new G4QHadron(PDG1,r1*t);
        theQHadrons.push_back(H1);                // (delete equivalent)
        G4QHadron* H2 = new G4QHadron(PDG2,r2*t);
        FillHadronVector(H2);                  // (delete equivalent)
	  }
      else
	  {
        G4cerr<<"***G4Quasm::FillHadrVec: PDG="<<thePDG<<"("<<t.m()<<") < Mes="<<PDG1
              <<"("<<m1<<") + ResA="<<PDG2<<"("<<m2<<"), d="<<fragMas-m1-m2<<G4endl;
	    G4Exception("G4Quasm::FillHadrVec: mass of decaying hadron is too small");
	  }
	}
    else if(abs(fragMas-GSMass)<.1)            // the Nucleus is too close the Ground State
    {
#ifdef pdebug
	  G4cout<<"G4Quasm::FillHadrVect: Ground state"<<G4endl;
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
            <<nN<<",nZ="<<nZ<<",nL="<<nS<<",totM="<<fragMas<<",n="<<fragMas-nResM-mNeut<<",p="
            <<fragMas-pResM-mProt<<",l="<<fragMas-lResM-mLamb<<G4endl;
#endif
      if(thePDG==90004004||bA>1&&(nN>0&&fragMas>nResM+mNeut||nZ>0&&fragMas>pResM+mProt
                                                           ||nS>0&&fragMas>lResM+mLamb))
	  {
        G4int barPDG = 90002002;
        G4int resPDG = 90002002;
        G4double barM= mAlph;
        G4double resM= mAlph;
		if     (fragMas>nResM+mNeut)           // Can radiate a neutron (priority 1)
		{
          barPDG = 90000001;
          resPDG = nResPDG;
          barM= mNeut;
          resM= nResM;
		}
		else if(fragMas>pResM+mProt)           // Can radiate a proton (priority 2)
		{
          barPDG=90001000;
          resPDG=pResPDG;
          barM  =mProt;
          resM  =pResM;
		}
		else if(fragMas>lResM+mLamb)           // Can radiate a Lambda (priority 3)
		{
          barPDG=91000000;
          resPDG=lResPDG;
          barM  =mLamb;
          resM  =lResM;
		}
        else if(thePDG!=90004004&&fragMas>GSMass)// If it's not Be8 decay in gamma
		{
          barPDG=22;
          resPDG=thePDG;
          barM  =0.;
          resM  =pResM;
		}
        else if(thePDG!=90004004)
		{
          G4cerr<<"***G4Q::FillHadrVect:PDG="<<thePDG<<",M="<<fragMas<<"< GSM="<<GSMass<<G4endl;
          G4Exception("***G4Quasmon::FillHadronVector: Below GSM but can not decay in p or n");
		}
        G4LorentzVector a4Mom(0.,0.,0.,barM);
        G4LorentzVector b4Mom(0.,0.,0.,resM);
        if(!qH->DecayIn2(a4Mom,b4Mom))
        {
          theQHadrons.push_back(qH);              // No decay (delete equivalent)
          G4cerr<<"***G4Q::FillHadronVector: Be8 decay did not succeed"<<G4endl;
	    }
        else
        {
          //qH->SetNFragments(2);                // Fill a#of fragments to decaying Hadron
          //theQHadrons.push_back(qH);              // Fill hadron with nf=2 (delete equivalent)
          // Instead
          delete qH;
          //
          G4QHadron* HadrB = new G4QHadron(barPDG,a4Mom);
          FillHadronVector(HadrB);             // Fill 1st Hadron (delete equivalent)
          G4QHadron* HadrR = new G4QHadron(resPDG,b4Mom);
          FillHadronVector(HadrR);             // Fill 2nd Hadron (delete equivalent)
        }
	  }
      else
      {
#ifdef pdebug
	    G4cout<<"G4Quasm::FillHadrVect: Leave as it is"<<G4endl;
#endif
        theQHadrons.push_back(qH);             // No decay  (delete equivalent)
	  }
    }
    else if (fragMas<GSMass)
	{
      G4cerr<<"***G4Quasm::FillHV: tM="<<fragMas<<" > GSM="<<GSMass<<", d="<<fragMas-GSMass<<G4endl;
      G4Exception("***G4Quasmon::FillHadronVector: Out Mass is below the Ground State value");
	}
    else if (bA==1&&fragMas>GSMass)
	{
      G4int gamPDG=22;
      G4double gamM=0.;
      if(fragMas>mPi0+GSMass)
	  {
        gamPDG=111;
        gamM=mPi0;
	  }
      G4LorentzVector a4Mom(0.,0.,0.,gamM);
      G4LorentzVector b4Mom(0.,0.,0.,GSMass);
      if(!qH->DecayIn2(a4Mom,b4Mom))
      {
        theQHadrons.push_back(qH);              // No decay (delete equivalent)
        G4cerr<<"***G4Q::FillHadronVector: N*->gamma/pi0+N decay did not succeed"<<G4endl;
	  }
      else
      {
        G4QHadron* HadrB = new G4QHadron(gamPDG,a4Mom);
        FillHadronVector(HadrB);             // Fill gamma/Pi0 Hadron (delete equivalent)
        qH->Set4Momentum(b4Mom);
        theQHadrons.push_back(qH);                // Fill corrected baryon in the HadronVector
      }
    }
    else                                       // ===> Evaporation of excited system
	{
#ifdef pdebug
      G4cout<<"G4Quasm::FillHadrVect:Evaporate "<<thePDG<<",tM="<<fragMas<<" > GS="<<GSMass
            <<qNuc.Get4Momentum()<<", m="<<qNuc.Get4Momentum().m()<<G4endl;
#endif
      G4QHadron* bHadron = new G4QHadron;
      G4QHadron* rHadron = new G4QHadron;
      if(!qNuc.EvaporateBaryon(bHadron,rHadron))
	  {
        G4cerr<<"***G4Quasmon::FillHadronVector:Evaporation, PDG="<<thePDG<<",tM="<<fragMas<<G4endl;
        delete bHadron;
        delete rHadron;
        theQHadrons.push_back(qH);                // Fill hadron in the HadronVector as it is
	  }
#ifdef pdebug
      G4cout<<"G4Quasm::FillHadrVec:Done b="<<bHadron->GetQPDG()<<",r="<<rHadron->GetQPDG()<<G4endl;
#endif
      //qH->SetNFragments(2);                    // Fill a#of fragments to decaying Hadron
      //theQHadrons.push_back(qH);                  // Fill hadron with nf=2 (delete equivalent)
      // Instead
      delete qH;
      //
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
} // End of "FillHadronVector"

// Calculate a momentum of quark-parton greater then minimum value kMin
//G4double G4Quasmon::GetQPartonMomentum(G4double mR2, G4double mC2)
G4double G4Quasmon::GetQPartonMomentum(G4double kMax, G4double mC2)
//       ==========================================================
{
  //gives k>kMin QParton Momentum for the current Quasmon
#ifdef pdebug
  //G4cout<<"G4Quasmon::GetQPartonMomentum: ***called*** mR="<<sqrt(mR2)<<",mC="<<sqrt(mC2)<<G4endl;
  G4cout<<"G4Quasmon::GetQPartonMomentum: ***called*** kMax="<<kMax<<",mC="<<sqrt(mC2)<<G4endl;
#endif
  G4double qMass = q4Mom.m();
  G4double kLim  = qMass/2.;               //Kinematikal limit for "k"
  G4double twM   = qMass+qMass;
  ////kMax  = kLim;                            //@@@@@@@@@
  if(kLim<kMax) kMax  = kLim;
  G4double kMin  = mC2/twM;
  //if(mR2!=0.)
  //{
  //  G4double qM2   = qMass*qMass;          // Squared Mass of Quasmon
  //  G4double frM   = twM+twM;              // 4*M_Quasmon
  //  G4double fM2m2 = frM*qMass*mC2;        // 4*MQ**2*mC2=0.
  //  G4double Mmum  = qM2+mC2-mR2;          // QM**2-mR2
  //  G4double Mmum2 = Mmum*Mmum;            // (QM**2-mR2)**2
  //  if(Mmum2<fM2m2)G4Exception("G4Quasmon::QPartMom: mR & mC are too big in comp. with mQ");
  //  G4double sqM   = sqrt(Mmum2-fM2m2);    // QM**2-mR2
  //  kMin  = (Mmum-sqM)/frM;                // kMin=0.
  //  kMax  = (Mmum+sqM)/frM;                // kMax=2*(QM**2-mR2)/4QM
  //}
  if (kMin<0 || kMax<0 || kMax>kLim || qMass<=0. || nOfQ<2)
  {
    G4cerr<<"***G4Q::GetQPMom: kMax="<<kMax<<", kMin="<<kMin<<", kLim="<<kLim<<", MQ="<<qMass
          <<", n="<<nOfQ<<G4endl;
	G4Exception("G4Quasmon::GetQPartonMomentum: Can not generate quark-parton");   
  }
#ifdef pdebug
  G4cout<<"G4Q::GetQPMom: kLim="<<kLim<<",kMin="<<kMin<<",kMax="<<kMax
        <<", nQ="<<nOfQ<<G4endl;
#endif
  if(kMin>kMax||nOfQ==2) return kMax;
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
      G4cout<<"G4Q::GetQPMom: f="<<f<<" is changer to 99"<<G4endl;
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
	G4cout<<"G4Q::GetQPMom:>>>First x="<<x<<", n="<<n<<", f="<<f<<", d/R(first)="<<od/vRndm<<G4endl;
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
	  G4cout<<"G4Q::GetQPMom: Iter#"<<it<<": (c="<<c<<" - R="<<vRndm<<")/R ="<<d/vRndm<<", x="<<x
            <<", f="<<f<<G4endl;
#endif
      od = d;
      ox = x;
      it++;
	};
    if(it>nitMax)
	{
#ifdef pdebug
      G4cout<<"G4Q::GetQPMom: a#of iterations > nitMax="<<nitMax<<G4endl;
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
} // End of "GetQPartonMomentum"

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
  // --- Exponent, Double Split, Poisson 1 ============
  //G4int b = valQ.GetBaryonNumber();
  //G4int mq= 3*b;
  //if (!b) mq=2;
  //G4double mean = ((1.+sqrt(1.+qMOverT*qMOverT))/2. - mq)/2.;
  //if(mean<0.) nOfQ=mq;
  // --- Uncomment up to here ================^^^^^^^^^
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
  G4cout<<"G4Q::Calc#ofQP:QM="<<q4Mom<<qMass<<",T="<<Temperature<<",QC="<<valQ<<",n="<<nOfQ<<G4endl;
#endif
  G4int absb = abs(valQ.GetBaryonNumber());
  G4int tabn = 0;
  if(absb)tabn=3*absb;      // Minimal QC for baryonic system fragmentation
  else if(tabn<4) tabn=4;   // Minimal QC for mesonic system fragmentation
  if (nOfQ<tabn) nOfQ=tabn;
  G4int nSeaPairs = (nOfQ-valc)/2;
  G4int stran = abs(valQ.GetS());
  G4int astra = abs(valQ.GetAS());
  if(astra>stran) stran=astra;
  G4int nMaxStrangeSea=static_cast<int>((qMass-stran*mK0)/(mK0+mK0));//KK is a minimum for s-sea
  if (absb) nMaxStrangeSea=static_cast<int>((qMass-absb)/672.); //LambdaK is a minimum for s-sea
#ifdef pdebug
  G4cout<<"G4Q::Calc#ofQP:"<<valQ<<",INtot="<<valc<<",nOfQ="<<nOfQ<<",SeaPairs="<<nSeaPairs<<G4endl;
#endif
  if (nSeaPairs)            // Add/subtract sea pairs to/from initial quark content
  {
    G4int morDec=0;
    if(nSeaPairs>0)valQ.IncQAQ(nSeaPairs,SSin2Gluons);
    else    morDec=valQ.DecQAQ(-nSeaPairs);
#ifdef pdebug
    if(morDec) G4cout<<"G4Q::Calc#ofQP: "<<morDec<<" pairs can be reduced more"<<G4endl;
#endif
    G4int sSea=valQ.GetS(); // Content of strange quarks
    G4int asSea=valQ.GetAS();
    if(asSea<sSea) sSea=asSea;
    if(sSea>nMaxStrangeSea) // @@@@@@@ Too many strange sea ??????
	{
#ifdef pdebug
      G4cout<<"G4Q::Calc#ofQP:**Reduce** S="<<sSea<<",aS="<<asSea<<",maxS="<<nMaxStrangeSea<<G4endl;
#endif
      sSea-=nMaxStrangeSea; // Strange sea excess
      valQ.DecS(sSea);      // Reduce strange sea to adoptable limit
      valQ.DecAS(sSea);
      valQ.IncQAQ(sSea,0.); // Add notstrange sea ????????
	}
  }
#ifdef pdebug
  G4cout<<"G4Quasmon::Calc#ofQP: *** RESULT IN*** nQ="<<nOfQ<<", FinalQC="<<valQ<<G4endl;
#endif
} // End of "CalculateNumberOfQPartons"

// Modify Candidate masses in nuclear matter and set possibilities
void G4Quasmon::ModifyInMatterCandidates()
//   ======================================
{
  ////////static const G4double mNeut= G4QPDGCode(2112).GetMass();
  G4double envM = theEnvironment.GetMass();      // Mass of the Current Environment
  G4QContent envQC=theEnvironment.GetQCZNS();    // QuarkContent of the current Nuclear Environ.
  G4int eP = theEnvironment.GetZ();              // A#of protons in the Current Environment
  G4int eN = theEnvironment.GetN();              // A#of neutrons in the Current Environment
  G4int eL = theEnvironment.GetS();              // A#of lambdas in the Current Environment
  G4QContent totQC=theEnvironment.GetQC()+valQ;  // Total Quark Comtent of the system
  G4int tP = totQC.GetP();                       // A#of protons in the Current Environment
  G4int tN = totQC.GetN();                       // A#of neutrons in the Current Environment
  G4int tL = totQC.GetL();                       // A#of lambdas in the Current Environment
  G4double totM=G4QNucleus(totQC).GetMZNS();     // Mass of total system
  for (unsigned ind=0; ind<theQCandidates.size(); ind++)
  {
    G4QCandidate* curCand=theQCandidates[ind];   // Pointer to the Candidate
    G4int  cPDG = curCand->GetPDGCode();         // PDGC of the Candidate
    G4bool poss = curCand->GetPossibility();     // Possibility for the Candidate
    G4QContent tmpTQ=totQC-curCand->GetQC();
    G4QNucleus tmpT(tmpTQ);                      // Nucleus for Total Resid. Nucleus for Fragment
    G4double   tmpTM=tmpT.GetMZNS();             // GSMass of Total Residual Nucleus for Fragment
    G4QPDGCode cQPDG(cPDG);                      // QPDG for the candidate
    G4double   frM=cQPDG.GetMass();              // Vacuum mass of the candidate
    if(cPDG>80000000&&cPDG!=90000000)            // Modify Fragments to take into account CurNuclEnv
	{
      if(totMass<tmpTM+frM)
      {
#ifdef sdebug
        G4cout<<"G4Q::MIM:C="<<cPDG<<tmpT<<tmpTM<<"+"<<frM<<"="<<tmpTM+frM<<">tM="<<totMass<<G4endl;
#endif
        curCand->SetPossibility(false);
	  }
      G4QNucleus cNuc(cPDG);                     // Fake nucleus for the Candidate
      G4int cP = cNuc.GetZ();                    // A#of protons in the Current Environment
      G4int cN = cNuc.GetN();                    // A#of neutrons in the Current Environment
      G4int cL = cNuc.GetS();                    // A#of lambdas in the Current Environment
      G4QPDGCode cQPDG(cPDG);                    // QPDG of the Current Cluster
#ifdef debug
      if(cPDG==90001000) G4cout<<"G4Q::MIM:>>>cPDG=90001000<<<, possibility="<<poss<<G4endl;
#endif
      if(eP>=cP&&eN>=cN&&eL>=cL&&poss)           // Cluster exists & possible
      {
        G4double clME = 0.;                      // Prototype of the BoundClusterMass in Environment
        G4double clMN = 0.;                      // Prototype of the BoundClusterMass in TotNucleus
        G4double renvM = 0;                      // Prototype of the residual Environ mass
        if(cP==eP&&cN==eN&&cL==eL)clME=cQPDG.GetMass();// The only not bound Cluster of Environment
        else                                     // Bound Cluster in the Environment
		{
          renvM = cQPDG.GetNuclMass(eP-cP,eN-cN,eL-cL); // Mass of residual for Environment
          clME  = envM-renvM;
        }
        if(cP==tP&&cN==tN&&cL==tL)clMN=cQPDG.GetMass(); // The only not bound Cluster of TotNucleus
        else                                     // Bound Cluster in Total Nucleus
		{
          renvM = cQPDG.GetNuclMass(tP-cP,tN-cN,tL-cL); // Mass of residual for Total Nucleus
          clMN   = totM-renvM;
        }
        curCand->SetParPossibility(true);
        curCand->SetEBMass(clME);
        curCand->SetNBMass(clMN);
#ifdef sdebug
        G4int envPDGC = theEnvironment.GetPDGCode();   // PDG Code of Current Environment
        G4cout<<"G4Q:ModMatC:C="<<cPDG<<cNuc<<clME<<","<<clMN<<",E="<<envPDGC<<",M="<<renvM<<G4endl;
#endif
	  }
      else curCand->SetParPossibility(false);
	} // @@ Modification of hadron masses in nuclear matter are not implemented yet
  }
} // End of "ModifyInMatterCandidates"

// Randomize the Resonance masses and calculate probabilities of hadronization for them
void G4Quasmon::CalculateHadronizationProbabilities(G4double E, G4double kVal, G4double kLS,
                                                    G4bool piF, G4bool gaF)
//   =======================================================================================
{
  G4int    vap = nOfQ-3;                         // Vacuum power
  G4double mQ2 = q4Mom.m2();                     // Squared Mass of the Quasmon
  G4double eQ  = q4Mom.e();                      // LS Energy of the Quasmon
  G4double mQ  = sqrt(mQ2);                      // Mass of the decaying Quasmon
  G4double dk  = kVal + kVal;                    // Double momentu of quark-parton in QCM
  G4double rQ2 = mQ2-dk*mQ;                      // Min Residual Colored Quasmon Squared Mass
  //////////////G4double rQ  = sqrt(rQ2);
  G4double mQk = mQ-dk;                          // For acceleration
  G4double vaf = theEnvironment.GetProbability()*mQk/kVal/vap; //@@ !! Vacuum factor
  //G4double vaf = mQk/kVal/vap;                   //@@ Vacuum factor (instead of in G4QNucleus)
  G4double accumulatedProbability = 0.;
  G4double secondAccumProbability = 0.;
  G4int qBar =valQ.GetBaryonNumber();            // BaryNum of Quasmon
  G4int absb = abs(qBar);                        // Abs BaryNum of Quasmon
  G4int nofU = valQ.GetU()- valQ.GetAU();        // A#of u-quarks
  G4int dofU = nofU+nofU;
  G4int nofD = valQ.GetD()- valQ.GetAD();        // A#of d-quarks
  G4int dofD = nofD+nofD;
  G4int qChg = valQ.GetCharge();
  G4int qIso = qBar-qChg-qChg;                   // Charge of Quasmon
  ////////////////G4int aQuI = abs(qIso);                        // Abs value for estimate
  ///////////////G4int oeQB = qBar%2;                           // odd(1)/even(0) QBaryn flag
  G4int maxC = theQCandidates.size();         // A#of candidates
  ///////////////G4double totZ = theEnvironment.GetZ() + valQ.GetCharge();       // Z of the Nucleus
  ///////////////G4double totA = theEnvironment.GetA() + valQ.GetBaryonNumber(); // A of the Nucleus
  G4double envM = theEnvironment.GetMass();      // Mass of the Current Environment
  G4int envPDGC = theEnvironment.GetPDGCode();   // PDG Code of Current Environment
  G4int envN    = theEnvironment.GetN();         // N of current Nuclear Environment
  G4int envZ    = theEnvironment.GetZ();         // Z of current Nuclear Environment
  ////////////////G4int envS    = theEnvironment.GetS();         // S of current Nuclear Environment
  G4int envA    = theEnvironment.GetA();         // A of current Nuclear Environment
  G4QContent envQC=theEnvironment.GetQCZNS();    // QuarkContent of the current Nuclear Environ.
  ///////G4double addPhoton=phot4M.e();                 // Photon energy for capture by quark-parton
#ifdef sdebug
  G4int maxB = theEnvironment.GetMaxClust();     // Maximum BaryNum for clusters
  G4cout<<"G4Q::CalculateHadronizationProbabilities:Q="<<mQ<<valQ<<",v="<<vaf<<",mC="<<maxB<<G4endl;
#endif
  // ================= Calculate probabilities for candidates
  for (unsigned index=0; index<theQCandidates.size(); index++)
  {
    G4QCandidate* curCand=theQCandidates[index];
    curCand->ClearParClustVector();             // Clear Parent Cluster Vector for the Fragment
    G4double probability = 0.;
    G4double secondProbab = 0.;
    G4int    resPDG=0;
  	G4double comb=0.;
    G4int cPDG = curCand->GetPDGCode();
    G4QContent candQC = curCand->GetQC();
    G4QContent tmpTQ=envQC+valQ-candQC;              // QC for Total ResidualNucleus for the Cluster
    G4QNucleus tmpT(tmpTQ);                          // Nucleus for Total ResidNucleus for Fragment
    G4double   tmpTM=tmpT.GetMZNS();                 // GSMass of Total ResidualNucleus for Fragment
    G4QPDGCode cQPDG(cPDG);                          // QPDG for the candidate
    G4double   frM=cQPDG.GetMass();                  // Vacuum mass of the candidate
    G4int cU=candQC.GetU()-candQC.GetAU();
    ////////////G4int dU=cU+cU;
    G4int cD=candQC.GetD()-candQC.GetAD();
    ////////////G4int dD=cD+cD;
    G4int dUD=abs(cU-cD);
    ////////////G4int cS=candQC.GetS()-candQC.GetAS();
    G4bool pos=curCand->GetPossibility()&&totMass>tmpTM+frM;
    //G4bool pos=curCand->GetPossibility();
#ifdef pdebug
    if(cPDG==111||cPDG==90001000||cPDG==90000001||cPDG==90000002||cPDG==90001001||cPDG==90002002)
	  G4cout<<"G4Q::CHP:==*****==>>>c="<<cPDG<<",dUD="<<dUD<<",pos="<<pos<<",eA="<<envA<<G4endl;
#endif
	if(pos&&(cPDG<80000000||(cPDG>80000000&&cPDG!=90000000&&dUD<2))) // 2 (OK) or 3
	{
      G4QContent curQ = valQ;                   // Make a current copy of the Quasmon QuarkCont
      G4int baryn= candQC.GetBaryonNumber();    // Baryon number of the Candidate
      G4int cC   = candQC.GetCharge();          // Charge of the Candidate
      G4double CB=0.;
      if(envA) CB=theEnvironment.CoulombBarrier(cC,baryn);
      /////////////G4int cI   = baryn-cC-cC;
      G4int cNQ= candQC.GetTot()-2;           // A#of quark-partons in the candidate - 2 (OK)
      //G4int cNQ= candQC.GetTot()+baryn-2;     // A#of q-partons + string-junktions - 2 (q-link)
      //G4int cNQ= candQC.GetTot()*2-4;         // A#of q-partons + string-junktions - 2 (bag-link?)
      G4double resM=0.;                         // Prototype for minMass of residual hadron
      if(cPDG>80000000&&cPDG!=90000000&&baryn<=envA) //===> "Nuclear Fragment" case (QUark EXchange)
      {
        G4int      pc=0;                        // Parents counter
        G4double   pcomb=0.;                    // Summed probability of parent clusters
        G4double   frM2=frM*frM;
        //////////G4double   qMax=frM+CB-kLS;
        //////////G4double   qM2=qMax+qMax;
        G4int iQmin=0;
        G4int iQmax=3;
        G4int oQmin=0;
        G4int oQmax=3;
		if     (dofU<=nofD) iQmax=1;            // d-quark (u-quark is forbiden)
        else if(dofD<=nofU) iQmin=1;            // u-quark (d-quark is forbiden)
        if(piF)                                 // pi- transfer it's charge to the quark
		{
          iQmin=0;
          iQmax=1;
		}
#ifdef pdebug
        if(baryn<5)
          G4cout<<"G4Q::CHP:***>>>c="<<cPDG<<",m="<<frM<<",imi="<<iQmin<<",ima="<<iQmax<<G4endl;
#endif
        if(oQmax>oQmin) for(int iq=iQmin; iq<iQmax; iq++) // Entering (from Quasmon) quark
		{
          G4double qFact=1.;
          //if(iq==1&&addPhoton>0.) qFact=4.;
          if(iq==1&&gaF)
          {
            qFact=4.;
#ifdef pdebug
            if(cPDG==90001000) G4cout<<"G4Q::CHP: proton gaF enhanced"<<G4endl;
#endif
		  }
          G4double nqInQ=0.;                    // A#of quarks of this sort in a Quasmon
          if     (!iq)   nqInQ=valQ.GetD();
          else if(iq==1) nqInQ=valQ.GetU();
          else if(iq==2) nqInQ=valQ.GetS();
          comb=0.;                              // Local summ for the i-quark of the Quasmon
#ifdef sdebug
          G4cout<<"G4Q::CHP:i="<<iq<<",cU="<<cU<<",cD="<<cD<<",omi="<<oQmin<<",oma="<<oQmax<<G4endl;
#endif
          if(oQmax>oQmin) for(int oq=oQmin; oq<oQmax; oq++) // Exiting (to Quasmon) quark
		  {
            G4int shift= cQPDG.GetRelCrossIndex(iq, oq);
            G4QContent ioQC=cQPDG.GetExQContent(iq, oq);
            G4QContent resQC=valQ+ioQC;         // Quark Content of the residual Quasmon
#ifdef sdebug
			G4cout<<"G4Q::CHP:iq="<<iq<<",oq="<<oq<<",QC="<<ioQC<<",rQC="<<resQC<<G4endl;
#endif
            G4QPDGCode resQPDG(resQC);          // QPDG of the residual Quasmon
            resPDG=resQPDG.GetPDGCode();        // PDG Code of the residual Quasmon
            G4int resQ=resQPDG.GetQCode();      // Q Code of the residual Quasmon
#ifdef pdebug
            if(baryn<5) G4cout<<"G4Q::CHP:i="<<iq<<",o="<<oq<<ioQC<<",shift="<<shift<<",cQPDG="
                              <<cQPDG<<", residQC="<<resQC<<resQPDG<<G4endl;
#endif
            G4int resD=resQC.GetD()-resQC.GetAD();
            G4int resU=resQC.GetU()-resQC.GetAU();
            G4int resS=resQC.GetS()-resQC.GetAS();
            G4int resA=resQC.GetBaryonNumber();
            G4bool rI= resA>0 && resU>=0 && resD>=0 && (resU+resS>resD+resD || resD+resS>resU+resU);
            //if(resQ>-2&&resPDG&&resPDG!=10&&!rI)// The Residual Quasmon is possible
            //if(resQ>-2&&resPDG&&resPDG!=10&&!rI&&(!piF||cPDG==90000001)&&(!gaF||baryn==1))
			if(resQ>-2&&resPDG&&resPDG!=10&&!rI&&(!piF||cPDG==90000001))
			//if(resQ>-2&&resPDG&&resPDG!=10&&!rI&&(!piF||cPDG==90000001||cPDG==90001002))
			//if(resQ>-2&&resPDG&&resPDG!=10&&!rI&&(!gaF||baryn==1))
			{
              G4int is=index+shift;
              if(shift!=7&&is<maxC)                      // This quark exchange is possible
			  {
                G4QCandidate* parCand=theQCandidates[is];// Pointer to ParentCluster of Cand.
                G4QContent parQC = parCand->GetQC();     // Quark Content of the Parent Cluster
                G4int barot = parQC.GetBaryonNumber();   // Baryon Number of the Parent Cluster
                G4int charge= parQC.GetCharge();         // Charge of the Parent Cluster
                G4int possib=parCand->GetParPossibility();
#ifdef pdebug
                if(baryn<5) G4cout<<"G4Q::CHP:parentPossibility="<<possib<<",pZ="<<charge<<",pN="
                                  <<barot-charge<<", cPDG="<<cPDG<<G4endl;
#endif
                if(possib&&charge<=envZ&&barot-charge<=envN)
                {
                  G4QContent rQQC = valQ+ioQC;           // Quark Content of Residual Quasmon
                  ///////////G4int rQU=rQQC.GetU()-rQQC.GetAU();
                  ///////////G4int rQD=rQQC.GetD()-rQQC.GetAD();
                  G4int isos  = barot-charge-charge;     // Isospin of the Parent Cluster
                  G4double pUD= pow(2.,abs(isos));
                  if(barot!=baryn) G4cerr<<"**G4Quasmon::CalcHadrProb:cand="<<candQC<<",par="<<parQC
                                         <<",shift="<<shift<<",ind="<<index<<",is="<<is<<G4endl;
                  G4int    dI=qIso-isos;                 // IsotopicShiftDifference for ParC & Quasm
                  G4int    dC=cC-charge;                 // ChargeDifference for outFragm & ParentCl
                  G4int    dS=dI+dC;                     // Isotop Symmetry Parameter
#ifdef pdebug
                  if(baryn<5) G4cout<<"G4Q::CHP: dI="<<dI<<", dC="<<dC<<", I="<<qIso<<",i="<<isos
						            <<", C="<<cC<<",c="<<charge<<G4endl;
#endif
                  //********* ISOTOPIC  FOCUSING  *******************
                  // ==== Old (First Publication) Complicated rule ====
				  //if(
				  //  //zZ<3 &&
                  //  //(
				  //  abs(dI)<1 ||
				  //  (barot==1 && (
                  //     abs(dI)<2&&abs(cC-charge)<2 ||
				  //     (dI>=2&&cC<charge)   || (dI<=-2&&cC>charge)
                  //     //dI==2&&cC<=charge || dI==-2&&cC>=charge ||
                  //     //dI>2&&cC<charge   || dI<-2&&cC>charge
				  //   )) || 
				  //   (barot>1&&barot<3
				  //     && (
				  //     abs(dI)<2&&abs(cC-charge)<2 ||
                  //     //dI>=2&&cC<charge   || dI<=-2&&cC>charge
				  //     dI<=2&&cC<=charge || dI==-2&&cC>=charge ||
				  //     dI>2&&cC<charge   || dI<-2&&cC>charge
				  //   )) ||
				  //   (barot>2&&barot<4
				  //    && (
				  //     abs(dI)<2&&abs(cC-charge)<2 ||
                  //     //dI>=2&&cC<charge   || dI<=-2&&cC>charge
				  //     dI<=2&&cC<=charge || dI==-3&&cC>=charge ||
				  //     dI>2&&cC<charge   || dI<-3&&cC>charge
				  //   )) ||
				  //   (barot>3
				  //    && (
				  //     abs(dI)<2&&abs(cC-charge)<2 ||
                  //     dI>=2&&cC<charge   || dI<=-2&&cC>charge
				  //     //dI<=2&&cC<=charge || dI==-3&&cC>=charge ||
				  //     //dI>2&&cC<charge   || dI<-3&&cC>charge
				  //    )
                  //   )
				  //  )
                  // ==== Just a coridor =======
				  if(abs(dS)<3 || (qIso>0&&dC<0||qIso<0&&dC>0) && baryn==1)
				  //if(abs(dS)<3 || qIso>0&&dC<0 || qIso<0&&dC>0)
				  //if(abs(dS)<4 || qIso>0&&dC<0 || qIso<0&&dC>0)
                  //if(3>2)                                   //***>>> ***NO*** Isotope Focusing ***
                  {
                    G4double pPP=parCand->GetPreProbability();// Probability of Parent Cluster
                    //G4double pPP=parCand->GetDenseProbability();// Probability of Parent Cluster
                    G4int    parPDG=parCand->GetPDGCode();    // PDGCode of the Parent Clucter
#ifdef pdebug
                    if(baryn<5) G4cout<<"G4Q::CHP:pPDG="<<parPDG<<",pPP="<<pPP<<G4endl;
#endif
                    G4double boundM=parCand->GetEBMass();     // EnvironBound mass of Parent Cluster
                    G4double nucBM =parCand->GetNBMass();     // TotNuclBound mass of Parent Cluster
#ifdef pdebug
                    if(baryn<5) G4cout<<"G4Q::CHP:c="<<cPDG<<",p="<<parPDG<<",bM="<<boundM
                                      <<",i="<<is<<",adr="<<parCand<<G4endl;
#endif
                    // Kinematical analysis of decay possibility
                    G4double   minM  =0.;                     // Prototype of minM of ResidQuasm
                    if (resPDG==10)minM=G4QChipolino(resQC).GetMass();// ResidualQuasmonChipo
                    else if(resPDG)minM=G4QPDGCode(resPDG).GetMass(); // ResidualQuasmonHadron
                    G4double bNM2=nucBM*nucBM;
                    G4double nDelta=(frM2-bNM2)/(nucBM+nucBM);
#ifdef pdebug
                    G4int    iniPDG =valQ.GetSPDGCode();
                    G4double iniQM = G4QPDGCode(iniPDG).GetMass(); // Not bounded quasmon GS mass
                    G4double freeE = (mQ-iniQM)*iniQM;
                    G4double kCut=boundM/2.+freeE/(iniQM+boundM);
                    if(baryn<5)G4cout<<"G4Q::CHP:rPDG="<<resPDG<<",miM="<<minM<<",kLS="<<kLS<<"<kM="
                                     <<kCut<<",E="<<E<<">nDel="<<nDelta<<",parPDG="<<parPDG<<G4endl;
#endif
					if(resPDG&&minM>0.) // Kinematical analysis of hadronization
                    {
#ifdef debug
                      if(baryn<5)G4cout<<"G4Q::CHP: fM="<<frM<<",bM="<<boundM<<",rM="<<tmpTM<<",tM="
                                       <<totMass<<G4endl;
#endif
                      //G4double pmk=rMo*boundM/kLS;
                      G4double pmk=rMo*nucBM/kLS;
                      G4double bM2=boundM*boundM;
                      G4double eDelta=(frM2-bM2)/(boundM+boundM);
                      G4double ked =kLS-eDelta;
                      G4double dked=ked+ked;
                      G4double dkedC=dked-CB-CB;
                      G4double kd =kLS-nDelta; // For the total Nucleus (including Quasmon)
                      G4double dkd=kd+kd;
                      //////////G4double dkdC=dkd-CB-CB;
                      G4double dkLS=kLS+kLS;
                      //////////G4double Em=(E-eDelta)*(1.-frM/totMass);
                      //G4double Em=(E-nDelta)*(1.-frM/totMass);
                      //G4double Em=(E-nDelta-CB)*(1.-frM/totMass);
                      // *** START LIMITS ***
                      //G4double ne=1.-dked/(boundM+dkLS); // q_min=DEFOULT=bndM*(k-delta)/(bndM+2k)
                      //G4double ne=1.-dkd/(nucBM+dkLS); // q_min=DEFOULT=bndM*(k-delta)/(bndM+2k)
                      G4double ne=1.-dkedC/(boundM+dkLS); // q_min=DEFOULT=bM*(k-delta-CB)/(bM+2k)
                      //G4double ne=1.-dkdC/(nucBM+dkLS); // q_min=DEFOULT=bndM*(k-del-CB)/(bndM+2k)
#ifdef pdebug
                      if(baryn<5) G4cout<<"G4Q::CHP:<qi_DEF>="<<ne<<",dked="<<dked<<",dkLS="<<dkLS
                                        <<G4endl;
#endif
                      G4double kf=1.;
                      if(ne>0.&&ne<1.)kf=pow(ne,cNQ);
                      // == Prepare for the residual nucleus restriction ==
                      // *** LIM ***
                      G4QContent rtQC=valQ+ioQC; // Total Residual Quark Content
                      if(envA-barot<=bEn) rtQC+=envQC-parQC; // Total Residual Quark Content
                      else rtQC+=bEnQC;
                      G4QNucleus rtN(rtQC);               // Create PseudoNucleus for totResid.
                      G4double rtM=rtN.GetGSMass();       // MinMass of residual Q+(Env-ParC) system
                      G4double rtEP=rEP;                  // E+Mom of the real TotColouredResidSys
                      // *** LIM ***
                      if(envA-barot<=bEn) rtEP+=envM-boundM;
                      else rtEP+=mbEn;
                      G4double rtE=rtEP-rMo;              // Energy of the real TotColouredResidSys
                      ////////////G4double rtEMP=rtE-rMo;
#ifdef pdebug
                      G4QContent tmpEQ=envQC-parQC;       // Quark content for Residual Environment
                      G4QNucleus tmpN(tmpEQ);             // Pseudo nucleus for Residual Environment
                      G4double   tmpNM=tmpN.GetMZNS();    // Mass of Residual Environment
                      if(baryn<5)G4cout<<"G4Q::CHP:RN="<<tmpNM<<"="<<envM-boundM<<"=eM="<<envM
                               <<"-bM="<<boundM<<",E="<<rtE<<",eQC="<<envQC<<",pQC="<<parQC<<G4endl;
#endif
                      G4double mintM2=rtM*rtM+.1;         // Mass of the Mini TotColouredResidualSys
                      G4double rtQ2=rtE*rtE-rMo*rMo;      // SquaredMinMass of Residual Q+(Env-ParC)
                      G4double minM2=minM*minM+.1;
                      // ***VRQ***
                      G4double minBM=minM;
                      if(envM>boundM)
                      {
                        minBM=rtM;
                        // *** LIM ***
                        if(envA-barot<=bEn) minBM-=envM-boundM; // MinResidualBoundedQuasmonMass
                        else minBM-=mbEn;
					  }
                      G4double minBM2=minBM*minBM+.1;
#ifdef pdebug
                      if(baryn<5)G4cout<<"G4Q::CHP:M2="<<minM2<<",R2="<<rQ2<<",N2="<<mintM2<<",RN2="
                             <<rtQ2<<",q="<<(minM2-rQ2)/rEP/2<<",qN="<<(mintM2-rtQ2)/rtEP/2<<G4endl;
#endif
					  G4double newh=1.;
                      // == (@@) Historical additional cuts for q_min ===
                      //G4double nc=1.-(dkLS-E-E)/boundM;   // q_min=k-E
                      ////G4double nc=1.-(dkLS-E-E+CB+CB)/boundM;   // q_min=k-E+CB
#ifdef pdebug
					  G4double newl=0.;
                      G4double ph=kf;                     // Just for printing
                      //if(baryn<5)
	                  //G4cout<<"G4Q::CHP:qi_k-E="<<nc<<",k="<<kLS<<",E="<<E<<",M="<<boundM<<G4endl;
#endif
                      //if(nc>0.&&nc<1.&&nc<ne)
                      //{
                      //  ne=nc;
                      //  newh=pow(nc,cNQ);
                      //  if(newh<kf) kf=newh;
                      //}
                      //else if(nc<=0.)kf=0.;

                      //G4double nk=1.-(dkd-Em-Em)/boundM;  // q_min=(k-delta)-E*(M-m)/M
#ifdef pdebug
                      //if(baryn<5)
	                  //G4cout<<"G4Q::CHP:qi_R="<<nk<<",kd="<<kd<<",E="<<Em<<",M="<<totMass<<G4endl;
#endif
                      //if(nk>0.&&nk<1.&&nk<ne)
                      //{
                      //  ne=nk;
                      //  newh=pow(nk,cNQ);
                      //  if(newh<kf) kf=newh;
                      //}
                      //else if(nk<=0.)kf=0.;

                      //G4double mex=frM+Em;
                      //G4double sr=sqrt(mex*mex-frM2);     // q_min=k-sqrt((m+E*(M-m)/M)^2-m^2)
                      //G4double np=1.-(dkLS-sr-sr)/boundM;
#ifdef pdebug
                      //if(baryn<5)
	                  //G4cout<<"G4Q::CHP:qi_k-sr="<<np<<",sr="<<sr<<",m="<<mex<<",M="<<frM<<G4endl;
#endif
                      //if(np>0.&&np<1.&&np<ne)
                      //{
                      //  ne=np;
                      //  newh=pow(np,cNQ);
                      //  if(newh<kf) kf=newh;
                      //}
                      //else if(np<=0.)kf=0.;

                      ////G4double mix=boundM+E;              // The same don't change ***
                      //G4double mix=nucBM+E;
                      ////G4double mix=boundM+E-CB;
                      ////G4double mix=nucBM+E-CB;
                      //G4double st=sqrt(mix*mix-frM2);
                      //G4double nq=1.-(dkLS-st-st)/boundM; // q_min=k-sqrt((m+E*(M-m)/M)^2-m^2)
#ifdef pdebug
                      //if(baryn<5)
	                  //G4cout<<"G4Q::CHP:qi_k-st="<<nq<<",st="<<st<<",m="<<mix<<",M="<<frM<<G4endl;
#endif
                      //if(nq>0.&&nq<1.&&nq<ne)
                      //{
                      //  ne=nq;
                      //  newh=pow(nq,cNQ);
                      //  if(newh<kf) kf=newh;
                      //}
                      //else if(nq<=0.)kf=0.;
                      // == This is the best for the residual nucleus cut (@@ can be improved)
                      G4bool atrest=(eQ-mQ)/mQ<.001;        // Condition of Q at rest (pion capture)
                      // ***VTN*** CHECK for the Total Nucleus
					  if(mintM2>rtQ2)            // ==> Check of residual Total Nucleus ***Always***
					  //if(2>3)
                      {
                        G4double nz=0.;
						//nz=1.-(mintM2-rtQ2)/(boundM*rtEP);
						//if(atrest) nz=1.-(mintM2-rtQ2+pmk*dked)/(boundM*(rtEP+pmk));
						nz=1.-(mintM2-rtQ2)/(nucBM*rtEP);
						if(atrest) nz=1.-(mintM2-rtQ2+pmk*dkd)/(nucBM*(rtEP+pmk));
#ifdef pdebug
                        if(baryn<5)
	                    G4cout<<"G4Q::CHP:q="<<nz<<",a="<<atrest<<",M2="<<mintM2<<">"<<rtQ2<<G4endl;
#endif
                        if(nz>0.&&nz<1.&&nz<ne)
                        {
                          ne=nz;
                          newh=pow(nz,cNQ);
                          if(newh<kf) kf=newh;
                        }
                        else if(nz<=0.)kf=0.;
                      }
                      // *** VRQ *** CHECK for the virtual Residual Quazmon
					  //if(minBM2>rQ2&&!piF&&!gaF&&baryn>3) // ==> Check of Residual Virtual Quasmon
					  if(minBM2>rQ2&&!piF&&!gaF&&baryn>2) // ==> Check of Residual Virtual Quasmon
					  //if(minBM2>rQ2&&baryn>2) // ==> Check of Residual Virtual Quasmon
					  //if(minBM2>rQ2&&!piF&&baryn>3) // ==> Check of Residual Virtual Quasmon
					  //if(minBM2>rQ2&&!piF&&baryn>2) // ==> Check of Residual Virtual Quasmon
					  //if(minBM2>rQ2&&!piF&&!gaF)    // ==> Check of Residual Virtual Quasmon
					  //if(minBM2>rQ2&&!piF)          // ==> Check of Residual Virtual Quasmon
					  //if(minBM2>rQ2)                // ==> Check of Residual (Virtual?) Quasmon
					  //if(2>3)
                      {
                        G4double nz=0;
                        nz=1.-(minBM2-rQ2)/(boundM*rEP);
                        if(atrest) nz=1.-(minBM2-rQ2+pmk*dked)/(boundM*(rEP+pmk));
                        //nz=1.-(minBM2-rQ2)/(nucBM*rEP);
                        //if(atrest) nz=1.-(minBM2-rQ2+pmk*dkd)/(nucBM*(rEP+pmk));
#ifdef pdebug		  
                        if(baryn<5)
	                    G4cout<<"G4Q::CHP:q="<<nz<<",a="<<atrest<<",VQM2="<<minM2<<">"<<rQ2<<G4endl;
#endif				  
                        if(nz>0.&&nz<1.&&nz<ne)
                        {
                          ne=nz;
                          newh=pow(nz,cNQ);
                          if(newh<kf) kf=newh;
                        }
                        else if(nz<=0.)kf=0.;
                      }
					  //if(minM2>rQ2&&baryn>3)       // ==> Check of Residual Quasmon
					  if(minM2>rQ2&&!piF&&!gaF&&baryn>3)       // ==> Check of Residual Quasmon
					  //if(minM2>rQ2&&!piF&&baryn>2)       // ==> Check of Residual Quasmon
					  //if(minM2>rQ2&&!piF&&baryn>3)       // ==> Check of Residual Quasmon
					  //if(minM2>rQ2&&!piF&&!gaF)       // ==> Check of Residual Quasmon
					  //if(minM2>rQ2&&!piF)       // ==> Check of Residual Quasmon
					  //if(minM2>rQ2)            // ==> Check of Residual (Virtual?) Quasmon
					  //if(2>3)
                      {
                        G4double nz=0;
                        nz=1.-(minM2-rQ2)/(boundM*rEP);
                        if(atrest) nz=1.-(minM2-rQ2+pmk*dked)/(boundM*(rEP+pmk));
                        //nz=1.-(minM2-rQ2)/(nucBM*rEP);
                        //if(atrest) nz=1.-(minM2-rQ2+pmk*dkd)/(nucBM*(rEP+pmk));
#ifdef pdebug
                        if(baryn<5)
	                    G4cout<<"G4Q::CHP:q="<<nz<<",a="<<atrest<<",VQM2="<<minM2<<">"<<rQ2<<G4endl;
#endif
                        if(nz>0.&&nz<1.&&nz<ne)
                        {
                          ne=nz;
                          newh=pow(nz,cNQ);
                          if(newh<kf) kf=newh;
                        }
                        else if(nz<=0.)kf=0.;
                      }
                      if(kf<0.)kf=0.;
                      if(kf>1.)kf=1.;
                      G4double high = kf;
#ifdef pdebug
                      if(baryn<5)
						G4cout<<"G4Q::CHP:kf="<<kf<<",mM2="<<minM2<<",rQ2="<<rQ2<<G4endl;
#endif
                      //G4double lz=1.-dked/boundM;            // q_max=DEFAULT=k-delta
                      //G4double lz=1.-dkd/nucBM;              // q_max=DEFAULT=k-delta
                      G4double lz=1.-dkedC/boundM;           // q_max=DEFAULT=k-delta-CB
                      //G4double lz=1.-dkdC/nucBM;             // q_max=DEFAULT=k-delta-CB
#ifdef pdebug
                      if(baryn<5)
	                  G4cout<<"G4Q::CHP:<qa_DEF>="<<lz<<", eDel="<<eDelta<<",nDel="<<nDelta<<G4endl;
#endif
                      G4double low=0.;
                      if(lz>0.&&lz<1.)low=pow(lz,cNQ);
                      else if(lz>=1.)low=1.;
                      // == (@@) Historical additional cuts for q_max ===
                      //G4double tms=kLS+nDelta+Em;
                      ////G4double tms=kLS+eDelta+Em;           // The same don't change ***
                      //G4double le=1.-(tms+tms)/boundM;      // q_max=k+delta+E*(M-m)/M
#ifdef pdebug
                      G4double pl=low;                      // Just for printing
                      //if(baryn<5)
	                  //G4cout<<"G4Q::CHP:qa_t="<<le<<",k="<<kLS<<",E="<<Em<<",bM="<<boundM<<G4endl;
#endif
                      //if(le>0.&&le<1.&&le>lz)
                      //{
                      //  lz=le;
                      //  newl=pow(le,cNQ);
                      //  if(newl>low) low=newl;
                      //}
                      //else if(le>=1.)low=1.;

                      //G4double lk=1.-(dkLS+E+E)/boundM;    // q_max=k+E
                      ////G4double lk=1.-(dkLS+E+E-CB-CB)/boundM; //q_max=k+E-CB (surface condition)
#ifdef pdebug
                      //if(baryn<5)
	                  //G4cout<<"G4Q::CHP:qa_k+E="<<lk<<",k="<<kLS<<",E="<<E<<",M="<<boundM<<G4endl;
#endif
                      //if(lk>0.&&lk<1.&&lk>lz)
                      //{
                      //  lz=lk;
                      //  newl=pow(lk,cNQ);
                      //  if(newl>low) low=newl;
                      //}
                      //else if(lk>=1.)low=1.;

                      //G4double lq=1.-(dkLS+st+st)/boundM;  // q_max=k+sqrt((E*(M-m)/M)^2-m^2)
#ifdef pdebug
                      //if(baryn<5)
	                  //G4cout<<"G4Q::CHP:qa_k+st="<<lq<<",st="<<st<<",m="<<mix<<",M="<<frM<<G4endl;
#endif
                      //if(lq>0.&&lq<1.&&lq>lz)
                      //{
                      //  lz=lq;
                      //  newl=pow(lq,cNQ);
                      //  if(newl>low) low=newl;
                      //}
                      //else if(lq>=1.)low=1.;

                      //G4double lp=1.-(dkLS+sr+sr)/boundM;  // q_max=k+sqrt((E*(M-m)/M)^2-m^2)
#ifdef pdebug
                      //if(baryn<5)
	                  //G4cout<<"G4Q::CHP:qa_k+sr="<<lp<<",sr="<<sr<<",m="<<mex<<",M="<<frM<<G4endl;
#endif
                      //if(lp>0.&&lp<1.&&lp>lz)
                      //{
                      //  lz=lp;
                      //  newl=pow(lp,cNQ);
                      //  if(newl>low) low=newl;
                      //}
                      //else if(lp>=1.)low=1.;
                      // It's a specific Coulomb Barrier limit for charged particles (can be skiped)
                      //if(totZ>cC)                         // ==> Check of Coulomb Barrier
                      //{
                      //  G4double qmaCB=boundM-qMax;
                      //  //G4double qmaCB=nucBM-qMax;
                      //  G4double nz=1.-(qmaCB+qmaCB)/boundM;// q_max=Mb-Mf-CB+kLS, qMax=Mf+CB-kLS
#ifdef pdebug
                      //  if(baryn<5)
	                  //  G4cout<<"G4Q::CHP:<qa_CB>="<<nz<<",qMax="<<qmaCB<<",CB="<<CB<<G4endl;
#endif
                      //  if(nz>0.&&nz>lz)
                      //  {
                      //    newl=pow(nz,cNQ);
                      //    if(newl>low) low=newl;
                      //  }
                      //  else if(nz>1.) low=10.;
                      //}
                      kf-=low;
#ifdef pdebug
                      if(baryn<5)
                      G4cout<<"G4Q::CHP:>>"<<cPDG<<",l="<<low<<",h="<<high<<",ol="<<pl<<",oh="
                            <<ph<<",nl="<<newl<<",nh="<<newh<<",kf="<<kf<<",d="<<eDelta<<G4endl;
#endif
                      G4double probab=0.;
                      if(kf>0)
                      {
                        kf*=boundM/kLS/cNQ;      // Final value of kinematical (i,o) factor
                        G4int noc=cQPDG.GetNumOfComb(iq, oq);
                        probab=qFact*kf*nqInQ*pPP*noc/pUD; // The main with wing suppresion
                        // qFact - squared charde for photons & u-quark, for others it is 1
                        // kf    - the phase space integral
                        // nqInQ - a#of i-quarks in the Quasmon
                        // pPP   - probability to find (a#of) the Parent Cluster
                        // noc   - a#of o-quarks in the Parent Cluster
#ifdef sdebug
                        G4cout<<"G4Q::CHP:p="<<probab<<",qF="<<qFact<<",iq="<<iq<<",oq="<<oq<<",Ph="
                              <<phot4M<<G4endl;
#endif
                        if(probab<0.) probab=0.;
                      }
                      pcomb += probab;           // Update integrated probability for ParntClust
                      G4QParentCluster* curParC = new G4QParentCluster(parPDG,pcomb);
                      curParC->SetTransQC(ioQC); // Keep Quark Content of the Exchange Meson
                      curParC->SetLow(low);      // Keep the Low limit of randomization
                      curParC->SetHigh(high);    // Keep the High limit of randomization
                      curParC->SetEBMass(boundM);// Keep Environ BoundedMass for future calculations
                      curParC->SetNBMass(nucBM); // Keep totNucl BoundedMass for future calculations
                      curParC->SetEBind(eDelta); // Keep EnvBindingEnerergy for future calculations
                      curParC->SetNBind(nDelta); // Keep NucBindingEnerergy for future calculations
                      curParC->SetNQPart2(cNQ);  // Keep a#of quark-partons in the fragment
#ifdef sdebug
		 	          G4cout<<"G4Q::CalcHP: FillParentClaster="<<*curParC<<G4endl;
#endif
                      curCand->FillPClustVec(curParC);//Fill possible ParentClust to FragmVector
                      delete curParC;
                      comb += probab;
#ifdef sdebug
		 	          G4cout<<"G4Q::CHP:i="<<index<<",cC="<<cPDG<<",pc"<<pc<<parQC<<",E="
						    <<theEnvironment<<comb<<",p="<<parCand->GetParPossibility()<<G4endl;
#endif
                      pc++;
					}
				  }
#ifdef sdebug
				  else G4cout<<"***G4Q::CHP:dI="<<dI<<",cC="<<cC<<G4endl;
#endif
				}
			  }                                // >>>> End of if of Quark Exchange possibility
			}                                  // +++> End of if of existinr residual Q Code
            probability+=comb;                 // Collect the probability for the fragment
#ifdef sdebug
	 	    G4cout<<"G4Q::CHPr:pr="<<probability<<"(+"<<comb<<"),iq="<<iq<<",oq="<<oq<<G4endl;
#endif
		  }                                    // ...> End of Quark Exchange "oq" Test LOOP
		}                                      // ...> End of Quark Exchange "iq" Test LOOP
	  }                                        // ---> End of Nuclear Case of fragmentation
      else if(cPDG<80000000)                   // ===> "Hadron" case (QUark FUsion mechanism)
      {
        comb = valQ.NOfCombinations(candQC);
        if(cPDG==111||cPDG==221||cPDG==331||cPDG==113||cPDG==223||cPDG==333||cPDG==115||cPDG==225||
           cPDG==335||cPDG==117||cPDG==227||cPDG==337||cPDG==110||cPDG==220||cPDG==330)
        {
          G4QContent tQCd(1,0,0,1,0,0);
          G4QContent tQCu(0,1,0,0,1,0);
          G4QContent tQCs(0,0,1,0,0,1);
          G4double cmd=valQ.NOfCombinations(tQCd);
          G4double cmu=valQ.NOfCombinations(tQCu);
          G4double cms=valQ.NOfCombinations(tQCs);
          if(cPDG!=333&&cPDG!=335&&cPDG!=337) comb=(cmd+cmu)/2.;
          if(cPDG==331||cPDG==221) comb =(comb + cms)/2.; //eta,eta'
          // Hadronization in eta/eta' suppression (conversion to rho/omega)
          if(cPDG==221||cPDG==331) comb*=EtaEtaprime;     // eta/eta': -2*(1-e)
          if(cPDG==113||cPDG==223) comb*=2.-EtaEtaprime;  // rho/omega: 2*(1-e)
#ifdef pdebug
          if(cPDG==111) G4cout<<"G4Q::CHP:comb="<<comb<<",d="<<cmd<<",u="<<cmu<<",s="<<cms<<G4endl;
#endif
        }
        curQ -= candQC;                           // This is a quark content of residual quasmon
        resPDG = curQ.GetSPDGCode();              // PDG for the lowest residual hadronic state
        G4QContent resTQC = curQ+envQC;           // Total nuclear Residual Quark Content
        G4double resTM=G4QPDGCode(resTQC.GetSPDGCode()).GetMass();
#ifdef pdebug
        G4int aPDG = abs(cPDG);
        if(aPDG<10000&&aPDG%10<3)
		G4cout<<"G4Q::CHP:***>>>PDG="<<cPDG<<",cQC="<<candQC<<",comb="<<comb<<",rQC="<<curQ<<G4endl;
#endif
        if (resPDG==221 && (mQ<683.|| absb)) resPDG=111;// pi0 minimum residual instead of eta
        if (resPDG==111 && (mQ>683.&&!absb)) resPDG=221;// eta minimum residual instead of pi0
#ifdef pdebug
        if(aPDG<10000&&aPDG%10<3)
		G4cout<<"G4Q::CHP:cPDG="<<cPDG<<",comb="<<comb<<",rPDG="<<resPDG<<curQ<<", tM="
                <<totMass<<">"<<frM-CB+resTM<<"=fM="<<frM<<"+RM="<<resTM<<"-CB="<<CB<<G4endl;
#endif
        if(comb&&resPDG && totMass>frM-CB+resTM &&(resPDG>80000000&&resPDG!=90000000||resPDG<10000))
	    {
#ifdef pdebug
          if(aPDG<10000&&aPDG%10<3)
		  G4cout<<"G4Q::CHP:ind="<<index<<",Q="<<valQ<<mQ<<",c="<<cPDG<<",r="<<resPDG<<curQ<<G4endl;
#endif
          if(resPDG!=10)resM=G4QPDGCode(resPDG).GetMass();// PDG mean mass for the resid. hadron
          else resM=G4QChipolino(curQ).GetMass();  // Chipolino mass for the residual hadron
          G4int resQCode=G4QPDGCode(curQ).GetQCode();
#ifdef pdebug
          if(aPDG<10000&&aPDG%10<3)
          G4cout<<"G4Q::CHP: RQMass/QC="<<resM<<curQ<<",ePDG="<<envPDGC<<",rQC="<<resQCode<<G4endl;
#endif
          //if(envPDGC>80000000 && envPDGC!=90000000 && resM>0. && aPDG>1000 && // @@??
          if(envPDGC>80000000 && envPDGC!=90000000 && resM>0. &&
             resPDG!=10 && resPDG!=1114 && resPDG!=2224) // => Take into account the Environment
          { 
            G4QContent rtQC=curQ+envQC;            // Total Residual Quark Content
            G4QNucleus rtN(rtQC);                  // Create a pseudo-nucleus for residual
            G4double rtM =rtN.GetMZNS();           // Min Mass of total residual Nucleus
            G4double bnRQ=rtM-envM;                // Bound mass of residual Quasmon
#ifdef pdebug
            if(aPDG<10000&&aPDG%10<3)
            G4cout<<"G4Q::CHP: **Recalculate** RQMass="<<bnRQ<<",envM="<<envM<<",rtM="<<rtM<<G4endl;
#endif
            // ***VRQ***
            if(bnRQ<resM) resM=bnRQ;
          }
          if(resM>0. && resQCode>-2)
	      {
            G4double limM=mQ-resM;
            G4double rndM=GetRandomMass(cPDG,limM);// Candidate's Mass randomization
#ifdef pdebug
            G4double cMass=G4QPDGCode(cPDG).GetMass();
            if(aPDG<10000&&aPDG%10<3)
			G4cout<<"G4Q::CHP:rndM="<<rndM<<",limM="<<limM<<" > cM="<<cMass<<" ,rM+fM="<<resM+rndM
                    <<" < mQ="<<mQ<<G4endl;
#endif
            // --- Kinematical Factors ---
            if(rndM>0. && resM+rndM<mQ)
	        {
              curCand->SetEBMass(rndM);            // Set a Randomized Mass value of the Candidate
              curCand->SetNBMass(rndM);            // Set a Randomized Mass value of the Candidate
              G4double mH2 = rndM*rndM;            // Squared mass of the candidate (Mu2)
              G4double rHk = mH2/dk;
              G4double zMax = 1.-rHk/mQ;           // z_max
              G4double mR2 = resM*resM;            // SquaredBoundedMass of the ResidualQuasmon
              G4double zMin=0.;
              if(qBar) zMin= mR2/mQ/(mQ-dk);       // z_min for Quasmon-Baryon
			  G4double possibility=zMax-zMin;
#ifdef pdebug
              if(aPDG<10000&&aPDG%10<3)
				G4cout<<"G4Q::CHP:M="<<rndM<<",pos="<<possibility<<",za="<<zMax<<",rH="<<rHk<<",mQ="
                      <<mQ<<","<<dk<<",zi="<<zMin<<",mR2="<<mR2<<","<<resM<<";"<<mQ*(mQ-dk)<<G4endl;
#endif
              if (resPDG==10)                      // Chipolino case - check minimum
		      {
                G4double rM2 = mQk*(mQ-rHk);
                if(rM2<resM*resM) possibility = 0.;
              }
              if (possibility>0. && vap>0 && zMax>zMin)
	          {
                probability = vaf*(pow(zMax, vap)-pow(zMin, vap));
#ifdef pdebug
                if(aPDG<10000&&aPDG%10<3)
                G4cout<<"G4Q::CHP:#"<<index<<",m2="<<mH2<<",n="<<nOfQ<<",p="<<probability
                      <<",vaf="<<vaf<<",vap="<<vap<<",zMax="<<zMax<<",zMin="<<zMin<<G4endl;
#endif
                
                if(qBar>1&&baryn>0)                //---> High baryon number ("nuclear") case
                {
                  G4QContent rtQC=curQ+envQC;            // Total Residual Quark Content
                  G4QNucleus rtN(rtQC);                  // Create a pseudo-nucleus for residual
                  //////G4double rtM =rtN.GetMZNS();           // Min Mass of total residual Nucleus
                  /////////G4double bnRQ=rtM-envM;                // Bound mass of residual Quasmon
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
#ifdef pdebug
              if(aPDG<10000&&aPDG%10<3)
              G4cout<<"G4Q::CHP:cM=0["<<cPDG<<"],mQ="<<mQ<<valQ<<",rM="<<resM<<curQ<<G4endl;
#endif
            }
	      }
          else
	      {
#ifdef pdebug
            if(aPDG<10000&&aPDG%10<3)
            G4cout<<"***G4Q::CHP:M=0,#"<<index<<valQ<<",cP="<<cPDG<<"+rP="<<resPDG<<curQ<<G4endl;
#endif
	      }
        }
        else
	    {
          probability=0.;
#ifdef pdebug
          if(aPDG<10000&&aPDG%10<3)
          G4cout<<"G4Q::CHP:"<<index<<valQ<<",c="<<cPDG<<"+r="<<resPDG<<curQ<<":comb=0 || tM="
                <<totMass<<"<"<<frM-CB+resTM<<" = fM="<<frM<<"+RM="<<resTM<<"-CB="<<CB<< G4endl;
#endif
	    }
        //if(cPDG==111) secondProbab = 1.;
      }// ---> End of Hadronic Case of fragmentation
      else probability=0.;
#ifdef pdebug
      G4int aPDG = abs(cPDG);
      if(cPDG>90000000&&baryn<5||aPDG<10000&&aPDG%10<3)
      G4cout<<"G4Q::CHP:^^^cPDG="<<cPDG<<",pos="<<pos<<",rPDG="<<resPDG<<curQ<<resM<<", p="
            <<probability<<", s="<<accumulatedProbability<<",sp="<<secondProbab<<G4endl;
#endif
    }                                                 // ===> End of possibility check
    curCand->SetRelProbability(probability);
    accumulatedProbability += probability;
    curCand->SetIntegProbability(accumulatedProbability);
    curCand->SetSecondRelProb(secondProbab);
    secondAccumProbability += secondProbab;
    curCand->SetSecondIntProb(secondAccumProbability);
  }                                                   // ***> End of LOOP over candidates}
} // End of "CalculateHadronizationProbabilities"

// Check that it's still possible to decay the Total Residual Nucleus in Quasmon + Environ & correct
G4bool G4Quasmon::CheckGroundState(G4bool corFlag) // Correction is forbidden by default
{ //   ===========================================
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  /////////static const G4double mNeut= G4QPDGCode(2112).GetMass();
  /////////static const G4double mProt= G4QPDGCode(2212).GetMass();
  /////////static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4QNucleus vacuum= G4QNucleus(G4LorentzVector(0.,0.,0.,0.),90000000);
  ///@@@///
  /////////corFlag=true;
  ///@@@///
  G4int    resQPDG=valQ.GetSPDGCode();           // Reachable in a member function
  G4double resQMa=G4QPDGCode(resQPDG).GetMass(); // GS Mass of the Residual Quasmon
  G4double resEMa=0.;                            // GS Mass of the Empty Residual Environment
  G4bool   bsCond=false;                         // BaryonSeparetionCondition for Quasmon in vacuum
  G4LorentzVector enva4M=G4LorentzVector(0.,0.,0.,0.);
  G4LorentzVector reTLV=q4Mom;                   // Prototyoe of the 4-Mom of the Residual Nucleus
  G4double resSMa=resQMa;                        // Prototype of MinimalSplitMass of ResidualNucleus
  if(theEnvironment.GetPDG()!=90000000)          // "Existing Environment" case
  {
    resEMa=theEnvironment.GetMZNS();             // GSMass of the Residual Environment
    enva4M=G4LorentzVector(0.,0.,0.,resEMa);     // 4-Mom of the Residual Environment
    reTLV+=enva4M;                               // 4-Mom of Residual Nucleus
    resSMa+=resEMa;                              // Minimal Split Mass of Residual Nucleus
  }
  else                                           //Calculate BaryonSeparetionCondition for vacQuasm.
  {
    //G4QNucleus tmpQN(valQ,reTLV);                // TemporaryNucleus for the ResidualVacuumQuasmon
    G4QNucleus tmpQN(valQ);                      // Temporary nucleus for the Residual vacuumQuasmon
    bsCond = tmpQN.SplitBaryon();                // Possibility to split a baryon from the ResidVacQ
  }
  G4QContent envaQC = theEnvironment.GetQCZNS(); // Quark Content of the Environment
  G4double resTMa=reTLV.m();                     // CM Mass of the Residual Nucleus (Quasm+Environ)
  G4int nOfOUT = theQHadrons.size();          // Total #of QHadrons at this point
#ifdef pdebug
  G4cout<<"G4Q::CheckGS: (totM="<<resTMa<<" < rQM+rEM="<<resSMa<<" || rEM="<<resEMa<<"=0 && bs="
        <<bsCond<<"=0) && n="<<nOfOUT<<" >0"<<G4endl;
#endif
  if((resTMa<resSMa || !resEMa&&!bsCond) && nOfOUT>0 && corFlag) // *** CORRECTION ***
  {
    G4QHadron*  theLast = theQHadrons[nOfOUT-1];
    if(!(theLast->GetNFragments()) && theLast->GetPDGCode()!=22) // Not decayed hadron && not gamma
	{
      G4LorentzVector hadr4M=theLast->Get4Momentum();
      G4double  hadrMa=hadr4M.m();                // Mass of the Last hadron (==GSMass)
      G4LorentzVector tmpTLV=reTLV+hadr4M;        // Tot (ResidNucl+LastHadron) 4-Mom
#ifdef pdebug
	  G4cout<<"G4Q::CheckGS:YES,tLV="<<tmpTLV<<tmpTLV.m()<<">rSM+hM="<<resSMa+hadrMa<<G4endl;
#endif
      if(tmpTLV.m()>resSMa+hadrMa)                // resMa contains 2 Hadrons: resQ and Environ
      {
        if(resEMa)                                // => "The not vacuum Environment exists" case
        {
          G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resQMa); // GS Mass of Quasmon
          if(!G4QHadron(tmpTLV).DecayIn3(hadr4M,quas4M,enva4M))
          {
            G4cerr<<"***G4Quasm::CheckGS: Decay in Fragm+ResQ+ResEnv did not succeeded"<<G4endl;
            return false;
          }
          else
          {
            //@@CHECK CoulBar (only for ResQuasm in respect to ResEnv) and may be evaporate instead
            //theEnvironment = G4QNucleus(envaQC,enva4M);
            G4QHadron* envH = new G4QHadron(envaQC,enva4M); //@@ Moving Environment !
            FillHadronVector(envH);               // Fill Moving Environment (delete equivalent)
            theEnvironment = vacuum;
            G4QHadron* quasH = new G4QHadron(valQ, quas4M);
            //quasH->Set4Momentum(quas4M);        // @@
            FillHadronVector(quasH);              // Fill ResidQuasm Hadron (delete equivalent)
            theLast->Set4Momentum(hadr4M);
		  }
	    }
        else                                      // => "The Environ is vacuum" case (DecayIn2)
        {
          G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resQMa); // GS Mass of Quasmon
          G4QHadron* quasH = new G4QHadron(valQ, quas4M);
          if(!G4QHadron(tmpTLV).DecayIn2(hadr4M,quas4M))
          {
            delete quasH;                         // Delete "new Quasmon Hadron"
            G4cerr<<"***G4Quasm::CheckGS: Decay in Fragm+ResQ did not succeeded"<<G4endl;
            return false;
          }
          else
          {
            //@@CHECK CoulBar (only for ResQuasm in respect to ResEnv) and may be evaporate instead
            theLast->Set4Momentum(hadr4M);
            quasH->Set4Momentum(quas4M);
            FillHadronVector(quasH);              // Fill ResidQuasmHadron (delete equivalent)
		  }
	    }
      }
      else                                        // "CORRECTION" !!!!
      {
        if(nOfOUT>1 && corFlag)
	    {
          G4QHadron*  thePrev = theQHadrons[nOfOUT-2];    // Get pointer to prev-before-last hadron
          if(thePrev->GetNFragments()||thePrev->GetPDGCode()==22) return false; // Decayed or gamma
          G4LorentzVector prev4M=thePrev->Get4Momentum(); // 4-mom of the Previous but Last Hadron
          G4double  prevMa=prev4M.m();                    // Mass of previous hadron (==GSMass)
          tmpTLV+=prev4M;                                 // Increment Total 4-Mom of TotalResidNucl
          G4QContent totQC=valQ+envaQC;                   // Quark Cont of the ResidNucl=ResidQ+Env
          G4int      totPDG=totQC.GetSPDGCode();          // PDG Code of Total Residual Nucleus 
          G4double   totQMa=G4QPDGCode(totPDG).GetMass(); // GS Mass of the Residual Nucleus
          G4double   totNMa=tmpTLV.m();                   // Real Mass of the Total Residual nucleus
#ifdef pdebug
	      G4cout<<"G4Q::CheckGS:NO, 4M/M="<<tmpTLV<<totNMa<<">r+p+l="<<totQMa+hadrMa+prevMa<<G4endl;
#endif
          if(totNMa>totQMa+hadrMa+prevMa)
          {
            G4LorentzVector nuc4M = G4LorentzVector(0.,0.,0.,totQMa); // 4-Mom of ResidNucl at rest
            if(!G4QHadron(tmpTLV).DecayIn3(hadr4M,prev4M,nuc4M))
            {
              G4cerr<<"***G4Quasm::CheckGS:DecayIn3 ResidNuc+LastH+PrevH did not succeeded"<<G4endl;
              return false;
            }
            else
            {
              //theEnvironment = G4QNucleus(totQC, nuc4M);
              G4QHadron* envH = new G4QHadron(totQC,nuc4M); //@@ Moving Environment !
              FillHadronVector(envH);                    // Fill Moving Environment (delete equiv.)
              theEnvironment = vacuum;                   //@@ Instead of Moving Environment !
              theLast->Set4Momentum(hadr4M);
              thePrev->Set4Momentum(prev4M);
#ifdef pdebug
	          G4cout<<"G4Q::CheckGS: Yes, Check D4M="<<tmpTLV-hadr4M-prev4M-nuc4M<<G4endl;
#endif
	        }
		  }
          else                                           // Try to decay in one hadron, absorb other
		  {
            G4QContent tmpLNQ=totQC+thePrev->GetQC();
            G4int    resLPDG =tmpLNQ.GetSPDGCode();
            G4double resLastM=G4QPDGCode(resLPDG).GetMass(); // GSMass of ResidNucl for the LastH
            G4QContent tmpPNQ=totQC+theLast->GetQC();
            G4int    resPPDG =tmpPNQ.GetSPDGCode();
            G4double resPrevM=G4QPDGCode(resPPDG).GetMass(); // GSMass of ResidNucl for the PrevH
            //////G4bool    which = true;                // The LastH is absorbed, PrevH is radiated
#ifdef pdebug
	        G4cout<<"G4Quasm::CheckGS: NO, tM="<<totNMa<<" > rp+l="<<resLastM+hadrMa<<" || > rl+p="
                  <<resPrevM+prevMa<<G4endl;
#endif
            if      (totNMa>resLastM+hadrMa)             // "Just exclude the Last" case
		    {
              theQHadrons.pop_back();                  // the LastH is excluded from OUTPUT HV
              theQHadrons.pop_back();                  // the PrevH is excluded from OUTPUT HV
              theQHadrons.push_back(theLast);               // LastH substitute Prev H in OUTPUT HV
              G4QHadron* destrP=thePrev;                 // destruction Pointer for the QHadron
              delete     destrP;                         // the Last QHadron is destructed
              thePrev=theLast;
              resPPDG=resLPDG;
              resPrevM=resLastM;
              prev4M  = hadr4M;
		    }
            else if (totNMa>resPrevM+prevMa)              // "Just exclude the Last" case
		    {
              theQHadrons.pop_back();
              delete theLast;
		    }
            else return false;
			//return false; // @@@@@@@@@@@@@@@@@@@@@ temporary
            G4LorentzVector nuc4M = G4LorentzVector(0.,0.,0.,resPrevM);// 4-M of ResidNucl to PrevH
            if(!G4QHadron(tmpTLV).DecayIn2(prev4M,nuc4M))
            {
              G4cerr<<"***G4Quasm::CheckGS:DecayIn2 (ResNucl+LastH)+PrevH didn't succeeded"<<G4endl;
              return false;
            }
            else
            {
              //theEnvironment = G4QNucleus(nuc4M, resPPDG);
              G4QHadron* envH = new G4QHadron(resPPDG,nuc4M); //@@ Moving Environment !
              FillHadronVector(envH);              // Fill Moving Environment (delete equivalent)
              theEnvironment = vacuum;
              thePrev->Set4Momentum(prev4M);
#ifdef pdebug
	          G4cout<<"G4Q::CheckGS:Yes, Check D4M="<<tmpTLV-prev4M-nuc4M<<G4endl;
#endif
	        }
		  }
	    }
        else return false;
      }
	}
    else return false;
  }
  else return false;
  return true;
} // End of "CheckGroundState"

// Decay of outgoing hadron
G4bool G4Quasmon::DecayOutHadron(G4QHadron* qHadron, G4int DFlag)
{ //   ==========================================================
  G4QPDGCode theQPDG  = qHadron->GetQPDG();
  G4int      thePDG   = theQPDG.GetPDGCode();         // Get the PDG code of decaying hadron
  G4int        pap = 0;                               // --- particle
  if(thePDG<0) pap = 1;                               // --- anti-particle
  G4LorentzVector t = qHadron->Get4Momentum();        // Get 4-momentum of decaying hadron
  G4double m = t.m();                                 // Get the mass value of decaying Hadron
  // --- Randomize a channel of decay
  G4QDecayChanVector decV = theWorld.GetQParticle(theQPDG)->GetDecayVector();
  G4int nChan = decV.size();
#ifdef tdebug
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
#ifdef tdebug
        G4cout<<"G4Quasmon::DecayOutHadron: i="<<i<<",r="<<rnd<<" < dl="<<dC->GetDecayChanLimit()
              <<", mm="<<dC->GetMinMass()<<G4endl;
#endif
        if(rnd<dC->GetDecayChanLimit() && m>dC->GetMinMass()) break;
	  }
	}
    G4QPDGCodeVector cV=decV[i]->GetVecOfSecHadrons();// PDGVector of the selected Decay Channel
    G4int nPart=cV.size();                         // A#of particles to decay in
#ifdef tdebug
	G4cout<<"G4Quasmon::DecayOutHadron: resi="<<i<<",nP="<<nPart<<":"<<cV[0]->GetPDGCode()
          <<","<<cV[1]->GetPDGCode();
    if(nPart>2) G4cout<<","<<cV[2]->GetPDGCode();
    G4cout<<G4endl;
#endif
    if(nPart<2||nPart>3)
    {
      G4cerr<<"***G4Q::DecayOutH: n="<<nPart<<", ch#"<<i<<",PDG="<<thePDG<<G4endl;
      return false;
	}
#ifdef tdebug
    G4cout<<"G4Quasmon::DecayOutHadron: Fill PDG= "<<thePDG<<t<<m<<"(nP="<<nPart<<")>>>>>"<<G4endl;
#endif
    if(nPart==2)
	{
      G4QHadron* fHadr;
      G4QHadron* sHadr;
      G4int fPDG=cV[0]->GetPDGCode();
      G4int sPDG=cV[1]->GetPDGCode();
      if(fPDG!=22 && sPDG!=22 || abs(thePDG)%10>2 || DFlag)
      {
        if(cV[0]->GetWidth()==0.)
	    { // Randomize only the second Hardon or none
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
          G4double fm=fHadr->GetMass();                   // Mass of the first hadron
          G4double sdm = m - fm;                          // MaxMassLimit for the Second Hadron
          sHadr = new G4QHadron(sPart,sdm);               // the 2-nd Hadron is initialized *2*
          if(sPDG<0) sHadr->MakeAntiHadron();
#ifdef tdebug
          G4cout<<"G4Q::DOH:M="<<m<<",mi="<<mim<<",fd="<<fdm<<",fm="<<fm<<",sd="<<sdm
                <<",sm="<<sHadr->GetMass()<<G4endl;
#endif
        }
#ifdef tdebug
        G4cout<<"G4Q::DOH:(Decay in 2) m1="<<fHadr->GetMass()<<", m2="<<sHadr->GetMass()<<G4endl;
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
          G4cerr<<"***G4Q::DecayOutHadron:in2, PDGC="<<thePDG<<", ch#"<<i<<G4endl;
          G4Exception("***Ex***G4Q::Failed to decay in 2");
          theQHadrons.push_back(qHadron);                   // Fill as it is (delete equivalent)
          return false;
	    }
        else
        {
          //qHadron->SetNFragments(2);
          //theQHadrons.push_back(qHadron);                   // Fill with NFr=2 (delete equivalent)
          // Instead
          delete qHadron;
          //
          fHadr->Set4Momentum(f4Mom);                    // Put the randomized 4Mom to 1-st Hadron
          FillHadronVector(fHadr);                       // Fill 1st Hadron (delete equivalent)
          sHadr->Set4Momentum(s4Mom);                    // Put the randomized 4Mom to 2-nd Hadron
          FillHadronVector(sHadr);                       // Fill 2nd Hadron (delete equivalent)
        }
	  }
      else theQHadrons.push_back(qHadron);                  // Fill hadron as it is (delete equivalent)
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
        tHadr = new G4QHadron(tPart,tdm);                  // the Third Hadron is created *3*
        if(tPDG<0) tHadr->MakeAntiHadron();
      }     
#ifdef tdebug
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
        G4cerr<<"***G4Q::DecayOutHadron:in3, PDGC="<<thePDG<<", ch#"<<i<<G4endl;
        theQHadrons.push_back(qHadron);                   // Fill as it is (delete equivalent)
        return false;
	  }
      else
      {
        //qHadron->SetNFragments(3);
        //theQHadrons.push_back(qHadron);                   // Fill with NFr=3 (delete equivalent)
        // Instead
        delete qHadron;
        //
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
    G4cout<<"G4Quasmon::DecayOutHadron: Fill PDG= "<<thePDG<<t<<m<<" as it is ***0***>>>>"<<G4endl;
#endif
    theQHadrons.push_back(qHadron);                       // Fill as it is (delete equivalent)
  }
  return true;                                         // FAKE (impossible to reach)
} // End of "DecayOutHadron"

// Random integer value for the Poiasson Distribution with meanValue
G4int G4Quasmon::RandomPoisson(G4double meanValue)
//               =================================
{
  if (meanValue<=0.)
  {
    G4cerr<<"***G4Quasmon::RandomPoisson: negative or zero Mean Value="<<meanValue<<G4endl;
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
// End of "RandomPoisson"

//The public Hadronisation routine with delete responsibility of User (!)
G4QHadronVector* G4Quasmon::Fragment(G4QNucleus& nucEnviron, G4int nQ)
{//              =====================================================
#ifdef pdebug
  G4cout<<"G4Quasmon::Fragment is called theEnviron="<<nucEnviron<<G4endl;
#endif
  HadronizeQuasmon(nucEnviron,nQ);
  G4int nHadrs=theQHadrons.size();
#ifdef pdebug
  G4cout<<"G4Quasmon::Fragment after HadronizeQuasmon nH="<<nHadrs<<G4endl;
#endif
  G4QHadronVector* theFragments = new G4QHadronVector;
  if(nHadrs) for (int hadron=0; hadron<nHadrs; hadron++)
  {
    G4QHadron* curHadr = new G4QHadron(theQHadrons[hadron]);
    theFragments->push_back(curHadr);         // (delete equivalent - user)
  }
#ifdef pdebug
  else G4cerr<<"*******G4Quasmon::Fragment *** Nothing is in the output ***"<<G4endl;
#endif
  return theFragments;
} // End of "Fragment"
