// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QNucleus.cc,v 1.2 2000-08-23 11:29:33 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QNucleus ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Nuclei/Nuclear Environment used by CHIPS Model
// ------------------------------------------------------------------

//#define debug
//#define pdebug
//#define ppdebug

#include "G4QNucleus.hh"

G4QNucleus::G4QNucleus() : Z(0),N(0),S(0),maxClust(0) {};

G4QNucleus::G4QNucleus(G4int z, G4int n, G4int s) :
  Z(z),N(n),S(s),maxClust(0)
{
  G4int zns=z+n+s;
  G4QContent nQC(n+zns,z+zns,s,0,0,0);
  SetQC(nQC);
  G4QPDGCode nPDG(90000000+s*1000000+z*1000+n);
  SetQPDG(nPDG);
  G4LorentzVector p(nPDG.GetMass(),0.,0.,0.);
  Set4Momentum(p);
  SetNFragments(0);
}

G4QNucleus::G4QNucleus(G4QContent nucQC):maxClust(0)
{
  G4int u=nucQC.GetU()-nucQC.GetAU();
  G4int d=nucQC.GetD()-nucQC.GetAD();
  S = nucQC.GetS()-nucQC.GetAS();
  G4int du= d-u;      // isotopic shift
  G4int b =(d+u+S)/3; // baryon number
  Z = (b-S-du)/2;     // protons
  N = Z+du;           // neutrons
  SetQC(nucQC);
  G4QPDGCode nPDG(90000000+S*1000000+Z*1000+N);
  SetQPDG(nPDG);
  G4LorentzVector p(nPDG.GetMass(),0.,0.,0.);
  Set4Momentum(p);
  SetNFragments(0);
}

G4QNucleus::G4QNucleus(G4int nucPDG):maxClust(0) {InitByPDG(nucPDG);}

G4QNucleus::G4QNucleus(G4LorentzVector p, G4int nucPDG):maxClust(0)
{
  InitByPDG(nucPDG);
  Set4Momentum(p);
}

G4QNucleus::G4QNucleus(G4int z, G4int n, G4int s, G4LorentzVector p) :
  Z(z),N(n),S(s),maxClust(0)
{
  Set4Momentum(p);
  SetNFragments(0);
  G4int ZNS=Z+N+S;
  G4QPDGCode nPDG(90000000+S*1000000+Z*1000+N);
  SetQPDG(nPDG);
  G4QContent nQC(N+ZNS,Z+ZNS,S,0,0,0);
  SetQC(nQC);
}

G4QNucleus::G4QNucleus(G4QContent nucQC, G4LorentzVector p):maxClust(0)
{
  Set4Momentum(p);
  G4int u=nucQC.GetU()-nucQC.GetAU();
  G4int d=nucQC.GetD()-nucQC.GetAD();
  S = nucQC.GetS()-nucQC.GetAS();
  G4int du= d-u;      // isotopic shift
  G4int b =(d+u+S)/3; // baryon number
  Z = (b-S-du)/2;     // protons
  N = Z+du;           // neutrons
  SetQC(nucQC);
  G4QPDGCode nPDG(90000000+S*1000000+Z*1000+N);
  SetQPDG(nPDG);
  SetNFragments(0);
}

G4QNucleus::G4QNucleus(const G4QNucleus &right):maxClust(0)
{
  Set4Momentum         (right.Get4Momentum());
  SetQPDG              (right.GetQPDG());
  SetQC                (right.GetQC());
  SetNFragments        (right.GetNFragments());
  Z = right.GetZ();
  N = right.GetN();
  S = right.GetS();
  maxClust=right.maxClust;
  for(int i=0; i<=maxClust; i++) probVect[i]=right.probVect[i];
}

G4double G4QNucleus::freeNuc=0.1;  
G4double G4QNucleus::freeDib=.05;  
G4double G4QNucleus::clustProb=4.;
// Fill the private parameters
void G4QNucleus::SetParameters(G4double fN, G4double fD, G4double cP)
{//  ================================================================
  freeNuc=fN; 
  freeDib=fD; 
  clustProb=cP;
}

// Init existing nucleus by new PDG Code
void G4QNucleus::InitByPDG(G4int nucPDG)
{//  ===================================
  static const G4int NUCPDG  = 90000000;
  if(nucPDG>=NUCPDG)
  {
    G4int zns=nucPDG-NUCPDG;
    if(zns)
    {
      G4int zs =zns/1000;
      N  = zns%1000;
      Z  = zs%1000;
      S  = zs/1000;
    }
    else
	{
      Z  = 0;
      N  = 0;
      S  = 0;
	}
    G4int ZNS=Z+N+S;
    G4QContent nQC(N+ZNS,Z+ZNS,S,0,0,0);
    SetQC(nQC);
    G4QPDGCode nPDG(nucPDG);
    SetQPDG(nPDG);
    G4LorentzVector p(0.,0.,0.,nPDG.GetMass());
    Set4Momentum(p);
    SetNFragments(0);
#ifdef debug
	cout<<"G4QNucleus::InitByPDG:->QPDG="<<nPDG<<": Z="<<Z<<",N="<<N<<",S="<<S<<",4M="<<p<<endl;
#endif
  }
  else cerr<<"***G4QNucleus::InitByPDG: Initialized by not nuclear PDGCode="<<nucPDG<<endl;
}

// Assignment operator
const G4QNucleus& G4QNucleus::operator=(const G4QNucleus &right)
{//               ==============================================
  Set4Momentum         (right.Get4Momentum());
  SetQPDG              (right.GetQPDG());
  SetQC                (right.GetQC());
  SetNFragments        (right.GetNFragments());
  Z = right.GetZ();
  N = right.GetN();
  S = right.GetS();
  maxClust=right.maxClust;
  for(int i=0; i<=maxClust; i++) probVect[i]=right.probVect[i];

  return *this;
}

G4QNucleus::~G4QNucleus() {}

// Standard output for QNucleus {Z - a#of protons, N - a#of neutrons, S - a#of lambdas}
ostream& operator<<(ostream& lhs, G4QNucleus& rhs)
{//      =========================================
  lhs<<"{Z="<<rhs.GetZ()<<",N="<<rhs.GetN()<<",S="<<rhs.GetS()<< "}";
  return lhs;
}

// Standard output for const QNucleus {Z - a#of protons, N - a#of neutrons, S - a#of lambdas}
ostream& operator<<(ostream& lhs, const G4QNucleus& rhs)
{//      ===============================================
  lhs<<"{Z="<<rhs.GetZ()<<",N="<<rhs.GetN()<<",S="<<rhs.GetS()<< "}";
  return lhs;
}

void G4QNucleus::UpdateClusters(G4int maxCls)
{//  ========================================
  //static const G4double r0 = 1.1;               // fm, for nuclear radius: r=r0*A^(1/3)
  //static const G4double del= .55;               // fm, for a difused surface of the nucleus
  //static const G4double rCl= 2.0;               // clusterization radius @@??
  //static const G4double freeibuc = 0.10;         // probab. of the quasi-free baryon on surface
  //static const G4double freeDib = 0.05;         // probab. of the quasi-free dibar. on surface
  //static const G4double clustProb = 4.0;        // clusterization probability in dense region
  static const G4double prQ = 1.0;                // relative probability for a Quasmon
  //static const G4double prQ = 0.;                //@@for pi@@relative probability for Quasmon
  probVect[0]=prQ;
  maxClust=maxCls;
  if(maxClust<0) maxClust=0;
#ifdef debug
  cout<<"G4QNucleus::UpdateClusters: for nucleus with Z="<<Z<<", N="<<N<<", S="<<S<<endl;
#endif
  G4int a = Z + N + S;                        // atomic number
  G4double A=a;
  if(maxCls<1 || A==0.)
  {
#ifdef debug
    cout<<"***G4QNucleus::UpdateClusters: no clusters demanded maxCls="<<maxCls<<",A="<<A<<endl;
#endif
    return;
  }
  G4double surf=freeNuc+freeDib;                 // surface relative population
  G4double surA=A*surf;                          // surface absolute population
  G4int sA=static_cast<G4int>(surA);
  if(surA>0&&surA+G4UniformRand()>sA+1)sA+=1;    // a#of nucleons on the surface of the nucleus
  G4int dA=a-sA;                                 // a#of nucleons in a dense part of the nucleus
#ifdef debug
  cout<<"G4QNucleus::UpdateClusters:dA="<<dA<<",totA="<<A<<",sf="<<surf<<",ints="<<sA<<endl;
#endif
  G4double pA=0.;
  G4double uA=0.;
  if(surf>0.)
  {
    pA=0.5*freeDib*sA/surf;                  //@@// a#of quasi-free Nucleon Pairs on the surface
    uA=sA-pA-pA;                                 // a#of quasi-free nucleons on the nuc. surface
  }
  if(dA<2)                                       // There is no dense phase at all
  {
    maxClust=1;
    probVect[1]= (uA+dA)/A;                      // a#of quasi-free nucleons (correct)
    probVect[254]= 0;                            // a#of dense nucleons (correct)
    if(A>1 && pA>0.)
    {
      probVect[2]= (pA+pA)/A/(A-1);              // a#of quasi-free "dibaryons" (correct)
      probVect[255]= 0;                          // a#of dense "dibaryons" (correct)
      maxClust=2;
	}
#ifdef debug
  cout<<"G4QNucleus::UpdateClusters: only quasi-free nucleons pV[1]="<<probVect[1]<<endl;
#endif
    dZ=0;
    dN=0;
    dS=0;
  }
  else
  {
    maxClust=dA;                                 // "dibaryon"-clusters can be found
    G4double wrd=clustProb/dA;                   // @@ can be "wrd=wcd/A" ... as a parameter
    G4double sud=pow(1.+wrd,dA-1);               // normalization factor for the dense region
    G4double rd= dA/sud;
    probVect[1]= (rd+uA)/A;                      // a#of quasi-free nucleons (correct)
    probVect[254]= rd/A;                         // a#of dense nucleons (correct)
    rd*=wrd*(dA-1.)/2;
    G4double comb=A*(A-1)/2.;
    probVect[2]= (rd+pA)/comb;                   // a#of quasi-free "dibaryons" (correct)
    probVect[255]= rd/comb;                      // a#of dense "dibaryons" (correct)
#ifdef debug
	cout<<"G4QNucleus::UpdateClusters:p1="<<probVect[1]<<", p2="<<probVect[2]<<",sA="<<sA
        <<",uA="<<uA<<",pA="<<pA<<",wrd="<<wrd<<",sud="<<sud<<endl;
#endif
    if(dA>2)
    {
      rd=wrd/sud;                                // Pseudo pair norm (in case of dA=A, pA=0)
      for (int i=3; i<=dA; i++)
      {
        rd*=wrd;
        probVect[i]=rd;                          // Combinations are included for N,Z,&S later
#ifdef debug
        cout<<"G4QNucleus::UpdateClusters: cluster of "<<i<<" baryons, pV="<<probVect[i]<<endl;
#endif
      }
	}
    dS = S;                                      // @@ Lambdas are always in the dense region
    dZ = static_cast<int>(static_cast<double>((dA-dS)*Z)/(Z+N) + 0.5);
    dN = dA - dZ;
  }
}


// Reduce nucleus by emitted cluster with PDG Code cPDG
void G4QNucleus::Reduce(G4int cPDG)
{//  ==============================
  static const G4int NUCPDG=90000000;
  if(cPDG>NUCPDG)
  {
    G4int curPDG=GetPDG();
    G4int newPDG=curPDG-cPDG+NUCPDG;             // PDG Code of Residual Nucleus
    if(newPDG==NUCPDG) InitByPDG(NUCPDG);        // Empty
    else
    {
      if(abs(newPDG)<NUCPDG)
	  {
        cerr<<"***G4QNucleus::Reduce:iPDG="<<curPDG<<" = newPDG="<<newPDG<<"+cPDG="<<cPDG<<endl;
        G4Exception("*E*:::G4QNucleus::Reduce: Abnormal Nuclear Reduction");
	  }
      InitByPDG(newPDG);                         // Reinit the Nucleus
	}
  }
  else if(cPDG!=NUCPDG) cout<<"***G4QNucl::Reduce:Subtracting not nuclear PDGCode="<<cPDG<<endl;
}

// Increase nucleus by cluster with PDG Code cPDG
void G4QNucleus::Increase(G4int cPDG)
{//  ================================
  static const G4int NUCPDG=90000000;
  if(cPDG>NUCPDG)
  {
    G4int newPDG=GetPDG()+cPDG-NUCPDG;        // PDG Code of the New Nucleus
    InitByPDG(newPDG);                        // Reinit the Nucleus
  }
  else cerr<<"***G4QNucleus::Increase: Adding not nuclear PDGCode="<<cPDG<<endl;
}

// Evaporate one Baryon (n,p,Lambda) (h1) from the Nucleus & get Residual Nucleus (h2)
G4bool G4QNucleus::EvaporateBaryon(G4QHadron* h1, G4QHadron* h2)
{//  ===========================================================
  static const G4double   uWell=1.7;              // Depth of potential well Bi
  static const G4double   gunA=20.;               // Switch parameter
  static const G4double   maSht=1.2;              // shift for maximal x
  static const G4double   coSht=.19;              // multiple for maximal x
  static const G4double   third=1./3.;            // multiple for maximal x
  static const G4int      gPDG =   22;            // PDGCode of gamma
  static const G4QPDGCode gQPDG(gPDG);            // QPDGCode of gamma
  static const G4int      nPDG = 2112;            // PDGCode of neutron
  static const G4QPDGCode nQPDG(nPDG);            // QPDGCode of neutron
  static const G4int      pPDG = 2212;            // PDGCode of proton
  static const G4QPDGCode pQPDG(pPDG);            // QPDGCode of proton
  static const G4int      lPDG = 3122;            // PDGCode of Lambda
  static const G4QPDGCode lQPDG(lPDG);            // QPDGCode of Lambda
  static const G4int      dPDG = 90001001;        // PDGCode of deutron
  static const G4QPDGCode NPQPDG(dPDG);           // QPDGCode of deutron
  static const G4QPDGCode NNQPDG(90000002);       // QPDGCode of dineutron
  static const G4QPDGCode PPQPDG(90002000);       // QPDGCode of diproton
  static const G4QPDGCode NLQPDG(91000001);       // QPDGCode of nL
  static const G4QPDGCode PLQPDG(91001000);       // QPDGCode of pL
  static const G4QPDGCode LLQPDG(92000000);       // QPDGCode of LL
  static const G4int      aPDG = 90002002;        // PDGCode of ALPHA
  static const G4QPDGCode aQPDG(aPDG);            // QPDGCode of ALPHA
  static const G4double   mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double   mNeut= G4QPDGCode(nPDG).GetMass();          // Mass of neutron
  static const G4double   mProt= G4QPDGCode(pPDG).GetMass();          // Mass of proton
  static const G4double   mLamb= G4QPDGCode(lPDG).GetMass();          // Mass of Lambda
  static const G4double   mDeut= G4QPDGCode(pPDG).GetNuclMass(1,1,0); // Mass of Deutron
  static const G4double   mN2  = mNeut*mNeut;     // Mass^2 of neutron
  static const G4double   mP2  = mProt*mProt;     // Mass^2 of proton
  static const G4double   mL2  = mLamb*mLamb;     // Mass^2 of Lambda
  static const G4double   mA2  = mAlph*mAlph;     // Mass^2 of Alpha
  static const G4double   mNP  = mNeut+mProt;     // proton and neutron mass
  //static const G4double   mNN  = mNeut+mNeut;     // 2 neutrons mass
  //static const G4double   mPP  = mProt+mProt;     // 2 protons mass
  //static const G4double   mNL  = mNeut+mLamb;     // neutron and Lambda mass
  //static const G4double   mPL  = mProt+mLamb;     // proton and Lambda mass
  //static const G4double   mLL  = mLamb+mLamb;     // 2 Lambdas mass
  G4double uW=uWell;
  G4int    a = GetA();
#ifdef ppdebug
  cout<<"G4QNucleus::EvaporateBaryon: Called with a="<<a<<GetThis()<<endl;
#endif
  G4double a1= a-1;
  G4double a2= a-2;
  G4double z = Z;
  G4double z1= Z-1;
  G4double zn= Z+N;
  G4double CBarr= z1/(pow(a1,third)+1.);          // Coulomb Barrier for proton
  G4double ABarr= 2*(z1-1.)/(pow(a2,third)+1.26); // Coulomb Barrier for alpha
  G4LorentzVector h1mom;
  G4LorentzVector h2mom;
  G4LorentzVector h3mom;
  G4double totMass= GetMass();                    // Total mass of the Nucleus
  if(a==2)
  {
    if     (Z==2)
    {
      h1mom=G4LorentzVector(0.,0.,0.,mProt);
      h2mom=h1mom;
      h1->SetQPDG(pQPDG);
      h2->SetQPDG(pQPDG);
	  if(!DecayIn2(h1mom,h2mom)) return false;
    }
    else if(N==2)
    {
      h1mom=G4LorentzVector(0.,0.,0.,mNeut);
      h2mom=h1mom;
      h1->SetQPDG(nQPDG);
      h2->SetQPDG(nQPDG);
	  if(!DecayIn2(h1mom,h2mom)) return false;
    }
    else if(N==1&&Z==1)
    {
      if(totMass<=mNP)
	  {
#ifdef ppdebug
      cout<<"G4QNucleus::EvaporateBaryon: Photon ### d+g ###, dM="<<totMass-mNP<<endl;
#endif
        h1mom=G4LorentzVector(0.,0.,0.,0.);
        h2mom=G4LorentzVector(0.,0.,0.,mDeut);
        h1->SetQPDG(gQPDG);
        h2->SetQPDG(NPQPDG);
	  }
      else
	  {
        h1mom=G4LorentzVector(0.,0.,0.,mProt);
        h2mom=G4LorentzVector(0.,0.,0.,mNeut);
        h1->SetQPDG(pQPDG);
        h2->SetQPDG(nQPDG);
      }
	  if(!DecayIn2(h1mom,h2mom)) return false;
    }
    else if(Z==1&&S==1)
    {
      h1mom=G4LorentzVector(0.,0.,0.,mProt);
      h2mom=G4LorentzVector(0.,0.,0.,mLamb);
      h1->SetQPDG(pQPDG);
      h2->SetQPDG(lQPDG);
	  if(!DecayIn2(h1mom,h2mom)) return false;
    }
    else
    {
      h1mom=G4LorentzVector(0.,0.,0.,mNeut);
      h2mom=G4LorentzVector(0.,0.,0.,mLamb);
      h1->SetQPDG(nQPDG);
      h2->SetQPDG(lQPDG);
	  if(!DecayIn2(h1mom,h2mom)) return false;
    }
    h1->Set4Momentum(h1mom);
    h2->Set4Momentum(h2mom);
	return true;
  }
  else if(a>2)
  {
    //G4double gunFact= true; // No Gun Flag
    G4bool gunFact  = G4UniformRand()>exp(-a/gunA); // No Gun Flag
    //gunFact=true;                              //@@
    G4bool nFlag    = false;                   // Flag of possibility to radiate neutron
    G4bool pFlag    = false;                   // Flag of possibility to radiate proton
    G4bool lFlag    = false;                   // Flag of possibility to radiate lambda
    G4bool aFlag    = false;                   // Flag of possibility to radiate alpha
    G4bool nnFlag   = false;                   // Flag of possibility to radiate 2 neutrons
    G4bool npFlag   = false;                   // Flag of possibility to radiate neutron+proton
    G4bool nlFlag   = false;                   // Flag of possibility to radiate neutron+lambda
    G4bool ppFlag   = false;                   // Flag of possibility to radiate 2 protons
    G4bool plFlag   = false;                   // Flag of possibility to radiate proton+lambda
    G4bool llFlag   = false;                   // Flag of possibility to radiate 2 lambdas
    G4bool nnnF     = false;
    G4bool nnpF     = false;
    G4bool nppF     = false;
    G4bool pppF     = false;
    G4bool nnlF     = false;
    G4bool nplF     = false;
    G4bool pplF     = false;
    G4bool nllF     = false;
    G4bool pllF     = false;
    G4bool lllF     = false;
    G4double GSMass = GetGSMass();             // Ground State mass of the Nucleus
    G4double GSResNN= GSMass;                  // Prototype of Residual Nuclear Mass for n+n
    G4double GSResNP= GSMass;                  // Prototype of Residual Nuclear Mass for n+p
    G4double GSResNL= GSMass;                  // Prototype of Residual Nuclear Mass for n+l
    G4double GSResPP= GSMass;                  // Prototype of Residual Nuclear Mass for p+p
    G4double GSResPL= GSMass;                  // Prototype of Residual Nuclear Mass for p+l
    G4double GSResLL= GSMass;                  // Prototype of Residual Nuclear Mass for l+l
    G4double GSResAL= GSMass;                  // Prototype of Residual Nuclear Mass for alpha
    G4double GSReNNN= GSMass;                  // Prototype of Residual Nuclear Mass for n+n+n
    G4double GSReNNP= GSMass;                  // Prototype of Residual Nuclear Mass for n+n+p
    G4double GSReNPP= GSMass;                  // Prototype of Residual Nuclear Mass for n+p+p
    G4double GSRePPP= GSMass;                  // Prototype of Residual Nuclear Mass for p+p+p
    G4double GSReNNL= GSMass;                  // Prototype of Residual Nuclear Mass for n+n+l
    G4double GSReNPL= GSMass;                  // Prototype of Residual Nuclear Mass for n+p+l
    G4double GSRePPL= GSMass;                  // Prototype of Residual Nuclear Mass for p+p+l
    G4double GSReNLL= GSMass;                  // Prototype of Residual Nuclear Mass for n+l+l
    G4double GSRePLL= GSMass;                  // Prototype of Residual Nuclear Mass for p+l+l
    G4double GSReLLL= GSMass;                  // Prototype of Residual Nuclear Mass for l+l+l
    G4QPDGCode nnQPDG(22);                     // Prototype of nn-dibaryon QPDG
    G4QPDGCode npQPDG(22);                     // Prototype of np-dibaryon QPDG
    G4QPDGCode nlQPDG(22);                     // Prototype of nl-dibaryon QPDG
    G4QPDGCode ppQPDG(22);                     // Prototype of pp-dibaryon QPDG
    G4QPDGCode plQPDG(22);                     // Prototype of pl-dibaryon QPDG
    G4QPDGCode llQPDG(22);                     // Prototype of ll-dibaryon QPDG
    G4QPDGCode dbQPDG(22);                     // Prototype of chosen dibaryon QPDG
    G4QPDGCode fQPDG(22);                      // Prototype of QPDG of the Second Baryon
    G4double rMass  = 0.;                      // Prototype of mass of Residual Nucleus
    G4double eMass  = 0.;                      // Prototype of mass of Evaporated Baryon
    G4double fMass  = 0.;                      // Prototype of mass of the Second Baryon
#ifdef ppdebug
    cout<<"G4QNucleus::EvaporateBaryon: a>2, totM="<<totMass<<" > GSMass="<<GSMass<<" ?"<<endl;
#endif
    G4double tM2    = totMass*totMass;
    G4double qtM2   = 4*tM2;
    G4double GSResNp= GSMass;                  // Prototype of Residual Nuclear Mass for proton
    G4double pExcess= 0.;                      // Prototype of excess energy for proton
    G4double aExcess= 0.;                      // Prototype of excess energy for alpha
    G4double pp2m   = 0.;                      // Prototype of max square momentum for proton
    G4double ap2m   = 0.;                      // Prototype of max square momentum for proton
    G4double pBnd   = 0.;                      // Binding energy for proton
    G4double aBnd   = 0.;                      // Binding energy for proton
    G4bool three=false;                        // Prototype of the Flag of b+b+ResNuc decay
    if(Z>0)
    {
      GSResNp=G4QNucleus(Z-1,N,S).GetGSMass();
      G4double mpls=GSResNp+mProt;
      G4double mmin=GSResNp-mProt;
      pp2m=(tM2-mpls*mpls)*(tM2-mmin*mmin)/qtM2;
      if(pp2m>=0.000001)
	  {
        pFlag=true;
        pBnd=mProt-GSMass+GSResNp;             // Binding energy for proton
        G4double eMax=sqrt(mP2+pp2m);
#ifdef ppdebug
	    cout<<"G4QNucleus::EvapBaryon:pm="<<eMax+sqrt(pp2m+GSResNp*GSResNp)<<" = M="<<totMass
            <<", sm="<<GSResNp+mProt+CBarr<<",pp2="<<pp2m<<",pB="<<pBnd<<endl;
#endif
        pExcess=eMax-mProt+pBnd;               // Max Kin Energy from bottom
	  }
      else pExcess=pBnd;
      if(Z>1)
	  {
        GSResPP=G4QNucleus(Z-2,N,S).GetGSMass();
        ppQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-2)+N);
        if(GSResPP+mProt+CBarr+mProt+CBarr<totMass) ppFlag=true;
        if(Z>2)
        {
          GSRePPP=G4QNucleus(Z-3,N,S).GetGSMass();
          if(GSRePPP+(mProt+CBarr)*3<totMass) pppF=true;
        }
        if(N>0)
        {
          GSReNPP=G4QNucleus(Z-2,N-1,S).GetGSMass();
          if(GSReNPP+(mProt+CBarr)*2+mNeut<totMass) nppF=true;
		}
        if(S>0)
        {
          GSRePPL=G4QNucleus(Z-2,N,S-1).GetGSMass();
          if(GSRePPL+(mProt+CBarr)*2+mLamb<totMass) pplF=true;
		}
        if(N>1)
        {
          GSResAL=G4QNucleus(Z-2,N-2,S).GetGSMass();
          mpls=GSResAL+mAlph;
          mmin=GSResAL-mAlph;
          ap2m=(tM2-mpls*mpls)*(tM2-mmin*mmin)/qtM2;
          if(ap2m>=0.000001)
	      {
            aFlag=true;
            aBnd=mAlph-GSMass+GSResAL;           // Binding energy for ALPHA
            G4double eMax=sqrt(mA2+ap2m);
#ifdef ppdebug
	        cout<<"G4QNucleus::EvapBaryon:am="<<eMax+sqrt(ap2m+GSResAL*GSResAL)<<" = M="
                <<totMass<<", sm="<<GSResNp+mProt+CBarr<<",pp2="<<pp2m<<",pB="<<pBnd<<endl;
#endif
            aExcess=eMax-mAlph+aBnd;             // Max Kin Energy from bottom
	      }
          else aExcess=pBnd;
		}
	  }
      if(N>0)
	  {
        GSResNP=G4QNucleus(Z-1,N-1,S).GetGSMass();
        npQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-1)+N-1);
        if(GSResNP+mNeut+mProt+CBarr<totMass) npFlag=true;
        if(N>1)
        {
          GSReNNP=G4QNucleus(Z-1,N-2,S).GetGSMass();
          if(GSReNNP+mProt+CBarr+mNeut+mNeut<totMass) nnpF=true;
		}
        if(S>0)
        {
          GSReNPL=G4QNucleus(Z-1,N-1,S-1).GetGSMass();
          if(GSReNPL+mProt+CBarr+mNeut+mLamb<totMass) nplF=true;
		}
	  }
      if(S>0)
	  {
        GSResPL=G4QNucleus(Z-1,N,S-1).GetGSMass();
        plQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z-1)+N);
        if(GSResPL+mProt+CBarr+mLamb<totMass) plFlag=true;
        if(S>1)
        {
          GSRePLL=G4QNucleus(Z-1,N,S-2).GetGSMass();
          if(GSRePLL+mProt+CBarr+mLamb+mLamb<totMass) pllF=true;
		}
	  }
	}
    G4double GSResNn= GSMass;                  // Prototype of Residual Nuclear Mass for neutron
    G4double nExcess= 0.;                      // Prototype of excess energy for neutron
    G4double np2m   = 0.;                      // Prototype of max square momentum for neutron
    G4double nBnd   = 0.;                      // Binding energy for neutron
    if(N>0)
    {
      GSResNn=G4QNucleus(Z,N-1,S).GetGSMass();
      G4double mpls=GSResNn+mNeut;
      G4double mmin=GSResNn-mNeut;
      np2m=(tM2-mpls*mpls)*(tM2-mmin*mmin)/qtM2;
      if(np2m>=0.000001)
      {
        nFlag=true;
        nBnd=mNeut-GSMass+GSResNn;               // Binding energy for neutron
        G4double eMax=sqrt(mN2+np2m);
#ifdef ppdebug
	    cout<<"G4QNucleus::EvapBaryon:nm="<<eMax+sqrt(np2m+GSResNn*GSResNn)<<" = M="<<totMass
            <<", sm="<<GSResNn+mNeut<<",np2="<<np2m<<",nB="<<nBnd<<endl;
#endif
        nExcess=eMax-mNeut+nBnd;
	  }
      else nExcess=nBnd;
      if(N>1)
	  {
        GSResNN=G4QNucleus(Z,N-2,S).GetGSMass();
        nnQPDG=G4QPDGCode(90000000+1000*(1000*S+Z)+N-2);
        if(GSResNN+mNeut+mNeut<totMass) nnFlag=true;
        if(N>2)
        {
          GSReNNN=G4QNucleus(Z,N-3,S).GetGSMass();
          if(GSReNNN+mNeut*3<totMass) nnnF=true;
		}
        if(S>0)
        {
          GSReNNL=G4QNucleus(Z,N-2,S-1).GetGSMass();
          if(GSReNNL+mNeut+mNeut+mLamb<totMass) nnlF=true;
		}
	  }
      if(S>0)
	  {
        GSResNL=G4QNucleus(Z,N-1,S-1).GetGSMass();
        nlQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z)+N-1);
        if(GSResNL+mProt+CBarr+mLamb<totMass) nlFlag=true;
        if(S>1)
        {
          GSReNLL=G4QNucleus(Z,N-1,S-2).GetGSMass();
          if(GSReNLL+mNeut+mNeut+mLamb<totMass) nllF=true;
		}
	  }
	}
    G4double GSResNl= GSMass;                  // Prototype of Residual Nuclear Mass for Lambda
    G4double lExcess= 0.;                      // Prototype of excess energy for Lambda
    G4double lp2m   = 0.;                      // Prototype of max square momentum for lambda
    G4double lBnd   = 0.;                      // Binding energy for lambda
    if(S>0)
    {
      GSResNl=G4QNucleus(Z,N,S-1).GetGSMass();
      G4double mpls=GSResNl+mLamb;
      G4double mmin=GSResNl-mLamb;
      lp2m=(tM2-mpls*mpls)*(tM2-mmin*mmin)/qtM2;
      if(lp2m>=0.000001)
      {
        lFlag=true;
        lBnd=mLamb-GSMass+GSResNl;               // Binding energy for lambda
        G4double eMax=sqrt(mL2+lp2m);
#ifdef ppdebug
	    cout<<"G4QNucleus::EvapBaryon:lm="<<eMax+sqrt(lp2m+GSResNl*GSResNl)<<" = M="<<totMass
            <<", sm="<<GSResNl+mLamb<<",lp2="<<lp2m<<",lB="<<lBnd<<endl;
#endif
        lExcess=eMax-mLamb+lBnd;
	  }
      else lExcess=lBnd;
      if(S>1)
	  {
        GSResLL=G4QNucleus(Z,N,S-2).GetGSMass();
        ppQPDG=G4QPDGCode(90000000+1000*(1000*(S-2)+Z)+N);
        if(GSResNL+mLamb+mLamb<totMass) llFlag=true;
        if(S>2)
        {
          GSReLLL=G4QNucleus(Z,N,S-3).GetGSMass();
          if(GSRePPP+mLamb*3<totMass) lllF=true;
		}
	  }
	}
    G4bool nSecF = nnFlag || npFlag || nlFlag; // Possibility of second baryon radiation after n
    G4bool pSecF = npFlag || ppFlag || plFlag; // Possibility of second baryon radiation after p
    G4bool lSecF = nlFlag || plFlag || llFlag; // Possibility of second baryon radiation after l
    G4bool nTrF=nnnF||nnpF||nppF||nnlF||nplF||nllF;//Possib of third baryon radiation after n
    G4bool pTrF=nnpF||nppF||pppF||nplF||pplF||pllF;//Possib of third baryon radiation after p
    G4bool lTrF=nnlF||nplF||pplF||nllF||pllF||lllF;//Possib of third baryon radiation after l
    G4bool secB  = nSecF || pSecF || lSecF;    // Possibility to decay in two baryons
    G4bool thdB  = nTrF || pTrF || lTrF;       // Possibility to decay in three baryons
#ifdef ppdebug
	cout<<"G4QNucl::EvapBar:nS="<<nSecF<<",pS="<<pSecF<<",lS="<<lSecF<<",secB="<<secB<<",nnF="
        <<nnFlag<<",npF="<<npFlag<<",nlF="<<nlFlag<<",ppF="<<ppFlag<<",plF="<<plFlag<<endl;
#endif
    G4QPDGCode bQPDG;
    G4QPDGCode rQPDG;
    if(thdB&&gunFact)                            // Decay in three baryons is possible
    {
      //if(!nSecF) nFlag=false;
      //if(!pSecF) pFlag=false;
      //if(!lSecF) lFlag=false;
#ifdef ppdebug
	cout<<"G4QNucl::EvapBar:nF="<<nFlag<<",pF="<<pFlag<<",lF="<<lFlag<<endl;
#endif
      G4double maxE=0.;                          // Prototype for maximum energy
      if(nFlag&&nExcess>maxE) maxE=nExcess;
      if(pFlag&&pExcess>maxE) maxE=pExcess;
      if(lFlag&&lExcess>maxE) maxE=lExcess;
      if(lFlag&&aExcess>maxE) maxE=aExcess;
      G4double pMin=pBnd;                        // Binding energy for proton
      if(pFlag)pMin+= CBarr;                     // Add Coulomb Barrier for protons
      G4double nMin=nBnd;                        // Binding energy for neutron
      G4double lMin=lBnd;                        // Binding energy for Lambda
      G4double aMin=aBnd;                        // Binding energy for alpha
      if(aFlag)aMin+= ABarr;                     // Add Coulomb Barrier for alpha
      G4double minE=GSMass;                      // Prototype for mimimum energy
      if(nFlag&&nMin<minE) minE=nMin;
      if(pFlag&&pMin<minE) minE=pMin;
      if(lFlag&&lMin<minE) minE=lMin;
      if(aFlag&&aMin<minE) minE=aMin;
#ifdef ppdebug
      cout<<"G4QNucl::EvBar:nEx="<<nExcess<<",pEx="<<pExcess<<",sEx="<<lExcess<<",aEx="<<aExcess
          <<",minE="<<minE<<" < maxE="<<maxE<<", ni="<<nMin<<", pi="<<pMin<<", li="<<lMin<<endl;
#endif
      // @@ Here one can put a condition for the Baryon Gun
      G4int    cntr= 0;
      G4int    cntm= 49;
      if((pFlag&&pExcess>pMin ||nFlag&&nExcess>nMin ||lFlag&&lExcess>lMin ||aFlag&&aExcess>aMin)
         && minE<maxE)
      {
        G4double mi=uWell+minE;                  // Minimum Kinetic Energy for minimal nucleon
        G4double mm=uWell+maxE;                  // Personal maximum
        G4double ma=uWell*a+maxE;                // Total Kinetic Energy of baryons
        if(mi<0.)
	    {
          uW-=mi;
          mm-=mi;
          mi=0.;
	    }
        G4bool good=true;
        if(ma<mm)
        {
          ma=mm;
          good=false;
        }
#ifdef ppdebug
	    cout<<"G4QNuc::EvapB:iE="<<minE<<",aE="<<maxE<<",mi="<<mi<<",mm="<<mm<<",ma="<<ma<<endl;
#endif
        G4double xMi=mi/ma;                     // Minimal value of x
        G4double xMm=mm/ma;                     // Personal maximum x
        //G4double xCa=maSht-coSht*log(a);        // Maximal value of x (approximation)
        //G4double xMa=xCa;                       // Maximal value of x
        //if(xMm<xMa) xMa=xMm;
        G4double xMa=xMm;
        if(xMa>1.)xMa=1.;
        if(xMi<0.)xMi=0.;
        if(xMi>xMa)
        {
	      cerr<<"***G4QNucleus::EvapBaryon: M="<<mm/ma<<",xi="<<xMi<<",xa="<<xMa<<endl;
          return false;
        }
        xMi=sqrt(xMi);
        xMa=sqrt(xMa);
#ifdef ppdebug
	    cout<<"G4QNucleus:EvaporateBaryon:mi="<<mi<<",ma="<<ma<<", xi="<<xMi<<",xa="<<xMa<<endl;
#endif
        G4double powr=1.5*a1;                   // Power for low & up limits
        G4double revP=1./powr;                  // Reversed power for randomization
#ifdef ppdebug
        cout<<"G4QNucleus::EvaporateBaryon: power="<<powr<<",rev.power="<<revP<<endl;
#endif
        G4double minR=pow(1.-xMa*xMa,powr);
        G4double maxR=pow(1.-xMi*xMi,powr);
#ifdef ppdebug
        cout<<"G4QNucleus::EvaporateBaryon: miR="<<minR<<", maR="<<maxR<<endl;
#endif
        G4bool   cond=true;
        G4int    PDG = 0;
        G4double tk  = 0.;                      // Kinetic energy over the well
        while(cond&&cntr<cntm)
        {
          G4double R = minR+(maxR-minR)*G4UniformRand();
          //if(!good)R = maxR;
          G4double x2= 1.-pow(R,revP);
          G4double x = sqrt(x2);
          if(x<xMi||x>xMa)
		  {
#ifdef ppdebug
            cerr<<"G4QNucleus::EvaporateBary:R="<<R<<",xi="<<xMi<<" < "<<x<<" < xa="<<xMa<<endl;
#endif
            if(x<xMi) x=xMi;
            else      x=xMa;
		  }
          G4double rn=G4UniformRand();
          //G4double al= Z*z1*N*(N-1)/4.;
          G4double al=0.;
          G4double r = (a+al)*G4UniformRand();
          //if(rn<x/xMa||!good)
          if(rn<x/xMa)
          {
            tk= ma*x2-uW;
            if     (lFlag&&r>a&&tk>lMin)
		    {
              PDG=aPDG;
              cond=false;
		    }
            else if(lFlag&&r>zn&&r<=a&&tk>lMin)
		    {
              PDG=lPDG;
              cond=false;
		    }
            else if(nFlag&&r>=z&&r<=zn&&tk>nMin)
	        {
              PDG=nPDG;
              cond=false;
		    }
            else if(pFlag&&r<z&&tk>pMin)
	        {
              PDG=pPDG;
              cond=false;
   	        }
		  }
#ifdef ppdebug
          cout<<"G4QNuc::EvapBar:c="<<cond<<",x="<<x<<",cnt="<<cntr<<",R="<<R<<",ma="<<ma<<",r="
              <<r<<",rn="<<rn<<"<rx="<<x/xMa<<",tk="<<tk<<",ni="<<nMin<<",pi="<<pMin<<endl;
#endif
          cntr++;
	    }
        if(cntr<cntm)
	    {
          G4double p2=0.;
          if     (PDG==aPDG)
          {
            tk-=aBnd-mAlph;
            p2=tk*tk-mA2;
            if(p2>ap2m)
            {
              p2=ap2m;
              tk=sqrt(p2+mA2);
		    }
            eMass=mAlph;
            bQPDG=aQPDG;
            rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-2)+N-2);
	      }
          else if(PDG==pPDG)
          {
            tk-=pBnd-mProt;
            p2=tk*tk-mP2;
            if(p2>pp2m)
            {
              p2=pp2m;
              tk=sqrt(p2+mP2);
		    }
            eMass=mProt;
            bQPDG=pQPDG;
            rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-1)+N);
	      }
          else if(PDG==nPDG)
          {
            tk-=nBnd-mNeut;
            p2=tk*tk-mN2;
#ifdef ppdebug
            cout<<"G4QNucleus::EvaporateBaryon:np2="<<p2<<",np2m="<<np2m<<endl;
#endif
            if(p2>np2m)
            {
              p2=np2m;
              tk=sqrt(p2+mN2);
		    }
            eMass=mNeut;
            bQPDG=nQPDG;
            rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z)+N-1);
	      }
          else if(PDG==lPDG)
          {
            tk-=lBnd-mLamb;
            p2=tk*tk-mL2;
            if(p2>lp2m)
            {
              p2=lp2m;
              tk=sqrt(p2+mL2);
		    }
            eMass=mLamb;
            bQPDG=lQPDG;
            rQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z)+N);
	      }
          else cerr<<"***G4QNucleus::EvaporateBaryon: PDG="<<PDG<<endl;
          G4double rEn=totMass-tk;
          rMass=sqrt(rEn*rEn-p2);                  // Mass of Residual Nucleus
          // Now one needs to define if it's below the second baryon decay limit
          G4bool nnCond=GSResNN+mNeut>rMass;
          G4bool npCond=GSResNP+mProt+CBarr>rMass;
          G4bool nlCond=GSResNL+mLamb>rMass;
          G4bool pnCond=GSResNP+mNeut>rMass;
          G4bool ppCond=GSResPP+mProt+CBarr>rMass;
          G4bool plCond=GSResPL+mLamb>rMass;
          G4bool lnCond=GSResNL+mNeut>rMass;
          G4bool lpCond=GSResPL+mProt+CBarr>rMass;
          G4bool llCond=GSResLL+mLamb>rMass;
#ifdef ppdebug
		  cout<<"G4QNucl::EvapBary:p:rM="<<rMass<<",NN="<<GSResNN+mNeut<<",NP="
              <<GSResNP+mProt+CBarr<<",NL="<<GSResNL+mLamb<<endl;
#endif
          three=true;                               // Flag of b+b+ResNuc decay
          if     (PDG==pPDG&&pnCond&&ppCond&&plCond)// p+b+RN decay can happen
		  {
            fMass=mProt;
            fQPDG=pQPDG;
            // Now select second baryon
            G4double nLim=0.;
            if(GSResNP!=GSMass&&fMass+CBarr+mNeut+GSResNP<totMass)nLim+=N;
            G4double zLim=nLim;
            if(GSResPP!=GSMass&&fMass+CBarr+mProt+CBarr+GSResPP<totMass)zLim+=Z;
            G4double sLim=zLim;
            if(GSResPL!=GSMass&&fMass+CBarr+mLamb+GSResPL<totMass)sLim+=S;
            G4double r = sLim*G4UniformRand();
#ifdef ppdebug
			cout<<"G4QNucl::EvapBary:p:r="<<r<<",nL="<<nLim<<",zL="<<zLim<<endl;
#endif
            if     (r>zLim)
		    {
              eMass = mLamb;
              dbQPDG= PLQPDG;
              rMass = GSResPL;
              rQPDG = plQPDG;
#ifdef ppdebug
			  cout<<"G4QNucleus::EvaporateBary: P+L"<<endl;
#endif
		    }
            else if(r>=nLim&&r<=zLim&&nLim<zLim)
	        {
              eMass = mProt;
              dbQPDG= PPQPDG;
              rMass = GSResPP;
              rQPDG = ppQPDG;
#ifdef ppdebug
			  cout<<"G4QNucleus::EvaporateBary: P+P"<<endl;
#endif
		    }
            else if(r<nLim)
	        {
              eMass = mNeut;
              dbQPDG= NPQPDG;
              rMass = GSResNP;
              rQPDG = npQPDG;
#ifdef ppdebug
			  cout<<"G4QNucleus::EvaporateBary: P+N"<<endl;
#endif
   	        }
            else three=false;
		  }
          else if(PDG==nPDG&&nnCond&&npCond&&nlCond)// n+b+RN decay can happen
		  {
            fMass=mNeut;
            fQPDG=nQPDG;
            G4double nLim=0.;
            if(GSResNN!=GSMass&&fMass+mNeut+GSResNN<totMass)nLim+=N;
            G4double zLim=nLim;
            if(GSResNP!=GSMass&&fMass+mProt+CBarr+GSResNP<totMass)zLim+=Z;
            G4double sLim=zLim;
            if(GSResNL!=GSMass&&fMass+mLamb+GSResNL<totMass)sLim+=S;
            G4double r = sLim*G4UniformRand();
#ifdef ppdebug
			cout<<"G4QNucl::EvapBary:n:r="<<r<<",nL="<<nLim<<",zL="<<zLim<<endl;
#endif
            if     (r>zLim)
		    {
              eMass = mLamb;
              dbQPDG= NLQPDG;
              rMass = GSResNL;
              rQPDG = nlQPDG;
#ifdef ppdebug
			  cout<<"G4QNucleus::EvaporateBary: N+L"<<endl;
#endif
		    }
            else if(r>=nLim&&r<=zLim&&nLim<zLim)
	        {
              eMass = mProt;
              dbQPDG= NPQPDG;
              rMass = GSResNP;
              rQPDG = npQPDG;
#ifdef ppdebug
			  cout<<"G4QNucleus::EvaporateBary: N+P"<<endl;
#endif
		    }
            else if(r<nLim)
	        {
              eMass = mNeut;
              dbQPDG= NNQPDG;
              rMass = GSResNN;
              rQPDG = nnQPDG;
#ifdef ppdebug
			  cout<<"G4QNucleus::EvaporateBary: N+N"<<endl;
#endif
   	        }     
            else three=false;
		  }
          else if(PDG==lPDG&&nlCond&&plCond&&llCond)// l+b+RN decay can happen
		  {
            fMass=mLamb;
            fQPDG=lQPDG;
            G4double nLim=0.;
            if(GSResNL!=GSMass&&fMass+mNeut+GSResNL<totMass)nLim+=N;
            G4double zLim=nLim;
            if(GSResPL!=GSMass&&fMass+mProt+CBarr+GSResPL<totMass)zLim+=Z;
            G4double sLim=zLim;
            if(GSResLL!=GSMass&&fMass+mLamb+GSResLL<totMass)sLim+=S;
            G4double r = sLim*G4UniformRand();
#ifdef ppdebug
			cout<<"G4QNucl::EvapBary:l:r="<<r<<",nL="<<nLim<<",zL="<<zLim<<endl;
#endif
            if     (r>zLim)
		    {
              eMass = mLamb;
              dbQPDG= LLQPDG;
              rMass = GSResLL;
              rQPDG = llQPDG;
#ifdef ppdebug
			  cout<<"G4QNucleus::EvaporateBary: L+L"<<endl;
#endif
		    }
            else if(r>=nLim&&r<=zLim&&nLim<zLim)
	        {
              eMass = mProt;
              dbQPDG= PLQPDG;
              rMass = GSResPL;
              rQPDG = plQPDG;
#ifdef ppdebug
			  cout<<"G4QNucleus::EvaporateBary: L+P"<<endl;
#endif
		    }
            else if(r<nLim)
	        {
              eMass = mNeut;
              dbQPDG= NLQPDG;
              rMass = GSResNL;
              rQPDG = nlQPDG;
#ifdef ppdebug
			  cout<<"G4QNucleus::EvaporateBary: L+N"<<endl;
#endif
   	        }
            else three=false;
		  }
          else three=false;
#ifdef ppdebug
          cout<<"G4QNucleus::EvaporateBary:evaBar="<<eMass<<bQPDG<<",resN="<<rMass<<rQPDG
              <<",secB="<<fMass<<",three="<<three<<endl;
#endif
	    }
	  }
      else // ==============> Just decay in baryon and residual (to avoid gamma-decay)
	  {
        G4double nLim=0.;
        if(nFlag&&mNeut+GSResNn<totMass)nLim+=N;
        G4double zLim=nLim;
        if(pFlag&&mProt+CBarr+GSResNp<totMass)zLim+=Z;
        G4double sLim=zLim;
        if(lFlag&&mLamb+GSResNl<totMass)sLim+=S;
        G4double r = sLim*G4UniformRand();
#ifdef ppdebug
        cout<<"G4QNucl::EvaporateBaryon:2Decay r="<<r<<",nLim="<<nLim<<",zLim="<<zLim
            <<",sLim="<<sLim<<",nF="<<nFlag<<",pF="<<pFlag<<",lF="<<lFlag<<endl;
#endif
        if     (lFlag&&r>zLim)
		{
          bQPDG=lQPDG;
          eMass=mLamb;
          rQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z)+N);
          rMass=GSResNl;
		}
        else if(nFlag&&r>=nLim&&r<=zLim&&nLim<zLim)
	    {
          bQPDG=pQPDG;
          eMass=mProt;
          rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-1)+N);
          rMass=GSResNp;
		}
        else if(pFlag&&r<nLim)
	    {
          bQPDG=nQPDG;
          eMass=mNeut;
          rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z)+N-1);
          rMass=GSResNn;
   	    }
        else
	    {
#ifdef ppdebug
          cout<<"G4QNucleus::EvaporateBaryon: Photon # 2-B #, dM="<<totMass-GSMass<<endl;
#endif
          bQPDG=gQPDG;
          rQPDG=GetQPDG();
          eMass=0.;
          rMass=GSMass;
	    }
#ifdef ppdebug
        cout<<"G4QNucl::EvaporateBaryon: b="<<eMass<<bQPDG<<",r="<<rMass<<rQPDG<<endl;
#endif
      }

      if(three)           // Decay in two baryons + Residual Nucleus
	  {
        h1mom=G4LorentzVector(0.,0.,0.,eMass);
        h2mom=G4LorentzVector(0.,0.,0.,rMass);
        h3mom=G4LorentzVector(0.,0.,0.,fMass);
	    if(!DecayIn3(h1mom,h2mom,h3mom))
        {
#ifdef ppdebug
          cout<<"*G4QNucl::EvaporateBaryon:Decay M="<<totMass<<",b="<<eMass<<bQPDG
		      <<",f="<<fMass<<fQPDG<<",r="<<rMass<<rQPDG<<endl;
#endif
          return false;
	    }
        h1mom+=h3mom;
        bQPDG=dbQPDG;
	  }
      else
	  {
        if(eMass+rMass<totMass&&cntr<cntm)
        {
          h1mom=G4LorentzVector(0.,0.,0.,eMass);
          h2mom=G4LorentzVector(0.,0.,0.,rMass);
        }
        else if(totMass>GSMass)               // Photon if 2-Decay failed
	    {
#ifdef ppdebug
		  cout<<"G4QNucleus::EvaporateBaryon: Photon ### 2 ###, dM="<<totMass-GSMass<<endl;
#endif
          bQPDG=gQPDG;
          rQPDG=GetQPDG();
          h1mom=G4LorentzVector(0.,0.,0.,0.);
          h2mom=G4LorentzVector(0.,0.,0.,GSMass);      
	    }
        else
        {
	      cerr<<"***G4QNucleus::EvaporateBaryon: Cann't evaporate even gamma (1)"<<endl;
          return false;
        }
      }
	}
    else // ==> Decay in 3 Baryons + Residual is impossible at this point
	{
      if(secB&&gunFact)                        // Decay in 2 Baryons+Residual is possible
	  {
#ifdef ppdebug
        cout<<"G4QNucleus::EvaporateBaryon: Decay in 2 baryons"<<endl;
#endif
        G4double  nnLim=0.;
        if(nnFlag)nnLim+=N*(N-1)/2;
        G4double  nzLim=nnLim;
        if(npFlag)nzLim+=N*Z;
        G4double  zzLim=nzLim;
        if(ppFlag)zzLim+=Z*(Z-1)/2;
        G4double  nlLim=zzLim;
        if(nlFlag)nlLim+=N*S;
        G4double  zlLim=nlLim;
        if(plFlag)zlLim+=Z*S;
        G4double  llLim=zlLim;
        if(llFlag)llLim+=S*(S-1)/2;
        G4double r = llLim*G4UniformRand();
        if     (llFlag&&r>zlLim)
	    {
          dbQPDG= LLQPDG;
          eMass=mLamb;
          fMass=mLamb;
          rQPDG=G4QPDGCode(90000000+1000*(1000*(S-2)+Z)+N);
          rMass=GSResLL;
	    }
        else if(plFlag&&nlLim<r&&r<zlLim)
	    {
          dbQPDG= PLQPDG;
          eMass=mLamb;
          fMass=mProt;
          rQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z-1)+N);
          rMass=GSResPL;
	    }
        else if(nlFlag&&zzLim<r&&r<nlLim)
	    {
          dbQPDG= NLQPDG;
          eMass=mLamb;
          fMass=mNeut;
          rQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z)+N-1);
          rMass=GSResNL;
	    }
        else if(ppFlag&&nzLim<r&&r<zzLim)
	    {
          dbQPDG= PPQPDG;
          eMass=mProt;
          fMass=mProt;
          rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-2)+N);
          rMass=GSResPP;
	    }
        else if(npFlag&&nnLim<r&&r<nzLim)
	    {
          dbQPDG= NPQPDG;
          eMass=mNeut;
          fMass=mProt;
          rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-1)+N-1);
          rMass=GSResNP;
	    }
        else if(nnFlag&&r<nnLim)
	    {
          dbQPDG= NNQPDG;
          eMass=mNeut;
          fMass=mNeut;
          rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z)+N-2);
          rMass=GSResNN;
   	    }
        else
	    {
          if     (nFlag)
		  {
            bQPDG=nQPDG;
            rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z)+N-1);
            eMass=mNeut;
            rMass=GSResNn;
		  }
          else if(pFlag)
		  {
            bQPDG=pQPDG;
            rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-1)+N);
            eMass=mProt;
            rMass=GSResNp;
		  }
          else if(pFlag)
		  {
            bQPDG=lQPDG;
            rQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z)+N);
            eMass=mLamb;
            rMass=GSResNl;
		  }
          else
		  {
#ifdef ppdebug
            cout<<"G4QNucleus::EvaporateBaryon:Photon ### 3-Big ###,dM="<<totMass-GSMass<<endl;
#endif
            bQPDG=gQPDG;
            rQPDG=GetQPDG();
            eMass=0.;
            rMass=GSMass;
		  }
	    }

        h1mom=G4LorentzVector(0.,0.,0.,eMass);
        h2mom=G4LorentzVector(0.,0.,0.,rMass);
        h3mom=G4LorentzVector(0.,0.,0.,fMass);
	    if(!DecayIn3(h1mom,h2mom,h3mom))
        {
#ifdef ppdebug
          cout<<"*G4QNucl::EvaporateBaryon:Decay M="<<totMass<<",b="<<eMass<<bQPDG
		      <<",f="<<fMass<<fQPDG<<",r="<<rMass<<rQPDG<<endl;
#endif
          return false;
	    }
        h1mom+=h3mom;
        bQPDG=dbQPDG;
#ifdef ppdebug
        G4double sma=h1mom.m();
        G4double dma=sma-eMass-fMass;
		cout<<"G4QNucleus::EvaporBary:S="<<sma<<",em="<<eMass<<",fM="<<fMass<<",d="<<dma<<endl;
#endif
	  }
      else                                     // Only decay in Baryon+Residual is possible
      {
#ifdef ppdebug
        cout<<"G4QNucleus::EvaporateBaryon: Decay in Baryon+Resid"<<endl;
#endif
        if(!gunFact&&nFlag&&pFlag)
		//if(2>3)
		{
          if(Z>N)nFlag=false;
          else   pFlag=false;
		}
        G4double nLim=0.;
        if(nFlag&&mNeut+GSResNn<totMass)nLim+=N;
        G4double zLim=nLim;
        if(pFlag&&mProt+CBarr+GSResNp<totMass)zLim+=Z;
        G4double sLim=zLim;
        if(lFlag&&mLamb+GSResNl<totMass)sLim+=S;
        G4double r = sLim*G4UniformRand();
#ifdef ppdebug
        cout<<"G4QNucl::EvaporateBaryon:2Decay#2# r="<<r<<",nLim="<<nLim<<",zLim="<<zLim
            <<",sLim="<<sLim<<",nF="<<nFlag<<",pF="<<pFlag<<",lF="<<lFlag<<endl;
#endif
        if     (lFlag&&r>zLim)
	    {
          bQPDG=lQPDG;
          eMass=mLamb;
          rQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z)+N);
          rMass=GSResNl;
	    }
        else if(pFlag&&r>nLim&&r<zLim)
	    {
          bQPDG=pQPDG;
          eMass=mProt;
          rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-1)+N);
          rMass=GSResNp;
	    }
        else if(nFlag&&r<nLim)
	    {
          bQPDG=nQPDG;
          eMass=mNeut;
          rQPDG=G4QPDGCode(90000000+1000*(1000*S+Z)+N-1);
          rMass=GSResNn;
   	    }
        else
	    {
#ifdef ppdebug
          cout<<"G4QNucleus::EvaporateBaryon: Photon ### 4-Big ###, dM="<<totMass-GSMass<<endl;
#endif
          bQPDG=gQPDG;
          rQPDG=GetQPDG();
          eMass=0.;
          rMass=GSMass;
	    }
        h1mom=G4LorentzVector(0.,0.,0.,eMass);
        h2mom=G4LorentzVector(0.,0.,0.,rMass);      
      }
	}
	if(!DecayIn2(h1mom,h2mom))
    {
#ifdef ppdebug
      cout<<"*G4QNucleus::EvaporateBaryon: Decay M="<<totMass<<",b="<<bQPDG<<h1->GetQC()
		  <<",r="<<rQPDG<<h2->GetQC()<<endl;
#endif
      return false;
	}
    h1->SetQPDG(bQPDG);
    h2->SetQPDG(rQPDG);
    h1->Set4Momentum(h1mom);
    h2->Set4Momentum(h2mom);
#ifdef ppdebug
    cout<<"G4QNucleus::EvaporateBaryon: Evaporation is done for b="<<bQPDG<<h1->GetQC()
        <<",r="<<rQPDG<<h2->GetQC()<<endl;
#endif
	return true;
  }
  else if(a==1)
  {
#ifdef ppdebug
    cerr<<"***G4QNucleus::EvaporateBaryon: ??? A=1"<<endl;
#endif
    h1->SetQPDG(gQPDG);
    h1->Set4Momentum(G4LorentzVector(0.,0.,0.,0.));
    if     (Z>0)
    {
      h2->SetQPDG(pQPDG);
      h2->Set4Momentum(G4LorentzVector(0.,0.,0.,mProt));
    }
    else if(N>0)
    {
      h2->SetQPDG(nQPDG);
      h2->Set4Momentum(G4LorentzVector(0.,0.,0.,mNeut));
    }
    else if(S>0)
    {
      h2->SetQPDG(lQPDG);
      h2->Set4Momentum(G4LorentzVector(0.,0.,0.,mLamb));
    }
    else return false;
	return true;
  }
  if(a<1) cerr<<"***G4QNucleus::EvaporateBaryon: A="<<a<<endl;
  return false;
}











