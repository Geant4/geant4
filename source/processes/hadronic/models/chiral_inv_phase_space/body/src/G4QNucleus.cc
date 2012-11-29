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
// $Id$
//
//      ---------------- G4QNucleus ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Nuclei/Nuclear Environment used by CHIPS Model
// ---------------------------------------------------------------------
// Short description: a class describing properties of nuclei, which
// are necessary for the CHIPS Model.
// ---------------------------------------------------------------------


//#define debug
//#define pdebug 
//#define cldebug
//#define qdebug
//#define cldebug
//#define pardeb
//#define ppdebug

#include <algorithm>
#include <cmath>
#include <vector>

#include "G4QNucleus.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

// Static parameters definition
G4double G4QNucleus::freeNuc=0.1;       // probability to find quasiFreeBaryon on Surface
G4double G4QNucleus::freeDib=.05;       // probability to find quasiFreeDiBaryon on Surface
G4double G4QNucleus::clustProb=4.;      // clusterization probability in dense region
G4double G4QNucleus::mediRatio=1.;      // relative vacuum hadronization probability
G4double G4QNucleus::nucleonDistance=.8*fermi; // Distance between nucleons (0.8 fm) (Body)
G4double G4QNucleus::WoodSaxonSurf=.545*fermi; // WoodSaxon Surface Param (0.545 fm) (Body)

G4QNucleus::G4QNucleus(): G4QHadron(), Z(0), N(0), S(0), dZ(0), dN(0), dS(0), maxClust(0),
                          theNucleons(),currentNucleon(-1),
                          rho0(1.), radius(1.), Tb(), TbActive(false), RhoActive(false)
{
  probVect[0]=mediRatio;
  for(G4int i=1; i<256; i++) {probVect[i] = 0.;}
#ifdef pardeb
  G4cout<<"G4QNucleus::Constructor:(1) N="<<freeNuc<<", D="<<freeDib<<", W="<<clustProb
        <<", R="<<mediRatio<<G4endl;
#endif
}

G4QNucleus::G4QNucleus(G4int z, G4int n, G4int s_value) :
  G4QHadron(90000000+s_value*1000000+z*1000+n), Z(z),N(n),S(s_value), dZ(0),dN(0),dS(0), maxClust(0),
  theNucleons(), currentNucleon(-1), rho0(1.), radius(1.),
  Tb(), TbActive(false), RhoActive(false)
{
  probVect[0]=mediRatio;
  for(G4int i=1; i<256; i++) {probVect[i] = 0.;}
#ifdef debug
  G4cout<<"G4QNucleus::Construction By Z="<<z<<",N="<<n<<",S="<<s_value<<G4endl;
#endif
  SetZNSQC(z,n,s_value);
  G4QPDGCode nQPDG(90000000+S*1000000+Z*1000+N); // Not necessary (? look above)
#ifdef debug
  G4cout<<"G4QNucleus::ConstructionByZNS: nQPDG="<<nQPDG<<G4endl;
#endif
  G4double mass=nQPDG.GetNuclMass(Z,N,S);
#ifdef debug
  G4cout<<"G4QNucleus::ConstructionByZNS: mass="<<mass<<G4endl;
#endif
  SetQPDG(nQPDG);                                // Not necessary (? look above)
#ifdef debug
  G4cout<<"G4QNucleus::ConstructionByZNS: nQPDG set"<<G4endl;
#endif
  G4LorentzVector p(0.,0.,0.,mass);
  Set4Momentum(p);
  SetNFragments(0);
#ifdef debug
  G4cout<<"G4QNucleus::Constructor:(2) N="<<freeNuc<<", D="<<freeDib<<", W="<<clustProb
        <<", R="<<mediRatio<<G4endl;
#endif
}

G4QNucleus::G4QNucleus(G4int nucPDG):
  G4QHadron(nucPDG), maxClust(0), theNucleons(),
  currentNucleon(-1), rho0(1.), radius(1.), Tb(), TbActive(false), RhoActive(false)
{
  if(nucPDG==22) nucPDG=90000000;
  InitByPDG(nucPDG);
  G4LorentzVector p(0.,0.,0.,GetGSMass());
  Set4Momentum(p);
#ifdef pardeb
  G4cout<<"G4QNucleus::Constructor:(3) N="<<freeNuc<<", D="<<freeDib<<", W="<<clustProb
        <<", R="<<mediRatio<<", 4M="<<p<<G4endl;
#endif
}

G4QNucleus::G4QNucleus(G4LorentzVector p, G4int nucPDG):
  G4QHadron(nucPDG, p), maxClust(0), theNucleons(),
  currentNucleon(-1), rho0(1.), radius(1.), Tb(), TbActive(false), RhoActive(false)
{
  InitByPDG(nucPDG);
  Set4Momentum(p);
#ifdef pardeb
  G4cout<<"G4QNucleus::Constructor:(4) N="<<freeNuc<<", D="<<freeDib<<", W="<<clustProb
        <<", R="<<mediRatio<<", 4M="<<p<<G4endl;
#endif
}

G4QNucleus::G4QNucleus(G4int z, G4int n, G4int s_value, G4LorentzVector p) :
  G4QHadron(90000000+s_value*1000000+z*1000+n,p), Z(z),N(n),S(s_value), dZ(0),dN(0),dS(0), maxClust(0),
  theNucleons(), currentNucleon(-1), rho0(1.),radius(1.),
  Tb(), TbActive(false), RhoActive(false)
{
  probVect[0]=mediRatio;
  for(G4int i=1; i<256; i++) {probVect[i] = 0.;}
  Set4Momentum(p);
  SetNFragments(0);
  G4int ZNS=Z+N+S;
  G4QPDGCode nPDG(90000000+S*1000000+Z*1000+N);
  SetQPDG(nPDG);
  G4QContent nQC(N+ZNS,Z+ZNS,S,0,0,0);
  SetZNSQC(z,n,s_value);
#ifdef pardeb
  G4cout<<"G4QNucleus::Constructor:(5) N="<<freeNuc<<", D="<<freeDib<<", W="<<clustProb
        <<", R="<<mediRatio<<G4endl;
#endif
}

G4QNucleus::G4QNucleus(G4QContent nucQC):
  G4QHadron(nucQC), dZ(0), dN(0), dS(0), maxClust(0), theNucleons(), currentNucleon(-1),
  rho0(1.), radius(1.), Tb(), TbActive(false), RhoActive(false)
{
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
#ifdef debug
  G4cout<<"G4QNucleus::Construction By QC="<<nucQC<<G4endl;
#endif
  probVect[0]=mediRatio;
  for(G4int i=1; i<256; i++) {probVect[i] = 0.;}
  G4int u=nucQC.GetU()-nucQC.GetAU();
  G4int d=nucQC.GetD()-nucQC.GetAD();
  S = nucQC.GetS()-nucQC.GetAS();     // a#of LAMBDA's in the nucleus
  G4int du= d-u;                      // isotopic shift
  G4int b =(d+u+S)/3;                 // baryon number
  Z = (b-S-du)/2;                     // protons
  N = Z+du;                           // neutrons
  SetQC(nucQC);
#ifdef debug
  G4cout<<"G4QNucleus::ConstructionByQC: N="<<N<<",Z="<<Z<<",S="<<S<<G4endl;
#endif
  G4int nucPDG=90000000+S*1000000+Z*1000+N;
  G4QPDGCode nQPDG(nucPDG);
#ifdef debug
  G4cout<<"G4QNucleus::ConstructionByQC: nQPDG="<<nQPDG<<G4endl;
#endif
  G4double mass=nQPDG.GetNuclMass(Z,N,S);
  if(nucPDG==90000000)
  {
    if(nucQC.GetTot()) mass=mPi0;
    else               mass=0.;
  }
#ifdef debug
  G4cout<<"G4QNucleus::ConstructionByQC: mass="<<mass<<G4endl;
#endif
  SetQPDG(nQPDG);
#ifdef debug
  G4cout<<"G4QNucleus::ConstructionByQC: nQPDG set"<<G4endl;
#endif
  G4LorentzVector p(0.,0.,0.,mass);
  Set4Momentum(p);
  SetNFragments(0);
#ifdef pardeb
  G4cout<<"G4QNucleus::Constructor:(6) N="<<freeNuc<<", D="<<freeDib<<", W="<<clustProb
        <<", R="<<mediRatio<<G4endl;
#endif
}

G4QNucleus::G4QNucleus(G4QContent nucQC, G4LorentzVector p):
  G4QHadron(nucQC,p), dZ(0), dN(0), dS(0), maxClust(0), theNucleons(), currentNucleon(-1),
  rho0(1.), radius(1.), Tb(), TbActive(false), RhoActive(false)
{
#ifdef debug
  G4cout<<"G4QNucleus::(LV)Construction By QC="<<nucQC<<G4endl;
#endif
  probVect[0]=mediRatio;
  for(G4int i=1; i<256; i++) {probVect[i] = 0.;}
  Set4Momentum(p);
  G4int u=nucQC.GetU()-nucQC.GetAU();
  G4int d=nucQC.GetD()-nucQC.GetAD();
  S = nucQC.GetS()-nucQC.GetAS();     // a#of LAMBDA's in the nucleus
  G4int du= d-u;                      // isotopic shift
  G4int b =(d+u+S)/3;                 // baryon number
  Z = (b-S-du)/2;                     // protons
  N = Z+du;                           // neutrons
  SetQC(nucQC);
#ifdef debug
  G4cout<<"G4QNucleus::(LV)ConstructionByQC: N="<<N<<",Z="<<Z<<",S="<<S<<G4endl;
#endif
  G4QPDGCode nPDG(90000000+S*1000000+Z*1000+N);
  SetQPDG(nPDG);
  SetNFragments(0);
#ifdef pardeb
  G4cout<<"G4QNucleus::Constructor:(7) N="<<freeNuc<<", D="<<freeDib<<", W="<<clustProb
        <<", R="<<mediRatio<<G4endl;
#endif
}

G4QNucleus::G4QNucleus(G4QNucleus* right, G4bool cop3D) : currentNucleon(-1)
{
  Z             = right->Z;
  N             = right->N;
  S             = right->S;
  dZ            = right->dZ;
  dN            = right->dN;
  dS            = right->dS;
  maxClust      = right->maxClust;
  for(G4int i=0; i<=maxClust; i++) probVect[i] = right->probVect[i];
  probVect[254] = right->probVect[254];
  probVect[255] = right->probVect[255];
  Tb            = right->Tb;
  TbActive      = right->TbActive;
  RhoActive     = right->RhoActive;
  Set4Momentum   (right->Get4Momentum());
  SetQPDG        (right->GetQPDG());
  SetQC          (right->GetQC());
  SetNFragments  (right->GetNFragments());
  rho0        = right->rho0;
  radius      = right->radius;
  if(cop3D)
  {
    G4int nn=right->theNucleons.size();
    for(G4int i=0; i<nn; ++i)
    {
      G4QHadron* nucleon = new G4QHadron(right->theNucleons[i]);
      theNucleons.push_back(nucleon);
    }
  }
#ifdef pardeb
  G4cout<<"G4QNucleus::Constructor:(8) N="<<freeNuc<<", D="<<freeDib<<", W="<<clustProb
        <<", R="<<mediRatio<<G4endl;
#endif
}

G4QNucleus::G4QNucleus(const G4QNucleus &right, G4bool cop3D):
  G4QHadron(), currentNucleon(-1)
{
  Z             = right.Z;
  N             = right.N;
  S             = right.S;
  dZ            = right.dZ;
  dN            = right.dN;
  dS            = right.dS;
  maxClust      = right.maxClust;
  for(G4int i=0; i<=maxClust; i++) probVect[i] = right.probVect[i];
  probVect[254] = right.probVect[254];
  probVect[255] = right.probVect[255];
  Tb            = right.Tb;
  TbActive      = right.TbActive;
  RhoActive     = right.RhoActive;
  Set4Momentum   (right.Get4Momentum());
  SetQPDG        (right.GetQPDG());
  SetQC          (right.GetQC());
  SetNFragments  (right.GetNFragments());
  rho0        = right.rho0;
  radius      = right.radius;
  if(cop3D)
  {
    G4int nn=right.theNucleons.size();
    for(G4int i=0; i<nn; ++i)
    {
      G4QHadron* nucleon = new G4QHadron(right.theNucleons[i]);
      theNucleons.push_back(nucleon);
    }
  }
#ifdef pardeb
  G4cout<<"G4QNucleus::Constructor:(9) N="<<freeNuc<<", D="<<freeDib<<", W="<<clustProb
        <<", R="<<mediRatio<<G4endl;
#endif
}

// Assignment operator
const G4QNucleus& G4QNucleus::operator=(const G4QNucleus& right)
{
  if(this != &right)                          // Beware of self assignment
  {
    currentNucleon= -1;
    TbActive      = right.TbActive;
    Tb            = right.Tb;
    RhoActive     = right.RhoActive;
    rho0          = right.rho0;
    radius        = right.radius;
    G4int nn      = right.theNucleons.size();
    for(G4int i=0; i < nn; ++i)
    {
      G4QHadron* nucleon = new G4QHadron(right.theNucleons[i]);
      theNucleons.push_back(nucleon);
    }
    Set4Momentum   (right.Get4Momentum());
    SetQPDG        (right.GetQPDG());
    SetQC          (right.GetQC());
    SetNFragments  (right.GetNFragments());
    Z             = right.Z;
    N             = right.N;
    S             = right.S;
    dZ            = right.dZ;
    dN            = right.dN;
    dS            = right.dS;
    maxClust      = right.maxClust;
    for(G4int i=0; i<=maxClust; i++) probVect[i] = right.probVect[i];
    probVect[254] = right.probVect[254];
    probVect[255] = right.probVect[255];
  }
  return *this;
}

G4QNucleus::~G4QNucleus()
{
  for_each(theNucleons.begin(),theNucleons.end(),DeleteQHadron());
}

// Fill the private parameters
void G4QNucleus::SetParameters(G4double a,G4double b, G4double c, G4double d, G4double e)
{
  freeNuc=a; 
  freeDib=b; 
  clustProb=c;
  mediRatio=d;
  nucleonDistance=e;
}

// Standard output for QNucleus {Z - a#of protons, N - a#of neutrons, S - a#of lambdas}
ostream& operator<<(ostream& lhs, G4QNucleus& rhs)
{
 lhs<<"{Z="<<rhs.GetZ()<<",N="<<rhs.GetN()<<",S="<<rhs.GetS()<<",M="<<rhs.GetGSMass()<<"}";
 return lhs;
}

// Standard output for QNucleus {Z - a#of protons, N - a#of neutrons, S - a#of lambdas}
ostream& operator<<(ostream& lhs, const G4QNucleus& rhs)
{
  lhs<<"{Z="<<rhs.GetZ()<<",N="<<rhs.GetN()<<",S="<<rhs.GetS()<< "}";
  return lhs;
}

// Init existing nucleus by new PDG Code
void G4QNucleus::InitByPDG(G4int nucPDG)
{
  static const G4int NUCPDG  = 90000000;
#ifdef debug
  G4cout<<"G4QNucleus::InitByPDG: >Called< PDG="<<nucPDG<<G4endl;
#endif
  dZ=0;
  dN=0;
  dS=0;
  probVect[0]=mediRatio;                        // init Vacuum/Medium probability
  for(G4int i=1; i<256; i++) {probVect[i] = 0.;}
  //std::uninitialized_fill( probVect+1, probVect+256, 0.0 ); // Worse in performance!
  if(nucPDG<80000000) nucPDG=HadrToNucPDG(nucPDG); // Convert HadrPDGCode to NucPDGCode
  G4int s_value=0;
  G4int z=0;
  G4int n=0;
  if(nucPDG>80000000 && nucPDG<100000000) // Try to convert the NUCCoding to PDGCoding
  {
    G4QPDGCode(22).ConvertPDGToZNS(nucPDG, z, n, s_value);
    Z  =z;
    N  =n;
    S  =s_value;
#ifdef debug
    G4cout<<"G4QNucleus::InitByPDG:Z="<<Z<<",N="<<N<<",S="<<S<<G4endl;
#endif
    SetZNSQC(Z,N,S);  // @@ ??
    G4QPDGCode nPDG(nucPDG);
    G4double PDGMass=0.;
    if(nucPDG!=NUCPDG) PDGMass=nPDG.GetMass();
    SetQPDG(nPDG);
    G4LorentzVector p(0.,0.,0.,PDGMass);
    Set4Momentum(p);
    SetNFragments(0);
#ifdef debug
    G4cout<<"G4QNucleus::InitByPDG:->QPDG="<<nPDG<<": 4M="<<p<<G4endl;
#endif
  }
  else
  {
    G4cerr<<"***G4QNucleus::InitByPDG:Initialized by not nuclear PDGCode="<<nucPDG<<G4endl;
    //throw G4QException("G4QNucleus::InitByPDG:PDGCode can't be converted to NucPDGCode");
  }
}
// End of "InitByPDG"

// Calculate probabilities of clusters and return the maximum baryon number of clusters
G4int G4QNucleus::UpdateClusters(G4bool din) // din true means use only dense nuclear part
{
  //static const G4double r0 = 1.1;          // fm, for nuclear radius: r=r0*A^(1/3)
  //static const G4double del= .55;          // fm, for a difused surface of the nucleus
  //static const G4double rCl= 2.0;          // clusterization radius @@??
  //static const G4double freeibuc = 0.10;   // probab. of the quasi-free baryon on surface
  //static const G4double freeDib = 0.05;    // probab. of the quasi-free dibar. on surface
  //static const G4double clustProb = 4.0;   // clusterization probability in dense region
  //static const G4double prQ = 1.0;         // relative probability for a Quasmon
  //static const G4double prQ = 0.;          //@@for pi@@relative probability for Quasmon
  //G4double probSInt[254];                    // integratedStaticProbabilities @@ not used
  probVect[0]=mediRatio;
  for (G4int in=1; in<256; in++) probVect[in]=0.; // Make preinit to avoid the postinit
  //probSInt[0]=0;                             // integrated static probabilities
  dZ=0;
  dN=0;
  dS=0;
  G4int a = Z + N + S;                       // atomic number
#ifdef debug
  G4cout<<"G4QN::UpdateCl:A="<<a<<"(Z="<<Z<<",N="<<N<<",S="<<S<<"),mR="<<mediRatio<<G4endl;
#endif
  G4double A=a;
  if(A<=0.)
  {
#ifdef debug
    G4cout<<"***G4QNucleus::UpdateClusters:No clusters can be calculated as A="<<A<<G4endl;
#endif
    return 0;
  }
  G4double surf=freeNuc+freeDib;             // surface relative population
  G4double surA=A*surf;                      // surface absolute population
  G4int sA=static_cast<G4int>(surA);
  if(surf>0.||surf<1.)sA=RandomizeBinom(surf,a); // randomize SurfaceNucleons by Binomial
#ifdef debug
  G4cout<<"G4QN::UpdateCl:surf="<<surf<<"= N="<<freeNuc<<"+D="<<freeDib<<",A="<<sA<<G4endl;
#endif
  G4int dA=a-sA;                             // a#of nucleons in dense part of the nucleus
  if (din && dA<2 && a>2)
  {
    dA=2;
    sA=a-2;
  }
#ifdef debug
  G4cout<<"G4QN::UpdtC:dA="<<dA<<",A="<<A<<",s="<<surf<<",S="<<sA<<",C="<<maxClust<<G4endl;
#endif
  G4int maxi=1;                              // A#of elements filled by the progran
  G4double pA=0.;
  G4double uA=0.;
  if(surf>0.)
  {
    pA=0.5*freeDib*sA/surf; //@@Randomize(?)// a#of quasi-free Nucleon Pairs on the surface
    uA=sA-pA-pA;                             // a#of quasi-free nucleons on Nuclear Surface
  }
  uA=uA/A;                                   // Normalization of probability
  pA=pA/A;
  G4double sum =0.;
  if(dA<2)                                   // There is no dense phase at all
  {
    //probVect[1]= dA/A;                       // a#of quasi-free nucleons (only dense)
    //probVect[1]= (uA+dA)/A;                  // a#of quasi-free nucleons (different norm)
    probVect[1]= uA+dA/A;                    // a#of quasi-free nucleons (correct)
    sum = probVect[1];
    //probSInt[1]=sum;                         // integrated static probabilities
    maxi=2;
    probVect[254]= 0;                        // a#of dense nucleons (correct)
    if(A>1 && pA>0.)
    {
      //probVect[2]= (pA+pA)/A/(A-1);          // a#of quasi-free "dibaryons" (correct)
      probVect[2]= pA;                       // a#of quasi-free "dibaryons" (correct)
      //probVect[2]= 0;                        // a#of quasi-free "dibaryons" (only dense)
      sum+= probVect[2]+probVect[2];
      //probSInt[2]=sum;                       // integrated static probabilities
      maxi=3;
      probVect[255]= 0;                      // a#of dense "dibaryons" (correct)
    }
#ifdef debug
    G4cout<<"G4QNucleus::UpdateClust:Only quasi-free nucleons pV[1]="<<probVect[1]<<G4endl;
#endif
  }
  else
  {
    G4double wrd=clustProb/dA;               // relative volume of clusterization (omega)
    G4double sud=pow(1.+wrd,dA-1);           // normalization factor for the dense region
    // dA=C*Sum_k=1-A[n*C^A_k*wrd^(k-1)]=C*dA*(1+wrd)^(dA-1) => C=1/sud, sud=(1+wrd)^(dA-1)
    // =1
    G4double rd= dA/sud/A;
    //G4double comb=A;
    //G4double prb=rd;                              // (only dense)
    G4double prb=rd+uA;
    sum =prb;
#ifdef debug
   G4cout<<"G4QNucl::UpdateCl:sud="<<sud<<",v[1]=s="<<sum<<",dA="<<dA<<",uA="<<uA<<G4endl;
#endif
    //probVect[1]= prb/comb;                   // a#of quasi-free nucleons (correct)
    //probVect[254]= rd/comb;                  // a#of dense nucleons (correct)
    probVect[1]= prb;                        // a#of quasi-free nucleons (correct)
    probVect[254]= rd;                       // a#of dense nucleons (correct)
    //probSInt[1]=sum;                         // integrated static probabilities
    // =2
    rd*=wrd*(dA-1.)/2;
    //comb*=(A-1.)/2;
    //prb=rd;                                  // (only dense)
    prb=rd+pA;
    sum+=prb+prb;
#ifdef debug
    G4cout<<"G4QNucl::UpdateCl:sud="<<sud<<",v[2]="<<prb<<",s="<<sum<<",pA="<<pA<<G4endl;
#endif
    //probVect[2]= prb/comb;                   // a#of quasi-free "dibaryons" (correct)
    //probVect[255]= rd/comb;                  // a#of dense "dibaryons" (correct)
    probVect[2]= prb;                        // a#of quasi-free "dibaryons" (correct)
    probVect[255]= rd;                       // a#of dense "dibaryons" (correct)
    //probSInt[2]=sum;                         // integrated static probabilities
    // >2
    maxi=3;
#ifdef debug
    G4cout<<"G4QNucleus::UpdateClusters:p1="<<probVect[1]<<", p2="<<probVect[2]<<",sA="<<sA
          <<",uA="<<uA<<",pA="<<pA<<",wrd="<<wrd<<",sud="<<sud<<G4endl;
#endif
    if(dA>2)
    {
      ///////////G4double itA=A+1.; 
      G4double idA=dA+1.; 
      G4int dLim=dA;
      if(maxClust<dA) dLim=maxClust;
      for (int i=3; i<=dLim; i++)
      {
        rd*=wrd*(idA-i)/i;
        sum+=rd*i;
#ifdef debug
     G4cout<<"G4QNucleus::UpdateCl:sud="<<sud<<", v["<<i<<"]="<<rd<<", s="<<sum<<G4endl;
#endif
        //comb*=(itA-i)/i;
        //probVect[i]=rd/comb;                 // Divide by sum of combinations for N+Z+S
        probVect[i]=rd;                      // Comb's for N,Z,S are canceled later(G4QNuc)
        //probSInt[i]=sum;                     // integrated static probabilities
        maxi=i+1;
#ifdef debug
        G4cout<<"G4QNucleus::UpdateCl:Cluster of "<<i<<" baryons,pV="<<probVect[i]<<G4endl;
#endif
      }
    }
    dS = S;                                  // @@ Lambdas are always in the dense region
    dZ = static_cast<int>(static_cast<double>((dA-dS)*Z)/(Z+N) + 0.5);
    dN = dA - dZ;
  }
#ifdef debug
  G4cout<<"G4QNucleus::UpdateClusters: Sum of weighted probabilities s="<<sum<<G4endl;
#endif
  maxClust=maxi-1;
  //for (G4int j=maxi; j<255; j++) probVect[j]=0.;//Make the rest to be 0 [preinited above]
  // =----------------= From here probability randomization starts =---------------=
  //  G4int rA=a;                              // Residual number of nucleons
  //#ifdef debug
  //G4cout<<"G4QNuc::UpdateClust:A="<<A<<",M="<<k<<",P1="<<probVect[1]<<",P2="<<probVect[2]
  //      <<G4endl;
  //#endif
  //if (k>1) for (j=k; j>1; j--)               // nucleons are not randomized
  //{
  //  G4int jmax=rA/j;                         // Max number of this kind of clusters
  //  if (jmax)
  //  {
  //    G4double prob=probVect[j]/probSInt[j]; // Probab of the cluster in the dest nucleus
  //#ifdef debug
  //    G4cout<<"G4QNucl::UpdateClusters: j="<<j<<",sP="<<probVect[j]<<",iP="<<probSInt[j]
  //          <<G4endl;
  //#endif
  // G4int m=RandomizeBinom(prob,jmax);     // A#of clusters of this type
  //    if(m)
  //    {
  //      probVect[j]=m;
  //      rA-=m*j;
  //    }
  //    else
  //    {
  //      probVect[j]=0.;
  //      if(j==maxClust) maxClust--;
  //    }
  //#ifdef debug
  //    G4cout<<"G4QNucl::UpdateClust:p="<<prob<<",r="<<rA<<",m="<<jmax<<",P="<<probVect[j]
  //          <<G4endl;
  //#endif
  //  }
  //  else
  //  {
  //    probVect[j]=0.;
  //    if(j==maxClust) maxClust--;
  //  }
  //}
  //probVect[1]=rA;
  // =------------------= From here probability randomization starts =-------------------=
  return maxClust;
}
// End of "UpdateClusters"

// Reduce the 3D Nucleus by the used nucleon + update nucleon's energies (in LS!)
void G4QNucleus::SubtractNucleon(G4QHadron* uNuc)
{
  G4int NotFound=true;                              // Not found flag
  G4QHadronVector::iterator u;                      // iterator of the used nucleon
  for(u=theNucleons.begin(); u!=theNucleons.end(); u++)
  {
#ifdef debug
    G4cout<<"G4QNucleus::SubtractNucleon: LOOP 4M="<<(*u)->Get4Momentum()<<G4endl;
#endif
    if (uNuc==*u)                                   // Find uNuceon-pointer
    {
      NotFound=false;
      break;
    }
  }
//  if(NotFound) throw G4QException("G4QNucleus::SubtractNucleon: The nucleon isn't found");
  if (NotFound) G4Exception("G4QNucleus::SubtractNucleon()", "HAD_CHPS_0000",
                            FatalException, "The nucleon isn't found");
  else
  {
    G4int tPDG=GetPDGCode();                    // Nucleus PDG before the subtraction
    G4LorentzVector t4M=Get4Momentum();         // Nucleus 4-mom before the subtraction
#ifdef debug
    G4cout<<"G4QNucleus::SubtractNucleon: InitialNucleus 4M="<<t4M<<", PDG="<<tPDG<<", nN="
          <<theNucleons.size()<<G4endl;
#endif
    G4int uPDG=(*u)->GetPDGCode();              // PDG code of the subtracted nucleon
    G4LorentzVector u4M=(*u)->Get4Momentum();   // 4-momentum of the subtracted nucleon
#ifdef debug
    G4cout<<"G4QNucleus::SubtractNucleon: subtractNucleon 4M="<<u4M<<",PDG="<<uPDG<<G4endl;
#endif
    delete *u;                                  // Delete the nucleon as an object
    theNucleons.erase(u);                       // exclude the nucleon pointer from the HV
    --currentNucleon;                           // Continue selection from theSame position
    t4M-=u4M;                                   // Update the nucleus 4-momentum VALUE
    if     (uPDG==2212) tPDG-=1000;             // Reduce the nucleus PDG Code by a proton
    else if(uPDG==2112) tPDG--;                 // Reduce the nucleus PDG Code by a neutron
    else
    {
      // G4cerr<<"***G4QNucleus::SubtractNucleon: Unexpected Nucleon PDGCode ="<<uPDG<<G4endl;
      // throw G4QException("G4QNucleus::SubtractNucleon: Impossible nucleon PDG Code");
      G4ExceptionDescription ed;
      ed << "Impossible nucleon PDG Code: Unexpected Nucleon PDGCode ="
         << uPDG << G4endl;
      G4Exception("G4QNucleus::SubtractNucleon()", "HAD_CHPS_0001",
                  FatalException, ed);
    }
#ifdef debug
    G4cout<<"G4QNucleus::SubtractNucleon: theResidualNucleus PDG="<<tPDG<<", 4M="<<t4M
          <<", nN="<<theNucleons.size()<<G4endl;
#endif
    InitByPDG(tPDG);                            // Reinitialize the nucleus, not 3D nucleus
    theMomentum=t4M;                            // Fill the residual 4-momentum
    //#ifdef debug
    G4double mR2=sqr(GetGSMass());              // Real squared residual nucleus mass 
    G4double tM2=t4M.m2();                      // Squared residual nucleus mass from 4M 
#ifdef debug
    G4cout<<"G4QNucleus::SubtractNucleon: rAm2="<<mR2<<" =? 4Mm2="<<tM2<<G4endl;
    G4int cnt=0;                                // Counter of nucleons for print
#endif
    if(std::fabs(mR2-tM2)>.01)G4cout<<"*G4QNucleus::SubNucleon:rM="<<mR2<<"#"<<tM2<<G4endl;
    //#endif
    G4double tE=t4M.e();                        // Energy of the residual nucleus (in CM!)
    G4double m2p=sqr(G4QNucleus(tPDG-1000).GetGSMass()); // subResid. nuclearM2 for protons
    G4double m2n=sqr(G4QNucleus(tPDG-1).GetGSMass()); // subResidual nuclearM2 for neutrons
    for(u=theNucleons.begin(); u!=theNucleons.end(); u++) // Correct the nucleon's energies
    {
      G4LorentzVector n4M=(*u)->Get4Momentum(); // 4-mom of the current nucleon
      G4double srP2=(t4M-n4M).vect().mag2();    // p2 of the subResNucleus
      G4double m2_value=m2n;                    // default subResNucleusM2 (for neutrons) 
      if((*u)->GetPDGCode()==2212) m2_value=m2p;// change it to subResNucleusM2 for protons
      G4double srE=std::sqrt(srP2+m2_value);    // Energy of the subResNucleus
#ifdef debug
      G4cout<<"G4QNucleus::SubtractNucleon:#"<<cnt++<<", correctedEnergy="<<tE-srE<<G4endl;
#endif
      n4M.setE(tE-srE);                         // Update the energy of the nucleon
      (*u)->Set4Momentum(n4M);                  // Update the 4-momentum of the nucleon
    }
  }
#ifdef debug
  G4cout<<"G4QNucleus::SubtractNucleon:ResNuc4M="<<theMomentum<<",Z="<<Z<<",N="<<N<<G4endl;
#endif
}

// Delete all residual nucleons
void G4QNucleus::DeleteNucleons()
{
  G4QHadronVector::iterator u;                      // iterator for the nucleons
  for(u=theNucleons.begin(); u!=theNucleons.end(); u++) delete *u;
  theMomentum=G4LorentzVector(0.,0.,0.,0.);
}

// Reduce nucleus by emitted cluster with PDG Code cPDG
void G4QNucleus::Reduce(G4int cPDG)
{
  static const G4int NUCPDG=90000000;
  if(cPDG>80000000&&cPDG!=NUCPDG)
  {
    G4int curPDG=GetPDG();
    G4int newPDG=curPDG-cPDG+NUCPDG;             // PDG Code of Residual Nucleus
    if(newPDG==NUCPDG) InitByPDG(NUCPDG);        // Empty
    else
    {
      //if(abs(newPDG)<NUCPDG)
      //{
      //  G4cerr<<"***G4QNucleus::Reduce:iPDG="<<curPDG<<"=newPDG="<<newPDG<<"+cPDG="<<cPDG
      //        <<G4endl;
      //  throw G4QException("*E*:::G4QNucleus::Reduce: Abnormal Nuclear Reduction");
      //}
      InitByPDG(newPDG);                         // Reinit the Nucleus
    }
  }
  else if(cPDG!=NUCPDG) G4cerr<<"***G4QN::Reduce:Subtract not nuclear PDGC="<<cPDG<<G4endl;
  // in case of cPDG=90000000 - subtract nothing
}

// Increase nucleus by cluster with PDG Code cPDG (4-mom is optional)
void G4QNucleus::Increase(G4int cPDG, G4LorentzVector c4M)
{
  static const G4int NUCPDG=90000000;
  if(cPDG>80000000&&cPDG!=NUCPDG)
  {
    G4int newPDG=GetPDG()+cPDG-NUCPDG;        // PDG Code of the New Nucleus
    InitByPDG(newPDG);                        // Reinit the Nucleus
    if (c4M!=G4LorentzVector(0.,0.,0.,0.))
    {
      G4LorentzVector t4M = Get4Momentum();   // 4Mom of the nucleus
      t4M +=c4M;
      Set4Momentum(t4M);
    }
  }
  else G4cerr<<"***G4QNucleus::Increase: PDGCode="<<cPDG<<",4M="<<c4M<<G4endl;
}

// Increase nucleus by Quasmon with Quark Content qQC (4-mom is optional)
void G4QNucleus::Increase(G4QContent qQC, G4LorentzVector q4M)
{
    G4LorentzVector t4M = Get4Momentum();     // 4Mom of the old nucleus
    G4QContent  newQC   = GetQC()+qQC;        // Quark Content of the New Nucleus
    InitByQC(newQC);                          // Reinit the Nucleus
    t4M +=q4M;
    Set4Momentum(t4M);                        // 4Mom of the new nucleus
}

// Set Quark Content, using Z,N,S of nucleus
void G4QNucleus::SetZNSQC(G4int z, G4int n, G4int s_value)
{
  G4int zns=z+n+s_value;
  G4int Dq=n+zns;
  G4int Uq=z+zns;
  G4int Sq=s_value;
  if      (Dq<0&&Uq<0&&Sq<0)SetQC(G4QContent( 0, 0, 0,-Dq,-Uq,-Sq));
  else if (Uq<0&&Sq<0)      SetQC(G4QContent(Dq, 0, 0,  0,-Uq,-Sq));
  else if (Dq<0&&Sq<0)      SetQC(G4QContent( 0,Uq, 0,-Dq,  0,-Sq));
  else if (Dq<0&&Uq<0)      SetQC(G4QContent( 0, 0,Sq,-Dq,-Uq,  0));
  else if (Uq<0)            SetQC(G4QContent(Dq, 0,Sq,  0,-Uq,  0));
  else if (Sq<0)            SetQC(G4QContent(Dq,Uq, 0,  0,  0,-Sq));
  else if (Dq<0)            SetQC(G4QContent(0 ,Uq,Sq,-Dq,  0,  0));
  else                      SetQC(G4QContent(Dq,Uq,Sq,  0,  0,  0));
}
  
// Tests if it is possible to split one Baryon (n,p,Lambda) or alpha from the Nucleus
G4int G4QNucleus::SplitBaryon()
{
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QContent deutQC(3,3,0,0,0,0);
  static const G4QContent alphQC(6,6,0,0,0,0);
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  G4int     baryn=GetA();                         // Baryon Number of the Nucleus
  if(baryn<2) return 0;
  //G4double   totM=GetGSMass();                    // GS Mass value of the Nucleus
  G4double   totM=Get4Momentum().m();             // Real Mass value of the Nucleus
  G4QContent valQC=GetQCZNS();                    // Quark Content of the Nucleus
#ifdef debug
  G4cout<<"G4QNucleus::SplitBaryon: B="<<baryn<<", M="<<totM<<valQC<<G4endl;
#endif
  G4int NQ=valQC.GetN();
  if(NQ)                                          // ===> "Can try to split a neutron" case
  {
    G4QContent resQC=valQC-neutQC;                // QC of Residual for the Neutron
    G4int    resPDG=resQC.GetSPDGCode();          // PDG of Residual for the Neutron
    G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
    G4double sM=resMas+mNeut;
#ifdef debug
    G4cout<<"G4QNucleus::SplitBaryon: (neutron),sM="<<sM<<",d="<<totM-sM<<G4endl;
#endif
    if(sM<totM+.001) return 2112;
  }
  G4int PQ=valQC.GetP();
  if(PQ)                                          // ===> "Can try to split a proton" case
  {
    G4QContent resQC=valQC-protQC;                // QC of Residual for the Proton
    G4int    resPDG=resQC.GetSPDGCode();          // PDG of Residual for the Proton
    G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
    G4double CB=CoulombBarrier(1,1);              // Coulomb Barrier for the proton
    G4double sM=resMas+mProt+CB;
    /////////G4double sM=resMas+mProt;
#ifdef debug
    G4cout<<"G4QNucleus::SplitBaryon: (proton),sM="<<sM<<",d="<<totM-sM<<G4endl;
#endif
    if(sM<totM+.001) return 2212;
  }
  G4int LQ=valQC.GetL();
  if(LQ)                                          // ===> "Can try to split a lambda" case
  {
    G4QContent resQC=valQC-lambQC;                // QC of Residual for the Lambda
    G4int    resPDG=resQC.GetSPDGCode();          // PDG of Residual for the Lambda
    G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
    G4double sM=resMas+mLamb;
#ifdef debug
    G4cout<<"G4QNucleus::SplitBaryon: (lambda),sM="<<sM<<",d="<<totM-sM<<G4endl;
#endif
    if(sM<totM+.001) return 3122;
  }
  G4int AQ=NQ+PQ+LQ;
  if(NQ>0&&PQ>0&&AQ>2)                            // ===> "Can try to split deuteron" case
  {
    G4QContent resQC=valQC-deutQC;                // QC of Residual for the Deuteron
    G4int    resPDG=resQC.GetSPDGCode();          // PDG of Residual for the Deuteron
    G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
    G4double CB=CoulombBarrier(1,2);              // Coulomb Barrier for the Deuteron
    G4double sM=resMas+mDeut+CB;
    //G4double sM=resMas+mDeut;
#ifdef debug
    G4cout<<"G4QNucleus::SplitBaryon: (deuteron),sM="<<sM<<",d="<<totM-sM<<G4endl;
#endif
    if(sM<totM+.001) return 90001001;
  }
  if(NQ>1&&PQ>1&&AQ>4)                            // ===> "Can try to split an alpha" case
  {
    G4QContent resQC=valQC-alphQC;                // QC of Residual for the Alpha
    G4int    resPDG=resQC.GetSPDGCode();          // PDG of Residual for the Alpha
    G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
    G4double CB=CoulombBarrier(2,4);              // Coulomb Barrier for the Alpha
    G4double sM=resMas+mAlph;
    if(NQ!=4||PQ!=4) sM+=CB;
#ifdef debug
    G4cout<<"G4QNucleus::SplitBaryon: (alpha),sM="<<sM<<",d="<<totM-sM<<G4endl;
#endif
    if(sM<totM+.001) return 90002002;
  }
  return 0;
}
  
// Tests if it is possible to split two Baryons (nn,np,pp,Ln,Lp,LL) from the Nucleus
G4bool G4QNucleus::Split2Baryons()
{
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  G4int     baryn=GetA();                        // Baryon Number of the Nucleus
  if(baryn<3) return false;
  G4double   totM=theMomentum.m();            // Real Mass value of the Nucleus
  G4QContent valQC=GetQCZNS();                   // Quark Content of the Nucleus
#ifdef debug
  G4cout<<"G4QNucleus::Split2Baryons: B="<<baryn<<", M="<<totM<<valQC<<G4endl;
#endif
  G4int NQ=valQC.GetN();
  if(NQ>1)                                       // ===> "Can try to split 2 neutrons" case
  {
    G4QContent resQC=valQC-neutQC-neutQC;        // QC of ResidNucleus for the Two Neutrons
    G4int    resPDG=resQC.GetSPDGCode();         // PDG of ResidNucleus for 2 Neutrons
    G4double resMas=G4QPDGCode(resPDG).GetMass();// GS Mass of the Residual Nucleus
    G4double sM=resMas+mNeut+mNeut;
#ifdef debug
    G4cout<<"G4QNucleus::Split2Baryons: (2 neutrons), sM="<<sM<<", d="<<totM-sM<<G4endl;
#endif
    if(sM<totM) return true;
  }
  G4int PQ=valQC.GetP();
  if(PQ>1)                                       // ===> "Can try to split 2 protons" case
  {
    G4QContent resQC=valQC-protQC-protQC;        // QC of ResidualNucleus for 2 Protons
    G4int    resPDG=resQC.GetSPDGCode();         // PDG of Residual Nucleus for 2 Proton
    G4double resMas=G4QPDGCode(resPDG).GetMass();// GS Mass of the Residual Nucleus
    G4double sM=resMas+mProt+mProt;
#ifdef debug
    G4cout<<"G4QNucleus::Split2Baryons: (2 protons), sM="<<sM<<", d="<<totM-sM<<G4endl;
#endif
    if(sM<totM) return true;
  }
  if(PQ&&NQ)                                     // ===> "Can try to split proton+neutron"
  {
    G4QContent resQC=valQC-protQC-neutQC;        // QC of ResidNucleus for Proton+Neutron
    G4int    resPDG=resQC.GetSPDGCode();         // PDG of ResidNucleus for Proton+Neutron
    G4double resMas=G4QPDGCode(resPDG).GetMass();// GS Mass of the Residual Nucleus
    G4double sM=resMas+mProt+mNeut;
#ifdef debug
    G4cout<<"G4QNucleus::Split2Baryons:(proton+neutron), sM="<<sM<<", d="<<totM-sM<<G4endl;
#endif
    if(sM<totM) return true;
  }
  G4int LQ=valQC.GetL();
  if(LQ&&NQ)                                     // ===> "Can try to split lambda+neutron"
  {
    G4QContent resQC=valQC-lambQC-neutQC;        // QC of ResidNucleus for Lambda+Neutron
    G4int    resPDG=resQC.GetSPDGCode();         // PDG of ResidNucleus for Lambda+Neutron
    G4double resMas=G4QPDGCode(resPDG).GetMass();// GS Mass of the Residual Nucleus
    G4double sM=resMas+mLamb+mNeut;
#ifdef debug
    G4cout<<"G4QNucleus::Split2Baryons:(lambda+neutron), sM="<<sM<<", d="<<totM-sM<<G4endl;
#endif
    if(sM<totM) return true;
  }
  if(LQ&&PQ)                                     // ===> "Can try to split lambda+proton"
  {
    G4QContent resQC=valQC-protQC-lambQC;        // QC of ResidNucleus for Proton+Lambda
    G4int    resPDG=resQC.GetSPDGCode();         // PDG of ResidNucleus for Proton+Lambda
    G4double resMas=G4QPDGCode(resPDG).GetMass();// GS Mass of the Residual Nucleus
    G4double sM=resMas+mProt+mLamb;
#ifdef debug
    G4cout<<"G4QNucleus::Split2Baryons: (proton+lambda), sM="<<sM<<", d="<<totM-sM<<G4endl;
#endif
    if(sM<totM) return true;
  }
  if(LQ>1)                                       // ===> "Can try to split 2 lambdas" case
  {
    G4QContent resQC=valQC-lambQC-lambQC;        // QC of ResidNucleus for the Two Lambdas
    G4int    resPDG=resQC.GetSPDGCode();         // PDG of ResidNucleus for the Two Lambdas
    G4double resMas=G4QPDGCode(resPDG).GetMass();// GS Mass of the Residual Nucleus
    G4double sM=resMas+mLamb+mLamb;
#ifdef debug
    G4cout<<"G4QNucleus::Split2Baryons: (two lambdas), sM="<<sM<<", d="<<totM-sM<<G4endl;
#endif
    if(sM<totM) return true;
  }
  return false;
}
  
// Evaporate one Baryon (n,p,Lambda) (h1) from the Nucleus & get Residual Nucleus (h2)
G4bool G4QNucleus::EvaporateBaryon(G4QHadron* h1, G4QHadron* h2)
{
  //static const G4double   uWell=2.7;              // EffectiveDepth of potential well B
  //static const G4double   uWell=7.;               // EffectiveDepth of potential well B
  static const G4double   uWell=1.7;              // EffectiveDepth of potential well B
  //static const G4double   uWell=0.0;              // EffectiveDepth of potential well B
  //////////static const G4double   gunA=80.;       // Switch A-parameter for BaryonGun
  //static const G4double   gunB=exp(1)/gunA;
  ///////////////////static const G4double   gunB=exp(2)/4/gunA/gunA;
  //////////////static const G4double   gunP2=200000.; // Switch P2-parameter for BaryonGun
  //////////////static const G4double   maSht=1.2;  // shift for maximal x approximation
  ///////////static const G4double   coSht=.19;     // multiple for maximal x approximation
  //////////////static const G4double   third=1./3.;// power for maximal x approximation
  static const G4int      gPDG =   22;            // PDGCode of gamma
  static const G4QPDGCode gQPDG(gPDG);            // QPDGCode of gamma
  static const G4int      nPDG = 2112;            // PDGCode of neutron
  static const G4QPDGCode nQPDG(nPDG);            // QPDGCode of neutron
  static const G4QPDGCode anQPDG(-nPDG);          // QPDGCode of anti-neutron
  static const G4int      pPDG = 2212;            // PDGCode of proton
  static const G4QPDGCode pQPDG(pPDG);            // QPDGCode of proton
  static const G4QPDGCode apQPDG(-pPDG);          // QPDGCode of anti-proton
  static const G4int      lPDG = 3122;            // PDGCode of Lambda
  static const G4QPDGCode lQPDG(lPDG);            // QPDGCode of Lambda
  static const G4QPDGCode aDppQPDG(-2224);        // QPDGCode of anti-Delta++
  static const G4QPDGCode aDmQPDG(-1114);         // QPDGCode of anti-Delta-
  static const G4QPDGCode alQPDG(-lPDG);          // QPDGCode of anti-Lambda
  static const G4int      dPDG = 90001001;        // PDGCode of deutron
  static const G4int      aPDG = 90002002;        // PDGCode of ALPHA
  static const G4QPDGCode aQPDG(aPDG);            // QPDGCode of ALPHA
  static const G4QPDGCode NPQPDG(dPDG);           // QPDGCode of deutron
  static const G4QPDGCode NNQPDG(90000002);       // QPDGCode of n+n
  static const G4QPDGCode PPQPDG(90002000);       // QPDGCode of p+p
  static const G4QPDGCode NLQPDG(91000001);       // QPDGCode of n+L
  static const G4QPDGCode PLQPDG(91001000);       // QPDGCode of p+L
  static const G4QPDGCode LLQPDG(92000000);       // QPDGCode of L+L
  static const G4QPDGCode NAQPDG(90002003);       // QPDGCode of N+ALPHA
  static const G4QPDGCode PAQPDG(90003002);       // QPDGCode of L+ALPHA
  static const G4QPDGCode LAQPDG(91002002);       // QPDGCode of L+ALPHA
  static const G4QPDGCode AAQPDG(90004004);       // QPDGCode of ALPHA+ALPHA
  static const G4QPDGCode PIPQPDG(211);           // QPDGCode of PI+
  static const G4QPDGCode PIMQPDG(-211);          // QPDGCode of PI+
  static const G4double   mNeut= G4QPDGCode(nPDG).GetMass(); // Mass of neutron
  static const G4double   mProt= G4QPDGCode(pPDG).GetMass(); // Mass of proton
  static const G4double   mLamb= G4QPDGCode(lPDG).GetMass(); // Mass of Lambda
  static const G4double   mDeut= G4QPDGCode(nPDG).GetNuclMass(1,1,0);// Mass of deutr
  static const G4double   mAlph= G4QPDGCode(nPDG).GetNuclMass(2,2,0);// Mass of alpha
  static const G4double   mPi  = G4QPDGCode(211).GetMass();  // Mass of charged pion
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
  G4bool barf=true;                               // Take into account CB in limits
  G4double uW=uWell;
  G4int    a = GetA();
  G4double  evalph=0.1;                            // Probability for alpha to evaporate
  //if(a>4.5) evalph=2.7/sqrt(a-4.);                // Probability for alpha to evaporate
  //G4double evalph=clustProb*clustProb*clustProb;
#ifdef debug
  G4cout<<"G4QNucleus::EvaporBaryon: *Called*, a="<<a<<GetThis()<<",alph="<<evalph<<G4endl;
#endif
  G4double a1= a-1;
  //////////G4double z = Z;
  //////////G4double zn= Z+N;
  G4double PBarr= CoulombBarrier(1,1);            // CoulombBarrier for proton
  G4double PPBarr= CoulombBarrier(1,1,1,1);       // CoulombBarrier for proton (after prot)
  G4double PABarr= CoulombBarrier(1,1,2,4);       // CoulombBarrier for proton (after alph)
  G4double APBarr= CoulombBarrier(2,4,1,1);       // CoulombBarrier for alpha (after prot)
  G4double ABarr= CoulombBarrier(2,4);            // CoulombBarrier for alpha
  G4double AABarr= CoulombBarrier(2,4,2,4);       // CoulombBarrier for alpha (after alpha)
  //G4double PPPBarr= CoulombBarrier(1,1,2,2);    // CoulombBarrier for proton (after 2 pr)
  //G4double AAABarr= CoulombBarrier(2,4,4,8);    // CoulombBarrier for alpha (after 2alph)
  //////G4double APABarr= CoulombBarrier(2,4,3,5);// CoulombBarrier for alpha (after p+al)
  //G4double PPABarr= CoulombBarrier(1,1,3,5);    // CoulombBarrier for proton (after p+al)
  G4double SPPBarr=PBarr+PPBarr;                  // SummedCoulombBarrier for p+p pair
  G4double SAABarr=ABarr+AABarr;                  // SummedCoulombBarrier for 2 alpha pair
  //G4double SPPPBarr=SPPBarr+PPPBarr;            // SummedCoulombBarrier for 3 protons
  //G4double SAAABarr=SAABarr+AAABarr;            // SummedCoulombBarrier for 3 alphas
  G4double SAPBarr=PABarr+ABarr;                  // SummedCoulombBarrier for alpha+p pair
  G4double DAPBarr=APBarr+PBarr;                  // Other SummedCoulombBarrier for alph+2p
  if(DAPBarr>SAPBarr)SAPBarr=DAPBarr;             // Get max to make possible BothSequences
  ///////G4double SAPABarr=APABarr+SAPBarr;       // Summed Coulomb Barrier for alph+p+alph
  //G4double SPPABarr=PPABarr+SAPBarr;            // Summed Coulomb Barrier for p+p+alpha
  G4LorentzVector h1mom;
  G4LorentzVector h2mom;
  G4LorentzVector h3mom;
  G4double totMass= GetMass();                    // Total mass of the Nucleus
#ifdef debug
  G4cout<<"G4QN::EB:pB="<<PBarr<<",aB="<<ABarr<<",ppB="<<PPBarr<<",paB="<<PABarr<<G4endl;
#endif
  if(a==-2)
  {
    if(Z==1 || N==1)
    {
      G4int  nucPDG  = -2112;
      G4int  piPDG   =  211;
      G4double nucM  = mNeut;
      G4QPDGCode del = aDmQPDG;
      G4QPDGCode nuc = anQPDG;
      if(N>0)
      {
        nucPDG = -2212;
        piPDG  = -211;
        nucM   = mProt;
        del    = aDppQPDG;
        nuc    = apQPDG;
      }
      if(totMass > mPi+nucM+nucM)
      {
        G4LorentzVector n14M(0.,0.,0.,nucM);
        G4LorentzVector n24M(0.,0.,0.,nucM);
        G4LorentzVector pi4M(0.,0.,0.,mPi);
        if(!DecayIn3(n14M, n24M, pi4M))
        {
          G4cerr<<"***G4QNucl::EvapBary: (anti) tM="<<totMass<<"-> 2N="<<nucPDG<<"(M="
                <<nucM<<") + pi="<<piPDG<<"(M="<<mPi<<")"<<G4endl;
          //throw G4QException("G4QNucl::EvapBary:ISO-dibaryon DecayIn3 did not succeed");
          return false;
        }
        n14M+=pi4M;
        h1->SetQPDG(del);
        h2->SetQPDG(nuc);
        h1->Set4Momentum(n14M);
        h2->Set4Momentum(n24M);
        return true;
      }
      else
      {
        G4cerr<<"***G4QNucleus::EvaporateBaryon: M="<<totMass
              <<", M="<<totMass<<" < M_2N+Pi, d="<<totMass-2*nucM-mPi<<G4endl;
        //throw G4QException("***G4QNucl::EvaporateBaryon: ISO-dibaryon under Mass Shell");
        return false;
      }      
    }
    else if(Z==2 || N==2)
    {
      G4int  nucPDG  = -2112;
      G4int  piPDG   =  211;
      G4double nucM  = mNeut;
      G4QPDGCode del = aDmQPDG;
      if(N==2)
      {
        nucPDG = -2212;
        piPDG  = -211;
        nucM   = mProt;
        del    = aDppQPDG;
      }
      if(totMass > mPi+mPi+nucM+nucM)
      {
        G4LorentzVector n14M(0.,0.,0.,nucM);
        G4LorentzVector n24M(0.,0.,0.,nucM);
        G4LorentzVector pi4M(0.,0.,0.,mPi+mPi);
        if(!DecayIn3(n14M, n24M, pi4M))
        {
          G4cerr<<"***G4QNucl::EvapBary: (anti) tM="<<totMass<<"-> 2N="<<nucPDG<<"(M="
                <<nucM<<") + 2pi="<<piPDG<<"(M="<<mPi<<")"<<G4endl;
          //throw G4QException("G4QNucl::EvapBary:ISO-dibaryon DecayIn3 did not succeed");
          return false;
        }
        G4LorentzVector hpi4M=pi4M/2.;
        n14M+=hpi4M;
        n24M+=hpi4M;
        h1->SetQPDG(del);
        h2->SetQPDG(del);
        h1->Set4Momentum(n14M);
        h2->Set4Momentum(n24M);
        return true;
      }
      else
      {
        G4cerr<<"***G4QNucleus::EvaporateBaryon: M="<<totMass
              <<", M="<<totMass<<" < M_2N+Pi, d="<<totMass-2*nucM-mPi<<G4endl;
        //throw G4QException("***G4QNucl::EvaporateBaryon: ISO-dibaryon under Mass Shell");
        return false;
      }      
    }
    else if(Z==-2)
    {
      h1mom=G4LorentzVector(0.,0.,0.,mProt);
      h2mom=h1mom;
      h1->SetQPDG(apQPDG);
      h2->SetQPDG(apQPDG);
      if(!DecayIn2(h1mom,h2mom)) return false;
    }
    else if(N==-2)
    {
      h1mom=G4LorentzVector(0.,0.,0.,mNeut);
      h2mom=h1mom;
      h1->SetQPDG(anQPDG);
      h2->SetQPDG(anQPDG);
      if(!DecayIn2(h1mom,h2mom)) return false;
    }
    else if(N==-1 && Z==-1)
    {
      h1mom=G4LorentzVector(0.,0.,0.,mProt);
      h2mom=G4LorentzVector(0.,0.,0.,mNeut);
      h1->SetQPDG(apQPDG);
      h2->SetQPDG(anQPDG);
      if(!DecayIn2(h1mom,h2mom)) return false;
    }
    else if(Z==-1 && S==-1)
    {
      h1mom=G4LorentzVector(0.,0.,0.,mProt);
      h2mom=G4LorentzVector(0.,0.,0.,mLamb);
      h1->SetQPDG(apQPDG);
      h2->SetQPDG(alQPDG);
      if(!DecayIn2(h1mom,h2mom)) return false;
    }
    else
    {
      h1mom=G4LorentzVector(0.,0.,0.,mNeut);
      h2mom=G4LorentzVector(0.,0.,0.,mLamb);
      h1->SetQPDG(anQPDG);
      h2->SetQPDG(alQPDG);
      if(!DecayIn2(h1mom,h2mom)) return false;
    }
    h1->Set4Momentum(h1mom);
    h2->Set4Momentum(h2mom);
    return true;
  }
  else if(a==2)
  {
    if(Z<0||N<0)
    {
      G4int  nucPDG = 2112;
      G4double nucM = mNeut;
      G4int   piPDG = -211;
      G4QPDGCode db = NNQPDG;
      G4QPDGCode pi_value = PIMQPDG;
      if(N<0)
      {
        nucPDG = 2212;
        nucM   = mProt;
        piPDG  = 211;
        db     = PPQPDG;
        pi_value = PIPQPDG;
      }
      if(totMass>mPi+nucM+nucM)
      {
        G4LorentzVector n14M(0.,0.,0.,nucM);
        G4LorentzVector n24M(0.,0.,0.,nucM);
        G4LorentzVector pi4M(0.,0.,0.,mPi);
        if(!DecayIn3(n14M,n24M,pi4M))
        {
          G4cerr<<"***G4QNucl::EvapBary: tM="<<totMass<<"-> 2N="<<nucPDG<<"(M="
          <<nucM<<") + pi="<<piPDG<<"(M="<<mPi<<")"<<G4endl;
          //throw G4QException("G4QNucl::EvapBary:ISO-dibaryon DecayIn3 did not succeed");
          return false;
        }
        n14M+=n24M;
        h1->SetQPDG(db);
        h2->SetQPDG(pi_value);
        h1->Set4Momentum(n14M);
        h2->Set4Momentum(pi4M);
        return true;
      }
      else
      {
        G4cerr<<"***G4QNucleus::EvaporateBaryon: M="<<totMass
              <<", M="<<totMass<<" < M_2N+Pi, d="<<totMass-2*nucM-mPi<<G4endl;
        //throw G4QException("***G4QNucl::EvaporateBaryon: ISO-dibaryon under Mass Shell");
        return false;
      }      
    }
    else if(Z==2)
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
#ifdef debug
        G4cout<<"G4QNucl::EvaporateBaryon: Photon ### d+g ###, dM="<<totMass-mNP<<G4endl;
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
    G4bool nFlag    = false;               // Flag of possibility to radiate neutron
    G4bool pFlag    = false;               // Flag of possibility to radiate proton
    G4bool lFlag    = false;               // Flag of possibility to radiate lambda
    G4bool aFlag    = false;               // Flag of possibility to radiate alpha
    G4bool nnFlag   = false;               // Flag of possibility to radiate 2 neutrons
    G4bool npFlag   = false;               // Flag of possibility to radiate neutron+proton
    G4bool nlFlag   = false;               // Flag of possibility to radiate neutron+lambda
    G4bool ppFlag   = false;               // Flag of possibility to radiate 2 protons
    G4bool plFlag   = false;               // Flag of possibility to radiate proton+lambda
    G4bool llFlag   = false;               // Flag of possibility to radiate 2 lambdas
    G4bool paFlag   = false;               // Flag of possibility to radiate proton+alpha
    G4bool naFlag   = false;               // Flag of possibility to radiate neutron+alpha
    G4bool laFlag   = false;               // Flag of possibility to radiate lambda+alpha
    G4bool aaFlag   = false;               // Flag of possibility to radiate alpha+alpha
    //G4bool nnnF     = false;             // Evaporation brunch is closed
    //G4bool nnpF     = false;
    //G4bool nppF     = false;
    //G4bool pppF     = false;
    //G4bool nnlF     = false;
    //G4bool nplF     = false;
    //G4bool pplF     = false;
    //G4bool nllF     = false;
    //G4bool pllF     = false;
    //G4bool lllF     = false;
    //G4bool nnaF     = false;
    //G4bool npaF     = false;
    //G4bool ppaF     = false;
    //G4bool nlaF     = false;
    //G4bool plaF     = false;
    //G4bool llaF     = false;
    //G4bool paaF     = false;
    //G4bool naaF     = false;
    //G4bool laaF     = false;
    //G4bool aaaF     = false;
    G4double GSMass = GetGSMass(); // Ground State mass of the Nucleus
    G4double GSResNN= GSMass;      // Prototype of Residual Nuclear Mass for n+n
    G4double GSResNP= GSMass;      // Prototype of Residual Nuclear Mass for n+p
    G4double GSResNL= GSMass;      // Prototype of Residual Nuclear Mass for n+l
    G4double GSResPP= GSMass;      // Prototype of Residual Nuclear Mass for p+p
    G4double GSResPL= GSMass;      // Prototype of Residual Nuclear Mass for p+l
    G4double GSResLL= GSMass;      // Prototype of Residual Nuclear Mass for l+l
    G4double GSResNA= GSMass;      // Prototype of Residual Nuclear Mass for n+alpha
    G4double GSResPA= GSMass;      // Prototype of Residual Nuclear Mass for p+alpha
    G4double GSResLA= GSMass;      // Prototype of Residual Nuclear Mass for l+alpha
    G4double GSResAA= GSMass;      // Prototype of Residual Nuclear Mass for alpha+alpha
    G4double GSResNa= GSMass;      // Prototype of Residual Nuclear Mass for alpha
    /*
    // DHW 16 June 2011 : these variables set but not used.  Comment out to fix
    //                    compiler warnings 
    G4double GSReNNN= GSMass;      // Prototype of Residual Nuclear Mass for n+n+n
    G4double GSReNNP= GSMass;      // Prototype of Residual Nuclear Mass for n+n+p
    G4double GSReNPP= GSMass;      // Prototype of Residual Nuclear Mass for n+p+p
    G4double GSRePPP= GSMass;      // Prototype of Residual Nuclear Mass for p+p+p
    G4double GSReNNL= GSMass;      // Prototype of Residual Nuclear Mass for n+n+l
    G4double GSReNPL= GSMass;      // Prototype of Residual Nuclear Mass for n+p+l
    G4double GSRePPL= GSMass;      // Prototype of Residual Nuclear Mass for p+p+l
    G4double GSReNLL= GSMass;      // Prototype of Residual Nuclear Mass for n+l+l
    G4double GSRePLL= GSMass;      // Prototype of Residual Nuclear Mass for p+l+l
    G4double GSReLLL= GSMass;      // Prototype of Residual Nuclear Mass for l+l+l
    G4double GSReNNA= GSMass;      // Prototype of Residual Nuclear Mass for n+n+a
    G4double GSReNPA= GSMass;      // Prototype of Residual Nuclear Mass for n+p+a
    G4double GSRePPA= GSMass;      // Prototype of Residual Nuclear Mass for p+p+a
    G4double GSReNLA= GSMass;      // Prototype of Residual Nuclear Mass for n+l+a
    G4double GSRePLA= GSMass;      // Prototype of Residual Nuclear Mass for p+l+a
    G4double GSReLLA= GSMass;      // Prototype of Residual Nuclear Mass for l+l+a
    G4double GSRePAA= GSMass;      // Prototype of Residual Nuclear Mass for p+a+a
    G4double GSReNAA= GSMass;      // Prototype of Residual Nuclear Mass for n+a+a
    G4double GSReLAA= GSMass;      // Prototype of Residual Nuclear Mass for l+a+a
    G4double GSReAAA= GSMass;      // Prototype of Residual Nuclear Mass for a+a+a
    */
    G4QPDGCode PQPDG(22);          // Prototype of QPDG for ResidualNucleus to proton
    G4QPDGCode NQPDG(22);          // Prototype of QPDG for ResidualNucleus to neutron
    G4QPDGCode LQPDG(22);          // Prototype of QPDG for ResidualNucleus to lambda
    G4QPDGCode AQPDG(22);          // Prototype of QPDG for ResidualNucleus to alpha
    G4QPDGCode nnQPDG(22);         // Prototype of QPDG for ResidualNucleus to nn-dibar.
    G4QPDGCode npQPDG(22);         // Prototype of QPDG for ResidualNucleus to np-dibar.
    G4QPDGCode nlQPDG(22);         // Prototype of QPDG for ResidualNucleus to nl-dibar.
    G4QPDGCode ppQPDG(22);         // Prototype of QPDG for ResidualNucleus to pp-dibar.
    G4QPDGCode plQPDG(22);         // Prototype of QPDG for ResidualNucleus to pl-dibar.
    G4QPDGCode llQPDG(22);         // Prototype of QPDG for ResidualNucleus to ll-dibar.
    G4QPDGCode naQPDG(22);         // Prototype of QPDG for ResidualNucleus to n+alpha
    G4QPDGCode paQPDG(22);         // Prototype of QPDG for ResidualNucleus to p+alpha
    G4QPDGCode laQPDG(22);         // Prototype of QPDG for ResidualNucleus to l+alpha
    G4QPDGCode aaQPDG(22);         // Prototype of QPDG for ResidualNucleus to alph+alph
    G4QPDGCode dbQPDG(22);         // Prototype of chosen dibaryon QPDG
    G4QPDGCode fQPDG(22);          // Prototype of QPDG of the Second Baryon
    G4double rMass  = 0.;          // Prototype of mass of Residual Nucleus
    G4double eMass  = 0.;          // Prototype of mass of Evaporated Baryon
    G4double fMass  = 0.;          // Prototype of mass of the Second Baryon
#ifdef debug
    G4cout<<"G4QNuc::EvaB:a>2, totM="<<totMass<<" > GSMass="<<GSMass<<",d="<<totMass-GSMass
          <<G4endl;
#endif
    G4double tM2    = totMass*totMass;
    G4double qtM2   = 4*tM2;
    G4double GSResNp= GSMass;     // Prototype of Residual Nuclear Mass for proton
    G4double pExcess= 0.;         // Prototype of excess energy for proton
    G4double aExcess= 0.;         // Prototype of excess energy for alpha
    G4double pp2m   = 0.;         // Prototype of max square momentum for proton
    G4double ap2m   = 0.;         // Prototype of max square momentum for proton
    G4double pBnd   = 0.;         // Binding energy for proton
    G4double aBnd   = 0.;         // Binding energy for proton
    G4bool three=false;           // Prototype of the Flag of b+b+ResNuc decay
    if(Z>0)
    {
      PQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-1)+N);
      GSResNp=PQPDG.GetMass();
      G4double mpls=GSResNp+mProt;
      G4double mmin=GSResNp-mProt;
      pp2m=(tM2-mpls*mpls)*(tM2-mmin*mmin)/qtM2;
      if(pp2m>=0.000001)
      {
        pFlag=true;
        pBnd=mProt-GSMass+GSResNp;             // Binding energy for proton
        G4double eMax=sqrt(mP2+pp2m);
#ifdef debug
        G4cout<<"G4QNuc::EvapBaryon:pm="<<eMax+sqrt(pp2m+GSResNp*GSResNp)<<" = M="<<totMass
              <<", sm="<<GSResNp+mProt+PBarr<<",pp2="<<pp2m<<",pB="<<pBnd<<G4endl;
#endif
        pExcess=eMax-mProt+pBnd;               // Max Kin Energy from bottom
      }
      else pExcess=pBnd;
      if(Z>1)
      {
        ppQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-2)+N);
        GSResPP=ppQPDG.GetMass();
#ifdef debug
        G4double sm=GSResPP+mProt+mProt+SPPBarr;
        G4cout<<"G4QNucl::EvapBaryon: ppM="<<GSResPP<<",T="<<sm-GSMass<<",E="<<totMass-sm
              <<",C="<<PBarr<<G4endl;
#endif
        if(GSResPP+mProt+mProt+SPPBarr<totMass) ppFlag=true;
        if(Z>2&&a>3)
        {
          /*
          // DHW 16 June 2011: variable set but not used.  See note at line 1317.
          GSRePPP=G4QPDGCode().GetNuclMass(Z-3,N,S);
          */
          //if(GSRePPP+mProt+mProt+mProt+SPPPBarr<totMass) pppF=true;
          if(N>1&&a>5)
          {
            paQPDG =G4QPDGCode(90000000+1000*(1000*S+Z-3)+N-2);
            GSResPA=paQPDG.GetMass();
#ifdef debug
            G4double s_value=GSResPA+mAlph+mProt+SAPBarr;
            G4cout<<"G4QN::EB:paM="<<GSResPA<<",T="<<s_value-GSMass<<",E="<<totMass-s_value<<G4endl;
#endif
            if(GSResPA+mProt+SAPBarr+mAlph<totMass) paFlag=true;
          }
        }
        if(N>0&&a>3)
        {
          /*
          // DHW  16 June 2011: variable set but not used.  See note at line 1317.
          GSReNPP=G4QPDGCode().GetNuclMass(Z-2,N-1,S);
          */
          //if(GSReNPP+mProt+mProt+SPPBarr+mNeut<totMass) nppF=true;
        }
        if(S>0&&a>3)
        {
          /*
          // DHW  16 June 2011: variable set but not used.  See note at line 1317.
          GSRePPL=G4QPDGCode().GetNuclMass(Z-2,N,S-1);
          */
          //if(GSRePPL+mProt+mProt+SPPBarr+mLamb<totMass) pplF=true;
        }
        if(N>1&&a>4)
        {
          if(a>6)
          {
            if(S>1)
            {
              /*
              // DHW  16 June 2011: variable set but not used.  See note at line 1317.
              GSReLLA=G4QPDGCode().GetNuclMass(Z-2,N-2,S-2);
              */
              //if(GSReLLA+mAlph+ABarr+mLamb+mLamb<totMass) llaF=true;
            }
            if(N>2&&S>0)
            {
              /*
              // DHW 16 June 2011:  variable set but not used.  See note at line 1317.
              GSReNLA=G4QPDGCode().GetNuclMass(Z-2,N-3,S-1);
              */
              //if(GSReNLA+mAlph+ABarr+mNeut+mLamb<totMass) nlaF=true;
            }
            if(Z>2&&S>0)
            {
              /*
              // DHW 16 June 2011:  variable set but not used.  See note at line 1317.
              GSRePLA=G4QPDGCode().GetNuclMass(Z-3,N-2,S-1);
              */
              //if(GSRePLA+mAlph+SAPBarr+mProt+mLamb<totMass) plaF=true;
            }
            if(N>3)
            {
              /*
              // DHW 16 June 2011:  variable set but not used.  See note at line 1317.
              GSReNNA=G4QPDGCode().GetNuclMass(Z-2,N-4,S);
              */
              //if(GSReNNA+mAlph+ABarr+mNeut+mNeut<totMass) nnaF=true;
            }
            if(Z>2&&N>2)
            {
              /*
              // DHW 16 June 2011:  variable set but not used.  See note at line 1317.
              GSReNPA=G4QPDGCode().GetNuclMass(Z-3,N-3,S);
              */
              //if(GSReNPA+mAlph+SAPBarr+mProt+mNeut<totMass) npaF=true;
            }
            if(N>3)
            {
              /*
              // DHW 16 June 2011:  variable set but not used.  See note at line 1317.
              GSRePPA=G4QPDGCode().GetNuclMass(Z-4,N-2,S);
              */
              //if(GSRePPA+mAlph+SPPABarr+mProt+mProt<totMass) ppaF=true;
            }
            if(a>9)
            {
              if(Z>3&&N>3&&S>0)
              {
                /*
                // DHW 16 June 2011:  variable set but not used.  See note at line 1317.
                GSReLAA=G4QPDGCode().GetNuclMass(Z-4,N-4,S-1);
                */
                //if(GSReLAA+mLamb+mAlph+mAlph+SAABarr<totMass) laaF=true;
              }
              if(Z>3&&N>4)
              {
                /*
                // DHW 16 June 2011:  variable set but not used.  See note at line 1317.
                GSReNAA=G4QPDGCode().GetNuclMass(Z-4,N-5,S);
                */
                //if(GSReNAA+mNeut+mAlph+mAlph+SAABarr<totMass) naaF=true;
              }
              if(Z>4&&N>3)
              {
                /*
                // DHW 16 June 2011:  variable set but not used.  See note at line 1317.
                GSRePAA=G4QPDGCode().GetNuclMass(Z-5,N-4,S);
                */
                //if(GSRePAA+mProt+mAlph+mAlph+SAABarr<totMass) paaF=true;
              }
              if(a>12&&N>5&&Z>5)
              {
                /*
                // DHW 16 June 2011:  variable set but not used.  See note at line 1317.
                GSReAAA=G4QPDGCode().GetNuclMass(Z-6,N-6,S);
                */
                //if(GSReAAA+mAlph+mAlph+mAlph+SAAABarr<totMass) aaaF=true;
              }
            }
          }
          if(N>3&&Z>3&&a>8)
          {
            aaQPDG =G4QPDGCode(90000000+1000*(1000*S+Z-4)+N-4);
            GSResAA=aaQPDG.GetMass();
#ifdef debug
            G4double s_value=GSResAA+mAlph+mAlph+SAABarr;
            G4cout<<"G4QNucl::EvapBaryon: a="<<GSResNP<<",T="<<s_value-GSMass<<",E="<<totMass-s_value
                  <<",A="<<SAABarr<<G4endl;
#endif
            if(GSResAA+mAlph+mAlph+SAABarr<totMass) aaFlag=true;
          }
          if(N>2&&a>5)
          {
            naQPDG =G4QPDGCode(90000000+1000*(1000*S+Z-2)+N-3);
            GSResNA=naQPDG.GetMass();
#ifdef debug
            G4double s_value=GSResNA+mAlph+mNeut;
            G4cout<<"G4QNucl::EvapBary: M="<<GSResNA<<",T="<<s_value-GSMass<<",E="<<totMass-s_value
                  <<",C="<<ABarr<<G4endl;
#endif
            if(GSResNA+mNeut+mAlph+ABarr<totMass) naFlag=true;
          }
          if(S>0&&a>5)
          {
            laQPDG =G4QPDGCode(90000000+1000*(1000*S+Z-1002)+N-2);
            GSResLA=laQPDG.GetMass();
            if(GSResLA+mLamb+mAlph+ABarr<totMass) laFlag=true;
          }
          AQPDG =G4QPDGCode(90000000+1000*(1000*S+Z-2)+N-2);
          GSResNa=AQPDG.GetMass();
          mpls=GSResNa+mAlph;
          mmin=GSResNa-mAlph;
          ap2m=(tM2-mpls*mpls)*(tM2-mmin*mmin)/qtM2;
          if(ap2m>=0.000001)
          {
            aFlag=true;
            aBnd=mAlph-GSMass+GSResNa;           // Binding energy for ALPHA
            G4double eMax=sqrt(mA2+ap2m);
#ifdef debug
            G4cout<<"G4QNuc::EvapBar:m="<<eMax+sqrt(ap2m+GSResNa*GSResNa)<<" = M="<<totMass
                  <<", sm="<<GSResNp+mProt+PBarr<<",pp2="<<pp2m<<",pB="<<pBnd<<G4endl;
#endif
            aExcess=eMax-mAlph+aBnd;             // Max Kin Energy from bottom
          }
          else aExcess=pBnd;
        }
      }
      if(N>0)
      {
        if(Z>0)
        {
          npQPDG=G4QPDGCode(90000000+1000*(1000*S+Z-1)+N-1);
          GSResNP=npQPDG.GetMass();
#ifdef debug
          G4double s_value=GSResNP+mNeut+mProt;
          G4cout<<"G4QNucl::EvapBaryon: npM="<<GSResNP<<",T="<<s_value-GSMass<<",E="<<totMass-s_value
                <<",C="<<PBarr<<G4endl;
#endif
          if(GSResNP+mNeut+mProt+PBarr<totMass) npFlag=true;
        }
        if(N>1)
        {
          /* 
          // DHW 16 June 2011: variable set but not used.  See note at line 1317;
          GSReNNP=G4QPDGCode().GetNuclMass(Z-1,N-2,S);
          */
          //if(GSReNNP+mProt+PBarr+mNeut+mNeut<totMass) nnpF=true;
        }
        if(S>0)
        {
          /*
          // DHW 16 June 2011: variable set but not used.  See note at line 1317;
          GSReNPL=G4QPDGCode().GetNuclMass(Z-1,N-1,S-1);
          */
          //if(GSReNPL+mProt+PBarr+mNeut+mLamb<totMass) nplF=true;
        }
      }
      if(S>0)
      {
        if(Z>0)
        {
          plQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z-1)+N);
          GSResPL=plQPDG.GetMass();
          if(GSResPL+mProt+PBarr+mLamb<totMass) plFlag=true;
        }
        if(S>1)
        {
          /*
          // DHW 16 June 2011: variable set but not used.  See note at line 1317;
          GSRePLL=G4QPDGCode().GetNuclMass(Z-1,N,S-2);
          */
          //if(GSRePLL+mProt+PBarr+mLamb+mLamb<totMass) pllF=true;
        }
      }
    }
    G4double GSResNn= GSMass;         // Prototype of Residual Nuclear Mass for neutron
    G4double nExcess= 0.;             // Prototype of excess energy for neutron
    G4double np2m   = 0.;             // Prototype of max square momentum for neutron
    G4double nBnd   = 0.;             // Binding energy for neutron
    if(N>0)
    {
      NQPDG=G4QPDGCode(90000000+1000*(1000*S+Z)+N-1);
      GSResNn=NQPDG.GetMass();
#ifdef debug
      G4cout<<"G4QNucleus::EvapBaryon: M(A-N)="<<GSResNn<<",Z="<<Z
            <<",N="<<N<<",S="<<S<<G4endl;
#endif
      G4double mpls=GSResNn+mNeut;
      G4double mmin=GSResNn-mNeut;
      np2m=(tM2-mpls*mpls)*(tM2-mmin*mmin)/qtM2;
      if(np2m>=0.000001)
      {
        nFlag=true;
        nBnd=mNeut-GSMass+GSResNn;    // Binding energy for neutron
        G4double eMax=sqrt(mN2+np2m);
#ifdef debug
        G4cout<<"G4QNuc::EvapBaryon:nm="<<eMax+sqrt(np2m+GSResNn*GSResNn)<<" = M="<<totMass
              <<", sm="<<GSResNn+mNeut<<",np2="<<np2m<<",nB="<<nBnd<<G4endl;
#endif
        nExcess=eMax-mNeut+nBnd;
      }
      else nExcess=nBnd;
      if(N>1)
      {
        nnQPDG=G4QPDGCode(90000000+1000*(1000*S+Z)+N-2);
        GSResNN=nnQPDG.GetMass();
        if(GSResNN+mNeut+mNeut<totMass) nnFlag=true;
        if(N>2)
        {
          /*
          // DHW 16 June 2011: variable set but not used.  See note at line 1317.
          GSReNNN=G4QPDGCode().GetNuclMass(Z,N-3,S);
          */
          //if(GSReNNN+mNeut*3<totMass) nnnF=true;
        }
        if(S>0)
        {
          /*
          // DHW 16 June 2011: variable set but not used.  See note at line 1317.
          GSReNNL=G4QPDGCode().GetNuclMass(Z,N-2,S-1);
          */
          //if(GSReNNL+mNeut+mNeut+mLamb<totMass) nnlF=true;
        }
      }
      if(S>0)
      {
        nlQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z)+N-1);
        GSResNL=nlQPDG.GetMass();
        if(GSResNL+mNeut+mLamb<totMass) nlFlag=true;
        if(S>1)
        {
          /*
          // DHW 16 June 2011: variable set but not used.  See note at line 1317.
          GSReNLL=G4QPDGCode().GetNuclMass(Z,N-1,S-2);
          */
          //if(GSReNLL+mNeut+mLamb+mLamb<totMass) nllF=true;
        }
      }
    }
    G4double GSResNl= GSMass;         // Prototype of Residual Nuclear Mass for Lambda
    G4double lExcess= 0.;             // Prototype of excess energy for Lambda
    G4double lp2m   = 0.;             // Prototype of max square momentum for lambda
    G4double lBnd   = 0.;             // Binding energy for lambda
    if(S>0)
    {
      LQPDG=G4QPDGCode(90000000+1000*(1000*(S-1)+Z)+N);
      GSResNl=LQPDG.GetMass();
      G4double mpls=GSResNl+mLamb;
      G4double mmin=GSResNl-mLamb;
      lp2m=(tM2-mpls*mpls)*(tM2-mmin*mmin)/qtM2;
      if(lp2m>=0.000001)
      {
        lFlag=true;
        lBnd=mLamb-GSMass+GSResNl;    // Binding energy for lambda
        G4double eMax=sqrt(mL2+lp2m);
#ifdef debug
        G4cout<<"G4QNuc::EvapBaryon:lm="<<eMax+sqrt(lp2m+GSResNl*GSResNl)<<" = M="<<totMass
              <<", sm="<<GSResNl+mLamb<<",lp2="<<lp2m<<",lB="<<lBnd<<G4endl;
#endif
        lExcess=eMax-mLamb+lBnd;
      }
      else lExcess=lBnd;
      if(S>1)
      {
        llQPDG=G4QPDGCode(90000000+1000*(1000*(S-2)+Z)+N);
        GSResLL=llQPDG.GetMass();
        if(GSResLL+mLamb+mLamb<totMass) llFlag=true;
        if(S>2)
        {
          /*
          // DHW 16 June 2011: variable set but not used.  See note at line 1317.
          GSReLLL=G4QPDGCode().GetNuclMass(Z,N,S-3);
          */
          //if(GSReLLL+mLamb*3<totMass) lllF=true;
        }
      }
    }
    G4bool nSecF = nnFlag||npFlag||nlFlag||naFlag; // Pos of second radiation after neutron
    G4bool pSecF = npFlag||ppFlag||plFlag||paFlag; // Pos of second radiation after proton
    G4bool lSecF = nlFlag||plFlag||llFlag||laFlag; // Pos of second radiation after lambda
    G4bool aSecF = naFlag||paFlag||laFlag||aaFlag; // Pos of second radiation after alpha
    //G4bool nTrF=nnnF||nnpF||nppF||nnlF||nplF||nllF; //Pos of 3-d baryon radiation after n
    //G4bool pTrF=nnpF||nppF||pppF||nplF||pplF||pllF; //Pos of 3-d baryon radiation after p
    //G4bool lTrF=nnlF||nplF||pplF||nllF||pllF||lllF; //Pos of 3-d baryon radiation after l
    //G4bool aTrF=nnaF||npaF||ppaF||nlaF||plaF||llaF; //Pos of 3-d baryon radiation after a
    G4bool secB  = nSecF||pSecF||lSecF||aSecF; // Possibili to decay in TwoBaryons (Alphas)
    //G4bool thdB  = nTrF||pTrF||lTrF||aTrF||naaF||paaF||laaF||aaaF;// Pos to radiate three
#ifdef debug
    G4cout<<"G4QNucl::EvapBary:n="<<nSecF<<",p="<<pSecF<<",l="<<lSecF<<",a="<<aSecF<<",nn="
          <<nnFlag<<", np="<<npFlag<<",pp="<<ppFlag<<",pa="<<paFlag<<",na="<<naFlag<<",aa="
          <<aaFlag<<G4endl;
#endif
    G4QPDGCode bQPDG;
    G4QPDGCode rQPDG;
    if(secB)                            // Decay in two baryons is possible
    //if(thdB)                            //@@CHECK@@ Decay in three baryons is possible
    {
      if(!nSecF) nFlag=false;
      if(!pSecF) pFlag=false;
      if(!lSecF) lFlag=false;
      if(!aSecF) aFlag=false;
#ifdef debug
      G4cout<<"G4QNuc::EB:nF="<<nFlag<<",pF="<<pFlag<<",lF="<<lFlag<<",aF="<<aFlag<<G4endl;
#endif
      G4double maxE=0.;                          // Prototype for maximum energy
      if(nFlag&&nExcess>maxE) maxE=nExcess;
      if(pFlag&&pExcess>maxE) maxE=pExcess;
      if(lFlag&&lExcess>maxE) maxE=lExcess;
      if(lFlag&&aExcess>maxE) maxE=aExcess;
      G4double pMin=pBnd;                        // Binding energy for proton
      if(pFlag)pMin+= PBarr;                     // Add Coulomb Barrier for protons
      G4double nMin=nBnd;                        // Binding energy for neutron
      G4double lMin=lBnd;                        // Binding energy for Lambda
      G4double aMin=aBnd;                        // Binding energy for alpha
      if(aFlag)aMin+= ABarr;                     // Add Coulomb Barrier for alpha
      G4double minE=GSMass;                      // Prototype for mimimum energy
      if(nFlag&&nMin<minE) minE=nMin;
      if(pFlag&&pMin<minE) minE=pMin;
      if(lFlag&&lMin<minE) minE=lMin;
      if(evalph&&aFlag&&aMin<minE) minE=aMin;

#ifdef debug
      G4cout<<"G4QNucleus::EvapBaryon: nE="<<nExcess<<">"<<nMin<<",pE="<<pExcess<<">"<<pMin
            <<",sE="<<lExcess<<">"<<lMin<<",E="<<aExcess<<">"<<aMin<<",miE="<<minE<<"<maE="
            <<maxE<<G4endl;
#endif
      // @@ Here one can put a condition for the Baryon Gun
      G4int    cntr= 0;
      //G4int    cntm= 27;
      //G4int    cntm= 72;               // Important difference !!DOn't change
      //G4int    cntm= 80;               // Important difference !!DOn'tChange"IsoNuclei"
      //G4int    cntm= 90;               // Important difference !!DOn'tChange "Lept/Hyper"
      G4int    cntm= 53;       // @@ NonClusters in CHIPSWorld (cntm=nQHM in G4QPDGCode.hh)
      if( ( (pFlag && pExcess > pMin) || 
            (nFlag && nExcess > nMin) || 
            (lFlag && lExcess > lMin) ||
            (aFlag && aExcess > aMin) ) && minE<maxE )
      {
        G4double mi=uWell+minE;          // Minimum Kinetic Energy for minimal nucleon
        G4double mm_value=uWell+maxE;    // Personal maximum for Kinetic Energy
        G4double ma=uWell*a+maxE;        // Total Kinetic Energy of baryons (@@alphas?)
        if(mi<0.)
        {
          uW-=mi;
          mm_value-=mi;
          mi=0.;
        }
        //G4bool good=true;
        if(ma<mm_value)
        {
          ma=mm_value;
          //good=false;
        }
#ifdef debug
        G4cout<<"G4QNuc::EvapBary:iE="<<minE<<",aE="<<maxE<<",mi="<<mi<<",mm="<<mm_value<<",ma="
              <<ma<<G4endl;
#endif
        G4double xMi=mi/ma;                       // Minimal value of x
        G4double xMm=mm_value/ma;                 // Personal maximum x
        //G4double xCa=maSht-coSht*log(a);        // Maximal value of x (approximation)
        //G4double xMa=xCa;                       // Maximal value of x
        //if(xMm<xMa) xMa=xMm;
        G4double xMa=xMm;
        if(xMa>1.)xMa=1.;
        if(xMi<0.)xMi=0.;
        if(xMi>xMa)
        {
          G4cerr<<"***G4QNucleus::EvapBaryon: M="<<mm_value/ma<<",xi="<<xMi<<",xa="<<xMa<<G4endl;
          return false;
        }
        xMi=sqrt(xMi);                          // @@ ?
        xMa=sqrt(xMa);                          // @@ ?
#ifdef debug
        G4cout<<"G4QNuc:EvapBaryon:mi="<<mi<<",ma="<<ma<<", xi="<<xMi<<",xa="<<xMa<<G4endl;
#endif
        G4double powr=1.5*a1;                   // Power for low & up limits
        G4double revP=1./powr;                  // Reversed power for randomization
#ifdef debug
        G4cout<<"G4QNucleus::EvaporateBaryon: Power="<<powr<<",RevPower="<<revP<<G4endl;
#endif
        G4double minR=pow(1.-xMa*xMa,powr);    // Look on @@ ? (up)
        G4double maxR=pow(1.-xMi*xMi,powr);
#ifdef debug
        G4cout<<"G4QNucleus::EvaporateBaryon: miR="<<minR<<", maR="<<maxR<<G4endl;
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
#ifdef debug
            G4cerr<<"**G4QNucl::EvapB:R="<<R<<",xi="<<xMi<<" < "<<x<<" < xa="<<xMa<<G4endl;
#endif
            if(x<xMi) x=xMi;
            else      x=xMa;
            x2 = x*x;
          }
          G4double rn=G4UniformRand();
          //if(rn<x/xMa||!good)
          if(rn<x/xMa)                                 // Randomization cut
          {
            tk= ma*x2-uW;                              // Kinetic energy of the fragment
            G4double psum =0.;
            G4double zCBPP=0.;                         // Probabylity for a proton
#ifdef debug
            G4cout<<"G4QNuc::EvapB:t="<<tk<<",pM="<<pMin<<",pB="<<pBnd<<",n="<<nMin<<",a="
                  <<aMin<<G4endl;
#endif
            if(pFlag&&tk>pMin)
            {
              G4double kin=tk-pBnd;
              //if(barf) kin-=PBarr; //@@ This is a mistake
#ifdef debug
              G4cout<<"G4QN::EB:Proton="<<kin<<",CB="<<PBarr<<",B="<<pBnd<<",M="<<pMin
                    <<",p="<<CoulBarPenProb(PBarr,kin,1,1)<<G4endl;
#endif
              zCBPP=Z*CoulBarPenProb(PBarr,kin,1,1)*sqrt(kin);
            }
            psum+=zCBPP;
            G4double nCBPP=0.;                       // Probability for a neutron (=> p+n)
            if(nFlag&&tk>nMin)
            {
              G4double kin=tk-nBnd;
#ifdef debug
              G4cout<<"G4QN::EB:Neutron="<<kin<<",p="<<CoulBarPenProb(0.,kin,0,1)<<G4endl;
#endif
              nCBPP=N*CoulBarPenProb(0.,kin,0,1)*sqrt(kin);
            }
            psum+=nCBPP;
            nCBPP+=zCBPP;
            G4double lCBPP=0.;                       // Probability for a lambda (=> p+n+l)
            if(lFlag&&tk>lMin)
            {
              G4double kin=tk-lBnd;
#ifdef debug
              G4cout<<"G4QN::EB:Lambda="<<kin<<",p="<<CoulBarPenProb(0,kin,0,1)<<G4endl;
#endif
              lCBPP=S*CoulBarPenProb(0.,kin,0,1)*sqrt(kin);
            }
            psum+=lCBPP;
            lCBPP+=nCBPP;
            if(evalph&&aFlag&&tk>aMin)
            {
              G4double kin=tk-aBnd;
              //if(barf) kin-=ABarr; //@@ This is a mistake
#ifdef debug
              G4cout<<"G4QN::EB:Alpha="<<kin<<",CB="<<ABarr<<",p="
                    <<CoulBarPenProb(ABarr,kin,2,4)<<G4endl;
#endif
              psum+=CoulBarPenProb(ABarr,kin,2,4)*sqrt(kin)*evalph*Z*(Z-1)*N*(N-1)
                                                 *6/a1/(a-2)/(a-3);
            }
            G4double r = psum*G4UniformRand();
#ifdef debug
            G4cout<<"G4QNuc::EvapB:"<<r<<",p="<<zCBPP<<",pn="<<nCBPP<<",pnl="<<lCBPP<<",t="
                  <<psum<<G4endl;
#endif
            cond = false;
            if     (r&&r>lCBPP)
            {
#ifdef debug
              G4cout<<"G4QNuc::EvaB:ALPHA is selected for evap, r="<<r<<">"<<lCBPP<<G4endl;
#endif
              PDG=aPDG;
            }
            else if(r&&r>nCBPP&&r<=lCBPP)
            {
#ifdef debug
              G4cout<<"G4QNuc::EvaB:LAMBDA is selected for evap,r="<<r<<"<"<<lCBPP<<G4endl;
#endif
              PDG=lPDG;
            }
            else if(r&&r>zCBPP&&r<=nCBPP)
            {
#ifdef debug
              G4cout<<"G4QNuc::EvaBar: N is selected for evapor,r="<<r<<"<"<<nCBPP<<G4endl;
#endif
              PDG=nPDG;
            }
            else if(r&&r<=zCBPP)
            {
#ifdef debug
              G4cout<<"G4QNuc::EvaBar: P is selected for evapor,r="<<r<<"<"<<zCBPP<<G4endl;
#endif
              PDG=pPDG;
            }
            else cond=true;
          }
#ifdef debug
          G4cout<<"G4QNuc::EvapBar:c="<<cond<<",x="<<x<<",cnt="<<cntr<<",R="<<R<<",ma="<<ma
                <<",rn="<<rn<<"<r="<<x/xMa<<",tk="<<tk<<",ni="<<nMin<<",pi="<<pMin<<G4endl;
#endif
          cntr++;
        }
        if(cntr<cntm)                       // => Succeeded to find the evaporation channel
        {
          G4double p2=0.;
          if     (PDG==aPDG)
          {
            tk-=aBnd-mAlph;                 // Pays for binding and convert to total energy
            p2=tk*tk-mA2;
            if(p2>ap2m)
            {
              p2=ap2m;
              tk=sqrt(p2+mA2);
            }
            eMass=mAlph;
            bQPDG=aQPDG;
            rQPDG=AQPDG;
          }
          else if(PDG==pPDG)
          {
            tk-=pBnd-mProt;                 // Pays for binding and convert to total energy
            p2=tk*tk-mP2;
            if(p2>pp2m)
            {
              p2=pp2m;
              tk=sqrt(p2+mP2);
            }
            eMass=mProt;
            bQPDG=pQPDG;
            rQPDG=PQPDG;
          }
          else if(PDG==nPDG)
          {
            tk-=nBnd-mNeut;                 // Pays for binding and convert to total energy
            p2=tk*tk-mN2;
#ifdef debug
            G4cout<<"G4QNucleus::EvaporateBaryon:np2="<<p2<<",np2m="<<np2m<<G4endl;
#endif
            if(p2>np2m)
            {
              p2=np2m;
              tk=sqrt(p2+mN2);
            }
            eMass=mNeut;
            bQPDG=nQPDG;
            rQPDG=NQPDG;
          }
          else if(PDG==lPDG)
          {
            tk-=lBnd-mLamb;                 // Pays for binding and convert to total energy
            p2=tk*tk-mL2;
            if(p2>lp2m)
            {
              p2=lp2m;
              tk=sqrt(p2+mL2);
            }
            eMass=mLamb;
            bQPDG=lQPDG;
            rQPDG=LQPDG;
          }
          else G4cerr<<"***G4QNucleus::EvaporateBaryon: PDG="<<PDG<<G4endl;
          G4double rEn=totMass-tk;
          G4double rEn2=rEn*rEn;
          if (rEn2 > p2) rMass=sqrt(rEn2-p2);      // Mass of the Residual Nucleus
          else           rMass=0.0;
          // Find out if the ResidualNucleus is below of the SecondBaryonDecayLimit
          //@@ Calculate it depending on PDG !!!!!!!
          G4bool nnCond = !nnFlag || (nnFlag && GSResNN+mNeut > rMass);
          G4bool npCond = !npFlag || (npFlag && GSResNP+mProt+PBarr > rMass);
          G4bool nlCond = !nlFlag || (nlFlag && GSResNL+mLamb > rMass);
          G4bool naCond = !naFlag || (naFlag && GSResNA+mAlph+ABarr > rMass);
          G4bool pnCond = !npFlag || (npFlag && GSResNP+mNeut > rMass);
          if(barf) pnCond = !npFlag || (npFlag && GSResNP+mNeut+PBarr > rMass);
          G4bool ppCond = !ppFlag || (ppFlag && GSResPP+mProt+PPBarr > rMass);
          if(barf) ppCond = !ppFlag || (ppFlag && GSResPP+mProt+SPPBarr > rMass);
          G4bool plCond = !plFlag || (plFlag && GSResPL+mLamb > rMass);
          if(barf) plCond = !plFlag || (plFlag && GSResPL+mLamb+PBarr > rMass);
          G4bool paCond = !paFlag || (paFlag && GSResPA+mAlph+APBarr > rMass);
          if(barf) paCond = !paFlag || (paFlag && GSResPA+mAlph+SAPBarr > rMass);
          G4bool lnCond = !nlFlag || (nlFlag && GSResNL+mNeut > rMass);
          G4bool lpCond = !plFlag || (plFlag && GSResPL+mProt+PBarr > rMass);
          G4bool llCond = !llFlag || (llFlag && GSResLL+mLamb > rMass);
          G4bool laCond = !laFlag || (laFlag && GSResLA+mAlph+ABarr > rMass);
          G4bool anCond = !naFlag || (naFlag && GSResNA+mNeut > rMass);
          if(barf) anCond = !naFlag || (naFlag && GSResNA+mNeut+ABarr > rMass);
          G4bool apCond = !paFlag || (paFlag && GSResPA+mProt+PABarr > rMass);
          if(barf) apCond = !paFlag || (paFlag && GSResPA+mProt+SAPBarr > rMass);
          G4bool alCond = !laFlag || (laFlag && GSResLA+mLamb > rMass);
          if(barf) alCond = !laFlag || (laFlag && GSResLA+mLamb+ABarr > rMass);
          G4bool aaCond = !aaFlag || (aaFlag && GSResAA+mAlph+AABarr > rMass);
          if(barf) aaCond = !aaFlag || (aaFlag && GSResAA+mAlph+SAABarr > rMass);
#ifdef debug
          G4cout<<"G4QNucl::EvaB:"<<PDG<<", E="<<tk<<", rM="<<rMass<<", ";
          if(PDG==pPDG)      G4cout<<"PN="<<GSResNP+mNeut<<"("<<pnCond<<"),PP="
                                   <<GSResPP+mProt+PPBarr<<"("<<ppCond<<"),PL="
                                   <<GSResPL+mLamb<<"("<<plCond<<"),PA="
                                   <<GSResPA+mAlph+APBarr<<"("<<paCond;
          else if(PDG==nPDG) G4cout<<"NN="<<GSResNN+mNeut<<"("<<nnCond<<"),NP="
                                   <<GSResNP+mProt+PBarr<<"("<<npCond<<"),NL="
                                   <<GSResNL+mLamb<<"("<<nlCond<<"),NA="
                                   <<GSResNA+mAlph+ABarr<<"("<<naCond;
          else if(PDG==lPDG) G4cout<<"LN="<<GSResNL+mNeut<<"("<<lnCond<<"),LP="
                                   <<GSResPL+mProt+PBarr<<"("<<lpCond<<"),LL="
                                   <<GSResLL+mLamb<<"("<<llCond<<"),LA="
                                   <<GSResLA+mAlph+ABarr<<"("<<laCond;
          else if(PDG==aPDG) G4cout<<"AN="<<GSResNA+mNeut<<"("<<anCond<<"),AP="
                                   <<GSResPA+mProt+PABarr<<"("<<apCond<<"),AL="
                                   <<GSResLA+mLamb<<"("<<alCond<<"),AA="
                                   <<GSResAA+mAlph+AABarr<<"("<<aaCond;
          G4cout<<")"<<G4endl;
#endif
          three=false;                               // Flag of b+b+ResNuc decay
          //if(3>2)three=false;                        // @@@@@@@@@@@@@@@@@@
          //else if(PDG==pPDG&&(pnCond&&ppCond&&plCond&&paCond)) // @@@@@@@@@@@@@@@@@@@
          if(PDG==pPDG&&(pnCond&&ppCond&&plCond&&paCond))//p+RN decay, p+b+RN dec is closed
          {
#ifdef debug
            G4cout<<"G4QN::EB:*p*: n="<<pnCond<<",p="<<ppCond<<",l="<<plCond<<",a="<<paCond
                  <<G4endl;
#endif
            fMass=mProt;
            fQPDG=pQPDG;
            G4double nLim=0.;
            if(N&&GSResNP!=GSMass&&fMass+PBarr+mNeut+GSResNP<totMass)
            {
              if(barf) nLim+=(N+N)*pow(totMass-mNeut-mProt-PBarr-GSResNP,2);
              else     nLim+=(N+N)*pow(totMass-mNeut-mProt-GSResNP,2);
            }
            G4double zLim=nLim;
            if(Z>1&&GSResPP!=GSMass&&fMass+mProt+SPPBarr+GSResPP<totMass)
            {
              if(barf) zLim+=(Z-1)*pow(totMass-mProt-mProt-SPPBarr-GSResPP,2);
              else     zLim+=(Z-1)*pow(totMass-mProt-mProt-GSResPP,2);
            }
            G4double sLim=zLim;
            if(S&&GSResPL!=GSMass&&fMass+PBarr+mLamb+GSResPL<totMass)
            {
              if(barf) sLim+=(S+S)*pow(totMass-mProt-mLamb-PBarr-GSResPL,2);
              else     sLim+=(S+S)*pow(totMass-mProt-mLamb-GSResPL,2);
            }
            G4double aLim=sLim;
            if(evalph&&Z>2&&N>1&&a>4&&GSResPL!=GSMass&&fMass+SAPBarr+mAlph+GSResPA<totMass)
            {
              if(barf) aLim+=pow(totMass-mProt-mAlph-SAPBarr-GSResPA,2)*evalph*
                             (Z-1)*(Z-2)*N*(N-1)*12/(a-2)/(a-3)/(a-4);
              else     aLim+=pow(totMass-mProt-mAlph-GSResPA,2)*evalph*(Z-1)*(Z-2)*N*(N-1)
                             *12/(a-2)/(a-3)/(a-4);
            }
            G4double r = aLim*G4UniformRand();
#ifdef debug
            G4cout<<"G4QNuc::EvaB:p, r="<<r<<",n="<<nLim<<",z="<<zLim<<",s="<<sLim<<",a="
                  <<aLim<<G4endl;
#endif
            three=true;                               // Flag of b+b+ResNuc decay
            if(!aLim) three=false;
            else if(r>sLim)
            {
              eMass = mAlph;
              dbQPDG= PAQPDG;
              rMass = GSResPA;
              rQPDG = paQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: P+A"<<G4endl;
#endif
            }
            else if(zLim<sLim&&r>zLim&&r<=sLim)
            {
              eMass = mLamb;
              dbQPDG= PLQPDG;
              rMass = GSResPL;
              rQPDG = plQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: P+L"<<G4endl;
#endif
            }
            else if(nLim<zLim&&r>nLim&&r<=zLim)
            {
              eMass = mProt;
              dbQPDG= PPQPDG;
              rMass = GSResPP;
              rQPDG = ppQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: P+P"<<G4endl;
#endif
            }
            else if(r<=nLim)
            {
              eMass = mNeut;
              dbQPDG= NPQPDG;
              rMass = GSResNP;
              rQPDG = npQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: P+N"<<G4endl;
#endif
            }
            else three=false;
          }
          else if(PDG==nPDG&&(nnCond&&npCond&&nlCond&&naCond)) // n+b+RN decay can't happen
          { //@@ Take into account Coulomb Barier Penetration Probability
#ifdef debug
            G4cout<<"G4QN::EB:*n*: n="<<nnCond<<",p="<<npCond<<",l="<<nlCond<<",a="<<naCond
                  <<G4endl;
#endif
            fMass=mNeut;
            fQPDG=nQPDG;
            G4double nLim=0.;
            if(N>1&&GSResNN!=GSMass&&fMass+mNeut+GSResNN<totMass)
              nLim+=(N-1)*pow(totMass-mNeut-mNeut-GSResNN,2);
            G4double zLim=nLim;
            if(Z&&GSResNP!=GSMass&&fMass+mProt+PBarr+GSResNP<totMass)
            {
              if(barf) zLim+=(Z+Z)*pow(totMass-mNeut-mProt-PBarr-GSResNP,2);
              else     zLim+=(Z+Z)*pow(totMass-mNeut-mProt-GSResNP,2);
            }
            G4double sLim=zLim;
            if(S&&GSResNL!=GSMass&&fMass+mLamb+GSResNL<totMass)
              sLim+=(S+S)*pow(totMass-mNeut-mLamb-GSResNL,2);
            G4double aLim=sLim;
            if(evalph&&Z>1&&N>2&&a>4&&GSResNA!=GSMass&&fMass+mAlph+ABarr+GSResNA<totMass)
            {
              if(barf) aLim+=pow(totMass-mNeut-mAlph-ABarr-GSResNA,2)*
                             evalph*Z*(Z-1)*(N-1)*(N-2)*12/(a-2)/(a-3)/(a-4);
              else     aLim+=pow(totMass-mNeut-mAlph-GSResNA,2)*
                             evalph*Z*(Z-1)*(N-1)*(N-2)*12/(a-2)/(a-3)/(a-4);
            }
            G4double r = aLim*G4UniformRand();
#ifdef debug
            G4cout<<"G4QN::EB:n, r="<<r<<",n="<<nLim<<",z="<<zLim<<",s="<<sLim<<",a="<<aLim
                  <<G4endl;
#endif
            three=true;                               // Flag of b+b+ResNuc decay
            if(!aLim) three=false;
            else if(r>sLim)
            {
              eMass = mAlph;
              dbQPDG= NAQPDG;
              rMass = GSResNA;
              rQPDG = naQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: N+A"<<G4endl;
#endif
            }
            else if(zLim<sLim&&r>zLim&&r<=sLim)
            {
              eMass = mLamb;
              dbQPDG= NLQPDG;
              rMass = GSResNL;
              rQPDG = nlQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: N+L"<<G4endl;
#endif
            }
            else if(nLim<zLim&&r>nLim&&r<=zLim)
            {
              eMass = mProt;
              dbQPDG= NPQPDG;
              rMass = GSResNP;
              rQPDG = npQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: N+P"<<G4endl;
#endif
            }
            else if(r<=nLim)
            {
              eMass = mNeut;
              dbQPDG= NNQPDG;
              rMass = GSResNN;
              rQPDG = nnQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: N+N"<<G4endl;
#endif
            }     
            else three=false;
          }
          else if(PDG==lPDG&&(lnCond&&lpCond&&llCond&&laCond)) // l+b+RN decay can't happen
          { //@@ Take into account Coulomb Barier Penetration Probability
#ifdef debug
            G4cout<<"G4QN::EB:*l*: n="<<lnCond<<",p="<<lpCond<<",l="<<llCond<<",a="<<laCond
                  <<G4endl;
#endif
            fMass=mLamb;
            fQPDG=lQPDG;
            G4double nLim=0.;
            if(N&&GSResNL!=GSMass&&fMass+mNeut+GSResNL<totMass)
              nLim+=(N+N)*pow(totMass-mNeut-mLamb-GSResNL,2);
            G4double zLim=nLim;
            if(Z&&GSResPL!=GSMass&&fMass+mProt+PBarr+GSResPL<totMass)
            {
              if(barf) zLim+=(Z+Z)*pow(totMass-mProt-mLamb-PBarr-GSResPL,2);
              else zLim+=(Z+Z)*pow(totMass-mProt-mLamb-GSResPL,2);
            }
            G4double sLim=zLim;
            if(S>1&&GSResLL!=GSMass&&fMass+mLamb+GSResLL<totMass)
              sLim+=(S-1)*pow(totMass-mLamb-mLamb-GSResLL,2);
            G4double aLim=sLim;
            if(evalph&&Z>1&&N>1&&a>4&&GSResLA!=GSMass&&fMass+mAlph+ABarr+GSResLA<totMass)
            {
              if(barf) aLim+=pow(totMass-mLamb-mAlph-ABarr-GSResLA,2)*
                             evalph*Z*(Z-1)*N*(N-1)*12/(a-2)/(a-3)/(a-4);
              else     aLim+=pow(totMass-mLamb-mAlph-GSResLA,2)*
                             evalph*Z*(Z-1)*N*(N-1)*12/(a-2)/(a-3)/(a-4);
            }
            G4double r = aLim*G4UniformRand();
#ifdef debug
            G4cout<<"G4QN::EB:l, r="<<r<<",n="<<nLim<<",z="<<zLim<<",s="<<sLim<<",a="<<aLim
                  <<G4endl;
#endif
            three=true;                               // Flag of b+b+ResNuc decay
            if(!aLim) three=false;
            else if(r>sLim)
            {
              eMass = mAlph;
              dbQPDG= LAQPDG;
              rMass = GSResLA;
              rQPDG = laQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: L+A"<<G4endl;
#endif
            }
            else if(zLim<sLim&&r>zLim&&r<=sLim)
            {
              eMass = mLamb;
              dbQPDG= LLQPDG;
              rMass = GSResLL;
              rQPDG = llQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: L+L"<<G4endl;
#endif
            }
            else if(nLim<zLim&&r>nLim&&r<=zLim)
            {
              eMass = mProt;
              dbQPDG= PLQPDG;
              rMass = GSResPL;
              rQPDG = plQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: L+P"<<G4endl;
#endif
            }
            else if(r<=nLim)
            {
              eMass = mNeut;
              dbQPDG= NLQPDG;
              rMass = GSResNL;
              rQPDG = nlQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: L+N"<<G4endl;
#endif
            }
            else three=false;
          }
          else if(PDG==aPDG&&(anCond&&apCond&&alCond&&aaCond)) // a+b+RN decay can't happen
          { //@@ Take into account Coulomb Barier Penetration Probability
#ifdef debug
            G4cout<<"G4QN::EB:*a*: n="<<anCond<<",p="<<apCond<<",l="<<alCond<<",a="<<aaCond
                  <<G4endl;
#endif
            fMass=mAlph;
            fQPDG=aQPDG;
            G4double nLim=0.;
            if(N>2&&GSResNA!=GSMass&&fMass+mNeut+ABarr+GSResNA<totMass)
            {
              if(barf) nLim+=(N-2)*pow(totMass-mNeut-mAlph-ABarr-GSResNA,2);
              else     nLim+=(N-2)*pow(totMass-mNeut-mAlph-GSResNA,2);
            }
            G4double zLim=nLim;
            if(Z>2&&GSResPA!=GSMass&&fMass+mProt+SAPBarr+GSResPA<totMass)
            {
              if(barf) zLim+=(Z-2)*pow(totMass-mProt-mAlph-SAPBarr-GSResPA,2);
              else     zLim+=(Z-2)*pow(totMass-mProt-mAlph-GSResPA,2);
            }
            G4double sLim=zLim;
            if(S&&GSResLA!=GSMass&&fMass+mLamb+ABarr+GSResLA<totMass)
            {
              if(barf) sLim+=S*pow(totMass-mLamb-mAlph-ABarr-GSResLA,2);
              else     sLim+=S*pow(totMass-mLamb-mAlph-GSResLA,2);
            }
            G4double aLim=sLim;
            if(evalph&&Z>3&&N>3&&a>7&&GSResAA!=GSMass&&fMass+mAlph+SAABarr+GSResAA<totMass)
            {
              if(barf) aLim+=pow(totMass-mAlph-mAlph-SAABarr-GSResAA,2)*
                             evalph*(Z-2)*(Z-3)*(N-2)*(N-3)*12/(a-5)/(a-6)/(a-7);
              else     aLim+=pow(totMass-mAlph-mAlph-GSResAA,2)*
                             evalph*(Z-2)*(Z-3)*(N-2)*(N-3)*12/(a-5)/(a-6)/(a-7);
            }
            G4double r = aLim*G4UniformRand();
#ifdef debug
            G4cout<<"G4QN::EB:a, r="<<r<<",n="<<nLim<<",z="<<zLim<<",s="<<sLim<<",a="<<aLim
                  <<G4endl;
#endif
            three=true;                               // Flag of b+b+ResNuc decay
            if(!aLim) three=false;
            else if(r>sLim)
            {
              eMass = mAlph;
              dbQPDG= AAQPDG;
              rMass = GSResAA;
              rQPDG = aaQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: A+A"<<G4endl;
#endif
            }
            else if(zLim<sLim&&r>zLim&&r<=sLim)
            {
              eMass = mLamb;
              dbQPDG= LAQPDG;
              rMass = GSResLA;
              rQPDG = laQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: A+L"<<G4endl;
#endif
            }
            else if(nLim<zLim&&r>nLim&&r<=zLim)
            {
              eMass = mProt;
              dbQPDG= PAQPDG;
              rMass = GSResPA;
              rQPDG = paQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: A+P"<<G4endl;
#endif
            }
            else if(r<=nLim)
            {
              eMass = mNeut;
              dbQPDG= NAQPDG;
              rMass = GSResNA;
              rQPDG = naQPDG;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateBary: A+N"<<G4endl;
#endif
            }
            else three=false;
          }
          else three=false;
          if(rMass<1600.)
          {
            if     (rQPDG==pQPDG)rMass=mProt;
            else if(rQPDG==nQPDG)rMass=mNeut;
            else if(rQPDG==lQPDG)rMass=mLamb;
          }
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBary:evaBar="<<eMass<<bQPDG<<",resN="<<rMass<<rQPDG
                <<",secB="<<fMass<<",three="<<three<<G4endl;
#endif
        }
      }
      else // =-------------=> Just decay in a baryon and a residual (to avoid gamma-decay)
      { //@@ Take into account Coulomb Barier Penetration Probability (?? - Emergency)
        G4double nLim=0.;
        if(nFlag&&mNeut+GSResNn<totMass)
        {
          G4double ken=totMass-mNeut-GSResNn;
          nLim+=N*CoulBarPenProb(0.,ken,0,1)*sqrt(ken);
        }
        G4double zLim=nLim;
        if(pFlag&&mProt+PBarr+GSResNp<totMass)
        {
          G4double ken=totMass-mProt-GSResNp;
          if(barf) ken-=PBarr;
          zLim+=Z*CoulBarPenProb(PBarr,ken,1,1)*sqrt(ken);
        }
        G4double sLim=zLim;
        if(lFlag&&mLamb+GSResNl<totMass)
        {
          G4double ken=totMass-mLamb-GSResNl;
          sLim+=S*CoulBarPenProb(0.,ken,0,1)*sqrt(ken);
        }
        G4double aLim=sLim;
        if(evalph&&aFlag&&mAlph+GSResNa+ABarr<totMass)
        {
          G4double ken=totMass-mAlph-GSResNa;
          if(barf) ken-=ABarr;
          aLim+=CoulBarPenProb(ABarr,ken,2,4)*sqrt(ken)*evalph*Z*(Z-1)*N*(N-1)
                *6/a1/(a-2)/(a-3);
        }
        G4double r = aLim*G4UniformRand();
#ifdef debug
        G4cout<<"G4QNucl::EvapBar:2Decay r="<<r<<",nLim="<<nLim<<",zLim="<<zLim<<",sLim="
              <<sLim<<",nF="<<nFlag<<",pF="<<pFlag<<",lF="<<lFlag<<",aF="<<aFlag<<G4endl;
#endif
        if     (aFlag&&r>sLim)
        {
          bQPDG=aQPDG;
          eMass=mAlph;
          rQPDG=AQPDG;
          rMass=GSResNa;
        }
        else if(lFlag&&r>=zLim&&r<=sLim&&zLim<sLim)
        {
          bQPDG=lQPDG;
          eMass=mLamb;
          rQPDG=LQPDG;
          rMass=GSResNl;
        }
        else if(nFlag&&r>=nLim&&r<=zLim&&nLim<zLim)
        {
          bQPDG=pQPDG;
          eMass=mProt;
          rQPDG=PQPDG;
          rMass=GSResNp;
        }
        else if(pFlag&&r<nLim)
        {
          bQPDG=nQPDG;
          eMass=mNeut;
          rQPDG=NQPDG;
          rMass=GSResNn;
        }
        else
        {
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBaryon: Photon #2-B#, dM="<<totMass-GSMass<<G4endl;
#endif
          bQPDG=gQPDG;
          rQPDG=GetQPDG();
          eMass=0.;
          rMass=GSMass;
        }
#ifdef debug
        G4cout<<"G4QNucl::EvaporateBaryon: b="<<eMass<<bQPDG<<",r="<<rMass<<rQPDG<<G4endl;
#endif
      }
      if(three)           // Decay in two baryons + Residual Nucleus
      {
#ifdef debug
          G4cout<<"G4QNucl::EvaporateBaryon:Decay in 3 particles"<<G4endl;
#endif
        h1mom=G4LorentzVector(0.,0.,0.,eMass);
        h2mom=G4LorentzVector(0.,0.,0.,rMass);
        h3mom=G4LorentzVector(0.,0.,0.,fMass);
        if(!DecayIn3(h1mom,h2mom,h3mom))
        {
#ifdef debug
          G4cout<<"*G4QNucl::EvaporateBaryon:Decay M="<<totMass<<",b="<<eMass<<bQPDG
          <<",f="<<fMass<<fQPDG<<",r="<<rMass<<rQPDG<<G4endl;
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
#ifdef debug
          G4cout<<"G4QN::EvaB:eM="<<eMass<<"+rM="<<rMass<<" ="<<eMass+rMass<<" < "<<totMass
                <<",c="<<cntr<<" < cm="<<cntm<<", bPDG="<<bQPDG<<", rPDG="<<rQPDG<<G4endl;
#endif
          if(rMass<1600.)
          {
            if     (rQPDG==pQPDG)rMass=mProt;
            else if(rQPDG==nQPDG)rMass=mNeut;
            else if(rQPDG==lQPDG)rMass=mLamb;
          }
          h1mom=G4LorentzVector(0.,0.,0.,eMass);
          h2mom=G4LorentzVector(0.,0.,0.,rMass);
        }
        else if(totMass>mNeut+GSResNn)               // Neutron if 2-Decay failed
        {
#ifdef debug
          G4cout<<"G4QNucl::EvaporateBaryon: Neutron , dM="<<totMass-GSResNn-mNeut<<G4endl;
#endif
          bQPDG=nQPDG;
          rQPDG=NQPDG;
          h1mom=G4LorentzVector(0.,0.,0.,mNeut);
          h2mom=G4LorentzVector(0.,0.,0.,GSResNn);      
        }
        else if(totMass>mProt+PBarr+GSResNp)               // Proton if 2-Decay failed
        {
#ifdef debug
          G4cout<<"G4QNucl::EvaporateBaryon: Proton , dM="<<totMass-GSResNp-mProt<<G4endl;
#endif
          bQPDG=pQPDG;
          rQPDG=PQPDG;
          h1mom=G4LorentzVector(0.,0.,0.,mProt);
          h2mom=G4LorentzVector(0.,0.,0.,GSResNp);      
        }
        else if(totMass>mAlph+ABarr+GSResNa)               // Alpha if 2-Decay failed
        {
#ifdef debug
          G4cout<<"G4QNucl::EvaporateBaryon: Alpha , dM="<<totMass-GSResNa-mAlph<<G4endl;
#endif
          bQPDG=aQPDG;
          rQPDG=AQPDG;
          h1mom=G4LorentzVector(0.,0.,0.,mAlph);
          h2mom=G4LorentzVector(0.,0.,0.,GSResNa);      
        }
        else if(totMass>GSMass)               // Photon if 2-Decay failed
        {
#ifdef debug
          G4cout<<"G4QNucl::EvaporateBaryon:Photon ### 2 ###, dM="<<totMass-GSMass<<G4endl;
#endif
          bQPDG=gQPDG;
          rQPDG=GetQPDG();
          h1mom=G4LorentzVector(0.,0.,0.,0.);
          h2mom=G4LorentzVector(0.,0.,0.,GSMass);      
        }
        else
        {
          G4cerr<<"***G4QNucl::EvaporateBaryon: Cann't evaporate even gamma (1)"<<G4endl;
          return false;
        }
      }
    }
    else // ==> Decay in 3 Baryons + Residual is impossible at this point
    {
      if(secB)                        // Decay in 2Baryons(2a,a+bary)+ResidN is possible
      //if(2>3)
      {
#ifdef debug
        G4cout<<"G4QNucleus::EvaporateBaryon: Decay in 2 baryons"<<G4endl;
#endif
        G4bool tpd=true;
        //@@ Coulomb Barrier penetration can be added
        G4double alp=0.;
        if(aSecF)alp=evalph*Z*(Z-1)*N*(N-1)*10/(a-2)/(a-3)/(a-4);
        G4double  nnLim=0.;
        if(nnFlag&&totMass>mNeut+mNeut+GSResNN)
          nnLim+=N*(N-1)*pow(totMass-mNeut-mNeut-GSResNN,2);
        G4double  nzLim=nnLim;
        if(npFlag&&totMass>mNeut+mProt+PBarr+GSResNP)
        {
          if(barf) nzLim+=2*N*Z*pow(totMass-mNeut-mProt-PBarr-GSResNP,2);
          else     nzLim+=2*N*Z*pow(totMass-mNeut-mProt-GSResNP,2);
        }
        G4double  zzLim=nzLim;
        if(ppFlag&&totMass>mProt+mProt+SPPBarr+GSResPP)
        {
          if(barf) zzLim+=Z*(Z-1)*pow(totMass-mProt-mProt-SPPBarr-GSResPP,2);
          else     zzLim+=Z*(Z-1)*pow(totMass-mProt-mProt-GSResPP,2);
        }
        G4double  nlLim=zzLim;
        if(nlFlag&&totMass>mNeut+mLamb+GSResNL)
          nlLim+=2*N*S*pow(totMass-mNeut-mLamb-GSResNL,2);
        G4double  zlLim=nlLim;
        if(plFlag&&totMass>mProt+PBarr+mLamb+GSResPL)
        {
          if(barf) zlLim+=2*Z*S*pow(totMass-mProt-mLamb-PBarr-GSResPL,2);
          else     zlLim+=2*Z*S*pow(totMass-mProt-mLamb-GSResPL,2);
        }
        G4double  llLim=zlLim;
        if(llFlag&&totMass>mLamb+mLamb+GSResLL)
          llLim+=S*(S-1)*pow(totMass-mLamb-mLamb-GSResLL,2);
        G4double  naLim=llLim;
        if(naFlag&&totMass>mNeut+mAlph+ABarr+GSResNA)
        {
          if(barf) naLim+=(N-3)*alp*pow(totMass-mNeut-mAlph-ABarr-GSResNA,2);
          else     naLim+=(N-3)*alp*pow(totMass-mNeut-mAlph-GSResNA,2);
        }
        G4double  zaLim=naLim;
        if(paFlag&&totMass>mProt+SAPBarr+mAlph+GSResPA)
        {
          if(barf) zaLim+=(Z-3)*alp*pow(totMass-mProt-mAlph-SAPBarr-GSResPA,2);
          else     zaLim+=(Z-3)*alp*pow(totMass-mProt-mAlph-GSResPA,2);
        }
        G4double  laLim=zaLim;
        if(laFlag&&totMass>mLamb+mAlph+ABarr+GSResLA)
        {
          if(barf) laLim+=S*alp*pow(totMass-mLamb-mAlph-ABarr-GSResLA,2);
          else     laLim+=S*alp*pow(totMass-mLamb-mAlph-GSResLA,2);
        }
        G4double  aaLim=laLim;
        if(evalph&&aaFlag&&totMass>mAlph+mAlph+SAABarr+GSResAA)
        {
          if(barf) aaLim+=alp*pow(totMass-mAlph-mAlph-SAABarr-GSResAA,2)*
                          evalph*(Z-2)*(Z-3)*(N-2)*(N-3)*7/(a-5)/(a-6)/(a-7);
          else     aaLim+=alp*pow(totMass-mAlph-mAlph-GSResAA,2)*
                          evalph*(Z-2)*(Z-3)*(N-2)*(N-3)*7/(a-5)/(a-6)/(a-7);
        }
        G4double r = aaLim*G4UniformRand();
        if     (aaLim&&laLim<r)
        {
          dbQPDG= AAQPDG;
          eMass=mAlph;
          fMass=mAlph;
          rQPDG=aaQPDG;
          rMass=GSResAA;
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: A+A, e="<<eMass<<",f="<<fMass<<",r="<<rMass<<G4endl;
#endif
        }
        else if(aaLim&&zaLim<r&&r<laLim)
        {
          dbQPDG= LAQPDG;
          eMass=mAlph;
          fMass=mLamb;
          rQPDG=laQPDG;
          rMass=GSResLA;
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: A+L, e="<<eMass<<",f="<<fMass<<",r="<<rMass<<G4endl;
#endif
        }
        else if(aaLim&&naLim<r&&r<zaLim)
        {
          dbQPDG= PAQPDG;
          eMass=mAlph;
          fMass=mProt;
          rQPDG=paQPDG;
          rMass=GSResPA;
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: A+P, e="<<eMass<<",f="<<fMass<<",r="<<rMass<<G4endl;
#endif
        }
        else if(aaLim&&llLim<r&&r<naLim)
        {
          dbQPDG= NAQPDG;
          eMass=mAlph;
          fMass=mNeut;
          rQPDG=naQPDG;
          rMass=GSResNA;
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: A+N, e="<<eMass<<",f="<<fMass<<",r="<<rMass<<G4endl;
#endif
        }
        else if(aaLim&&zlLim<r&&r<llLim)
        {
          dbQPDG= LLQPDG;
          eMass=mLamb;
          fMass=mLamb;
          rQPDG=llQPDG;
          rMass=GSResLL;
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: L+L, e="<<eMass<<",f="<<fMass<<",r="<<rMass<<G4endl;
#endif
        }
        else if(aaLim&&nlLim<r&&r<zlLim)
        {
          dbQPDG= PLQPDG;
          eMass=mLamb;
          fMass=mProt;
          rQPDG=plQPDG;
          rMass=GSResPL;
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: L+p, e="<<eMass<<",f="<<fMass<<",r="<<rMass<<G4endl;
#endif
        }
        else if(aaLim&&zzLim<r&&r<nlLim)
        {
          dbQPDG= NLQPDG;
          eMass=mLamb;
          fMass=mNeut;
          rQPDG=nlQPDG;
          rMass=GSResNL;
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: L+n, e="<<eMass<<",f="<<fMass<<",r="<<rMass<<G4endl;
#endif

        }
        else if(aaLim&&nzLim<r&&r<zzLim)
        {
          dbQPDG= PPQPDG;
          eMass=mProt;
          fMass=mProt;
          rQPDG=ppQPDG;
          rMass=GSResPP;
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: p+p, e="<<eMass<<",f="<<fMass<<",r="<<rMass<<G4endl;
#endif
        }
        else if(aaLim&&nnLim<r&&r<nzLim)
        {
          dbQPDG= NPQPDG;
          eMass=mNeut;
          fMass=mProt;
          rQPDG=npQPDG;
          rMass=GSResNP;
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: n+p, e="<<eMass<<",f="<<fMass<<",r="<<rMass<<G4endl;
#endif
        }
        else if(aaLim&&r<nnLim)
        {
          dbQPDG= NNQPDG;
          eMass=mNeut;
          fMass=mNeut;
          rQPDG=nnQPDG;
          rMass=GSResNN;
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: n+n, e="<<eMass<<",f="<<fMass<<",r="<<rMass<<G4endl;
#endif
        }
        //Two particle decay only possible (not frequent event!)
        else if(nFlag)
        {
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBaryon:Photon ### Decay in neutron ###"<<G4endl;
#endif
          tpd=false;
          bQPDG=nQPDG;
          rQPDG=NQPDG;
          eMass=mNeut;
          rMass=GSResNn;
        }
        else if(pFlag)
        {
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBaryon:Photon ### Decay in proton ###"<<G4endl;
#endif
          tpd=false;
          bQPDG=pQPDG;
          rQPDG=PQPDG;
          eMass=mProt;
          rMass=GSResNp;
        }
        else if(aFlag)
        {
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBaryon:Photon ### Decay in alpha ###"<<G4endl;
#endif
          tpd=false;
          bQPDG=aQPDG;
          rQPDG=AQPDG;
          eMass=mAlph;
          rMass=GSResNa;
        }
        else if(lFlag)
        {
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBaryon:Photon ### Decay in Lambda ###"<<G4endl;
#endif
          tpd=false;
          bQPDG=lQPDG;
          rQPDG=LQPDG;
          eMass=mLamb;
          rMass=GSResNl;
        }
        else
        {
#ifdef debug
          G4cout<<"G4QNuc::EvaporBaryon: Photon ### 3-Big ###,dM="<<totMass-GSMass<<G4endl;
#endif
          tpd=false;
          bQPDG=gQPDG;
          rQPDG=GetQPDG();
          eMass=0.;
          rMass=GSMass;
        }
        if(tpd)
        {
          h1mom=G4LorentzVector(0.,0.,0.,eMass);
          h2mom=G4LorentzVector(0.,0.,0.,rMass);
          h3mom=G4LorentzVector(0.,0.,0.,fMass);
          if(!DecayIn3(h1mom,h2mom,h3mom))
          {
#ifdef debug
            G4cout<<"*G4QNucl::EvaporateBaryon:Decay M="<<totMass<<",b="<<eMass<<bQPDG
            <<",f="<<fMass<<fQPDG<<",r="<<rMass<<rQPDG<<G4endl;
#endif
            return false;
          }
          h1mom+=h3mom;
          bQPDG=dbQPDG;
#ifdef debug
          G4double sma=h1mom.m();
          G4double dma=sma-eMass-fMass;
          G4cout<<"G4QNuc::EvapBar:s="<<sma<<",e="<<eMass<<",f="<<fMass<<",d="<<dma<<",rM="
                <<rMass<<G4endl;
#endif
        }
        else
        {
          if(rMass<1600.)
          {
            if     (rQPDG==pQPDG)rMass=mProt;
            else if(rQPDG==nQPDG)rMass=mNeut;
            else if(rQPDG==lQPDG)rMass=mLamb;
          }
          h1mom=G4LorentzVector(0.,0.,0.,eMass);
          h2mom=G4LorentzVector(0.,0.,0.,rMass);      
          if(!DecayIn2(h1mom,h2mom))
          {
#ifdef debug
            G4cout<<"***G4QNucleus::EvaporateBaryon: Emergency Decay M="<<totMass<<",b="
                  <<bQPDG<<h1->GetQC()<<eMass<<",r="<<rQPDG<<h2->GetQC()<<rMass<<G4endl;
#endif
            return false;
          }
          h1->SetQPDG(bQPDG);
          h2->SetQPDG(rQPDG);
          h1->Set4Momentum(h1mom);
          h2->Set4Momentum(h2mom);
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: Emergency decay is done for b="<<bQPDG<<h1->GetQC()
                <<h1mom<<h1mom.m()<<", r="<<rQPDG<<h2->GetQC()<<h2mom<<h2mom.m()<<G4endl;
#endif
          return true;
        }
      }
      else                                     // Only decay in Baryon+Residual is possible
      {
#ifdef debug
        G4cout<<"G4QNucleus::EvaporateBaryon: Decay in Baryon+Resid"<<G4endl;
#endif
        //@@ Take into account Coulomb Barier Penetration Probability
        G4double nLim=0.;
        if(nFlag&&mNeut+GSResNn<totMass)
        {
          G4double ken=totMass-mNeut-GSResNn;
          nLim+=N*CoulBarPenProb(0.,ken,0,1)*sqrt(ken);
        }
        G4double zLim=nLim;
        if(pFlag&&mProt+PBarr+GSResNp<totMass)
        {
          G4double ken=totMass-mProt-GSResNp;
          if(barf) ken-=PBarr;
          zLim+=Z*CoulBarPenProb(PBarr,ken,1,1)*sqrt(ken);
        }
        G4double sLim=zLim;
        if(lFlag&&mLamb+GSResNl<totMass)
        {
          G4double ken=totMass-mLamb-GSResNl;
          sLim+=S*CoulBarPenProb(0.,ken,0,1)*sqrt(ken);
        }
        G4double aLim=sLim;
        if(aFlag&&mAlph+GSResNa+ABarr<totMass)
        {
          G4double ken=totMass-mAlph-GSResNa;
          if(barf) ken-=ABarr;
          aLim+=CoulBarPenProb(ABarr,ken,2,4)*sqrt(ken)*evalph*Z*(Z-1)*N*(N-1)
                *6/a1/(a-2)/(a-3);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBaryon:al="<<evalph<<",k="<<ken<<",P="
                <<CoulBarPenProb(ABarr,ken,2,4)<<G4endl;
#endif
        }
        G4double r = aLim*G4UniformRand();
#ifdef debug
        G4cout<<"G4QN::EB:DecIn2#2#r="<<r<<",nL="<<nLim<<",zL="<<zLim<<",sL="<<sLim<<",aL="
              <<aLim<<",nF="<<nFlag<<",pF="<<pFlag<<",lF="<<lFlag<<",aF="<<aFlag<<G4endl;
#endif
        if     (aFlag&&r>sLim)
        {
          bQPDG=aQPDG;
          eMass=mAlph;
          rQPDG=AQPDG;
          rMass=GSResNa;
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBaryon: Decay in A + rA="<<GSResNa+mAlph<<G4endl;
#endif
        }
        else if(lFlag&&r>zLim&&r<sLim)
        {
          bQPDG=lQPDG;
          eMass=mLamb;
          rQPDG=LQPDG;
          rMass=GSResNl;
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBaryon: Decay in L + rA="<<GSResNl+mLamb<<G4endl;
#endif
        }
        else if(pFlag&&r>nLim&&r<zLim)
        {
          bQPDG=pQPDG;
          eMass=mProt;
          rQPDG=PQPDG;
          rMass=GSResNp;
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBaryon: Decay in P + rA="<<GSResNp+mProt<<G4endl;
#endif
        }
        else if(nFlag&&r<nLim)
        {
          bQPDG=nQPDG;
          eMass=mNeut;
          rQPDG=NQPDG;
          rMass=GSResNn;
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateBaryon: Decay in N + rA="<<GSResNn+mNeut<<G4endl;
#endif
        }
        else if(mProt+GSResNp<totMass)
        {
#ifdef debug
          G4cout<<"G4QNucl::EvapBar: Emergency Proton, dM="<<totMass-GSResNp-mProt<<G4endl;
#endif
          bQPDG=pQPDG;
          rQPDG=PQPDG;
          eMass=mProt;
          rMass=GSResNp;
        }
        else if(mAlph+GSResNa<totMass)
        {
#ifdef debug
          G4cout<<"G4QNucl::EvapBar: Emergency Alpha, dM="<<totMass-GSResNa-mAlph<<G4endl;
#endif
          bQPDG=aQPDG;
          rQPDG=AQPDG;
          eMass=mAlph;
          rMass=GSResNa;
        }
        else
        {
#ifdef debug
          G4cout<<"G4QNuc::EvapBaryon: Photon ### 4-Big ###, dM="<<totMass-GSMass<<G4endl;
#endif
          bQPDG=gQPDG;
          rQPDG=GetQPDG();
          eMass=0.;
          rMass=GSMass;
        }
        if(rMass<1600.)
        {
          if     (rQPDG==pQPDG)rMass=mProt;
          else if(rQPDG==nQPDG)rMass=mNeut;
          else if(rQPDG==lQPDG)rMass=mLamb;
        }
        h1mom=G4LorentzVector(0.,0.,0.,eMass);
        h2mom=G4LorentzVector(0.,0.,0.,rMass);      
      }
    }
    if(!DecayIn2(h1mom,h2mom))
    {
#ifdef debug
      G4cout<<"*G4QNucleus::EvaporateBaryon: Decay M="<<totMass<<",b="<<bQPDG<<h1->GetQC()
      <<eMass<<",r="<<rQPDG<<h2->GetQC()<<rMass<<G4endl;
#endif
      return false;
    }
#ifdef debug
    G4cout<<"G4QN::EvaB: **RESULT** b="<<bQPDG<<h1mom<<", r="<<rQPDG<<h2mom<<G4endl;
#endif
    h1->SetQPDG(bQPDG);
    h2->SetQPDG(rQPDG);
    h1->Set4Momentum(h1mom);
    h2->Set4Momentum(h2mom);
#ifdef debug
    G4cout<<"G4QNucleus::EvaporateBaryon: Evaporation is done for b="<<bQPDG<<h1->GetQC()
       <<h1mom<<h1mom.m()<<", r="<<rQPDG<<h2->GetQC()<<h2mom<<h2mom.m()<<G4endl;
#endif
    return true;
  }
  else if(a==1)
  {
#ifdef debug
    G4cerr<<"***G4QNucleus::EvaporateBaryon: ??? A=1"<<G4endl;
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
  if(a<-2) G4cerr<<"***G4QNucleus::EvaporateBaryon: A="<<a<<G4endl;
  return false;
}
// End of "EvaporateBaryon"

// Randomize bimomial low 
G4int G4QNucleus::RandomizeBinom(G4double p,G4int aN)
{
  G4double r = G4UniformRand();
  G4double d = 1.-p;
  if(d<=0.) return 0;
  G4double v = pow(d,aN);
  G4double s_value = v;
  if(r<s_value) return 0;
  G4int i=0;
  G4int j=aN+1;
  G4double f=p/d;
  while(r>s_value && i<aN)
  {
    j--;
    i++;
    v*=j*f/i;
    s_value+=v;
  }
  return i;
}

//Initialize a Candidate vector for the instance of a Quasmon
void G4QNucleus::InitCandidateVector(G4QCandidateVector& theQCandidates,
                                    G4int maxMes, G4int maxBar, G4int maxClst)
{
  static const G4int nOfMesons =45; //a#of S=0,1,2,3,4 Mesons, => candidates to hadrons
  static const G4int nOfBaryons=72; //a#of 1/2,3/2,5/2,7/2 Baryons => candidates to hadrons
  // Scalar resonances   (0):           Eta,Pi0,Pi+,APi-,Ka0,Ka+,AKa0,AKa-,Eta*
  static G4int mesonPDG[nOfMesons]  =  {221,111,211,-211,311,321,-311,-321,331,    //  0- 8
  // Vector resonances   (1):           omega,Rh0,Rh+,Rho-,K0*,K+*,AK0*,AK-*,Phi
                                        223,113,213,-213,313,323,-313,-323,333,    //  9-18
  // Tensor D-resonances (2):           f2 ,a20,a2+, a2-,K20,K2+,AK20,AK2-,f2'
                                        225,115,215,-215,315,325,-315,-325,335,    // 19-27
  // Tensor F-resonances (3):           om3,ro3,r3+,rh3-,K30,K3+,AK30,AK3-,Phi3
                                        227,117,217,-217,317,327,-317,-327,337,    // 28-35
  // Tensor G-resonances (4):           f4 ,a40,a4+, a4-,K40,K4+,AK40,AK4-,f4'
                                        229,119,219,-219,319,329,-319,-329,339};   // 36-44
  // Baryon octet      (1/2):          n  , an  , p  , ap  ,lamb,alamb, sig-,asig-
  static G4int baryonPDG[nOfBaryons]={2112,-2112,2212,-2212,3122,-3122,3112,-3112, // 45-52
  // Hyperon octet     (1/2):         sig0,asig0,sig+,asig+,ksi-,aksi-,ksi0,aksi0
                                      3212,-3212,3222,-3222,3312,-3312,3322,-3322, // 53-60
  // Baryon decuplet   (3/2):   del-,adel-,del0,adel0,del+,adel+,dl++,adl++,sis-,asis-
                                1114,-1114,2114,-2114,2214,-2214,2224,-2224,3114,-3114,//70
  //                            sis0,asis0,sis+,asis+,kss-,akss-,kss0,akss0,omeg,aomeg
                                3214,-3214,3224,-3224,3314,-3314,3324,-3324,3334,-3334,//80
  // Baryon octet      (5/2):         n5/2,an5/2,p5/2,ap5/2,l5/2,al5/2,si5-,asi5-
                                      2116,-2116,2216,-2216,3126,-3126,3116,-3116, // 81-88
  //                                  si50,asi50,si5+,asi5+,ks5-,aks5-,ks50,aks50
                                      3216,-3216,3226,-3226,3316,-3316,3326,-3326, // 89-96
  // Baryon decuplet   (7/2):  dl5-,adl5-,dl50,adl50,dl5+,adl5+,d5++,ad5++,si5-,asi5-
                              1118,-1118,2118,-2118,2218,-2218,2228,-2228,3118,-3118, //106
  //                          si50,asi50,si5+,asi5+,ks5-,aks5-,ks50,aks50,ome5,aome5
                              3218,-3218,3228,-3228,3318,-3318,3328,-3328,3338,-3338};//116
  G4int i=0;
#ifdef debug
  G4int ind=0;
#endif
  G4int iQC = theQCandidates.size();
  if(iQC) for(G4int jq=0; jq<iQC; jq++) delete theQCandidates[jq];
  theQCandidates.clear();
  if(maxMes>nOfMesons) maxMes=nOfMesons;
  if(maxMes>=0) for (i=0; i<maxMes; i++) 
  {
    theQCandidates.push_back(new G4QCandidate(mesonPDG[i]));
#ifdef debug
    G4cout<<"G4QNucleus::InitCandidateVector: "<<ind++<<", Meson # "<<i<<" with code = "
          <<mesonPDG[i]<<", QC="<<theQCandidates[i]->GetQC()<<" is initialized"<<G4endl;
#endif
  }
  if(maxBar>nOfBaryons) maxBar=nOfBaryons;
  if(maxBar>=0) for (i=0; i<maxBar; i++) 
  {
#ifdef debug
    G4cout<<"G4QNucleus::InitCandidateVector: define PDG="<<baryonPDG[i]<<G4endl;
#endif
    G4QCandidate* curBar=new G4QCandidate(baryonPDG[i]);
#ifdef debug
    G4cout<<"G4QNucleus::InitCandidateVector: current baryon is defined"<<G4endl;
#endif
    theQCandidates.push_back(curBar); // delete equivalent
#ifdef debug
    G4cout<<"G4Nucleus::InitCandidateVector: "<<ind++<<", Baryon # "<<i<<" with code = "
          <<baryonPDG[i]<< ", QC="<<theQCandidates[i]->GetQC()<<" is initialized"<<G4endl;
#endif
  }
  if(maxClst>=0) for (i=0; i<maxClst; i++) 
  {
    G4int clustQCode = i+G4QPDGCode().GetNQHadr(); //Q-codes of cluster in the CHIPS world
    G4QPDGCode clustQPDG;
    clustQPDG.InitByQCode(clustQCode);
    G4int clusterPDG=clustQPDG.GetPDGCode();
    theQCandidates.push_back(new G4QCandidate(clusterPDG)); // delete equivalent
#ifdef debug
    G4cout<<"G4QNucleus::InitCandidateVector:"<<ind++<<", Cluster # "<<i<<" with code = "
          <<clusterPDG<<", QC="<<clustQPDG.GetQuarkContent()<<" is initialized"<<G4endl;
#endif
  }
} // End of "InitCandidateVector"

// Calculate a#of Z,N,L-clusters in the nucleus and fill candidate's probabilities
void G4QNucleus::PrepareCandidates(G4QCandidateVector& theQCandidates, G4bool piF,
                                   G4bool gaF, G4LorentzVector pLV)
{
  static const G4LorentzVector zeroLV(0.,0.,0.,0.);
  G4double ze  = Z;
  G4double ne  = N;
  G4double se  = S;
  G4double ae  = Z+N+S;
  G4double aea = ae*(ae-1);
  G4double ae2 = aea/2.;
  G4double ze1 = dZ + 1;
  G4double ne1 = dN + 1;
  G4double se1 = dS + 1;
  G4double ae0 = dZ + dN + dS;
  G4double ae1 = ae0 + 1;
  G4double pos = probVect[0];           // Value of Pre-Probability for VacuumHadronization
#ifdef cldebug
  G4int mac=6;                          // Maximum cluster # for fixed baryon number
#endif
  G4int cca=0;                          // Counter of clusters for the same baryon number
  G4int acm=0;                          // Threshold ac value
  G4int mCand=theQCandidates.size();    // Full set of candidates made in UpdateClusters
  G4double s_value=0.;                  // Prototype of summ for constant A (=ac>2)
  G4double comb=ae0*(ae0-1)/2;          // Product up to ac=2
  if(comb<=0.) comb=1.;
#ifdef cldebug
  G4double sZ=0.;                       // Percent of protons
  G4double sN=0.;                       // Percent of neutrons
  G4cout<<"G4QN::PC:C#"<<mCand<<",dZ="<<dZ<<",dN="<<dN<<",ZNS="<<Z<<","<<N<<","<<S<<G4endl;
#endif
  for (G4int index=0; index<mCand; index++)
  {
    G4QCandidate* curCand=theQCandidates[index];
    G4int cPDG  = curCand->GetPDGCode();
    G4int cBN   = curCand->GetBaryonNumber();
    G4int cST   = curCand->GetStrangeness();
    // ***********************************************************************************
    // These are first 117 candidates which are defined in G4QNucleus::InitCandidateVector
    // ***!!!*** if they are changed there the corresponding change must be done here
    //static const G4int nOfMesons =45;//a#of S=0,1,2,3,4 Mesons, => candidates to hadrons
    //static const G4int nOfBaryons=72;//a#of 1/2,3/2,5/2,7/2 Baryons => cand's to hadrons
    // Scalar resonances   (0):           Eta,Pi0,Pi+,APi-,Ka0,Ka+,AKa0,AKa-,Eta*
    //static G4int mesonPDG[45]  =  {221,111,211,-211,311,321,-311,-321,331,       //  0- 8
    // Vector resonances   (1):    omega,Rh0,Rh+,Rho-,K0*,K+*,AK0*,AK-*,Phi
    //                              223,113,213,-213,313,323,-313,-323,333,        //  9-18
    // Tensor D-resonances (2):     f2 ,a20,a2+, a2-,K20,K2+,AK20,AK2-,f2'
    //                              225,115,215,-215,315,325,-315,-325,335,        // 19-27
    // Tensor F-resonances (3):     om3,ro3,r3+,rh3-,K30,K3+,AK30,AK3-,Phi3
    //                              227,117,217,-217,317,327,-317,-327,337,        // 28-35
    // Tensor G-resonances (4):     f4 ,a40,a4+, a4-,K40,K4+,AK40,AK4-,f4'
    //                              229,119,219,-219,319,329,-319,-329,339};       // 36-44
    // Baryon octet       (1/2):    n  , an  , p  , ap  ,lamb,alamb, sig-,asig-
    //static G4int baryonPDG[72]={2112,-2112,2212,-2212,3122,-3122,3112,-3112,     // 45-52
    // Hyperon octet     (1/2): sig0,asig0,sig+,asig+,ksi-,aksi-,ksi0,aksi0
    //                            3212,-3212,3222,-3222,3312,-3312,3322,-3322,     // 53-60
    // Baryon decuplet   (3/2): del-,adel-,del0,adel0,del+,adel+,dl++,adl++,sis-,asis-
    //                          1114,-1114,2114,-2114,2214,-2214,2224,-2224,3114,-3114,//70
    //                          sis0,asis0,sis+,asis+,kss-,akss-,kss0,akss0,omeg,aomeg
    //                          3214,-3214,3224,-3224,3314,-3314,3324,-3324,3334,-3334,//80
    // Baryon octet      (5/2): n5/2,an5/2,p5/2,ap5/2,l5/2,al5/2,si5-,asi5-
    //                          2116,-2116,2216,-2216,3126,-3126,3116,-3116,       // 81-88
    //                          si50,asi50,si5+,asi5+,ks5-,aks5-,ks50,aks50
    //                          3216,-3216,3226,-3226,3316,-3316,3326,-3326,       // 89-96
    // Baryon decuplet  (7/2): dl5-,adl5-,dl50,adl50,dl5+,adl5+,d5++,ad5++,si5-,asi5-
    //                        1118,-1118,2118,-2118,2218,-2218,2228,-2228,3118,-3118, //106
    //                        si50,asi50,si5+,asi5+,ks5-,aks5-,ks50,aks50,ome5,aome5
    //                        3218,-3218,3228,-3228,3318,-3318,3328,-3328,3338,-3338};//116
    // One should take into account, that #of mesons & baryons can be cut in G4Quas::HadrQE
    //G4int nP= theWorld->GetQPEntries(); // A#of initialized particles in CHIPS World
    ////@@ Make special parametyer to cut high resonances for nuclear fragmentation !!
    //G4int          nMesons  = 45;
    //if     (nP<34) nMesons  =  9;
    //else if(nP<51) nMesons  = 18;
    //else if(nP<65) nMesons  = 27;
    //else if(nP<82) nMesons  = 36;
    //G4int          nBaryons = 72;
    //if     (nP<45) nBaryons = 16;
    //else if(nP<59) nBaryons = 36;
    //else if(nP<76) nBaryons = 52;
    // **********************************************************************************
    //G4int cS    = curCand->GetStrangeness();
    //if(piF&&gaF&&cPDG!=90000001&&cPDG!=90001000) // Both flags, in case of pi-first-int
    //if(piF&&gaF&&cBN!=1&&cBN!=3) // Both flags, which is in case of pi-first-int
    if(piF&&gaF&&cBN!=1)// @@ Should be both, which is in case of pi-first-interaction @@ ?
    //if(piF&&gaF&&cBN!=1&&cBN!=4) // Should be both, in case of pi-first-interaction
    {
      curCand->SetPreProbability(0.);  
      curCand->SetDenseProbability(0.); 
      curCand->SetPossibility(false);    
#ifdef cldebug
      if(cPDG==90001001) G4cout<<"G4QNuc::PrepCand: piF/gaF fragments are blocked"<<G4endl;
#endif
    }
    // @@ in case of the Ksi or Omega- capture it can disturb the simulation
    else if(cPDG<80000000&&(abs(cPDG)%10>4||cST>2))// @@ PreClosed HighSpin/HighStrange
    {
      curCand->SetPreProbability(0.);  
      curCand->SetDenseProbability(0.); 
      curCand->SetPossibility(false);    
    }
    else
    {
      G4double tnM=GetQPDG().GetMass();          // Total mass of this nucleus
      if(cPDG>80000000&&cPDG!=90000000)          // ===> Cluster case
      {
        G4int sc = cST;                          // "S" of the cluster
        G4int zc = curCand->GetCharge();         // "Z" of the cluster
        G4int ac = cBN;                          // "A" of the cluster
        G4int nc = ac-zc-sc;                     // "N" of the cluster
        G4double cM=tnM-G4QNucleus(Z-zc,N-nc,S-sc).GetGSMass(); // BoundMass of the cluster
        G4LorentzVector intLV=pLV+G4LorentzVector(0.,0.,0.,cM); // 4-mom of the proj+clust
        pos      = probVect[ac];                 // Cluster Probability NormalizationFactor
        if(ac<=maxClust&&pos>0.&&(pLV==zeroLV||intLV.m()>.00001+cM))
        {

#ifdef cldebug
          G4cout<<"G4QNucleus::PrepareCand: ac="<<ac<<", pV="<<pos<<G4endl;
#endif
          G4int dac=ac+ac;
          if(ac && (piF || gaF))                 // zc>=0
          {
            if     (piF&&!gaF&&zc+ac) pos*=(zc+ac)/ac;  // piF interaction (#of u-quarks)
            else if(gaF&&!piF&&zc+dac) pos*=(zc+dac)/ac; // gaF interaction (sum of Q_q^2)
          }
          G4double dense=1.;
          if     (ac==1&&pos>0.)dense=probVect[254]/pos;
          else if(ac==2&&pos>0.)dense=probVect[255]/pos;
#ifdef cldebug
          G4cout<<"G4QNucleus::PrepC: cPDG="<<cPDG<<",norm="<<pos<<",zc="<<zc<<",nc="<<nc
                <<",sc="<<sc<<",ac="<<ac<<",ze1="<<ze1<<",ne1="<<ne1<<",se1="<<se1<<G4endl;
          G4double mp=pos;
#endif
          if     (ac==1 && ae)                   // ae=0 protection (since here no /pos)
          {
            if     (zc) pos*=ze/ae;
            else if(nc) pos*=ne/ae;
            else if(sc) pos*=se/ae;
            //if     (zc) pos*=ze;
            //else if(nc) pos*=ne;
            //else if(sc) pos*=se;
            acm=1;
#ifdef cldebug
            if(pos)
            G4cout<<"G4QN::PrC:mp="<<mp<<",pos="<<pos<<",ae="<<ae
                  <<",Z="<<zc<<",N="<<nc<<",mac="<<mac<<G4endl;
            sZ+=pos*zc;
            sN+=pos*nc;
#endif
          }
          else if(ac==2)
          {
            if(ze<zc||ne<nc||se<sc) pos=0.;
            else if(aea)                         // Protection against aea=0.
            {
              if     (zc==2) pos*=ze*(ze-1)/aea;
              else if(nc==2) pos*=ne*(ne-1)/aea;
              else if(sc==2) pos*=se*(se-1)/aea;
              else if(zc==1&&nc==1) pos*=ze*ne/ae2;
              else if(zc==1&&sc==1) pos*=ze*se/ae2;
              else if(sc==1&&nc==1) pos*=se*ne/ae2;
              //if     (zc==2) pos*=ze*(ze-1)/2.;
              //else if(nc==2) pos*=ne*(ne-1)/2.;
              //else if(sc==2) pos*=se*(se-1)/2.;
              //else if(zc==1&&nc==1) pos*=ze*ne;
              //else if(zc==1&&sc==1) pos*=ze*se;
              //else if(sc==1&&nc==1) pos*=se*ne;
              else G4cout<<"***G4QNucl::PrepCand: z="<<zc<<", n="<<nc<<", s="<<sc<<G4endl;
              // Normalization for only not strange matter
            }
            acm=2;
#ifdef cldebug
            if(pos)
            G4cout<<"G4QN::PrC:mp="<<mp<<",p="<<pos<<",A=2,(Z="<<zc<<",N="<<nc<<"),m="
                  <<mac<<G4endl;
            sZ+=pos*zc;
            sN+=pos*nc;
#endif
          }
          else                                     // ac>2
          {
            if(acm<ac)                             // first time that big cluster
            {
              if(ac<ae1 && ac>0) comb*=(ae1-ac)/ac;
              acm=ac;
              s_value=0.;
              cca=0;
#ifdef cldebug
              if(ac%2) mac=7;                      // @@Change it if cluster set is changed
              else     mac=8;                      // @@ It is not yet automatic
              G4cout<<"G4QNuc::PrepCl:c="<<comb<<",ac="<<ac<<"("<<index<<"),m="<<mac<<",a="
                    <<ae0<<G4endl;
#endif
            }
            G4double prod=1.;
            if(ze1<=zc||ne1<=nc||se1<=sc) prod=0.;
            else
            {
              if(zc>0) for(int iz=1; iz<=zc; iz++) prod*=(ze1-iz)/iz;
              if(nc>0) for(int in=1; in<=nc; in++) prod*=(ne1-in)/in;
              if(sc>0) for(int is=1; is<=sc; is++) prod*=(se1-is)/is;
            }
            s_value+=prod;
            pos*=prod;
            pos/=comb;
#ifdef cldebug
            if(pos) G4cout<<"G4QN::PreC:c="<<cPDG<<",p="<<pos<<",i="<<index<<",m="<<mac
                          <<",pr="<<prod<<",c="<<cca<<G4endl;
            sZ+=pos*zc;
            sN+=pos*nc;
#endif
            cca++;
          }
          curCand->SetPreProbability(pos);
          curCand->SetDenseProbability(pos*dense);
#ifdef cldebug
          G4cout<<"G4QN::PrepC: ClusterPDG="<<cPDG<<",preProb="<<pos<<",d="<<dense<<G4endl;
#endif
        }
        else                                       // => "Cluster is too big" case
        {
          curCand->SetPreProbability(0.);     
          curCand->SetDenseProbability(0.);
          curCand->SetPossibility(false);          // This candidate is not possible 
        }
      }
      else
      {
#ifdef cldebug
        G4cout<<"G4QNucl::PrepCand:cPDG="<<cPDG<<",pos="<<pos<<G4endl;
#endif
        curCand->SetPreProbability(pos);           // ===> Hadronic case in Vacuum     
        curCand->SetDenseProbability(0.);          // ===> Hadronic case in Vacuum
      }
      curCand->SetPossibility(true);           // All candidates are possible at this point
    }
  } // End of the LOOP over Candidates
#ifdef cldebug
  G4cout<<"G4QNucl::PrepCand:covP="<<ae*sZ/ze<<",covN="<<ae*sN/ne<<",totP="<<sZ+sN<<G4endl;
  //throw G4QException("G4QNucleus::PrepareCandidate: Temporary stop");
#endif
}// End of PrepareCandidates

//Coulomb Barrier Calculation - a general (external call) function, @@ to be moved
G4double G4QNucleus::CoulombBarGen(const G4double& rZ, const G4double& rA,
                                   const G4double& cZ, const G4double& cA)
{
  static const G4double third=1./3.;
  G4double ca=cA; // @@ can be integer
  if(cA < 0.) ca=-cA;
#ifdef debug
  if(rA < 0.) G4cout<<"-Warning-G4QNucl::CoulombBarGen: NucleusA="<<rA<<", Z="<<rZ<<G4endl;
#endif
  G4double ra=rA; // @@ can be integer
  if(rA < 0.) ra=-rA;
  G4double zz=rZ*cZ;
  // Naitive CHIPS radius: CB={1.46=200(MeV)/137}*z*Z/{R=1.13}*((a*z)**1/3+A**1/3) (?)
  //G4double cb=1.29*zz/(pow(rA,third)+pow(cA,third));
  //double cb=zz/(pow(rA,third)+pow(ca,third)+.1);
  G4double cb=1.29*zz/(pow(ra,third)+pow(ca,third)+.1); //CHIPS like potential
  // Geant4 solution for protons is practically the same:
  // G4double cb=1.263*Z/(1.0 + pow(rA,third));
  // @@ --- Temporary "Lambda/Delta barrier for mesons"
  //if(!cA) cb+=40.;
  // --- End of Temporary ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifdef debug
  G4cout<<"G4QNucl::CoulBG:rA="<<cA<<",rZ="<<cZ<<",cA="<<cA<<",cZ="<<cZ<<",C="<<cb<<G4endl;
#endif
  return cb;
} // End of "CoulombBarier"

//Coulomb Barrier Calculation: (cZ,cA)= fragment from the nucleus (delZ,dA)= second pt red
G4double G4QNucleus::CoulombBarrier(const G4double& cZ, const G4double& cA, G4double delZ,
                                    G4double dA) // (delZ,dA) are 0 by default
{
  G4double rA=GetA()-cA;
  if (dA != 0.)  rA-=dA;
#ifdef debug
  if(rA<0.) G4cout<<"-Warning-G4QNucl::CoulombBarrier: NucleusA="<<rA<<", rZ="<<cZ<<G4endl;
#endif
  G4double rZ=Z-cZ;
  if(delZ != 0.) rZ-=delZ;
#ifdef debug
  if(rA<0.) G4cout<<"-Warning-G4QNucl::CoulombBarrier: NucleusA="<<rA<<", rZ="<<cZ<<G4endl;
#endif
  return CoulombBarGen(rZ, rA, cZ, cA);
} // End of "CoulombBarier"

// Fission Coulomb Barrier Calculation
G4double G4QNucleus::FissionCoulombBarrier(const G4double& cZ, const G4double& cA,
                                           G4double delZ, G4double dA)
{
  static const G4double third=1./3.;
  if(cZ<=0.) return 0.;
  G4double rA=GetA()-cA;
  if(dA) rA-=dA;                        // Reduce rA f CB is calculated for wounded nucleus
  G4double rZ=Z-cZ;
  if(delZ) rZ-=delZ;                    // Reduce rZ f CB is calculated for wounded nucleus
  G4double zz=rZ*cZ;                    // Product of charges
  G4double r=(pow(rA,third)+pow(cA,third))*(1.51+.00921*zz)/(1.+.009443*zz);
  return 1.44*zz/r;
} // End of "FissionCoulombBarier"

//Coulomb Binding Energy for the cluster
G4double G4QNucleus::BindingEnergy(const G4double& cZ, const G4double& cA, G4double delZ,
                                   G4double dA)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass(); // Mass of neutron
  static const G4double mProt= G4QPDGCode(2212).GetMass(); // Mass of proton
  if(!cZ && !cA) return Z*mProt+N*mNeut-GetGSMass();       // Total binding energy far all
  G4double GSM=GetGSMass();
  G4int iZ=static_cast<int>(cZ);
  G4int cN=static_cast<int>(cA-cZ);
  G4int sZ=iZ;
  G4int sN=cN;
  if(delZ||dA)
  {
    G4int dz=static_cast<int>(delZ);
    G4int dn=static_cast<int>(dA-delZ);
    GSM=G4QNucleus(Z-dz,N-dn,S).GetGSMass();
    sZ-=dz;
    sN-=dn;
  }
  return G4QNucleus(Z-sZ,N-sN,S).GetGSMass()+G4QNucleus(iZ,cN,0).GetGSMass()-GSM;
} // End of "BindingEnergy"

//Coulomb Barrier Reflection Probability (CB - Coulomb Barrier, E - Kinetic Energy) 
G4double G4QNucleus::CoulBarPenProb(const G4double& CB, const G4double& E,
                                    const G4int& C, const G4int& B)
{//                 = A.Lepretre et al, Nucl.Phys., A390 (1982) 221-239 =
  static const G4double mNeut= G4QPDGCode(2112).GetMass();          // Mass of neutron
  static const G4double dNeut= mNeut+mNeut;                         // DiMass of neutron
  static const G4double mProt= G4QPDGCode(2212).GetMass();          // Mass of proton
  static const G4double dProt= mProt+mProt;                         // DiMass of proton
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0); // Mass of deuteron
  //static const G4double mTrit= G4QPDGCode(2112).GetNuclMass(1,2,0); // Mass of tritium
  //static const G4double mHel3= G4QPDGCode(2112).GetNuclMass(2,1,0); // Mass of Helium 3
  //static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0); // Mass of alpha
  static const G4double wellDebth=40.;          //@@ Should be jus binding energy @@ done
  // @@ --- Temporary 1 ---> close the OverBarrierReflection for all
  //return 1.;
  // ^^^^^^^---> End of Themporary 1
  // @@ --- Temporary 2 ---> close the OverBarrierReflection for fragments and mesons
  //if(E<CB) return 0.;
  //if(B!=1) return 1.;
  if(B<1 || B>2) return 1.;
  if(C>B+1)
  {
#ifdef debug
    G4cout<<"-Warning-G4QN::CBPP:SubtractedChrg="<<C<<" >SubtractedBaryonNmbr="<<B<<G4endl;
#endif
    return 1.;
  }
  // ^^^^^^^---> End of Themporary 2
  //G4double nA=GetA();
  //G4double nA=GetA()-B;
  //if(nA==40) G4cout<<"G4QN::CBPP:Z="<<GetZ()<<",C="<<C<<",B="<<B<<G4endl;
  //else     return 1.;           // @@@@@ Over barrier reflection is closed @@@ !!! @@@
  //      Li6      C12           Al27
  //else if(nA<7||nA>8&&nA<12||nA>16&&nA<40) return 1.;// "OverBarrierReflection is closed"
  //else if(nA>8&&nA<12||nA>16&&nA<40) return 1.; // "OverBarrierReflection is closed" Cond
  //else if(nA<12||nA>16&&nA<40) return 1.; // "OverBarrierReflection is closed" Condition
  //else if(nA<12||nA>16) return 1.; // "OverBarrierReflection is closed" Condition
  //else if(nA<12) return 1.;    // @@@@@ Over barrier reflection is closed @@@ !!! @@@
  //if(B+B>Z+N+S) return 1.;
  //G4double wD=wellDebth*B;
  G4double wD=wellDebth;
  //G4double wD=0.;
  // @@ --- Temporary 3 ---> close the OverBarrierReflection for mesons
  //if(!B) wD=0.;
  // ^^^^^^^---> End of Themporary 3
  G4double GSM=GetGSMass();
#ifdef debug
  G4cout<<"G4QNucl::CBPenProb:GSM="<<GSM<<",Z="<<Z<<",N="<<N<<",C="<<C<<",B="<<B<<G4endl;
#endif
  if(2>3);
  // @@ Temporary "Mass Barrier for mesons" @@ __________________
  //else if(!B) wD=40.;
  // @@ End of Temporary^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ////else if(nA<7&&B>0)  wD=0.;    // Only Coulomb Barrier can reflect !!!
  ////else if((nA<12||nA>16)&&B>0) wD=0.;// Only CoulombB can reflect !!! O16 E-dep of gamA
  ////else if((nA<12||nA>27)&&B>0) wD=0.;// Only CoulombB can reflect !!! O16 E-dep of gamA
  ////else if(nA<9&&B>0) return 1.;// Only CoulombBarrier can reflect !!! O16 E-dep of gamA
  ////else if(B>0)  wD=0.;    // Only Coulomb Barrier can reflect !!!
  ////else if(B==1)  wD=0.;
  else if(B==1&&C==1) wD=G4QNucleus(Z-1,N,S).GetGSMass()+mProt-GSM;
  else if(B==1&&C==0) wD=G4QNucleus(Z,N-1,S).GetGSMass()+mNeut-GSM;
  ////else if(B==1&&C==0) wD=0.;
  ////else if(B>1)  return 1.;
  ////else if(B>1)  wD=0.;
  ////else if(B==2)  wD=0.;
  else if(B==2)
  {
    if       (!C) wD=G4QNucleus(Z,N-2,S).GetGSMass()+dNeut-GSM;
    else if(C==1) wD=G4QNucleus(Z-1,N-1,S).GetGSMass()+mDeut-GSM;
    else if(C==2) wD=G4QNucleus(Z-2,N,S).GetGSMass()+dProt-GSM;
    // @@ Temporary "Local B=2 Anti Virial factor" @@ __________________
    wD=80.; // 40 MeV per each nucleon
    //wD=wD/2;
    // @@ End of Temporary^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  }
  ////else if(B>2)  wD=0.;
  ////else if(B>2)  return 1.;
  ////else if(B==3)  wD=0.;
  //else if(B==3&&C==1) wD=G4QNucleus(Z-1,N-2,S).GetGSMass()+mTrit-GSM;
  ////else if(B==3&&C==1) wD=0.;
  //else if(B==3&&C==2) wD=G4QNucleus(Z-2,N-1,S).GetGSMass()+mHel3-GSM;
  ////else if(B>3)  wD=0.;
  ////else if(B==4)  wD=0.;
  //else if(B==4&&C==2) wD=G4QNucleus(Z-2,N-2,S).GetGSMass()+mAlph-GSM;
  ////else if(B>4)  wD=0.;
  ////else if(B>4)  return 1.;
  //else if(B>4)wD=G4QNucleus(Z-C,N-B+C,S).GetGSMass()+G4QNucleus(C,B-C,S).GetGSMass()-GSM;
  if(wD<0.) wD=0.;
#ifdef debug
  G4cout<<"G4QNucl::CBPenProb: wD="<<wD<<",E="<<E<<",CB="<<CB<<G4endl;
#endif
  // @@ Temporary "Virial factor" @@ __________________
  wD=wD+wD;
  // @@ End of Temporary^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  G4double sR=0.;
  G4double CBD=CB+wD;
  G4double ED=E+wD;
  if(CBD<0.) return 1.;
  if(ED<=0.)  return 0.;
  //if(nA<27) sR=sqrt(wD/ED);
  //else      sR=sqrt(CBD/ED);
  sR=sqrt(CBD/ED);
  //sR=sqrt(wD/ED);
#ifdef debug
  G4cout<<"G4QN::CBPP:s="<<sR<<",E="<<E<<",w="<<wD<<",CB="<<CB<<",B="<<B<<",C="<<C<<G4endl;
#endif
  if(sR>=1.) return 0.;
  return   1.-sR*sR*sR;
} // End of "CoulBarPenProb"

// Randomize 2D-vector b = 2D impact parameter
pair<G4double, G4double> G4QNucleus::ChooseImpactXandY(G4double maxImpact)
{
  G4double x=1.;
  G4double y=1.;
  do
  {
    x = G4UniformRand();
    x = x+x-1.;
    y = G4UniformRand();
    y = y+y-1.;
  }
  while(x*x+y*y > 1.);
  std::pair<G4double, G4double> theImpactParameter;
  theImpactParameter.first  = x*maxImpact;
  theImpactParameter.second = y*maxImpact;
  return theImpactParameter;
} // End of "ChooseImpactXandY"

// Initializes 3D Nucleons in the nucleus using random sequence
void G4QNucleus::ChooseNucleons()
{
#ifdef debug
  G4cout<<"G4QNucleus::ChooseNucleons: Nucleons search is started"<<rho0<<G4endl;
#endif
  G4int protons =0;
  G4int nucleons=0;
  G4int theA=GetA();
  while (nucleons < theA)
  {
    if(protons<Z && G4UniformRand() < G4double(Z-protons)/G4double(theA-nucleons) )
    {
      protons++;
      nucleons++;
      G4QHadron* proton = new G4QHadron(2212);
      theNucleons.push_back(proton);
    }
    else if ( (nucleons-protons) < N )
    {
      nucleons++;
      G4QHadron* neutron = new G4QHadron(2112);
      theNucleons.push_back(neutron);
    }
    else G4cout<<"G4QNucleus::ChooseNucleons not efficient"<<G4endl;
  }
  return;
} // End of ChooseNucleons

// Initializes positions of 3D nucleons (@@ in QGS only 2D impact par positions are needed)
void G4QNucleus::ChoosePositions()
{
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mProt2= mProt*mProt;
  static const G4double third= 1./3.;
#ifdef debug
  G4cout<<"G4QNucl::ChoosePositions: is called"<<G4endl;
#endif
  G4int          i=0;                         // nucleon index
  G4int          theA=GetA();                 // cashed value of A
  G4int          lastN=theA-1;                // cashed value of A-1 (theLastNucleon index)
  G4ThreeVector  aPos(0.,0.,0.);              // Prototype of the nucleon position
  G4double       rPos=0.;                     // Radius of the nucleon position
  G4ThreeVector  delta(0.,0.,0.);             // Prototype of the distance between nucleons
  G4ThreeVector* places= new G4ThreeVector[theA]; // Vector of 3D positions
  G4bool         freeplace= false;            // flag of free space available
  G4double nucDist2= nucleonDistance*nucleonDistance; // @@ can be a common static
  G4double maxR= GetRadius(0.01);             // there are cond no nucleons at this density
  G4ThreeVector  sumPos(0.,0.,0.);            // Vector of the current 3D sum
  G4ThreeVector  minPos(0.,0.,0.);            // Minimum nucleon 3D position
  G4double       mirPos=maxR;                 // Proto MinimumRadius of the nucleonPosition
  G4int failCNT=0;                            // Counter of bad attempts
  G4int maxCNT=27;                            // Limit for the bad attempts
  while( i < theA && maxR > 0.)               // LOOP over all nucleons
  {
    rPos = maxR*pow(G4UniformRand(),third);   // Get random radius of the nucleon position
    G4double density=rPos*rPos;               // Density at R (temporary squared radius)
    if(theA<17) density=GetRelOMDensity(density); // Oscilator model (M.K.?)
    else        density=GetRelWSDensity(rPos);    // Wood-Saxon model
#ifdef debug
    G4cout<<"G4QNucl::ChoosePositions: i="<<i<<", pos="<<aPos<<", dens="<<density<<G4endl;
#endif
    if(G4UniformRand()<density)               // Try this position with frequency ~Density
    {
      // @@ Gaussian oscilator distribution is good only up to He4 (s-wave). Above: p-wave
      // (1+k*(r^2/R^2)]*exp[-r^2/R^2]. A=s+p=4+3*4=16 (M.K.) So Li,Be,C,N,O are wrong
      if(i==lastN) aPos=-rPos*sumPos.unit();  // TheLast tries toCompensate CenterOfGravity
      else     aPos=rPos*G4RandomDirection(); // It uses the standard G4 function
      freeplace = true;                       // Imply that there is a free space
      for(G4int j=0; j<i && freeplace; j++)   // Check that there is no overlap with others
      {
        delta = places[j] - aPos;             // Distance to nucleon j
        freeplace= delta.mag2()>nucDist2;     // If false break the LOOP
      }
      //  protons must at least have binding energy of CoulombBarrier (@@ ? M.K.), so
      //  assuming Fermi Energy corresponds to Potential, we must place protons such
      //  that the Fermi Energy > CoulombBarrier (?)
      //  @@ M.K.: 1. CoulBar depends on aPos; 2. Makes Isotopic assymetry (!); 3. Perform.
      G4int nucPDG= theNucleons[i]->GetPDGCode();
#ifdef debug
      G4cout<<"G4QNucl::ChoosePositions: frpl="<<freeplace<<", nucPDG="<<nucPDG<<G4endl;
#endif      
      if(freeplace && nucPDG == 2212) // Free Space Protons
      {
        G4double pFermi=GetFermiMomentum(GetDensity(aPos));
        G4double eFermi= sqrt(pFermi*pFermi+mProt2)-mProt;  // Kinetic energy
        if (eFermi <= CoulombBarrier()) freeplace=false;
      }
      if(rPos<mirPos)
      {
        mirPos=rPos;
        minPos=aPos;
      }
      if( freeplace || failCNT > maxCNT )
      {
        if( failCNT > maxCNT ) aPos=minPos;
#ifdef debug
        G4cout<<"G4QNuc::ChoosePos:->> fill N["<<i<<"], R="<<aPos<<", f="<<failCNT<<G4endl;
#endif      
        places[i]=aPos;
        sumPos+=aPos;
        ++i;
        failCNT=0;
        mirPos=maxR;
      }
      else ++failCNT;
    }
  }
#ifdef debug
  G4cout<<"G4QNucl::ChoosePositions: Out of the positioning LOOP"<<G4endl;
#endif
  if(theA > 2) ReduceSum(places,sumPos);              // Reduce the CM shift (equal weights)
#ifdef debug
  G4cout<<"G4QNucl::ChoosePositions: The reduced summ is made"<<G4endl;
#endif
  for(i=0; i<theA; i++) theNucleons[i]->SetPosition(places[i]);
  delete [] places;
#ifdef debug
  G4cout << "G4QNucleus::ChoosePositions: The positions are defined for A="<<theA<<G4endl;
#endif
} // End of ChoosePositions

// Initializes density of 3D nucleons in the 3D nucleus
void G4QNucleus::InitDensity()
{
  static const G4double r0sq=0.8133*fermi*fermi;      // Base for A-dep of rel.mean.radius
  static const G4double third=1./3.;
  G4int    iA = GetA();
  G4double rA = iA;
  G4double At = pow(rA,third);
  G4double At2= At*At;
#ifdef debug
  G4cout<<"G4QNucleus::InitDensity: rA=iA=A="<<iA<<", A^1/3="<<At<<", A^2/3="<<At2<<G4endl;
#endif
  if(iA<17)                                           // Gaussian density distribution
  {
    radius = r0sq*At2;                                // R2 Mean Squared Radius (fm^2)
    if(radius<=0.)
    {
      G4cout<<"-Warning-G4QNucl::ChoosePositions:L,iA="<<iA<<",Radius(?)="<<radius<<G4endl;
      radius=1.;
    }
    rho0   = pow(2*pi*radius, -1.5);                  // Central Density (M.K. 2 is added)
    // V=4pi*R2*sqrt(pi*R2/2)=(sqrt(2*pi*R2))^3
  }
  else                                                // Wood-Saxon density distribution
  {
    G4double r0=1.16*(1.-1.16/At2)*fermi;             // Base for A-dependent radius
    radius = r0*At;                                   // Half Density Radius (fm)
    if(radius<=0.)
    {
      G4cout<<"-Warning-G4QNucl::ChoosePositions:H,iA="<<iA<<",Radius(?)="<<radius<<G4endl;
      radius=1.;
    }
    G4double rd=WoodSaxonSurf/radius;                 // Relative thickness of the surface
    if(!(rd<=0.1) && !(rd>-0.1))                      // NAN for rd
    {
      G4cout<<"-Warning-G4QNucl::ChoosePositions:H,NAN,iA="<<iA<<", rd="<<rd<<G4endl;
      rd=1.;
    }
    rho0=0.75/(pi*pow(radius,3)*(1.+rd*rd*pi2));      // Central Density
  }
  RhoActive=true;
} // End of InitDensity

// Calculates Derivity of the nuclear density
G4double G4QNucleus::GetDeriv(const G4ThreeVector& aPosition)
{
  if(radius==0.) InitDensity();
  G4double rPos=aPosition.mag();
  if(GetA()<17) return -GetDensity(aPosition)*(rPos+rPos)/radius; // Gaussian density
  // Wood-Saxon density distribution
  G4double dens=GetRelativeDensity(aPosition);
  return -exp((rPos-radius)/WoodSaxonSurf)*dens*dens*rho0/WoodSaxonSurf;
} // End of GetDeriv

// Radius of the deffinit % of nuclear density (very strange solution ?? @@ M.K.)
G4double G4QNucleus::GetRadius(const G4double maxRelDens)
{
  if(radius==0.) InitDensity();
  if(GetA()<17)                                     // Gaussian density distribution
    return (maxRelDens>0 && maxRelDens <= 1. ) ? sqrt(-radius*log(maxRelDens) ) : DBL_MAX;
  // Wood-Saxon density distribution
  return (maxRelDens>0 && maxRelDens <= 1. ) ? (radius + WoodSaxonSurf*
         log((1.-maxRelDens+exp(-radius/WoodSaxonSurf))/maxRelDens) ) : DBL_MAX;
} // End of GetRadius (check @@ radius=sqr0 (fm^2) for A<17,  r0 (fm) for A>16 (units)

// Calculates Densyty/rho0
G4double G4QNucleus::GetRelativeDensity(const G4ThreeVector& aPosition)
{
  if(radius==0.) InitDensity();
  if(GetA()<17) return GetRelOMDensity(aPosition.mag2());// Gaussian distribution (OscMod?)
  return GetRelWSDensity(aPosition.mag());               // Wood-Saxon density distribution
} // End of GetRelativeDensity

// Calculates modul of Fermy Momentum  depending on density
G4double G4QNucleus::GetFermiMomentum(G4double density)
{
  static const G4double third=1./3.;
  static const G4double constofpmax=hbarc*pow(3.*pi2,third);
  return constofpmax * pow(density*GetA(),third);
} // End of GetFermiMomentum

// Calculates 3D Fermy Momentuma for 3D nucleons in 3D nucleus
void G4QNucleus::ChooseFermiMomenta()
{
  static const G4double mProt= G4QPDGCode(2212).GetMass(); // Mass of proton
  static const G4double mProt2= mProt*mProt;
  //static const G4double mNeut= G4QPDGCode(2112).GetMass(); // Mass of neutron
  static const G4double third= 1./3.;
  G4int i=0;
  G4double density=0.;                                // Prototype of density for Loop calc
  G4int theA=GetA();                                  // Atomic weight of the nucleus
  G4int am1=theA-1;                                   // The last index in the Loop
  G4ThreeVector* momentum = new G4ThreeVector[theA];  // Temporary array for nucleon's moms
  G4ThreeVector  sumMom(0.,0.,0.);                    // Sum of all momenta for mom-conserv
#ifdef debug
  G4cout<<"G4QNucleus::ChooseFermiMomentum is called for Z="<<Z<<", N="<<N<<G4endl;
#endif
  for(i=0; i<theA; i++) // momenta for all, including the last, in case we swap nucleons
  {
    density=GetDensity(theNucleons[i]->GetPosition());// density around nucleon i
    G4double ferm = GetFermiMomentum(density);        // module of momentum for nucleon i
    G4ThreeVector mom(0.,0.,0.);                      // proto 3vector for nucleon momentum
    G4double rn3=pow(G4UniformRand(),third);          // Spherical randomization
    G4ThreeVector dir(0.,0.,0.);                      // proto 3vector for the momDirection
    if( i == am1) dir=-sumMom.unit();                 // try to compensate the mom noncons.
    else          dir=G4RandomDirection();            // free randomization for i < A-1
    if(theNucleons[i]->GetPDGCode() == 2212)          // the nucleon is a proton
    {
      G4double eMax = sqrt(ferm*ferm+mProt2)-CoulombBarrier();
      if(eMax>mProt) mom=sqrt(eMax*eMax - mProt2)*rn3*dir; // 3D proton momentum
#ifdef debug
      else G4cerr<<"-Warning-G4QNucleus::ChooseFermM: FailToGetProtonMomentum,p=0"<<G4endl;
#endif
    }
    else mom=ferm*rn3*dir;                            // 3-vector for the neutron momentum
    momentum[i]= mom;
    sumMom+= mom;
#ifdef debug
    G4cout<<"G4QNucleus::ChooseFermiMomentum: for i="<<i<<", candidate mom="<<mom<<G4endl;
#endif
  }
  if(theA > 2) SimpleSumReduction(momentum, sumMom);  // Reduse momentum nonconservation
  //G4double bindEn=BindingEnergy()/theA;
  G4int thisPDG=GetPDG();
  G4double rMp=G4QPDGCode(thisPDG-1000).GetMass();    // Residual for the proton
  G4double rMn=G4QPDGCode(thisPDG-1).GetMass();       // Residual for the neutron
  G4double rMp2=rMp*rMp;
  G4double rMn2=rMn*rMn;
  //G4double rM=rMn;
  G4double rM2=rMn2;
  G4double thisM=GetGSMass();
#ifdef debug
  G4LorentzVector sum(0.,0.,0.,0.);
#endif
  for(i=0; i< theA ; i++ )
  {
    if(theNucleons[i]->GetPDGCode() == 2212)
    {
      //rM=rMp;
      rM2=rMp2;
    }
    else
    {
      //rM=rMn;
      rM2=rMn2;
    }
    G4ThreeVector curMom = momentum[i];
    G4double energy = thisM-std::sqrt(rM2+curMom.mag2()); // @@ update after splitting
    G4LorentzVector tempV(curMom,energy);
#ifdef debug
    G4cout<<"G4QNucleus::ChooseFermiMomentum: FINALLY for i="<<i<<", 4mom="<<tempV<<G4endl;
    sum+=tempV;
#endif
    theNucleons[i]->Set4Momentum(tempV);
  }
#ifdef debug
    G4cout<<"G4QNucleus::ChooseFermiMomentum: FINALLY sum4M="<<sum<<G4endl;
#endif
  delete [] momentum;
} // End of ChooseFermiMomenta

// Reduce momentum nonconservation or center of mass shift (Changes the momena!)
void G4QNucleus::SimpleSumReduction(G4ThreeVector* vect, G4ThreeVector sum)
{
  G4int theA=GetA();                                // A#of nucleons
  sum/=theA;
  for(G4int i=0; i<theA; i++) vect[i]-=sum;        // Simple reduction
}

// Reduce momentum nonconservation or center of mass shift (Keep values of momena) @@Bug!@@
G4bool G4QNucleus::ReduceSum(G4ThreeVector* vect, G4ThreeVector sum)
{
  G4int theA=GetA();                                // A#of nucleons
  if(theA<3)                                        // Can not reduce for 1 or 2 nucleons
  {
    G4cout<<"-Warning-G4QNucleus::ReduceSum: *Failed* A="<<theA<<" < 3"<<G4endl;
    return false;
  }
  // The last vector must have the same direction as the SUM (do not take into account
  G4int am1=theA-1;                                 // A-1 elements, which canBeCorrected
  G4double  sum2=sum.mag2();                        // Initial squared sum
  G4double  hsum2=sum2/2;                           // Half squared sum
  G4double*  dp= new G4double[am1];                 // Displacements
  G4int     m_value=am1;                            // #0fVectors used for correction
  G4double  minS=DBL_MAX;                           // Min value of Fermi Momentum
  G4int     minI=0;                                 // Index of maximum Fermi Momentum
  for(G4int i=0; i<am1; i++) dp[i]=sum.dot(vect[i]);// Calculation of dot-products
  while(m_value)
  {
    m_value=0;
    for(G4int i=0; i<am1; i++) if(dp[i]>0 && dp[i]<sum2) // can be used for the reduction
    {
      m_value++;
      G4double shift=fabs(dp[i]-hsum2);
      if(shift < minS)
      {
        minS=shift;
        minI=i;
      }
    }
    if(m_value)                                      // There is a vector reducing the sum
    {
      G4ThreeVector x=(dp[minI]/hsum2)*sum;          // turn-reduction of the sum-vector
      vect[minI]-=x;                                 // turn the minI-th vector
      sum-=x;                                        // reduce the sum
      sum2=sum.mag2();                               // Current squared sum
      hsum2=sum2/2;                                  // Current half squared sum
    }
  }
  if(sum2 > 0.)
  {
    sum/=theA;
    for(G4int i=0; i<theA; i++) vect[i]-=sum;        // Final reduction
  }
  delete[] dp;
  return true;
} // End of ReduceSum

// Initializes 3D nucleus with (Pos,4M)-Nucleons. It automatically starts the LOOP
void G4QNucleus::Init3D()
{
#ifdef debug
  G4cout<<"G4QNucleus::Init3D: is called currentNucleon="<<currentNucleon<<G4endl;
#endif
  for_each(theNucleons.begin(),theNucleons.end(),DeleteQHadron());
  theNucleons.clear();
  G4int theA = GetA();
  ChooseNucleons();
#ifdef debug
  G4cout<<"G4QNucleus::Init3D: Nucleons are initialized, nN="<<theNucleons.size()<<G4endl;
#endif
  InitDensity();
#ifdef debug
  G4cout<<"G4QNucl::Init3D: DensityPars for A="<<theA<<":R="<<radius <<",r0="<<rho0<<G4endl;
#endif
  ChoosePositions(); // CMS positions! No Lorentz boost! Use properely!
#ifdef debug
  G4cout<<"G4QNucleus::Init3D: Nucleons are positioned in the coordinate space"<<G4endl;
#endif
  ChooseFermiMomenta(); // CMS Fermi Momenta! Must be transfered to the LS if not at rest!
  G4ThreeVector n3M=Get3Momentum();                   // Velocity of the nucleus in LS
  if(n3M.x() || n3M.y() || n3M.z())                   // Boost the nucleons to LS
  {
    n3M/=GetEnergy();                                 // Now this is the boost velocity
    DoLorentzBoost(n3M);                              // Now nucleons are in LS
  }
#ifdef debug
  G4cout<<"G4QNucleus::Init3D: Nucleons are positioned in the momentum space"<<G4endl;
#endif
  G4double Ebind= BindingEnergy()/theA;
  for (G4int i=0; i<theA; i++) theNucleons[i]->SetBindingEnergy(Ebind); // @@ ? M.K.
  currentNucleon=0;                                   // Automatically starts the LOOP
  return;
} // End of Init3D

// Get radius of the most far nucleon (+ nucleon radius)
G4double G4QNucleus::GetOuterRadius()
{
  G4double maxradius2=0;
  G4int theA=theNucleons.size();
  if(theA) for(G4int i=0; i<theA; i++)
  {
    G4double nucr2=theNucleons[i]->GetPosition().mag2();
    if(nucr2 > maxradius2) maxradius2=nucr2;
  }
  return sqrt(maxradius2)+nucleonDistance;
} // End of GetOuterRadius

//
void G4QNucleus::DoLorentzContraction(const G4ThreeVector& theBeta)
{
  G4double bet2=theBeta.mag2();
  G4double factor=(1.-sqrt(1.-bet2))/bet2;         // 1./(beta2*gamma2)
  G4int theA=theNucleons.size();
  if(theA) for (G4int i=0; i< theA; i++)
  {
    G4ThreeVector pos=theNucleons[i]->GetPosition(); 
    pos -= factor*(theBeta*pos)*theBeta;  
    theNucleons[i]->SetPosition(pos);
  }    
} // End of DoLorentzContraction(G4ThreeVector)

// Shift all nucleons of the 3D Nucleus (Used in GHAD-TFT)
void G4QNucleus::DoTranslation(const G4ThreeVector& theShift)
{
  G4int theA=theNucleons.size();
  if(theA) for(G4int i=0; i<theA; i++)
    theNucleons[i]->SetPosition(theNucleons[i]->GetPosition() + theShift);
} // End of DoTranslation

// Initializes currentNucleon=0 returns size of theNucleons vector
G4bool G4QNucleus::StartLoop()
{
  G4int theA=theNucleons.size();
  if(theA) currentNucleon=0;
  else G4cout<<"-Warning-G4QNucleus::StartLoop: LOOP starts for uninited nucleons"<<G4endl;
  return theA;
} // End of StartLoop

//Calculate T(b) with step .1 fm
void G4QNucleus::ActivateBThickness()
{
  static const G4double aT= .0008;          // pred exponent parameter
  static const G4double sT= .42;            // slope parameter
  static const G4double pT=-.26;            // power parameter
  static const G4double db= .1;             // step in b (fm)
  // @@ make better approximation for light nuclei
  G4double A = GetA();                      // atomic weight
  G4double B = aT*A*A;                      // predexponent (no units)
  G4double D = sT*std::pow(A,pT);           // b^2 slope (fm^-2)
  G4double C = A*D/pi/std::log(1.+B);       // Norm for plane density (fm^-2)
  G4double mT= C*B/(1+B);                   // Max (b=0) b-thickness
  G4double T = mT;                          // Current b-thickness
  mT/=1000.;                                // Min b-thickness (@@ make 1000 a parameter)
  G4double b = 0.;
  while(T>mT)
  {
    //Tb->push_back(T);                      // Fill the thickness vector starting with b=0
    Tb.push_back(T);                        // Fill the thickness vector starting with b=0
    b+=db;                                  // increment impact parameter
    G4double E=B*std::exp(-D*b*b);          // b-dependent factor
    T=C*E/(1.+E);                           // T(b) in fm^-2
  }
  TbActive=true;                            // Flag of activation
} // End of "ActivateBThickness"

// Calculate the integral of T(b)
G4double G4QNucleus::GetTbIntegral() // Calculate the integral of T(b)
{
  if(!TbActive) ActivateBThickness();
  G4int nt = Tb.size();
  G4double sum=0.;
  for(G4int i=0; i<nt; ++i) sum+=i*Tb[i];
  sum*=.02*pi;
#ifdef debug
  G4cout<<"G4QNucleus::GetTbIntegral:TI="<<sum<<", RI="<<4*pi*rho0*pow(radius,3)/3<<G4endl;
#endif
  return sum;
}

// Calculates T(b)
G4double G4QNucleus::GetBThickness(G4double b)
{
  static const G4double dfermi=fermi/10.;
  static const G4double sfermi=fermi*fermi;
  if(!TbActive) ActivateBThickness();
  G4double bf = b/dfermi;
  G4int nb = static_cast<int>(bf);
  G4int eb = nb+1;
  G4int nt = Tb.size();
  if(eb>=nt) return 0.;
  G4double nT=Tb[nb];
  G4double eT=Tb[eb];
  return (nT-(bf-nb)*(nT-eT))/sfermi; // Independent units
}

// Calculates T(b)/rho0
G4double G4QNucleus::GetThickness(G4double b)
{
  G4int tA=GetA();
  if(tA<1)
  {
    G4cout<<"-Warning-G4QNucleus::GetThickness: for A="<<tA<<", => return 0"<<G4endl;
    return 0.;
  }
  else if(tA==1) return 0.;
  if(!TbActive) ActivateBThickness();
  if(!RhoActive) InitDensity();
  return GetBThickness(b)/rho0/tA;
}

// Add Cluster
G4QNucleus G4QNucleus::operator+=(const G4QNucleus& rhs)
{
  Z+=rhs.Z;
  N+=rhs.N;
  S+=rhs.S;
  dZ+=rhs.dZ;
  dN+=rhs.dN;
  dS+=rhs.dS;
  // Atributes of aHadron
  G4int           newPDG= GetPDGCode() + rhs.GetPDGCode() - 90000000;
  SetQPDG        (newPDG);
  G4QContent      newQC = GetQC()      + rhs.GetQC();
  SetQC          (newQC);
  theMomentum += rhs.Get4Momentum();
  return *this;
} 

// Subtract Cluster
G4QNucleus G4QNucleus::operator-=(const G4QNucleus& rhs)
{
  Z-=rhs.Z;
  N-=rhs.N;
  S-=rhs.S;
  dZ-=rhs.dZ;
  dN-=rhs.dN;
  dS-=rhs.dS;
  // Atributes of aHadron
  G4int           newPDG= GetPDGCode()   - rhs.GetPDGCode() + 90000000;
  SetQPDG        (newPDG);
  G4QContent      newQC = GetQC()        - rhs.GetQC();
  SetQC          (newQC);
  theMomentum -= rhs.Get4Momentum();
  return *this;
} 

// Multiply Nucleus by integer value
G4QNucleus G4QNucleus::operator*=(const G4int& rhs)
{
  Z*=rhs;
  N*=rhs;
  S*=rhs;
  dZ*=rhs;
  dN*=rhs;
  dS*=rhs;
  // Atributes of aHadron
  G4int           newPDG= rhs*(GetPDGCode() - 90000000) + 90000000;
  SetQPDG        (newPDG);
  G4QContent      newQC = rhs*GetQC();
  SetQC          (newQC);
  theMomentum *= rhs;
  return *this;
} 

// Converts hadronic PDG Code to nuclear PDG Code
G4int G4QNucleus::HadrToNucPDG(G4int hPDG)
{
  G4int  nPDG=hPDG;
  if     (hPDG==2212) nPDG=90001000; // p
  else if(hPDG==2112) nPDG=90000001; // n
  else if(hPDG==3122||hPDG==3212) nPDG=91000000; // Lambda
  else if(hPDG== 211) nPDG=90000999; // pi+
  else if(hPDG==-211) nPDG=89999001; // pi-
  else if(hPDG== 311) nPDG=89000001; // K0 (anti-strange)
  else if(hPDG== 321) nPDG=89001000; // K+ (anti-strange)
  else if(hPDG==-311) nPDG=90999999; // anti-K0 (strange)
  else if(hPDG==-321) nPDG=90999000; // K-      (strange)
  else if(hPDG==1114) nPDG=89999002; // Delta-
  else if(hPDG==2224) nPDG=90001999; // Delta++
  else if(hPDG==3112) nPDG=90999000; // Sigma-
  else if(hPDG==3222) nPDG=91000999; // Sigma+
  else if(hPDG==3312) nPDG=91999000; // Ksi-
  else if(hPDG==3322) nPDG=91999999; // Ksi0
  else if(hPDG==3334) nPDG=92998999; // Omega-
  else if(hPDG==-2212) nPDG=8999000; // anti-p
  else if(hPDG==-2112) nPDG=8999999; // anti-n
  else if(hPDG==-3122||hPDG==3212) nPDG=89000000; //anti- Lambda
  else if(hPDG==-3112) nPDG=89000999; // anti-Sigma-
  else if(hPDG==-3222) nPDG=88999001; // anti-Sigma+
  else if(hPDG==-3312) nPDG=88001000; // anti-Ksi-
  else if(hPDG==-3322) nPDG=88000001; // anti-Ksi0
  else if(hPDG==-3334) nPDG=87001001; // anti-Omega-
  return nPDG;
} 

// Converts nuclear PDG Code to hadronic PDG Code
G4int G4QNucleus::NucToHadrPDG(G4int nPDG)
{
  G4int  hPDG=nPDG;
  if     (nPDG==90001000) hPDG=2212; // p
  else if(nPDG==90000001) hPDG=2112; // n
  else if(nPDG==91000000) hPDG=3122; // Lambda
  else if(nPDG==90000999) hPDG= 211; // pi+
  else if(nPDG==89999001) hPDG=-211; // pi-
  else if(nPDG==89001000) hPDG= 213; // K0 (anti-strange)
  else if(nPDG==89000001) hPDG= 213; // K+ (anti-strange)
  else if(nPDG==90999000) hPDG=-213; // anti-K0 (strange)
  else if(nPDG==90999999) hPDG=-213; // K-      (strange)
  else if(nPDG==90001999) hPDG=1114; // Delta-
  else if(nPDG==89999002) hPDG=2224; // Delta++
  else if(nPDG==91000999) hPDG=3112; // Sigma-
  else if(nPDG==90999001) hPDG=3222; // Sigma+
  else if(nPDG==91999999) hPDG=3312; // Ksi-
  else if(nPDG==91999000) hPDG=3322; // Ksi0
  else if(nPDG==92998999) hPDG=3334; // Omega-
  return hPDG;
} 

//Evaporate Residual Nucleus
void G4QNucleus::EvaporateNucleus(G4QHadron* qH, G4QHadronVector* evaHV)
{
  static const G4double mHel6 = G4QPDGCode(2112).GetNuclMass(2,4,0);
  static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mDeut = G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mNeut = G4QPDGCode(2112).GetMass();
  static const G4double mProt = G4QPDGCode(2212).GetMass();
  static const G4double mLamb = G4QPDGCode(3122).GetMass();
  static const G4double mPi   = G4QPDGCode(211).GetMass();
  static const G4double mPi0  = G4QPDGCode(111).GetMass();
  static const G4double mK    = G4QPDGCode(321).GetMass();
  static const G4double mK0   = G4QPDGCode(311).GetMass();
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QContent deutQC(3,3,0,0,0,0);
  static const G4QContent alphQC(6,6,0,0,0,0);
  G4int      thePDG = qH->GetPDGCode();     // Get PDG code of the Residual Nucleus
  G4QContent theQC  = qH->GetQC();          // Quark Content of the hadron
#ifdef pdebug
  G4cout<<"G4QNucleus::EvaporateNucleus:-Called-PDG="<<thePDG<<",QC="<<theQC<<G4endl;
#endif
  G4int      theBN  = qH->GetBaryonNumber();// Baryon number of the nucleus
#ifdef pdebug
  G4cout<<"G4QNucleus::EvaporateNucleus: theBN="<<theBN<<G4endl;
#endif
  if((thePDG || thePDG==10) && theQC.GetBaryonNumber()>0) thePDG=theQC.GetZNSPDGCode();
  G4LorentzVector q4M = qH->Get4Momentum(); // Get 4-momentum of theTotalNucleus
  if(!theBN || thePDG<80000000 || thePDG==90000000) // Hadron, anti-nucleous, or vacuum
  {
#ifdef debug
    G4cout<<"G4QNucleus::EvaporateNucleus: Nucleus="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
    if(thePDG==90000000)
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (1) qH=0"<<G4endl;
#endif
      delete qH;
      G4cout<<"-Warning-G4QNucleus::EvapNuc:vacuum,4Mom="<<q4M<<G4endl;
    }
    else   // Put input to the output (delete equivalent)
    {
      G4cout<<"-Warning-G4QNucleus::EvapNuc:vacuum, Called for Meson PDG="<<thePDG<<G4endl;
      evaHV->push_back(qH);
    }
    return;
  }
  /// @@@@@@@ *** TEMPORARY TO AVOID HYPERMUCLEI FOR GEANT4 *** @@@@@@@
  if(thePDG>91000000) //@@MadeForGeant4@@: If there is a Lambda, substitute it by A neutron
  {
    G4int SSS=(thePDG-90000000)/1000000;      // A # of Lambdas
    thePDG-=SSS*999999;                       // S Neutrons instead of S Lambdas
    qH->SetQPDG(G4QPDGCode(thePDG));
    theQC  = qH->GetQC();          // Quark Content of the hadron
#ifdef debug
    G4cout<<"=>Hyper Change=>G4QNucleus::EvaporateNuceus: NewNucPDG="<<thePDG<<G4endl;
#endif
  }
  /// @@@ *** ^^^ END OF TEMPORARY ^^^ *** @@@
  if(thePDG<80000000)
  {
#ifdef debug
    G4cout<<"G4QN::EvaporateNuc: FundamentalParticle="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
    evaHV->push_back(qH);     // TheFundamentalParticles must be FilledAsTheyAre (del.eq)
    return;
  }
  G4int    theC=theQC.GetCharge();         // P
  G4int    theS=theQC.GetStrangeness();    // S
  G4int    theN=theBN-theC-theS;           // N
  G4double totGSM = G4QNucleus(thePDG).GetGSMass();// TheGroundStateMass of theTotalNucleus
  G4double totMass = q4M.m();              // Get the Real(Excited?)Mass of theTotalNucleus
#ifdef debug
  G4cout<<"G4QNucleus::EvaporateNucleus(EVA):===IN==> PDG="<<thePDG<<",4Mom="<<q4M<<", B="
        <<theBN<<", Z="<<theC<<", N="<<theN<<", S="<<theS<<G4endl;
#endif
  if(theBN<-2)
  {
    G4cout<<"-Warning-G4QNuc::EvapNuc: Evapor of anti-nuclei is not implemented yet PDG="
          <<thePDG<<G4endl;
    evaHV->push_back(qH);
    return;
  }
  else if(thePDG==91000000||thePDG==90001000||thePDG==90000001) // Excited Lambda* or N*
  //else if(2>3)// One can easily close this decay as it will be done later (time of calc?)
  {
    G4double gsM=mNeut;
    if(thePDG==90001000)      gsM=mProt;
    else if(thePDG==91000000) gsM=mLamb;
    if(fabs(totMass-gsM)<.001)
    {
#ifdef debug
      G4cout<<"G4QNu::EvaporateNucl:GSM="<<gsM<<", H="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
      evaHV->push_back(qH); // (delete equivalent)
      return;
    }
    else if(totMass<gsM)
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (2) qH=0"<<G4endl;
#endif
      delete qH;
      // G4cerr<<"***G4QN::EvaNuc:Baryon "<<thePDG<<" is belowMassShell, M="<<totMass<<G4endl;
      // throw G4QException("G4QNucleus::EvaporateNucleus: Baryon is below Mass Shell");
      G4ExceptionDescription ed;
      ed << "Baryon is below Mass Shell: Baryon " << thePDG
         << " is belowMassShell, M=" << totMass << G4endl;
      G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0000",
                  FatalException, ed);
    }
    else                                 // Decay in gamma or charged pion (@@ neutral)
    {
      G4double d=totMass-gsM;
#ifdef debug
      G4cout<<"G4QN::EvaporNucl: PDG="<<thePDG<<",M="<<totMass<<">"<<gsM<<",d="<<d<<G4endl;
#endif
      G4int decPDG=22;
      G4double decM=0.;
      if(d>142.)                           // @@ to avoid more precise calculations
      {
        if(thePDG==90001000)               // p* -> n + pi+
        {
          gsM=mNeut;
          thePDG=90000001;
          decPDG=211;
          decM=mPi;
        }
        else if(thePDG==90000001)          // n* -> p + pi-
        {
          gsM=mProt;
          thePDG=90001000;
          decPDG=-211;
          decM=mPi;
        }
        else                               // decay in Pi0
        {
          decPDG=111;
          decM=mPi0;
        }
      }
      G4LorentzVector h4Mom(0.,0.,0.,gsM); // GSMass must be updated in previous while-LOOP
      G4LorentzVector g4Mom(0.,0.,0.,decM);
      if(!G4QHadron(q4M).DecayIn2(h4Mom, g4Mom))
      {
#ifdef qdebug
        if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (3) qH=0"<<G4endl;
#endif
        delete qH;
        // G4cerr<<"***G4QNuc::EvaNuc:h="<<thePDG<<"(GSM="<<gsM<<")+gam>tM="<<totMass<<G4endl;
        // throw G4QException("G4QNucleus::EvaporateNucleus:BaryonDecay In Baryon+Gam Error");
        G4ExceptionDescription ed;
        ed << "BaryonDecay In Baryon+Gam Error: h=" << thePDG << "(GSM="
           << gsM << ")+gam>tM=" << totMass << G4endl;
        G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0001",
                    FatalException, ed);
      }
#ifdef debug
      G4cout<<"G4QNucl::EvaNuc:"<<totMass<<q4M<<"->"<<thePDG<<h4Mom<<"+g="<<g4Mom<<",n="
            <<evaHV->size()<<G4endl;
#endif
      G4QHadron* curH = new G4QHadron(thePDG,h4Mom);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: Hadr="<<thePDG<<h4Mom<<G4endl;
#endif
      evaHV->push_back(curH);         // Fill Baryon (delete equiv.)
      G4QHadron* curG = new G4QHadron(decPDG,g4Mom);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: Gamma(pion)4M="<<g4Mom<<G4endl;
#endif
      evaHV->push_back(curG);         // Fill gamma/pion (delete equivalent)
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (4) qH=0"<<G4endl;
#endif
      delete qH;
    }
  }
  else if(thePDG==89000000||thePDG==89999000||thePDG==89999999) // anti-Lambda* or anti-N*
  //else if(2>3)// One can easily close this decay as it will be done later (time of calc?)
  {
    G4double gsM=mNeut;
    if(thePDG==89999000)      gsM=mProt;
    else if(thePDG==89000000) gsM=mLamb;
    if(fabs(totMass-gsM)<.001)
    {
#ifdef debug
      G4cout<<"G4QNu::EvaNucl:(aB*),GSM="<<gsM<<", H="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
      evaHV->push_back(qH); // (delete equivalent)
      return;
    }
    else if(totMass<gsM)
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (2a) qH=0"<<G4endl;
#endif
      delete qH;
      // G4cerr<<"*G4QN::EvaNuc:antiBaryon="<<thePDG<<" below MassShell, M="<<totMass<<G4endl;
      // throw G4QException("G4QNucleus::EvaporateNucleus: anti-Baryon is below Mass Shell");
      G4ExceptionDescription ed;
      ed << "anti-Baryon is below Mass Shell: antiBaryon=" << thePDG
         << " below MassShell, M=" << totMass << G4endl;
      G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0002",
                  FatalException, ed);
    }
    else                                 // Decay in gamma or charged pion (@@ neutral)
    {
      G4double d=totMass-gsM;
#ifdef debug
      G4cout<<"G4QN::EvaporNucl: PDG="<<thePDG<<",M="<<totMass<<">"<<gsM<<",d="<<d<<G4endl;
#endif
      G4int decPDG=22;
      G4double decM=0.;
      if(d>142.)                           // @@ to avoid more precise calculations
      {
        if(thePDG==89999000)               // anti (p* -> n + pi+)
        {
          gsM=mNeut;
          thePDG=89999999;
          decPDG=-211;
          decM=mPi;
        }
        else if(thePDG==89999999)          // anti (n* -> p + pi-)
        {
          gsM=mProt;
          thePDG=89999000;
          decPDG=211;
          decM=mPi;
        }
        else                               // decay in Pi0
        {
          decPDG=111;
          decM=mPi0;
        }
      }
      G4LorentzVector h4Mom(0.,0.,0.,gsM); // GSMass must be updated in previous while-LOOP
      G4LorentzVector g4Mom(0.,0.,0.,decM);
      if(!G4QHadron(q4M).DecayIn2(h4Mom, g4Mom))
      {
#ifdef qdebug
        if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (3a) qH=0"<<G4endl;
#endif
        delete qH;
        // G4cerr<<"**G4QNuc::EvaNuc:ah="<<thePDG<<"(GSM="<<gsM<<")+gam>tM="<<totMass<<G4endl;
        // throw G4QException("G4QNucleus::EvaporateNucleus:BaryonDecay In Baryon+Gam Error");
        G4ExceptionDescription ed;
        ed << "BaryonDecay In Baryon+Gam Error: ah=" << thePDG << "(GSM="
           << gsM << ")+gam>tM=" << totMass << G4endl;
        G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0003",
                    FatalException, ed);
      }
#ifdef debug
      G4cout<<"G4QNucl::EvaNuc:aM="<<totMass<<q4M<<"->"<<thePDG<<h4Mom<<"+g="<<g4Mom<<",n="
            <<evaHV->size()<<G4endl;
#endif
      G4QHadron* curH = new G4QHadron(thePDG, h4Mom);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: antiHadr="<<thePDG<<h4Mom<<G4endl;
#endif
      evaHV->push_back(curH);         // Fill Baryon (delete equiv.)
      G4QHadron* curMes = new G4QHadron(decPDG, g4Mom);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: (anti) Gamma(pion)4M="<<g4Mom<<G4endl;
#endif
      evaHV->push_back(curMes);         // Fill gamma/pion (delete equivalent)
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (4a) qH=0"<<G4endl;
#endif
      delete qH;
    }
  }
  else if((thePDG==90001999||thePDG==89999002) && totMass>1080.) // @@ to avoid threshold
  //else if(2>3)// One can easily close this decay as it will be done later (time of calc?)
  {
    G4double gsM=mNeut;
    G4int barPDG=2112;
    G4int mesPDG=-211;
    if(thePDG==90001999)
    {
      gsM=mProt;
      barPDG=2212;
      mesPDG=211;
    }
    if(fabs(totMass-gsM-mPi)<.001)
    {
#ifdef debug
      G4cout<<"G4QN::EvaporateNuc:(D)GSM="<<gsM<<",H="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
      G4LorentzVector h4Mom=q4M*(gsM/totMass);           // At rest in CM
      G4QHadron* curB = new G4QHadron(barPDG,h4Mom);
      evaHV->push_back(curB); // (delete equivalent)
      G4LorentzVector g4Mom=q4M*(mPi/totMass);
      G4QHadron* curM = new G4QHadron(mesPDG,g4Mom);
      evaHV->push_back(curM); // (delete equivalent)
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (5) qH=0"<<G4endl;
#endif
      delete qH;
      return;
    }
    else if(totMass<gsM+mPi)
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (6) qH=0"<<G4endl;
#endif
      delete qH;
      // G4cerr<<"***G4QN::EvaNuc:Delta "<<thePDG<<" is belowMassShell, M="<<totMass<<G4endl;
      // throw G4QException("G4QNucleus::EvaporateNucleus: Delta is below Mass Shell");
      G4ExceptionDescription ed;
      ed << "Delta is below Mass Shell: Delta " << thePDG
         << " is belowMassShell, M=" << totMass << G4endl;
      G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0004",
                  FatalException, ed);
    }
    else                                 // Decay in gamma or charged pion (@@ neutral)
    {
      G4LorentzVector h4Mom(0.,0.,0.,gsM); // GSMass must be updated in previous while-LOOP
      G4LorentzVector g4Mom(0.,0.,0.,mPi);
      if(!G4QHadron(q4M).DecayIn2(h4Mom, g4Mom))
      {
#ifdef qdebug
        if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (7) qH=0"<<G4endl;
#endif
        delete qH;
        // G4cerr<<"***G4QNuc::EvaNuc:Dh="<<thePDG<<"N+pi="<<gsM+mPi<<">tM="<<totMass<<G4endl;
        // throw G4QException("G4QNucleus::EvaporateNucleus: DeltaDecInBaryon+Pi Error");
        G4ExceptionDescription ed;
        ed << "DeltaDecInBaryon+Pi Error: Dh=" << thePDG << "N+pi=" << gsM+mPi
           << ">tM=" << totMass << G4endl;
        G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0005",
                  FatalException, ed);
      }
#ifdef debug
      G4cout<<"G4QNuc::EvaNuc:"<<totMass<<q4M<<"->"<<thePDG<<h4Mom<<"+pi="<<g4Mom<<", nH="
            <<evaHV->size()<<G4endl;
#endif
      G4QHadron* curH = new G4QHadron(barPDG,h4Mom);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucl: Nucleon="<<thePDG<<h4Mom<<G4endl;
#endif
      evaHV->push_back(curH);         // Fill the nucleon (delete equiv.)
      G4QHadron* curG = new G4QHadron(mesPDG,g4Mom);
#ifdef debug
      G4cout<<"G4QE::EvaporateR: Pion="<<g4Mom<<G4endl;
#endif
      evaHV->push_back(curG);         // Fill the pion (delete equivalent)
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (8) qH=0"<<G4endl;
#endif
      delete qH;
    }
  }
  else if((thePDG==89998001||thePDG==90000998) && totMass>1080.) // @@ to avoid threshold
  //else if(2>3)// One can easily close this decay as it will be done later (time of calc?)
  {
    G4double gsM=mNeut;
    G4int barPDG=-2112;
    G4int mesPDG=211;
    if(thePDG==89998001)
    {
      gsM=mProt;
      barPDG=-2212;
      mesPDG=-211;
    }
    if(fabs(totMass-gsM-mPi)<.001)
    {
#ifdef debug
      G4cout<<"G4QN::EvaporateNuc:(A)GSM="<<gsM<<",H="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
      G4LorentzVector h4Mom=q4M*(gsM/totMass);           // At rest in CM
      G4QHadron* curB = new G4QHadron(barPDG,h4Mom);
      evaHV->push_back(curB); // (delete equivalent)
      G4LorentzVector g4Mom=q4M*(mPi/totMass);
      G4QHadron* curM = new G4QHadron(mesPDG,g4Mom);
      evaHV->push_back(curM); // (delete equivalent)
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (5a) qH=0"<<G4endl;
#endif
      delete qH;
      return;
    }
    else if(totMass<gsM+mPi)
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (6a) qH=0"<<G4endl;
#endif
      delete qH;
      // G4cerr<<"***G4QN::EvaNuc:aDelta "<<thePDG<<" is belowMassShell, M="<<totMass<<G4endl;
      // throw G4QException("G4QNucleus::EvaporateNucleus: anti-Delta is below Mass Shell");
      G4ExceptionDescription ed;
      ed << "anti-Delta is below Mass Shell: aDelta " << thePDG
         << " is belowMassShell, M=" << totMass << G4endl;
      G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0006",
                  FatalException, ed);
    }
    else                                 // Decay in gamma or charged pion (@@ neutral)
    {
      G4LorentzVector h4Mom(0.,0.,0.,gsM); // GSMass must be updated in previous while-LOOP
      G4LorentzVector g4Mom(0.,0.,0.,mPi);
      if(!G4QHadron(q4M).DecayIn2(h4Mom, g4Mom))
      {
#ifdef qdebug
        if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (7a) qH=0"<<G4endl;
#endif
        delete qH;
        // G4cerr<<"***G4QNuc::EvaNuc:aD="<<thePDG<<"N+pi="<<gsM+mPi<<">tM="<<totMass<<G4endl;
        // throw G4QException("G4QNucleus::EvaporateNucleus:AntiDeltaDecayInBaryon+Pi Error");
        G4ExceptionDescription ed;
        ed << "AntiDeltaDecayInBaryon+Pi Error: aD=" << thePDG << "N+pi="
           << gsM+mPi << ">tM=" << totMass << G4endl;
        G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0007",
                    FatalException, ed);
      }
#ifdef debug
      G4cout<<"G4QNuc::EvaNuc:(aD) "<<totMass<<q4M<<"->"<<thePDG<<h4Mom<<" + pi="<<g4Mom
            <<", nH="<<evaHV->size()<<G4endl;
#endif
      G4QHadron* curH = new G4QHadron(barPDG,h4Mom);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucl: Nucleon="<<thePDG<<h4Mom<<G4endl;
#endif
      evaHV->push_back(curH);         // Fill the nucleon (delete equiv.)
      G4QHadron* curMes = new G4QHadron(mesPDG,g4Mom);
#ifdef debug
      G4cout<<"G4QE::EvaporateR: Pion="<<g4Mom<<G4endl;
#endif
      evaHV->push_back(curMes);         // Fill the pion (delete equivalent)
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (8a) qH=0"<<G4endl;
#endif
      delete qH;
    }
  }
  else if(theBN>0&&theS<0) DecayAntiStrange(qH,evaHV); // "AntyStrangeNucleus" (del.eq.)
  else if(theBN>0&&(theC<0||theC>theBN-theS))DecayIsonucleus(qH,evaHV);//"Isonucleus"(d.e.)
  else if((thePDG==89999003 || thePDG==90002999) && totMass>2020.) //=> "ISO-dibarion"
  {
    // @@ Check that it never comes here !!
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
#ifdef qdebug
        if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (9) qH=0"<<G4endl;
#endif
        delete qH;
        // G4cerr<<"***G4QNucleus::EvaporateNucleus: tM="<<totMass<<"-> 2N="<<nucPDG<<"(M="
        //       <<nucM<<") + pi="<<piPDG<<"(M="<<mPi<<")"<<G4endl;
        // throw G4QException("G4QNucleus::EvaporateNucleus: ISO-Dibaryon DecayIn3 error");
        G4ExceptionDescription ed;
        ed << " ISO-Dibaryon DecayIn3 error: tM=" << totMass << "-> 2N="
           << nucPDG << "(M=" << nucM << ") + pi=" << piPDG << "(M="
           << mPi << ")" << G4endl;
        G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0008",
                    FatalException, ed);
      }
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (10) qH=0"<<G4endl;
#endif
      delete qH;
      G4QHadron* h1H = new G4QHadron(nucPDG,n14M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: Bar1="<<nucPDG<<n14M<<G4endl;
#endif
      evaHV->push_back(h1H);                // (delete equivalent)
      G4QHadron* h2H = new G4QHadron(nucPDG,n24M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: Bar2="<<nucPDG<<n24M<<G4endl;
#endif
      evaHV->push_back(h2H);                // (delete equivalent)
      G4QHadron* piH = new G4QHadron(piPDG,pi4M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: Pi="<<piPDG<<pi4M<<G4endl;
#endif
      evaHV->push_back(piH);                // (delete equivalent)
    }
    else
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (11) qH=0"<<G4endl;
#endif
      delete qH;
      // G4cerr<<"***G4QNucleus::EvapNucleus: IdPDG="<<thePDG<<", q4M="<<q4M<<", M="<<totMass
      //       <<" < M_2N+Pi, d="<<totMass-2*nucM-mPi<<G4endl;
      // throw G4QException("G4QNucleus::EvaporateNucleus:ISO-DiBaryon is under MassShell");
      G4ExceptionDescription ed;
      ed << "ISO-DiBaryon is under MassShell: IdPDG=" << thePDG << ", q4M="
         << q4M << ", M=" << totMass << " < M_2N+Pi, d=" << totMass-2*nucM-mPi
         << G4endl;
      G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0009",
                  FatalException, ed);
    }
  }
  else if((thePDG==90000002||thePDG==90001001||thePDG==90002000)&&totMass>2020.) //=> IsoBP
  {
    // @@ Pi0 can be taken into account !
    G4int  n1PDG = 2212;
    G4int  n2PDG = 2112;
    G4int  piPDG = -211;
    G4double n1M = mProt;
    G4double n2M = mNeut;
    if     (thePDG==90002000) piPDG  =  211;      // pp -> np + pi-
    else if(thePDG==90000002) piPDG  = -211;      // nn -> np + pi-
    else                                          // np -> 50%(nnpi+) 50%(pppi-)
    {
      if(G4UniformRand()>.5)
      {
        n1PDG=2112;
        n1M=mNeut;
        piPDG  =  211;
      }
      else
      {
        n2PDG=2212;
        n2M=mProt;
        piPDG  = -211;
      }
    }
    if(totMass>mPi+n1M+n2M)
    {
      G4LorentzVector n14M(0.,0.,0.,n1M);
      G4LorentzVector n24M(0.,0.,0.,n2M);
      G4LorentzVector pi4M(0.,0.,0.,mPi);
      if(!G4QHadron(q4M).DecayIn3(n14M,n24M,pi4M))
      {
        // G4cerr<<"***G4QNucl::EvapNucleus:IsoDN, tM="<<totMass<<"-> N1="<<n1PDG<<"(M="<<n1M
        //       <<") + N2="<<n2PDG<<"(M="<<n2M<<") + pi="<<piPDG<<" (Mpi="<<mPi<<")"<<G4endl;
        // throw G4QException("G4QNucl::EvaporateNucleus:ISO-dibaryon excit. DecayIn3 error");
        G4ExceptionDescription ed;
        ed << "ISO-dibaryon excit. DecayIn3 error: IsoDN, tM=" << totMass
           << "-> N1=" << n1PDG << "(M=" << n1M << ") + N2=" << n2PDG
           << "(M=" << n2M << ") + pi=" << piPDG << " (Mpi=" << mPi << ")"
           << G4endl;
        G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0010",
                    FatalException, ed);
      }
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (12) qH=0"<<G4endl;
#endif
      delete qH;
      G4QHadron* h1H = new G4QHadron(n1PDG,n14M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: Bar1="<<n1PDG<<n14M<<G4endl;
#endif
      evaHV->push_back(h1H);                // (delete equivalent)
      G4QHadron* h2H = new G4QHadron(n2PDG,n24M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: Bar2="<<n2PDG<<n24M<<G4endl;
#endif
      evaHV->push_back(h2H);                // (delete equivalent)
      G4QHadron* piH = new G4QHadron(piPDG,pi4M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: Pi="<<piPDG<<pi4M<<G4endl;
#endif
      evaHV->push_back(piH);                // (delete equivalent)
    }
    else
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (13) qH=0"<<G4endl;
#endif
      delete qH;
      // G4cerr<<"***G4QNuc::EvaporateNucleus: IdPDG="<<thePDG<<", q4M="<<q4M<<", M="<<totMass
      //       <<" < M1+M2+Pi, d="<<totMass-n1M-n2M-mPi<<G4endl;
      // throw G4QException("G4QNucleus::EvaporateNucleus:IsoDiBarState is under MassShell");
      G4ExceptionDescription ed;
      ed << "IsoDiBarState is under MassShell:  IdPDG=" << thePDG << ", q4M="
         << q4M << ", M=" << totMass << " < M1+M2+Pi, d="
         << totMass-n1M-n2M-mPi << G4endl;
        G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0011",
                    FatalException, ed);
    }
  }
  else if(theBN==2) DecayDibaryon(qH, evaHV);       //=> "Dibaryon" case (del eq.)
  else if((thePDG==90000997 || thePDG==89997001) && totMass>2020.) //=> "anti-ISO-dibarion"
  {
    // @@ Check that it never comes here !!
    G4int  nucPDG = -2112;
    G4double nucM = mNeut;
    G4int   piPDG = 211;
    if(thePDG==90002999)
    {
      nucPDG = -2212;
      nucM   = mProt;
      piPDG  = -211;
    }
    if(totMass>mPi+nucM+nucM)
    {
      G4LorentzVector n14M(0.,0.,0.,nucM);
      G4LorentzVector n24M(0.,0.,0.,nucM);
      G4LorentzVector pi4M(0.,0.,0.,mPi);
      if(!G4QHadron(q4M).DecayIn3(n14M,n24M,pi4M))
      {
#ifdef qdebug
        if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (9a) qH=0"<<G4endl;
#endif
        delete qH;
        // G4cerr<<"***G4QNucleus::EvaporateNucleus:antiM="<<totMass<<"-> 2aN="<<nucPDG<<"(M="
        //       <<nucM<<") + pi="<<piPDG<<"(M="<<mPi<<")"<<G4endl;
        // throw G4QException("G4QNucleus::EvaporateNucleus:Anti-ISO-DibaryonDecayIn3 error");
        G4ExceptionDescription ed;
        ed << "Anti-ISO-DibaryonDecayIn3 error: antiM=" << totMass
           << "-> 2aN=" << nucPDG << "(M=" << nucM << ") + pi=" << piPDG
           << "(M=" << mPi << ")" << G4endl;
        G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0012",
                    FatalException, ed);
      }
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (10a) qH=0"<<G4endl;
#endif
      delete qH;
      G4QHadron* h1H = new G4QHadron(nucPDG,n14M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus:(I) antiBar1="<<nucPDG<<n14M<<G4endl;
#endif
      evaHV->push_back(h1H);                // (delete equivalent)
      G4QHadron* h2H = new G4QHadron(nucPDG,n24M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus:(I) antiBar2="<<nucPDG<<n24M<<G4endl;
#endif
      evaHV->push_back(h2H);                // (delete equivalent)
      G4QHadron* piH = new G4QHadron(piPDG,pi4M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus:(I) (anti) Pi="<<piPDG<<pi4M<<G4endl;
#endif
      evaHV->push_back(piH);                // (delete equivalent)
    }
    else
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (11a) qH=0"<<G4endl;
#endif
      delete qH;
      // G4cerr<<"***G4QNucleus::EvapNucleus: aIdPDG="<<thePDG<<", q4M="<<q4M<<", M="<<totMass
      //       <<" < M_2N+Pi, d="<<totMass-2*nucM-mPi<<G4endl;
      // throw G4QException("G4QNucleus::EvaporateNucleus:AntiISODiBaryon is underMassShell");
      G4ExceptionDescription ed;
      ed << "AntiISODiBaryon is underMassShell: aIdPDG=" << thePDG << ", q4M="
         << q4M << ", M=" << totMass << " < M_2N+Pi, d=" << totMass-2*nucM-mPi
         << G4endl;
       G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0013",
                   FatalException, ed);
    }
  }
  else if((thePDG==89999998||thePDG==89998999||thePDG==89998000)&&totMass>2020.)//=>AnIsoBP
  {
    // @@ Pi0 can be taken into account !
    G4int  n1PDG = -2212;
    G4int  n2PDG = -2112;
    G4int  piPDG = 211;                           // dummy initialization
    G4double n1M = mProt;
    G4double n2M = mNeut;
    if     (thePDG==89998000) piPDG  = -211;      // anti ( pp -> np + pi- )
    else if(thePDG==89999998) piPDG  =  211;      // anti ( nn -> np + pi- )
    else                                          // anti ( np -> 50%(nnpi+) 50%(pppi-) )
    {
      if(G4UniformRand()>.5)
      {
        n1PDG=-2112;
        n1M=mNeut;
        piPDG  = -211;
      }
      else
      {
        n2PDG=-2212;
        n2M=mProt;
        piPDG  =  211;
      }
    }
    if(totMass>mPi+n1M+n2M)
    {
      G4LorentzVector n14M(0.,0.,0.,n1M);
      G4LorentzVector n24M(0.,0.,0.,n2M);
      G4LorentzVector pi4M(0.,0.,0.,mPi);
      if(!G4QHadron(q4M).DecayIn3(n14M,n24M,pi4M))
      {
        // G4cerr<<"**G4QNucl::EvapNucleus:IsoDN,antiM="<<totMass<<"-> N1="<<n1PDG<<"(M="<<n1M
        //       <<") + N2="<<n2PDG<<"(M="<<n2M<<") + pi="<<piPDG<<" (Mpi="<<mPi<<")"<<G4endl;
        // throw G4QException("G4QNucl::EvaporateNucleus:AntiExcitedDibaryon DecayIn3 error");
        G4ExceptionDescription ed;
        ed << "AntiExcitedDibaryon DecayIn3 error: IsoDN,antiM=" << totMass
           << "-> N1=" << n1PDG << "(M=" << n1M << ") + N2=" << n2PDG << "(M="
           << n2M << ") + pi=" << piPDG << " (Mpi=" << mPi << ")" << G4endl;
        G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0014",
                    FatalException, ed);
      }
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (12a) qH=0"<<G4endl;
#endif
      delete qH;
      G4QHadron* h1H = new G4QHadron(n1PDG,n14M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: antiBar1="<<n1PDG<<n14M<<G4endl;
#endif
      evaHV->push_back(h1H);                // (delete equivalent)
      G4QHadron* h2H = new G4QHadron(n2PDG,n24M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: antiBar2="<<n2PDG<<n24M<<G4endl;
#endif
      evaHV->push_back(h2H);                // (delete equivalent)
      G4QHadron* piH = new G4QHadron(piPDG,pi4M);
#ifdef debug
      G4cout<<"G4QNucleus::EvaporateNucleus: (anti)Pi="<<piPDG<<pi4M<<G4endl;
#endif
      evaHV->push_back(piH);                // (delete equivalent)
    }
    else
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (13a) qH=0"<<G4endl;
#endif
      delete qH;
      // G4cerr<<"***G4QNuc::EvaporateNucleus:andPDG="<<thePDG<<", q4M="<<q4M<<", M="<<totMass
      //       <<" < M1+M2+Pi, d="<<totMass-n1M-n2M-mPi<<G4endl;
      // throw G4QException("G4QNucleus::EvaporateNucleus:AntiDiBarState is under MassShell");
      G4ExceptionDescription ed;
      ed << "AntiDiBarState is under MassShell: andPDG=" << thePDG << ", q4M="
         << q4M << ", M=" << totMass << " < M1+M2+Pi, d="
         << totMass-n1M-n2M-mPi << G4endl;
      G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0015",
                  FatalException, ed);
    }
  }
  else if(theBN==-2) DecayAntiDibaryon(qH,evaHV);       //=> "Anti-Dibaryon" case (del eq.)
  else if(fabs(totMass-totGSM)<.001)  // Fill as it is or decay Be8, He5, Li5 (@@ add more)
  {
#ifdef debug
    G4cout<<"G4QNucleus::EvaNuc:GS "<<qH->GetQC()<<qH->Get4Momentum()<<" FillAsIs"<<G4endl;
#endif
    if(thePDG==90004004)
    {
      DecayAlphaAlpha(qH,evaHV);
    } // "Alpha+Alpha Decay" case (del eq.)
    else if(thePDG==90004002)
    {
      DecayAlphaDiN(qH,evaHV);
    } // Decay alpha+2p(alpha+2n is stable)
    else if((theC==theBN||theN==theBN||theS==theBN)&&theBN>1)
    {
      DecayMultyBaryon(qH,evaHV);
    }
    else if(theBN==5)
    {
      DecayAlphaBar(qH,evaHV);
    }   // Decay unstable A5 system (del eq.)
    else
    {
      evaHV->push_back(qH);
    }      // Fill as it is (del eq.)
  }
  else if(theBN>1 && thePDG>88000000 && thePDG<89000000) //==> 2antiK in the nucleus
  {
    G4cout<<"---Warning---G4QNucl::EvaNuc:MustNotBeHere.PDG="<<thePDG<<",S="<<theS<<G4endl;
    G4int bZ=theQC.GetCharge();
    G4int bN=theBN-bZ;
    G4int k1PDG = 321;
    G4double mK1= mK;
    G4int k2PDG = 321;
    G4double mK2= mK;
    G4int  nucPDG = thePDG;
    if(bZ>=bN) nucPDG+=999000;
    else
    {
      nucPDG+=999999;
      k1PDG = 311;
      mK1= mK0;
    }
    if(bZ>bN) nucPDG+=999000;
    else
    {
      nucPDG+=999999;
      k2PDG = 311;
      mK2= mK0;
    }
    G4double nucM = G4QNucleus(nucPDG).GetGSMass();
    G4cout<<"-Warning-G4QN::EvN:M="<<nucM<<","<<nucPDG<<",1="<<k1PDG<<",2="<<k2PDG<<G4endl;
    G4LorentzVector n4M(0.,0.,0.,nucM);
    G4LorentzVector k14M(0.,0.,0.,mK1);
    G4LorentzVector k24M(0.,0.,0.,mK2);
    if(!G4QHadron(q4M).DecayIn3(n4M,k14M,k24M))
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (14) qH=0"<<G4endl;
#endif
      delete qH;
      // G4cout<<"***G4QNucleus::EvaNuc:tM="<<totMass<<"-> N="<<nucPDG<<"(M="<<nucM<<") + k1="
      //       <<k1PDG<<"(M="<<mK1<<") + k2="<<k2PDG<<"(M="<<mK2<<")"<<G4endl;
      // throw G4QException("G4QNucleus::EvaporateNucleus: Nucleus+2antiK DecayIn3 error");
      G4ExceptionDescription ed;
      ed << " Nucleus+2antiK DecayIn3 error: tM=" << totMass << "-> N="
         << nucPDG << "(M=" << nucM << ") + k1=" << k1PDG << "(M=" << mK1
         << ") + k2=" << k2PDG << "(M=" << mK2 << ")" << G4endl;
      G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0016",
                  FatalException, ed);
    }
#ifdef qdebug
    if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (15) qH=0"<<G4endl;
#endif
    delete qH;
    G4QHadron* k1H = new G4QHadron(k1PDG,k14M);
#ifdef debug
    G4cout<<"G4QNucleus::EvaporateNucleus: k1="<<k1PDG<<k14M<<G4endl;
#endif
    evaHV->push_back(k1H);                // (delete equivalent)
    G4QHadron* k2H = new G4QHadron(k2PDG,k24M);
#ifdef debug
    G4cout<<"G4QNucleus::EvaporateNucleus: k2="<<k2PDG<<k24M<<G4endl;
#endif
    evaHV->push_back(k2H);                // (delete equivalent)
    G4QHadron* nH = new G4QHadron(nucPDG,n4M);
#ifdef debug
    G4cout<<"G4QNucleus::EvaporateNucleus: Nuc="<<nucPDG<<n4M<<G4endl;
#endif
    evaHV->push_back(nH);                 // (delete equivalent)
  }
  // ***>> From here the EVA code starts (baryons/hyperons can be excited) <<***
  else if ( (thePDG > 80000000 && thePDG != 90000000) || 
     thePDG == 2112 || thePDG == 2212 || thePDG == 3122)
  { // @@ Improve for Sigma+, Sigma-, Ksi0 & Ksi- content in the Total Np/Nn Nuclei
    if(thePDG<80000000)                        // Switch from QHadronCode to QNuclearCode
    {
      if     (thePDG==2112) thePDG=90000001;   // n
      else if(thePDG==2212) thePDG=90001000;   // p
      else if(thePDG==3122) thePDG=91000000;   // lambda
    }
    G4QNucleus qNuc(q4M,thePDG);               // Make a Nucleus for theTotalResidNucleus
    G4double GSMass =qNuc.GetGSMass();         // GSMass of the Total Residual Nucleus
    G4QContent totQC=qNuc.GetQCZNS();          // QuarkCont of theTotalResidNucleus (theQC)
    G4int    bA     =qNuc.GetA();              // A#of baryons in Total Residual Nucleus
    G4int    bZ     =qNuc.GetZ();              // A#of protons in the Total ResidualNucleus
    G4int    bN     =qNuc.GetN();              // A#of neutrons in the TotalResidualNucleus
#ifdef ppdebug
    G4cout<<"G4QN::EvaNuc: theBN="<<theBN<<", bA="<<bA<<", bZ="<<bZ<<", bN="<<bN<<G4endl;
#endif
    G4int    bS     =qNuc.GetS();              // A#of lambdas in the Total ResidualNucleus
#ifdef debug
    if(bZ==2&&bN==5)G4cout<<"G4QNucleus::EvaNucl: (2,5) GSM="<<GSMass<<" > "
                          <<G4QPDGCode(2112).GetNuclMass(2,4,0)+mNeut<<G4endl;
    if(bZ==1&&bN==3)G4cout<<"G4QNucl::EvaNucl: (1,3) GSM="<<GSMass<<" > "
                          <<G4QPDGCode(2112).GetNuclMass(1,2,0)+mNeut<<G4endl;
    G4double dM=totMass-GSMass;
    G4cout<<"G4QNucl::EvaNuc:"<<qNuc<<",PDG="<<thePDG<<",M="<<totMass<<",dM="<<dM<<G4endl;
    ////////if(dM>7.) throw G4QException("G4QNucleus::EvaporateNucleus: Big Excitation");
#endif
    G4int   bsCond =qNuc.SplitBaryon();        // (Bary/Deut/Alph)SeparCond for TotResNucl
    G4bool  dbsCond=qNuc.Split2Baryons();      // (Two Baryons)SeparCond for TotResidNucl
#ifdef debug
    G4cout<<"G4QNucleus::EvaporateNuc: bs="<<bsCond<<", dbs="<<dbsCond<<", A="<<bA<<G4endl;
#endif
    if(fabs(totMass-GSMass)<.003&&!bsCond&&!dbsCond) // GS or can't split 1(2)B FillAsItIs
    {
#ifdef debug
      G4cout<<"G4QN::EvaNuc: GS direct "<<qH->GetQC()<<qH->Get4Momentum()<<" AsIs"<<G4endl;
#endif
      evaHV->push_back(qH);
      return;
    }
    else if ( ( bA == 1 || (!bsCond && !dbsCond) ) && totMass > GSMass+.003 )//=>Fuse&Decay
    //else if(2>3)                                // Close "Fuse&Decay Technology"***@@@***
    {
#ifdef debug
      G4cout<<"G4QN::EvaN:SplitBar, s="<<bsCond<<",M="<<totMass<<" > GSM="<<GSMass<<G4endl;
#endif
      G4int nOfOUT = evaHV->size();            // Total#of QHadrons in Vector at this point
      G4bool bnfound=true;                     // Cure "back fusion fragment not found"
      while(nOfOUT)                            // Try BackFusionDecays till something is in
      {
        G4QHadron*     theLast = (*evaHV)[nOfOUT-1];
        G4int          lastBN = theLast->GetBaryonNumber();
        G4int          nFrag  = theLast->GetNFragments();
        //////////////////G4int          gam    = theLast->GetPDGCode();
#ifdef debug
        G4cout<<"G4QN::EvaNuc:*BackFus*,BN="<<lastBN<<",nF="<<nFrag<<",n="<<nOfOUT<<G4endl;
#endif
        while(nFrag)                       // => "Delete Last Decayed Hadron" case
        {
          G4QHadron* thePrev = (*evaHV)[nOfOUT-2];
          nFrag  = thePrev->GetNFragments();
          G4int      prevBN = thePrev->GetBaryonNumber();
#ifdef debug
          G4int     prevPDG = thePrev->GetPDGCode();
          G4cout<<"G4QNucl::EvaNucl: DelTheLast, nFr="<<nFrag<<", pPDG="<<prevPDG<<G4endl;
#endif
          evaHV->pop_back();               // the prev QHadron is excluded from OUTPUT
          delete theLast;//!!When kill,DON'T forget to del. theLastQHadron as an instance!!
          theLast = thePrev;               // Update theLastPntr(Prev instead of Last)
          lastBN=prevBN;
          nOfOUT--;                        // Reduce the stack for the Last decayed hadron
        }
        if(nOfOUT)
        {
          if(lastBN<1&&nOfOUT>1)           // Try Once To Skip Mesons/Gammas & Antibaryons
          {
            G4QHadron* thePrev = (*evaHV)[nOfOUT-2];//***Exchange between Last & Prev***
            evaHV->pop_back();             // the last QHadron is excluded from OUTPUT
            evaHV->pop_back();             // the prev QHadron is excluded from OUTPUT
            evaHV->push_back(theLast);     // the Last becomes the Prev (1st part of exch)
            evaHV->push_back(thePrev);     // the Prev becomes the Last (2nd part of exch)
            theLast = thePrev;             // Update theLastPointer (Prev instead of Last)
          }
          G4LorentzVector last4M = theLast->Get4Momentum(); // Selected the last 4-Mom
          G4QContent  lastQC = theLast->GetQC();
          G4double lastM  = last4M.m();    // Mass of the Probable Last BackFused Fragment
          totQC+=lastQC;                   // Update (increase) the total QC (prototype)
          q4M+=last4M;                     // Update (increase) the total 4-momentum
          totMass=q4M.m();                 // Calculate new real total mass of the fused
          G4int totPDG=totQC.GetSPDGCode();// The updated PDG for the TotalResidualNucleus
          if(totPDG==10&&totQC.GetBaryonNumber()>0) totPDG=totQC.GetZNSPDGCode();
          G4int totBN=totQC.GetBaryonNumber();// BaryonNumber of the Total Residual Nucleus
          G4double dM=totMass-GSMass-lastM;
#ifdef debug
          G4cout<<"G4QN::EvaNuc: tM="<<totMass<<"-LM="<<lastM<<lastQC<<"-GSM="<<GSMass<<"="
                <<dM<<G4endl;
#endif
          if(dM>-0.001)                    // Decay in the qH and the last is impossible
          {
            G4QHadron* evH = new G4QHadron(totPDG,q4M); // Create QHadron for CompResidNuc
            if(dM<=0.)
            {
              evaHV->pop_back();    // lastQHadron is excluded from QHadrV asIs in TRN
              delete theLast; //When kill, DON'T forget to delete lastQHadron asAnInstance!
#ifdef qdebug
              if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (16) qH=0"<<G4endl;
#endif
              delete qH;
#ifdef debug
              G4cout<<"G4QNucleus::EvaporateNucl: EVH "<<totPDG<<q4M<<" fill AsIs"<<G4endl;
#endif
              if(totBN==2)DecayDibaryon(evH,evaHV); // Fill dibaryon (with decay products)
              else evaHV->push_back(evH);// Fill TRN to HVect asIs (delete equivalent)
            }
            else                        // Decay TotalResidualNucleus in GSM+Last
            {
              G4LorentzVector r4Mom(0.,0.,0.,GSMass);
              if(!G4QHadron(q4M).DecayIn2(last4M,r4Mom)) // Decay failed
              {
                evaHV->pop_back(); // lastQHadron is excluded from QHadrV as is in TRN
                delete theLast; //When kill,DON'T forget to delete lastQHadron asAnInstance
#ifdef qdebug
                if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (17) qH=0"<<G4endl;
#endif
                delete qH;
#ifdef debug
                G4cout<<"***G4QNucleus::EvaNucl: EVH "<<totPDG<<q4M<<" fill AsIs"<<G4endl;
#endif
                evaHV->push_back(evH);// Fill TRN to Vect as it is (delete equivalent)
#ifdef debug
                G4cout<<"***G4QN::EvaN:DecayIn L"<<lastQC<<"+RN"<<totQC<<" failed"<<G4endl;
#endif
              }
              else                        // Decay in GSM+theLast succeeded
              {
                delete evH;
#ifdef qdebug
                if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (18) qH=0"<<G4endl;
#endif
                delete qH;
                theLast->Set4Momentum(last4M);// Already exists:don't create&fill,->set4Mom
                G4QHadron* nucH = new G4QHadron(thePDG,r4Mom); // Create QHadron for qH-nuc
#ifdef debug
                G4cout<<"G4QNucleus::EvaNuc:fill NH "<<totPDG<<r4Mom<<G4endl;
#endif
                // @@ What about others, not DB possibilities?
                if(thePDG==92000000||thePDG==90002000||thePDG==90000002)
                                                DecayDibaryon(nucH,evaHV);//DekayDibarions
                else evaHV->push_back(nucH);// Fill the Residual Nucleus (del.eq.)
              }
            }
            bnfound=false;
            break;
          }
          thePDG=totPDG;                   // Make ResidualNucleus outOf theTotResidualNucl
          GSMass=G4QPDGCode(thePDG).GetMass();// Update the Total Residual Nucleus mass
          evaHV->pop_back();               // the last QHadron is excluded from OUTPUT
          delete theLast;//!!When kill,DON'T forget to delete theLastQHadron asAnInstance!!
          nOfOUT--;                        // Update the value of OUTPUT entries
        }
      }
      if(bnfound)
      {
        G4LorentzVector h4Mom(0.,0.,0.,GSMass);//GSMass must be updated inPreviousWhileLOOP
        G4LorentzVector g4Mom(0.,0.,0.,0.);
        if(!G4QHadron(q4M).DecayIn2(h4Mom, g4Mom))
        {
#ifdef qdebug
          if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (19) qH=0"<<G4endl;
#endif
          delete qH;
          // G4cerr<<"**G4QN::EvaNuc:h="<<thePDG<<"(GSM="<<GSMass<<")+g>tM="<<totMass<<G4endl;
          // throw G4QException("G4QNucleus::EvaporateNucleus: Decay in Gamma failed");
          G4ExceptionDescription ed;
          ed << " Decay in Gamma failed: h=" << thePDG << "(GSM=" << GSMass
             << ")+g>tM=" << totMass << G4endl;
          G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0017",
                      FatalException, ed);
        }
#ifdef debug
        G4cout<<"G4QNuc::EvaNuc: "<<q4M<<"->totResN="<<thePDG<<h4Mom<<"+g="<<g4Mom<<G4endl;
#endif
        G4QHadron* curH = new G4QHadron(thePDG,h4Mom);
#ifdef debug
        G4cout<<"G4QNucleus::EvaporateNucleus: Fill a Fragment="<<thePDG<<h4Mom<<G4endl;
#endif
        if(thePDG==92000000||thePDG==90002000||thePDG==90000002) DecayDibaryon(curH,evaHV);
        else evaHV->push_back(curH);  // Fill the TotalResidualNucleus (del.equiv.)
        G4QHadron* curG = new G4QHadron(22,g4Mom);
#ifdef debug
        G4cout<<"G4QNucleus::EvaporateNucleus: Fill  a Gamma="<<g4Mom<<G4endl;
#endif
        evaHV->push_back(curG);       // Fill the gamma (delete equivalent)
#ifdef qdebug
        if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (20) qH=0"<<G4endl;
#endif
        delete qH;
      }
    }
    else if(bA>0&&bS<0) DecayAntiStrange(qH,evaHV);// Decay nucleus with antistrangeness
    else if(bA==2) DecayDibaryon(qH,evaHV); // Decay the residual dibaryon (del.equivalent)
    else if(bA==-2) DecayAntiDibaryon(qH,evaHV);   // Decay residual anti-dibaryon (del.eq)
    else if(totMass<GSMass+.003&&(bsCond||dbsCond))//==>" M<GSM but decay is possible" case
    {
#ifdef pdebug
      G4cout<<"G4QN::EvN:2B="<<dbsCond<<",B="<<bsCond<<",M="<<totMass<<"<"<<GSMass<<G4endl;
#endif
      //G4double gResM  =1000000.;           // Prototype of mass of residual for a gamma
      G4int    gResPDG=0;                  // Prototype of residualPDGCode for a gamma
      if(bN==4&&bZ==2&&!bS)                // It's He6 nucleus
      {
        gResPDG= thePDG;                   // PDG of the Residual Nucleus
        //gResM  = mHel6;                    // min mass of the Residual Nucleus
      }
      G4double nResM  =1000000.;           // Prototype of mass of residual for a neutron
      G4int    nResPDG=0;                  // Prototype of ResidualPDGCode for a neutron
      if(bsCond==2112&&bN>0&&bA>1)         // There's aNeutr in theNucl, which can be split
      {
        G4QContent resQC=totQC-neutQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        nResPDG=resN.GetPDG();             // PDG of the Residual Nucleus
        if     (nResPDG==90000001) nResM=mNeut;
        else if(nResPDG==90001000) nResM=mProt;
        else if(nResPDG==91000000) nResM=mLamb;
        else nResM=resN.GetMZNS();         // min mass of the Residual Nucleus
      }
      G4double pResM  =1000000.;           // Prototype of mass of residual for a proton
      G4int    pResPDG=0;                  // Prototype of PDGCode of residual for a proton
      if(bsCond==2212&&bZ>0&&bA>1)         // There's aProton in Nucl, which can be split
      {
        G4QContent resQC=totQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        pResPDG=resN.GetPDG();             // PDG of the Residual Nucleus
        if     (pResPDG==90000001) pResM=mNeut;
        else if(pResPDG==90001000) pResM=mProt;
        else if(pResPDG==91000000) pResM=mLamb;
        else pResM  =resN.GetMZNS();       // min mass of the Residual Nucleus
      }
      G4double lResM  =1000000.;           // Prototype of mass of residual for a Lambda
      G4int    lResPDG=0;                  // Prototype of PDGCode of residual for a Lambda
      if(bsCond==3122&&bS>0&&bA>1)         // There's aLambd in theNucl, which can be split
      {
        G4QContent resQC=totQC-lambQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        lResPDG=resN.GetPDG();             // PDG of the Residual Nucleus
        if     (lResPDG==90000001) lResM=mNeut;
        else if(lResPDG==90001000) lResM=mProt;
        else if(lResPDG==91000000) lResM=mLamb;
        else lResM  =resN.GetMZNS();       // min mass of the Residual Nucleus
      }
      G4double dResM  =1000000.;           // Prototype of mass of residual for a Alpha
      G4int    dResPDG=0;                  // Prototype of PDGCode of residual for a Alpha
      if(bsCond==90001001&&bN>0&&bZ>0&&bA>2)// There's aDeuter in Nucl, which canBeRadiated
      {
        G4QContent resQC=totQC-deutQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        dResPDG=resN.GetPDG();             // PDG of the Residual Nucleus
        if     (dResPDG==90000001) dResM=mNeut;
        else if(dResPDG==90001000) dResM=mProt;
        else if(dResPDG==91000000) dResM=mLamb;
        else dResM  =resN.GetMZNS();       // minMass of the Residual Nucleus
      }
      G4double aResM  =1000000.;           // Prototype of mass of residual for a Alpha
      G4int    aResPDG=0;                  // Prototype of PDGCode of residual for a Alpha
      if(bsCond==90002002&&bN>1&&bZ>1&&bA>4)// There's Alpha in theNucl, which can be split
      {
        G4QContent resQC=totQC-alphQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        aResPDG=resN.GetPDG();             // PDG of the Residual Nucleus
        if     (aResPDG==90000001) aResM=mNeut;
        else if(aResPDG==90001000) aResM=mProt;
        else if(aResPDG==91000000) aResM=mLamb;
        else aResM  =resN.GetMZNS();       // minMass of the Residual Nucleus
      }
      G4double nnResM  =1000000.;          // Prototype of mass of residual for a dineutron
      G4int    nnResPDG=0;                 // Prototype of ResidualPDGCode for a dineutron
      if(dbsCond&&bN>1&&bA>2)              // It's nucleus and there is a dineutron
      {
        G4QContent resQC=totQC-neutQC-neutQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        nnResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (nnResPDG==90000001) nnResM=mNeut;
        else if(nnResPDG==90001000) nnResM=mProt;
        else if(nnResPDG==91000000) nnResM=mLamb;
        else nnResM  =resN.GetMZNS();      // min mass of the Residual Nucleus
      }
      G4double ppResM  =1000000.;          // Prototype of mass of residual for a diproton
      G4int    ppResPDG=0;                 // Prototype of ResidualPDGCode for a diproton
      if(dbsCond&&bZ>1&&bA>2)              // It's nucleus and there is a diproton
      {
        G4QContent resQC=totQC-protQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        ppResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (ppResPDG==90000001) ppResM=mNeut;
        else if(ppResPDG==90001000) ppResM=mProt;
        else if(ppResPDG==91000000) ppResM=mLamb;
        else ppResM  =resN.GetMZNS();      // min mass of the Residual Nucleus
      }
      G4double npResM  =1000000.;          // Prototype of ResidualMass for proton+neutron
      G4int    npResPDG=0;                 // Prototype of ResidualPDGCode for a prot+neut
      if(dbsCond&&bN>0&&bZ>0&&bA>2)        // It's nucleus and there is aProton & aNeutron
      {
        G4QContent resQC=totQC-neutQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        npResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (npResPDG==90000001) npResM=mNeut;
        else if(npResPDG==90001000) npResM=mProt;
        else if(npResPDG==91000000) npResM=mLamb;
        else npResM  =resN.GetMZNS();      // min mass of the Residual Nucleus
      }
      G4double lnResM  =1000000.;          // Prototype of residualMass for lambda+neutron
      G4int    lnResPDG=0;                 // Prototype of ResidualPDGCode for aLambda+Neut
      if(dbsCond&&bN>0&&bS>0&&bA>2)        // It's nucleus and there is aLambda & aNeutron
      {
        G4QContent resQC=totQC-lambQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        lnResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (lnResPDG==90000001) lnResM=mNeut;
        else if(lnResPDG==90001000) lnResM=mProt;
        else if(lnResPDG==91000000) lnResM=mLamb;
        else lnResM  =resN.GetMZNS();      // min mass of the Residual Nucleus
      }
      G4double lpResM  =1000000.;          // Prototype of residualMass for a proton+lambda
      G4int    lpResPDG=0;                 // Prototype of ResidualPDGCode for theProt+lamb
      if(dbsCond&&bS>0&&bZ>0&&bA>2)        // It's nucleus and there is aProton and aLambda
      {
        G4QContent resQC=totQC-neutQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        lpResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (lpResPDG==90000001) lpResM=mNeut;
        else if(lpResPDG==90001000) lpResM=mProt;
        else if(lpResPDG==91000000) lpResM=mLamb;
        else lpResM  =resN.GetMZNS();      // minMass of the Residual Nucleus
      }
      G4double llResM  =1000000.;          // Prototype of mass of residual for a di-lambda
      G4int    llResPDG=0;                 // Prototype of ResidPDGCode for the di-lambda
      if(dbsCond&&bS>1&&bA>2)              // It's nucleus and there is a di-lambda
      {
        G4QContent resQC=totQC-neutQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        llResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (llResPDG==90000001) llResM=mNeut;
        else if(llResPDG==90001000) llResM=mProt;
        else if(llResPDG==91000000) llResM=mLamb;
        else llResM  =resN.GetMZNS();      // min mass of the Residual Nucleus
      }
#ifdef debug
      G4cout<<"G4QNucleus::EvaNucl: rP="<<pResPDG<<",rN="<<nResPDG<<",rL="<<lResPDG<<",N="
            <<bN<<",Z="<<bZ<<",nL="<<bS<<",totM="<<totMass<<",n="<<totMass-nResM-mNeut
            <<",p="<<totMass-pResM-mProt<<",l="<<totMass-lResM-mLamb<<G4endl;
#endif
      if ( thePDG == 90004004                                                 || 
          (thePDG == 90002004 && totMass > mHel6+.003)                        ||
          (bA > 4 && bsCond && bN > 1 && bZ > 1 && totMass > aResM+mAlph)     ||
          (bA > 1 && bsCond && ( (bN > 0 && totMass > nResM+mNeut) || 
                                 (bZ > 0 && totMass > pResM+mProt) || 
                                 (bS > 0 && totMass > lResM+mLamb) ) )        ||
          (bA > 2 && 
           (( bN > 0 && bZ > 0 && 
              ( (bsCond && totMass > dResM+mDeut) || (dbsCond && totMass > dResM+mDeut) )
            ) || ( dbsCond && ( (bN > 1   && totMass > nnResM+mNeut+mNeut) ||
                                (bZ > 1   && totMass > ppResM+mProt+mProt) ||
                                (bS > 1   && totMass > llResM+mLamb+mLamb) ||
                                (bN && bS && totMass > lnResM+mLamb+mNeut) ||
                                (bZ && bS && totMass > lpResM+mLamb+mProt)
                              )
                 )
           )
          )
         )
      {
        G4int barPDG = 90002002;           // Just for the default case of Be8->alpha+alpha
        G4int resPDG = 90002002;
        G4int thdPDG = 0;
        G4double barM= mAlph;
        G4double resM= mAlph;
        G4double thdM= mNeut;              // This default value is used in the IF
        G4double tMC=totMass+.0002;
        if(gResPDG&&tMC>mHel6+.003)        // Can make radiative decay of He6 (priority 0)
        {
          barPDG=90002004;
          resPDG=22;
          barM  =mHel6;
          resM  =0.;
        }
        else if(nResPDG&&tMC>nResM+mNeut)  // Can radiate a neutron (priority 1)
        {
          barPDG=90000001;
          resPDG=nResPDG;
          barM  =mNeut;
          resM  =nResM;
        }
        else if(pResPDG&&totMass+.001>pResM+mProt)   // Can radiate a proton (priority 2)
        {
          barPDG=90001000;
          resPDG=pResPDG;
          barM  =mProt;
          resM  =pResM;
        }
        else if(lResPDG&&tMC>lResM+mLamb)  // Can radiate a Lambda (priority 3) @@ Sigma0
        {
          barPDG=91000000;
          resPDG=lResPDG;
          barM  =mLamb;
          resM  =lResM;
        }
        else if(thePDG!=90004004&&bN>1&&bZ>1&&bA>4&&tMC>aResM+mAlph)// Decay in alpha (p4)
        {
          barPDG=90002002;
          resPDG=aResPDG;
          barM  =mAlph;
          resM  =aResM;
        }
        else if(dResPDG&&tMC>dResM+mDeut)  // Can radiate a Deuteron (priority 5)
        {
          barPDG=90001001;
          resPDG=dResPDG;
          barM  =mDeut;
          resM  =dResM;
        }
        else if(ppResPDG&&tMC>ppResM+mProt+mProt)// Can radiate a DiProton (priority 6)
        {
          barPDG=90001000;
          resPDG=ppResPDG;
          thdPDG=90001000;
          barM  =mProt;
          resM  =ppResM;
          thdM  =mProt;
        }
        else if(nnResPDG&&tMC>nnResM+mNeut+mNeut)// Can radiate a DiNeutron (priority 7)
        {
          barPDG=90000001;
          resPDG=nnResPDG;
          thdPDG=90000001;
          barM  =mNeut;
          resM  =nnResM;
        }
        else if(npResPDG&&tMC>npResM+mProt+mNeut)// Can radiate a neutron+proton (prior 8)
        {
          barPDG=90001000;
          resPDG=npResPDG;
          thdPDG=90000001;
          barM  =mProt;
          resM  =npResM;
        }
        else if(lnResPDG&&tMC>lnResM+mLamb+mNeut)// Can radiate a Lambda+neutron (prior 9)
        {
          barPDG=91000000; // @@ Sigma0
          resPDG=lnResPDG;
          thdPDG=90000001;
          barM  =mLamb;
          resM  =lnResM;
        }
        else if(lpResPDG&&tMC>lpResM+mLamb+mProt)// Can radiate a Lambda+proton (prior 10)
        {
          barPDG=91000000; // @@ Sigma0
          resPDG=lpResPDG;
          thdPDG=90001000;
          barM  =mLamb;
          resM  =lpResM;
          thdM  =mProt;
        }
        else if(llResPDG&&tMC>llResM+mLamb+mLamb)// Can radiate a DiLambda (priority 11)
        {
          barPDG=91000000; // @@ Sigma0
          resPDG=llResPDG;
          thdPDG=91000000; // @@ Sigma0
          barM  =mLamb;
          resM  =llResM;
          thdM  =mLamb;
        }
        else if(thePDG!=90004004&&tMC>GSMass)// If it's not Be8 decay in gamma & GSM
        {
          barPDG=thePDG;
          resPDG=22;
          barM  =GSMass;
          resM  =0.;
        }
        else if(thePDG!=90004004)
        {
          // G4cerr<<"***G4QNuc::EvaN:PDG="<<thePDG<<",M="<<totMass<<"< GSM="<<GSMass<<G4endl;
          // throw G4QException("G4QNucleus::EvaporateNucleus: M<GSM & can't decayInPNL");
          G4ExceptionDescription ed;
          ed << "M<GSM & can't decayInPNL: PDG=" << thePDG << ",M=" << totMass
             << "< GSM=" << GSMass << G4endl;
          G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0018",
                      FatalException, ed);
        }
        G4LorentzVector a4Mom(0.,0.,0.,barM);
        G4LorentzVector b4Mom(0.,0.,0.,resM);
        if(!thdPDG)
        {
          if(!qH->DecayIn2(a4Mom,b4Mom))
          {
            evaHV->push_back(qH);     // Fill as it is (delete equivalent)
            G4cout<<"---Warning---G4QNucleus::EvaNuc:rP="<<pResPDG<<",rN="<<nResPDG<<",rL="
                  <<lResPDG<<",N="<<bN<<",Z="<<bZ<<",L="<<bS<<",totM="<<totMass<<",n="
                  <<totMass-nResM-mNeut<<",p="<<totMass-pResM-mProt<<",l="
                  <<totMass-lResM-mLamb<<G4endl;
            G4cout<<"---Warning---G4QN::EvN:DecIn2Error b="<<barPDG<<",r="<<resPDG<<G4endl;
            return;
          }
          else
          {
#ifdef qdebug
            if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (21) qH=0"<<G4endl;
#endif
            delete qH;
            G4QHadron* HadrB = new G4QHadron(barPDG,a4Mom);
#ifdef debug
            G4cout<<"G4QNucleus::EvaNucleus:(1) Baryon="<<barPDG<<a4Mom<<G4endl;
#endif
            evaHV->push_back(HadrB);       // Fill the baryon (delete equivalent)
            G4QHadron* HadrR = new G4QHadron(resPDG,b4Mom);
#ifdef debug
            G4cout<<"G4QNucleus::EvaNucleus:(1) Residual="<<resPDG<<b4Mom<<G4endl;
#endif
            // @@ Self-call !!
            if(HadrR->GetBaryonNumber()>1) EvaporateNucleus(HadrR,evaHV);//ContinueDecay
            else evaHV->push_back(HadrR);  // Fill ResidNucl=Baryon to OutHadronVector
          }
        }
        else
        {
          G4LorentzVector c4Mom(0.,0.,0.,thdM);
          if(!qH->DecayIn3(a4Mom,b4Mom,c4Mom))
          {
            evaHV->push_back(qH);    // Fill as it is (delete equivalent)
            G4cout<<"---Warning---G4QN::EvN:rNN="<<nnResPDG<<",rNP="<<npResPDG<<",rPP="
                  <<ppResPDG<<",N="<<bN<<",Z="<<bZ<<",L="<<bS<<",tM="<<totMass<<",nn="
                  <<totMass-nnResM-mNeut-mNeut<<",np="<<totMass-npResM-mProt-mNeut<<",pp="
                  <<totMass-ppResM-mProt-mProt<<G4endl;
            G4cout<<"---Warning---G4QN::EvN:DecIn2Error,b="<<barPDG<<",r="<<resPDG<<G4endl;
            return;
          }
          else
          {
#ifdef qdebug
            if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (22) qH=0"<<G4endl;
#endif
            delete qH;
            G4QHadron* HadrB = new G4QHadron(barPDG,a4Mom);
#ifdef debug
            G4cout<<"G4QNucleus::EvaporateNucleus:(2) Baryon1="<<barPDG<<a4Mom<<G4endl;
#endif
            evaHV->push_back(HadrB);       // Fill the first baryon (del.equiv.)
            G4QHadron* HadrC = new G4QHadron(thdPDG,c4Mom);
#ifdef debug
            G4cout<<"G4QNucleus::EvaporateNucleus:(2) Baryon2="<<thdPDG<<c4Mom<<G4endl;
#endif
            evaHV->push_back(HadrC);       // Fill the second baryon (del.equiv.)
            G4QHadron* HadrR = new G4QHadron(resPDG,b4Mom);
#ifdef debug
            G4cout<<"G4QNucleus::EvaporateNucleus:(2) Residual="<<resPDG<<b4Mom<<G4endl;
#endif
            // @@ Self-call !!
            if(HadrR->GetBaryonNumber()>1) EvaporateNucleus(HadrR,evaHV); // Continue decay
            else evaHV->push_back(HadrR); // Fill ResidNucl=Baryon to OutputHadrVector
          }
        }
      }
      else if (fabs(totMass-GSMass)<.003) // @@ Looks like a duplication of the prev. check
      {
#ifdef debug
        G4cout<<"*|*|*|*G4QNucleus::EvaporateNuc: fill AsIs. Should never be here"<<G4endl;
#endif
        evaHV->push_back(qH);  // FillAsItIs (del.eq.)
        return;
      }
      else                             // "System is below mass shell and can't decay" case
      {
#ifdef debug
        G4cout<<"***G4QNucl::EvaNuc: tM="<<totMass<<"("<<thePDG<<") < GSM="<<GSMass<<", d="
              <<totMass-GSMass<<", QC="<<qH->GetQC()<<qH->Get4Momentum()<<"*AsIs*"<<G4endl;
#endif
        evaHV->push_back(qH);                   // Correct or fill as it is
        return;
      }
    }
    else                                        // ===> Evaporation of the excited system
    {
#ifdef pdebug
      G4cout<<"G4QN::EvaNuc:***EVA***tPDG="<<thePDG<<",M="<<totMass<<">GSM="<<GSMass<<",d="
            <<totMass-GSMass<<", N="<<qNuc.Get4Momentum()<<qNuc.Get4Momentum().m()<<G4endl;
#endif
      G4LorentzVector b4M;
      G4LorentzVector r4M;
      G4bool evC=true;                  // @@ It makes only one attempt to be possible
      G4int bPDG=0;
      G4int rPDG=0;
      //G4double bM = 0.;               // Prototype of Real Mass of the EvaporatedDibaryon
      G4double rM = 0.;                 // Prototype of Real Mass of the residual nucleus
      G4int bB=0;                       // Proto of Baryon Number of the evaporated baryon
      G4int rB=0;                       // Proto of Baryon Number of the residual nucleus
      G4QHadron* bHadron = new G4QHadron;// Proto of the evaporated baryon @@where deleted?
      G4QHadron* rHadron = new G4QHadron;// Proto of the residual nucleus @@where deleted?
      G4int evcn=0;
      //G4int evcm=27;
      G4int evcm=9;                     // Max numder of attempts to evaporate
      // @@ Does not look like it can make two attempts @@ Improve, Simplify @@
      while(evC&&evcn<evcm)
      {
        evC=true;
        evcn++;
        if(!qNuc.EvaporateBaryon(bHadron,rHadron)) // Evaporation did not succeed
        {
#ifdef debug
          G4cout<<"***G4QNuc::EvaNuc:***EVA Failed***PDG="<<thePDG<<",M="<<totMass<<G4endl;
#endif
          delete bHadron;
          delete rHadron;
#ifdef debug
          G4cout<<"***G4QNucl::EvaNuc: Residual="<<qH->GetQC()<<qH->Get4Momentum()<<G4endl;
#endif
          evaHV->push_back(qH);               // fill AsItIs
          return;
        }
        evC=false;
        b4M=bHadron->Get4Momentum();
        r4M=rHadron->Get4Momentum();
        //bM   = b4M.m();                       // Real mass of the evaporated dibaryon
        rM   = r4M.m();                       // Real mass of the residual nucleus
        bB=bHadron->GetBaryonNumber();        // Baryon number of the evaporated baryon
        rB=rHadron->GetBaryonNumber();        // Baryon number of the residual nucleus
        bPDG=bHadron->GetPDGCode();
        rPDG=rHadron->GetPDGCode();
#ifdef debug
        G4int bC=bHadron->GetCharge();        // Baryon number of the evaporated baryon
        //G4int rC=rHadron->GetCharge();       // Baryon number of the residual nucleus
        G4double bCB=qNuc.CoulombBarrier(bC,bB);
        //G4double rCB=qNuc.CoulombBarrier(rC,rB);
        G4cout<<"G4QNucl::EvaNuc:Attempt #"<<evcn<<" > "<<evcm<<", rPDG="<<rPDG<<", bPDG="
              <<bPDG<<", bE="<<b4M.e()-b4M.m()<<" > bCB="<<bCB<<G4endl;
#endif
        //if(b4M.e()-b4M.m()<bCB&&evcn<evcm) evC=true;
      }  // End of while
#ifdef debug
      G4cout<<"G4QNucl::EvaNuc:*** EVA IS DONE *** F="<<bPDG<<b4M<<",bB="<<bB<<", ResNuc="
            <<rPDG<<r4M<<",rB="<<rB<<G4endl;
#endif
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (23) qH=0"<<G4endl;
#endif
      delete qH;
      if(bB<2) evaHV->push_back(bHadron);         // Fill EvaporatedBaryon (del.equivalent)
      else if(bB==2) DecayDibaryon(bHadron,evaHV);// => "Dibaryon" case needs decay
      else if(bB==4) evaHV->push_back(bHadron);   // "Alpha radiation" case (del.eq.)
      else if(bB==5) DecayAlphaBar(bHadron,evaHV);// "Alpha+Baryon Decay" case (del.equiv.)
      else if(bPDG==90004002) DecayAlphaDiN(bHadron,evaHV); // alph+2p(alph+2n is stable)
      else if(bPDG==90004004) DecayAlphaAlpha(bHadron,evaHV);// Alph+Alph Decay (del.eq.)
      else
      {
        delete bHadron;
        // G4cerr<<"***G4QNuc::EvaNuc:bB="<<bB<<">2 - unexpected evaporated fragment"<<G4endl;
        // throw G4QException("G4QNucleus::EvaporateNucleus: Wrong evaporation act");
        G4ExceptionDescription ed;
        ed << "Wrong evaporation act: EvaNuc:bB=" << bB
           << ">2 - unexpected evaporated fragment" << G4endl;
        G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0019",
                    FatalException, ed);
      }
      if(rB>2) EvaporateNucleus(rHadron,evaHV);    // Continue evaporation (@@ Self-call)
      else if(rB==2)                   // => "Dibaryon" case needs decay @@ DecayDibaryon
      {
        G4double rGSM = rHadron->GetQPDG().GetMass(); // Ground State mass of the dibaryon
#ifdef debug
        G4cout<<"G4QNuc::EvaNuc:ResidDibM="<<rM<<",GSM="<<rGSM<<",M-GSM="<<rM-rGSM<<G4endl;
#endif
        if(rM<=rGSM-0.01)
        {
          delete rHadron;
          // G4cerr<<"***G4QNucleus::EvaporNucl: <residual> M="<<rM<<" < GSM="<<rGSM<<G4endl;
          // throw G4QException("G4QNucleus::EvaporateNucleus: Evaporation below MassShell");
          G4ExceptionDescription ed;
          ed << "Evaporation below MassShell: <residual> M=" << rM << " < GSM="
             << rGSM << G4endl;
          G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0020",
                      FatalException, ed);
        }
        else if(fabs(rM-rGSM)<0.01&&rPDG==90001001) evaHV->push_back(rHadron); // (DE)
        else DecayDibaryon(rHadron,evaHV);   // => "Dibaryon Decay" case (del.equivalent)
      }
      else if(rB==5) DecayAlphaBar(rHadron,evaHV);// "Alpha+Baryon Decay" case (del.equiv.)
      else if(rPDG==90004002) DecayAlphaDiN(rHadron,evaHV);//alph+2p (alph+2n is stable)
      else if(rPDG==90004004) DecayAlphaAlpha(rHadron,evaHV);//Alpha+Alpha Decay (delEq)
      else evaHV->push_back(rHadron);        // Fill ResidNucl=Baryon to OutputHadronVector
    } // End of Evaporation of excited system
#ifdef debug
    G4cout<<"G4QNucleus::EvaporateNucleus: === End of the evaporation attempt"<<G4endl;
#endif
  }
  else                                          // => "Decay if impossible evaporate" case
  {
#ifdef debug
    G4cout<<"*G4QNucleus::EvaporateNucleus: InputHadron4M="<<q4M<<", PDG="<<thePDG<<G4endl;
#endif
    if(thePDG)
    {
      if(thePDG==10)                            // "Chipolino decay" case 
      {
        G4QContent totQC = qH->GetQC();         // Quark content of the hadron
        G4QChipolino resChip(totQC);            // define the Residual as a Chipolino
        G4QPDGCode h1=resChip.GetQPDG1();
        G4double m1  =h1.GetMass();             // Mass of the first hadron
        G4QPDGCode h2=resChip.GetQPDG2();
        G4double m2_value  =h2.GetMass();       // Mass of the second hadron
        if(totMass+.0001>m1+m2_value)
        {
#ifdef qdebug
          if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (24) qH=0"<<G4endl;
#endif
          delete qH;                            // Chipolino should not be in a sequence
          G4LorentzVector fq4M(0.,0.,0.,m1);
          G4LorentzVector qe4M(0.,0.,0.,m2_value);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
          {
            // G4cerr<<"***G4QNuc::EvaNuc:tM="<<totMass<<"-> h1M="<<m1<<" + h2M="<<m2_value<<G4endl;
            // throw G4QException("G4QNucleus::EvaporateNucleus: Chipol->h1+h2 DecIn2 error");
            G4ExceptionDescription ed;
            ed << "Chipol->h1+h2 DecIn2 error: tM=" << totMass << "-> h1M="
               << m1 <<" + h2M=" << m2_value << G4endl;
            G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0021",
                        FatalException, ed);
          }
          G4QHadron* H2 = new G4QHadron(h2.GetPDGCode(),qe4M);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateNucleus:(2) h2="<<h2.GetPDGCode()<<qe4M<<G4endl;
#endif
          evaHV->push_back(H2);            // (delete equivalent)
          G4QHadron* H1 = new G4QHadron(h1.GetPDGCode(),fq4M);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateNucleus:(2) h1="<<h1.GetPDGCode()<<fq4M<<G4endl;
#endif
          evaHV->push_back(H1);            // (delete equivalent)
        }
        else
        {
#ifdef qdebug
          if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (25) qH=0"<<G4endl;
#endif
          delete qH;
          // G4cerr<<"**G4QN::EN:M="<<totMass<<"<"<<m1<<"+"<<m2_value<<",d="<<m1+m2_value-totMass<<G4endl;
          // throw G4QException("G4QNucleus::EvaporateNucleus: Chipolino is under MassShell");
          G4ExceptionDescription ed;
          ed << "Chipolino is under MassShell: M=" << totMass << "<" << m1
             << "+" << m2_value << ",d=" << m1+m2_value-totMass << G4endl;
          G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0022",
                      FatalException, ed);
        }
      }
      else                                      // "Hadron" case
      {
        G4double totM=G4QPDGCode(thePDG).GetMass();
        if(fabs(totMass-totM)<0.001||abs(thePDG)-10*(abs(thePDG)/10)>2)
        {
#ifdef debug
          G4cout<<"**G4QNuc::EvaNuc:EmerFill(2) "<<qH->GetQC()<<qH->Get4Momentum()<<G4endl;
#endif
          evaHV->push_back(qH);
        }
        else if ((thePDG==221||thePDG==331)&&totMass>mPi+mPi) // "Decay in pipi" case
        {
#ifdef qdebug
          if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (26) qH=0"<<G4endl;
#endif
          delete qH;
          G4LorentzVector fq4M(0.,0.,0.,mPi);
          G4LorentzVector qe4M(0.,0.,0.,mPi);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
          {
            // G4cerr<<"***G4QNucleus::EvaporateNucleus:tM="<<totMass<<"-> pi+ + pi-"<<G4endl;
            // throw G4QException("G4QNucleus::EvaporateNucleus: H->Pi+Pi DecayIn2 error");
            G4ExceptionDescription ed;
            ed << "H->Pi+Pi DecayIn2 error: tM=" << totMass << "-> pi+ + pi-"
               << G4endl;
            G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0023",
                        FatalException, ed);
          }
          G4QHadron* H1 = new G4QHadron(211,fq4M);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateNucleus:(3) PiPlus="<<fq4M<<G4endl;
#endif
          evaHV->push_back(H1);            // (delete equivalent)
          G4QHadron* H2 = new G4QHadron(-211,qe4M);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateNucleus:(3) PiMinus="<<qe4M<<G4endl;
#endif
          evaHV->push_back(H2);            // (delete equivalent)
        }
        else if ((thePDG==221||thePDG==331)&&totMass>mPi0+mPi0) // "Decay in 2pi0" case
        {
#ifdef qdebug
          if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (27) qH=0"<<G4endl;
#endif
          delete qH;
          G4LorentzVector fq4M(0.,0.,0.,mPi0);
          G4LorentzVector qe4M(0.,0.,0.,mPi0);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
          {
            // G4cerr<<"***G4QNucleus::EvaporateNucleus:tM="<<totMass<<"-> pi0 + pi0"<<G4endl;
            // throw G4QException("G4QNucleus::EvaporateNucleus: H->Pi+Pi DecayIn2 error");
            G4ExceptionDescription ed;
            ed << "H->Pi+Pi DecayIn2 error: tM=" << totMass << "-> pi0 + pi0"
               << G4endl;
            G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0024",
                        FatalException, ed);
          }
          G4QHadron* H1 = new G4QHadron(111,fq4M);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateNucleus:(4) Pi01="<<fq4M<<G4endl;
#endif
          evaHV->push_back(H1);            // (delete equivalent)
          G4QHadron* H2 = new G4QHadron(111,qe4M);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateNucleus:(4) Pi02="<<qe4M<<G4endl;
#endif
          evaHV->push_back(H2);            // (delete equivalent)
        }
        else if (totMass>totM)                  // "Radiative Hadron decay" case
        {
#ifdef qdebug
          if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (28) qH=0"<<G4endl;
#endif
          delete qH;
          G4LorentzVector fq4M(0.,0.,0.,0.);
          G4LorentzVector qe4M(0.,0.,0.,totM);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
          {
            // G4cerr<<"***G4QNuc::EvaporateNuc:tM="<<totMass<<"->h1M="<<totM<<"+gam"<<G4endl;
            // throw G4QException("G4QNucleus::EvaporateNucleus: H*->H+gamma DecIn2 error");
            G4ExceptionDescription ed;
            ed << "H*->H+gamma DecIn2 error: tM=" << totMass << "->h1M="
               << totM << "+gam" << G4endl;
            G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0025",
                        FatalException, ed);
          }
          G4QHadron* H2 = new G4QHadron(thePDG,qe4M);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateNucleus:(5) tot="<<thePDG<<qe4M<<G4endl;
#endif
          evaHV->push_back(H2);            // (delete equivalent)
          G4QHadron* H1 = new G4QHadron(22,fq4M);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateNucleus:(5) GamFortot="<<fq4M<<G4endl;
#endif
          evaHV->push_back(H1);            // (delete equivalent)
        }
        else if (thePDG==111||thePDG==221||thePDG==331) // "Gamma+Gamma decay" case
        {
#ifdef qdebug
          if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (29) qH=0"<<G4endl;
#endif
          delete qH;
          G4LorentzVector fq4M(0.,0.,0.,0.);
          G4LorentzVector qe4M(0.,0.,0.,0.);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
          {
            // G4cerr<<"***G4QNucl::EvaporateNucleus:tM="<<totMass<<"->gamma + gamma"<<G4endl;
            // throw G4QException("G4QNucleus::EvaporateNucleus: pi/eta->g+g DecIn2 error");
            G4ExceptionDescription ed;
            ed << "pi/eta->g+g DecIn2 error: tM=" << totMass
               << "->gamma + gamma" << G4endl;
            G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0026",
                        FatalException, ed);
          }
          G4QHadron* H2 = new G4QHadron(22,qe4M);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateNucleus:(6) gam1="<<qe4M<<G4endl;
#endif
          evaHV->push_back(H2);            // (delete equivalent)
          G4QHadron* H1 = new G4QHadron(22,fq4M);
#ifdef debug
          G4cout<<"G4QNucleus::EvaporateNucleus:(6) gam2="<<fq4M<<G4endl;
#endif
          evaHV->push_back(H1);            // (delete equivalent)
        }
        else
        {
#ifdef qdebug
          if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (30) qH=0"<<G4endl;
#endif
          delete qH;
          // G4cerr<<"***G4QNucl::EvaNuc: Nuc="<<thePDG<<theQC<<", q4M="<<q4M<<", M="<<totMass
          //       <<" < GSM="<<totM<<", 2Pi="<<mPi+mPi<<", 2Pi0="<<mPi0+mPi0<<G4endl;
          // throw G4QException("G4QNucleus::EvaporateNucleus: Hadron is under MassShell");
          G4ExceptionDescription ed;
          ed << "Hadron is under MassShell: Nuc=" << thePDG << theQC
             << ", q4M=" << q4M << ", M=" << totMass <<" < GSM=" << totM
             <<", 2Pi=" << mPi+mPi << ", 2Pi0=" << mPi0+mPi0 << G4endl;
          G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0027",
                      FatalException, ed);
        }
      }
    }
    else
    {
#ifdef qdebug
      if(!qH) G4cout<<"-Warning-G4QNucleus::EvaporateNucleus: (31) qH=0"<<G4endl;
#endif
      delete qH;
      // G4cerr<<"**G4QNuc::EvaNuc:RN="<<thePDG<<theQC<<",q4M="<<q4M<<",qM="<<totMass<<G4endl;
      // throw G4QException("G4QNucleus::EvaporateNucleus: This is not aNucleus nor aHadron");
      G4ExceptionDescription ed;
      ed << "This is not aNucleus nor aHadron: RN=" << thePDG << theQC
         << ",q4M=" << q4M <<",qM=" << totMass << G4endl;
      G4Exception("G4QNucleus::EvaporateNucleus()", "HAD_CHPS_0028",
                  FatalException, ed);
    }
  }
#ifdef qdebug
  if (qH)
  {
    G4cout<<"G4QNucleus::EvaporateNucleus: deletedAtEnd, PDG="<<qH->GetPDGCode()<<G4endl;
    if(!qH) G4cout<<"G4QNucleus::EvaporateNucleus: (20) qH="<<G4endl;
    else delete qH;
  }
#endif
#ifdef debug
  G4cout<<"G4QNucleus::EvaporateNucleus: =---=>> End. "<<G4endl;
#endif
  return;
} // End of EvaporateNucleus

//Unavoidable decay of the Isonucleus in nP+(Pi+) or nN+(Pi-)
void G4QNucleus::DecayIsonucleus(G4QHadron* qH, G4QHadronVector* evaHV)
{
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  //static const G4double mSigZ= G4QPDGCode(3212).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double eps  = 0.003;
  G4LorentzVector q4M = qH->Get4Momentum();      // Get 4-momentum of the Isonucleus
  G4double qM=q4M.m();                           // Real mass of the Isonucleus
  G4QContent qQC = qH->GetQC();                  // Get QuarcContent of the Isonucleus
  G4int qBN=qQC.GetBaryonNumber();               // Baryon number of the IsoNucleus
  G4int qC=qQC.GetCharge();                      // Charge of the IsoNucleus
  G4int qS=qQC.GetStrangeness();                 // Strangness of the IsoNucleus
#ifdef debug
  G4cout<<"G4QNuc:DecayIson:QC="<<qQC<<",M="<<qM<<",B="<<qBN<<",S="<<qS<<",C="<<qC<<G4endl;
#endif
  if(qS<0||qS>qBN)                               // *** Should not be here ***
  {
    G4cout<<"--Warning(Upgrade)--G4QNuc::DecIsonuc:FillAsIs,4M="<<q4M<<",QC="<<qQC<<G4endl;
    evaHV->push_back(qH);                        // fill as it is (delete equivalent)
    return;
  }
  G4int          qPN=qC-qBN;                     // Number of pions in the Isonucleus
  G4int          fPDG = 2212;                    // Prototype for nP+(Pi+) case
  G4int          sPDG = 211;
  G4int          tPDG = 3122;                    // @@ Sigma0 (?)
  G4double       fMass= mProt;
  G4double       sMass= mPi;
  G4double       tMass= mLamb;                   // @@ Sigma0 (?)
  // =---------= Negative state =-----------=
  if(qC<0)                                       // =-------= Only Pi- can help
  {
    if(qS&&qBN==qS)                              // --- n*Lamb + k*(Pi-) State ---
    {
      sPDG = -211;
      if(-qC==qS && qS==1)                       // Only one Sigma- like (qBN=1)
      {
        if(fabs(qM-mSigM)<eps)
        {
          evaHV->push_back(qH);                  // Fill Sigma- as it is
          return;
        }
        else if(qM>mLamb+mPi)                    //(2) Sigma- => Lambda + Pi- decay
        {
          fPDG = 3122;
          fMass= mLamb;
        }
        else if(qM>mSigM)                        //(2) Sigma+ => Sigma+ + gamma decay
        {
          fPDG = 3112;
          fMass= mSigM;
          sPDG = 22;
          sMass= 0.;
        }
        else                                      //(2) Sigma- => Neutron + Pi- decay
        {
          fPDG = 2112;
          fMass= mNeut;
        }
        qPN  = 1;                                 // #of (Pi+ or gamma)'s = 1
      }
      else if(-qC==qS)                            //(2) a few Sigma- like
      {
        qPN  = 1;                                 // One separated Sigma-
        fPDG = 3112;
        sPDG = 3112;
        sMass= mSigM;
        qBN--;
        fMass= mSigM;
      }
      else if(-qC>qS)                             //(2) n*(Sigma-)+m*(Pi-)
      {
        qPN  = -qC-qS;                            // #of Pi-'s
        fPDG = 3112;
        fMass= mSigM;
      }
      else                                        //(2) n*(Sigma-)+m*Lambda (-qC<qS)
      {
        qBN += qC;                                // #of Lambda's
        fPDG = 3122;
        fMass= mLamb;
        qPN  = -qC;                               // #of Sigma+'s
        sPDG = 3112;
        sMass= mSigM;
      }
      qS   = 0;                                   // Only decays in two are above
    }
    else if(qS)                                   // ->n*Lamb+m*Neut+k*(Pi-) State (qS<qBN)
    {
      qBN -= qS;                                  // #of neutrons
      fPDG = 2112;
      fMass= mNeut;
      G4int nPin = -qC;                           // #of Pi-'s                    
      if(qS==nPin)                                //(2) m*Neut+n*Sigma-
      {
        qPN  = qS;                                // #of Sigma-
        sPDG = 3112;
        sMass= mSigM;
        qS   = 0;
      }
      else if(qS>nPin)                            //(3) m*P+n*(Sigma+)+k*Lambda
      {
        qS-=nPin;                                 // #of Lambdas
        qPN  = nPin;                              // #of Sigma+
        sPDG = 3112;
        sMass= mSigM;
      }
      else                                        //(3) m*N+n*(Sigma-)+k*(Pi-) (qS<nPin)
      {
        qPN  = nPin-qS;                           // #of Pi-
        sPDG = -211;
        tPDG = 3112;
        tMass= mSigM;
      }
    }
    else                                          //(2) n*N+m*(Pi-)   (qS=0)
    {
      sPDG = -211;
      qPN  = -qC;
      fPDG = 2112;
      fMass= mNeut;
    }
  }
  else if(!qC)                                   // *** Should not be here ***
  {
    if(qS && qS<qBN)                             //(2) n*Lamb+m*N ***Should not be here***
    {
      qPN  = qS;
      fPDG = 2112;                               // mN+nL case
      sPDG = 3122;
      sMass= mLamb;
      qBN -= qS;
      fMass= mNeut;
      qS   = 0;
    }
    else if(qS>1 && qBN==qS)                     //(2) m*Lamb(m>1) ***Should not be here***
    {
      qPN  = 1;
      fPDG = 3122;
      sPDG = 3122;
      sMass= mLamb;
      qBN--;
      fMass= mLamb;
    }
    else if(!qS && qBN>1)                        //(2) n*Neut(n>1) ***Should not be here***
    {
      qPN  = 1;
      fPDG = 2112;
      sPDG = 2112;
      sMass= mNeut;
      qBN--;
      fMass= mNeut;
    }
    else G4cout<<"*?*G4QNuc::DecayIsonucleus: (1) QC="<<qQC<<G4endl;
  }
  else if(qC>0)                                  // n*Lamb+(m*P)+(k*Pi+)
  {
    if(qS && qS+qC==qBN)                         //(2) n*Lamb+m*P ***Should not be here***
    {
      qPN  = qS;
      qS   = 0;
      fPDG = 2212;
      sPDG = 3122;
      sMass= mLamb;
      qBN  = qC;
      fMass= mProt;
    }
    else if(qS  && qC<qBN-qS)                     //(3)n*L+m*P+k*N ***Should not be here***
    {
      qPN  = qC;                                  // #of protons
      fPDG = 2112;                                // mP+nL case
      sPDG = 2212;
      sMass= mProt;
      qBN -= qS+qC;                               // #of neutrons
      fMass= mNeut;
    }
    else if(qS  && qBN==qS)                       // ---> n*L+m*Pi+ State
    {
      if(qC==qS && qS==1)                         // Only one Sigma+ like State
      {
        if(fabs(qM-mSigP)<eps)                    // Fill Sigma+ as it is
        {
          evaHV->push_back(qH);
          return;
        }
        else if(qM>mLamb+mPi)                     //(2) Sigma+ => Lambda + Pi+ decay
        {
          fPDG = 3122;
          fMass= mLamb;
        }
        else if(qM>mNeut+mPi)                     //(2) Sigma+ => Neutron + Pi+ decay
        {
          fPDG = 2112;
          fMass= mNeut;
        }
        else if(qM>mSigP)                         //(2) Sigma+ => Sigma+ + gamma decay
        {
          fPDG = 3222;
          fMass= mSigP;
          sPDG = 22;
          sMass= 0.;
        }
        else                                      //(2) Sigma+ => Proton + gamma decay
        {
          fPDG = 2212;
          fMass= mProt;
          sPDG = 22;
          sMass= 0.;
        }
        qPN  = 1;                                 // #of (Pi+ or gamma)'s = 1
      }
      else if(qC==qS)                             //(2) a few Sigma+ like hyperons
      {
        qPN  = 1;
        fPDG = 3222;
        sPDG = 3222;
        sMass= mSigP;
        qBN--;
        fMass= mSigP;
      }
      else if(qC>qS)                              //(2) n*(Sigma+)+m*(Pi+)
      {
        qPN  = qC-qS;                             // #of Pi+'s
        fPDG = 3222;
        qBN  = qS;                                // #of Sigma+'s
        fMass= mSigP;
      }
      else                                        //(2) n*(Sigma+)+m*Lambda
      {
        qBN -= qC;                                // #of Lambda's
        fPDG = 3122;
        fMass= mLamb;
        qPN  = qC;                                // #of Sigma+'s
        sPDG = 3222;
        sMass= mSigP;
      }
      qS   = 0;                                   // All above are decays in 2
    }
    else if(qS && qC>qBN-qS)                      // n*Lamb+m*P+k*Pi+
    {
      qBN -= qS;                                  // #of protons
      G4int nPip = qC-qBN;                        // #of Pi+'s                    
      if(qS==nPip)                                //(2) m*P+n*Sigma+
      {
        qPN  = qS;                                // #of Sigma+
        sPDG = 3222;
        sMass= mSigP;
        qS   = 0;
      }
      else if(qS>nPip)                            //(3) m*P+n*(Sigma+)+k*Lambda
      {
        qS  -= nPip;                              // #of Lambdas
        qPN  = nPip;                              // #of Sigma+
        sPDG = 3222;
        sMass= mSigP;
      }
      else                                        //(3) m*P+n*(Sigma+)+k*(Pi+)
      {
        qPN  = nPip-qS;                           // #of Pi+
        tPDG = 3222;
        tMass= mSigP;
      }
    }
    if(qC<qBN)                                    //(2) n*P+m*N ***Should not be here***
    {
      fPDG = 2112;
      fMass= mNeut;
      qPN  = qC;
      sPDG = 2212;
      sMass= mProt;
    }
    else if(qBN==qC && qC>1)                     //(2) m*Prot(m>1) ***Should not be here***
    {
      qPN  = 1;
      fPDG = 2212;
      sPDG = 2212;
      sMass= mProt;
      qBN--;
      fMass= mProt;
    }
    else if(qC<=qBN||!qBN) G4cout<<"*?*G4QNuc::DecayIsonucleus: (2) QC="<<qQC<<G4endl;
    // !qS && qC>qBN                             //(2) Default condition n*P+m*(Pi+)
  }
  G4double tfM=qBN*fMass;
  G4double tsM=qPN*sMass;
  G4double ttM=0.;
  if(qS) ttM=qS*tMass;
  G4LorentzVector f4Mom(0.,0.,0.,tfM);
  G4LorentzVector s4Mom(0.,0.,0.,tsM);
  G4LorentzVector t4Mom(0.,0.,0.,ttM);
  G4double sum=tfM+tsM+ttM;
  if(fabs(qM-sum)<eps)
  {
    f4Mom=q4M*(tfM/sum);
    s4Mom=q4M*(tsM/sum);
    if(qS) t4Mom=q4M*(ttM/sum);
  }
  else if(!qS && (qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom)))
  {
#ifdef debug
    G4cout<<"***G4QNuc::DecIsonuc:fPDG="<<fPDG<<"*"<<qBN<<"(fM="<<fMass<<")+sPDG="<<sPDG
          <<"*"<<qPN<<"(sM="<<sMass<<")"<<"="<<sum<<" > TotM="<<qM<<q4M<<qQC<<qS<<G4endl;
#endif
    evaHV->push_back(qH);                  // fill as it is (delete equivalent)
    return;
  }
  else if(qS && (qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom)))
  {
#ifdef debug
    G4cout<<"***G4QNuc::DecIsonuc: "<<fPDG<<"*"<<qBN<<"("<<fMass<<")+"<<sPDG<<"*"<<qPN<<"("
          <<sMass<<")+Lamb*"<<qS<<"="<<sum<<" > TotM="<<qM<<q4M<<qQC<<G4endl;
#endif
    evaHV->push_back(qH);                  // fill as it is (delete equivalent)
    return;
  }
#ifdef debug
  G4cout<<"G4QNuc::DecayIsonucleus: *DONE* n="<<qPN<<f4Mom<<fPDG<<", m="<<qPN<<s4Mom<<sPDG
        <<", l="<<qS<<t4Mom<<G4endl;
#endif
  delete qH;
  if(qBN)
  {
    f4Mom/=qBN;
    for(G4int ih=0; ih<qBN; ih++)
    {
      G4QHadron* Hi = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the hyperon
      evaHV->push_back(Hi);                 // Fill "Hi" (delete equivalent)
    }
  }
  if(qPN)
  {
    s4Mom/=qPN;
    for(G4int ip=0; ip<qPN; ip++)
    {
      G4QHadron* Hj = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the meson
      evaHV->push_back(Hj);                 // Fill "Hj" (delete equivalent)
    }
  }
  if(qS)
  {
    t4Mom/=qS;
    for(G4int il=0; il<qS; il++)
    {
      G4QHadron* Hk = new G4QHadron(tPDG,t4Mom); // Create a Hadron for the lambda
      evaHV->push_back(Hk);                 // Fill "Hk" (delete equivalent)
    }
  }
#ifdef qdebug
  if (qH)
  {
    G4cout << "G4QNucleus::DecayIsonucleus: deleted at end - PDG: "
           << qH->GetPDGCode() << G4endl;
    delete qH;
  }
#endif
} // End of DecayIsonucleus

//Decay of the excited dibaryon in two baryons
void G4QNucleus::DecayDibaryon(G4QHadron* qH, G4QHadronVector* evaHV)
{
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
  static const G4double mKsiM= G4QPDGCode(3312).GetMass();
  static const G4double mKsiZ= G4QPDGCode(3322).GetMass();
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mPiN = mPi+mNeut;
  static const G4double mPiP = mPi+mProt;
  static const G4double dmPiN= mPiN+mPiN;
  static const G4double dmPiP= mPiP+mPiP;
  static const G4double nnPi = mNeut+mPiN;
  static const G4double ppPi = mProt+mPiP;
  static const G4double lnPi = mLamb+mPiN;
  static const G4double lpPi = mLamb+mPiP;
  static const G4double dNeut= mNeut+mNeut;
  static const G4double dProt= mProt+mProt;
  static const G4double dLamb= mLamb+mLamb;
  static const G4double dLaNe= mLamb+mNeut;
  static const G4double dLaPr= mLamb+mProt;
  static const G4double dSiPr= mSigP+mProt;
  static const G4double dSiNe= mSigM+mNeut;
  static const G4double dKsPr= mKsiZ+mProt;
  static const G4double dKsNe= mKsiM+mNeut;
  static const G4double eps  = 0.003;
  static const G4QNucleus vacuum(90000000);
  G4bool four=false;                           // defFALSE for 4-particle decay of diDelta
  G4LorentzVector q4M = qH->Get4Momentum();    // Get 4-momentum of the Dibaryon
  G4int          qPDG = qH->GetPDGCode();      // PDG Code of the decaying dybaryon
  G4double         qM = q4M.m();
  G4double         rM = qM+eps;                // Just to avoid the computer accuracy
#ifdef debug
  G4cout<<"G4QNucl::DecayDibaryon: *Called* PDG="<<qPDG<<",4Mom="<<q4M<<", M="<<qM<<G4endl;
#endif
  // Select a chanel of the dibaryon decay (including Delta+Delta-> 4 particle decay
  G4int          fPDG = 2212;                  // Prototype for pp case
  G4int          sPDG = 2212;
  G4int          tPDG = 0;                     // Zero prototype to separate 3 from 2 
  G4double       fMass= mProt;
  G4double       sMass= mProt;
  G4double       tMass= mPi;
  if     (qPDG==90003998 && rM>=dmPiP)         // "diDelta++" case
  {
    sPDG = 211;
    sMass= mPi;
    four = true;
  }
  else if(qPDG==89998004 && rM>=dmPiN)         // "diDelta--" case
  {
    sPDG = -211;
    fPDG = 2112;
    sMass= mPi;
    fMass= mNeut;
    four = true;
  }
  else if(qPDG==90000002 && rM>=dNeut)         // "dineutron" case
  {
    fPDG = 2112;
    sPDG = 2112;
    fMass= mNeut;
    sMass= mNeut;    
  }
  else if(qPDG==90001001 && rM>=mDeut)         // "exited deutron" case
  {
    if(fabs(qM-mDeut)<eps)
    {
      evaHV->push_back(qH);               // Fill as it is (delete equivalent)
      return;
    }
    else if(mProt+mNeut<rM)
    {
      fPDG = 2112;
      fMass= mNeut;    
    }
    else
    {
      fPDG = 22;
      sPDG = 90001001;
      fMass= 0.;
      sMass= mDeut;    
    }
  }
  else if(qPDG==91000001 && rM>=dLaNe)         // "Lambda-neutron" case
  {
    fPDG = 2112;
    sPDG = 3122;
    fMass= mNeut;
    sMass= mLamb;    
  }
  else if(qPDG==91001000 && rM>=dLaPr)         // "Lambda-proton" case
  {
    sPDG = 3122;
    sMass= mLamb;    
  }
  else if(qPDG==89999003 && rM>=nnPi)         // "neutron/Delta-" case
  {
    fPDG = 2112;
    sPDG = 2112;
    tPDG = -211;
    fMass= mNeut;
    sMass= mNeut;    
  }
  else if(qPDG==90002999 && rM>=ppPi)         // "proton/Delta++" case
  {
    tPDG = 211;
  }
  else if(qPDG==90999002 && rM>=lnPi)         // "lambda/Delta-" case
  {
    fPDG = 2112;
    sPDG = 3122;
    tPDG = -211;
    fMass= mNeut;
    sMass= mLamb;    
  }
  else if(qPDG==91001999 && rM>=lpPi)         // "lambda/Delta+" case
  {
    sPDG = 3122;
    tPDG = 211;
    sMass= mLamb;    
  }
  else if(qPDG==90999002 && rM>=dSiNe)         // "Sigma-/neutron" case
  {
    fPDG = 2112;
    sPDG = 3112;
    fMass= mNeut;
    sMass= mSigM;    
  }
  else if(qPDG==91001999 && rM>=dSiPr)         // "Sigma+/proton" case
  {
    sPDG = 3222;
    sMass= mSigP;    
  }
  else if(qPDG==92000000 && rM>=dLamb)         // "diLambda" case
  {
    fPDG = 3122;
    sPDG = 3122;
    fMass= mLamb;
    sMass= mLamb;    
  }
  else if(qPDG==91999001 && rM>=dKsNe)         // "Ksi-/neutron" case
  {
    fPDG = 2112;
    sPDG = 3312;
    fMass= mNeut;
    sMass= mKsiM;    
  }
  else if(qPDG==92000999 && rM>=dKsPr)         // "Ksi0/proton" case
  {
    sPDG = 3322;
    sMass= mKsiZ;    
  }
  else if(qPDG!=90002000|| rM<dProt)           // Other possibilities (if not a default)
  {
    G4int qS = qH->GetStrangeness();
    G4int qB = qH->GetBaryonNumber();
    if(qB>0&&qS<0)                             // Antistrange diBarion
    {
      DecayAntiStrange(qH,evaHV);
      return;
    }
    else
    {
      delete qH;
      G4cerr<<"***G4QN::DecDiBar: badPDG="<<qPDG<<" or smallM="<<qM<<",2mP="<<dProt
            <<",2mN="<<dNeut<<G4endl;
      // @@ Nothing to do. Just 2 GeV disappears... Very rare! Just to avoid the exception.
      //throw G4QException("G4QNucleus::DecayDibar: Unknown PDG code or small Mass of DB");
    }
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,sMass);
  G4LorentzVector t4Mom(0.,0.,0.,tMass);
  if(!tPDG&&!four)
  {
    G4double sum=fMass+sMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cout<<"---Warning---G4QN::DecDib:fPDG="<<fPDG<<"(M="<<fMass<<")+sPDG="<<sPDG<<"(M="
            <<sMass<<")="<<sum<<" >? TotM="<<q4M.m()<<q4M<<G4endl;
      //G4cerr<<"***G4QNucl::DecayDiBar: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("***G4QNucleus::DecayDibaryon: DiBaryon DecayIn2 error");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNucleus::DecayDibaryon:(2) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
          <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                 // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                 // Fill "H2" (delete equivalent)
  }
  else if(four)
  {
    q4M=q4M/2.;                                // Divided in 2 !!!
    qM/=2.;                                    // Divide the mass in 2 !
    G4double sum=fMass+sMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cout<<"---Warning---G4QN::DecDib:fPDG="<<fPDG<<"(M="<<fMass<<")+sPDG="<<sPDG<<"(M="
            <<sMass<<")"<<"="<<sum<<">tM="<<q4M.m()<<q4M<<G4endl;
      //G4cerr<<"***G4QN::DecayDibaryon: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("***G4QNucleus::DecDibaryon: Dibaryon DecayIn2 error");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNucleus::DecayDibaryon:(3) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
          <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                 // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                 // Fill "H2" (delete equivalent)
    // Now the second pair mus be decayed
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(!G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      // Should not be here as sum was already compared with qM above for the first delta
      delete qH;
      // G4cerr<<"***G4QNucl::DecDibar:fPDG="<<fPDG<<"(fM="<<fMass<<") + sPDG="<<sPDG<<"(sM="
      //       <<sMass<<")="<<sum<<" >? (DD2,Can't be here) TotM="<<q4M.m()<<q4M<<G4endl;
      // throw G4QException("G4QNucleus::DecayDibaryon: General DecayIn2 error");
      G4ExceptionDescription ed;
      ed << "General DecayIn2 error: fPDG=" << fPDG << "(fM=" << fMass
         << ") + sPDG=" << sPDG <<"(sM=" << sMass << ")=" << sum
         << " >? (DD2,Can't be here) TotM=" << q4M.m() << q4M << G4endl;
      G4Exception("G4QNucleus::DecayDibaryon()", "HAD_CHPS_0000",
                  FatalException, ed);
    }
#ifdef debug
    G4cout<<"G4QNucl::DecayDibaryon:(4) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
          <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    G4QHadron* H3 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    evaHV->push_back(H3);                 // Fill "H1" (delete equivalent)
    G4QHadron* H4 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H4);                 // Fill "H2" (delete equivalent)
    delete qH;
  }
  else
  {
    G4double sum=fMass+sMass+tMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
      t4Mom=q4M*(tMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
    {
      G4cout<<"---Warning---G4QN::DecDib:fPDG="<<fPDG<<"(M="<<fMass<<")+sPDG="<<sPDG<<"(M="
            <<sMass<<")+tPDG="<<tPDG<<"(tM="<<tMass<<")="<<sum<<">TotM="<<q4M.m()<<G4endl;
      //G4cerr<<"***G4QNuc::DecayDibaryon:qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("G4QNucleus::DecayDibaryon: diBar DecayIn3 error");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNuc::DecayDibaryon:(5) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG<<", s4M="<<s4Mom
          <<",sPDG="<<sPDG<<", t4M="<<t4Mom<<",tPDG="<<tPDG<<G4endl;
#endif
    //qH->SetNFragments(2);                    // Fill a#of fragments to decaying Dibaryon
    //evaHV->push_back(qH);               // Fill hadron with nf=2 (delete equivalent)
    // Instead
    delete qH;
    //
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                 // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                 // Fill "H2" (delete equivalent)
    G4QHadron* H3 = new G4QHadron(tPDG,t4Mom); // Create a Hadron for the meson
    evaHV->push_back(H3);                 // Fill "H3" (delete equivalent)
  }
#ifdef qdebug
  if (qH)
  {
    G4cout << "G4QNucleus::DecayDiBaryon: deleted at end - PDG: "
           << qH->GetPDGCode() << G4endl;
    delete qH;
  }
#endif
} // End of DecayDibaryon

//Decay of the excited anti-dibaryon in two anti-baryons
void G4QNucleus::DecayAntiDibaryon(G4QHadron* qH, G4QHadronVector* evaHV)
{
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
  static const G4double mKsiM= G4QPDGCode(3312).GetMass();
  static const G4double mKsiZ= G4QPDGCode(3322).GetMass();
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mPiN = mPi+mNeut;
  static const G4double mPiP = mPi+mProt;
  static const G4double dmPiN= mPiN+mPiN;
  static const G4double dmPiP= mPiP+mPiP;
  static const G4double nnPi = mNeut+mPiN;
  static const G4double ppPi = mProt+mPiP;
  static const G4double lnPi = mLamb+mPiN;
  static const G4double lpPi = mLamb+mPiP;
  static const G4double dNeut= mNeut+mNeut;
  static const G4double dProt= mProt+mProt;
  static const G4double dLamb= mLamb+mLamb;
  static const G4double dLaNe= mLamb+mNeut;
  static const G4double dLaPr= mLamb+mProt;
  static const G4double dSiPr= mSigP+mProt;
  static const G4double dSiNe= mSigM+mNeut;
  static const G4double dKsPr= mKsiZ+mProt;
  static const G4double dKsNe= mKsiM+mNeut;
  static const G4double eps  = 0.003;
  static const G4QNucleus vacuum(90000000);
  G4bool four=false;                           // defFALSE for 4-particle decay of diDelta
  G4LorentzVector q4M = qH->Get4Momentum();    // Get 4-momentum of the Dibaryon
  G4int          qPDG = qH->GetPDGCode();      // PDG Code of the decaying dybaryon
  G4double         qM = q4M.m();               // Mass of the decaying anti-di-baryon
  G4double         rM = qM+eps;                // Just to avoid the computer accuracy
#ifdef debug
  G4cout<<"G4QNucl::DecayAntiDibar:*Called* PDG="<<qPDG<<",4Mom="<<q4M<<", M="<<qM<<G4endl;
#endif
  // Select a chanel of the dibaryon decay (including Delta+Delta-> 4 particle decay
  G4int          fPDG = -2212;                 // Prototype for anti-pp case
  G4int          sPDG = -2212;
  G4int          tPDG = 0;                     // Zero prototype to separate 3 from 2 
  G4double       fMass= mProt;
  G4double       sMass= mProt;
  G4double       tMass= mPi;
  if     (qPDG==89996002 && rM>=dmPiP)         // "anti-diDelta++" case
  {
    sPDG = -211;
    sMass= mPi;
    four = true;
  }
  else if(qPDG==90001996 && rM>=dmPiN)         // "diDelta--" case
  {
    sPDG = 211;
    fPDG = -2112;
    sMass= mPi;
    fMass= mNeut;
    four = true;
  }
  else if(qPDG==89999998 && rM>=dNeut)         // "dineutron" case
  {
    fPDG = -2112;
    sPDG = -2112;
    fMass= mNeut;
    sMass= mNeut;    
  }
  else if(qPDG==89998999 && rM>=mDeut)         // "exited deutron" case
  {
    if(fabs(qM-mDeut)<eps)
    {
      evaHV->push_back(qH);                    // Fill as it is (delete equivalent)
      return;
    }
    else if(mProt+mNeut<rM)
    {
      fPDG = -2112;
      fMass= mNeut;    
    }
    else
    {
      fPDG = 22;
      sPDG = 89998999;                         // Anti-deuteron
      fMass= 0.;
      sMass= mDeut;    
      G4cout<<"--Warning--G4QNucl::DecayAntiDibar:ANTI-DEUTERON is created M="<<rM<<G4endl;
    }
  }
  else if(qPDG==88999999 && rM>=dLaNe)         // "Lambda-neutron" case
  {
    fPDG = -2112;
    sPDG = -3122;
    fMass= mNeut;
    sMass= mLamb;    
  }
  else if(qPDG==88999999 && rM>=dLaPr)         // "Lambda-proton" case
  {
    sPDG = -3122;
    sMass= mLamb;    
  }
  else if(qPDG==90000997 && rM>=nnPi)         // "neutron/Delta-" case
  {
    fPDG = -2112;
    sPDG = -2112;
    tPDG = 211;
    fMass= mNeut;
    sMass= mNeut;    
  }
  else if(qPDG==89997001 && rM>=ppPi)         // "proton/Delta++" case
  {
    tPDG = -211;
  }
  else if(qPDG==89000998 && rM>=lnPi)         // "lambda/Delta-" case
  {
    fPDG = -2112;
    sPDG = -3122;
    tPDG = 211;
    fMass= mNeut;
    sMass= mLamb;    
  }
  else if(qPDG==889998001 && rM>=lpPi)         // "lambda/Delta+" case
  {
    sPDG = -3122;
    tPDG = -211;
    sMass= mLamb;    
  }
  else if(qPDG==89000998 && rM>=dSiNe)         // "Sigma-/neutron" case
  {
    fPDG = -2112;
    sPDG = -3112;
    fMass= mNeut;
    sMass= mSigM;    
  }
  else if(qPDG==88998001 && rM>=dSiPr)         // "Sigma+/proton" case
  {
    sPDG = -3222;
    sMass= mSigP;    
  }
  else if(qPDG==88000000 && rM>=dLamb)         // "diLambda" case
  {
    fPDG = -3122;
    sPDG = -3122;
    fMass= mLamb;
    sMass= mLamb;    
  }
  else if(qPDG==88000999 && rM>=dKsNe)         // "Ksi-/neutron" case
  {
    fPDG = -2112;
    sPDG = -3312;
    fMass= mNeut;
    sMass= mKsiM;    
  }
  else if(qPDG==87999001 && rM>=dKsPr)         // "Ksi0/proton" case
  {
    sPDG = -3322;
    sMass= mKsiZ;    
  }
  else if(qPDG!=89998000|| rM<dProt)           // Other possibilities (if not a default)
  {
    G4int qS = qH->GetStrangeness();
    G4int qB = qH->GetBaryonNumber();
    if(qB>0&&qS<0)                             // Antistrange diBarion
    {
      DecayAntiStrange(qH,evaHV);
      return;
    }
    else
    {
      delete qH;
      G4cerr<<"**G4QNuc::DecayAntiDiBar: badPDG="<<qPDG<<" or smallM="<<qM<<", 2mP="<<dProt
            <<", 2mN="<<dNeut<<G4endl;
      // @@ Nothing to do. Just 2 GeV disappears... Very rare! Just to avoid the exception.
      //throw G4QException("G4QNucleus::DecayDibar: Unknown PDG code or small Mass of DB");
    }
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,sMass);
  G4LorentzVector t4Mom(0.,0.,0.,tMass);
  if(!tPDG&&!four)
  {
    G4double sum=fMass+sMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cout<<"---Warning---G4QN::DecAntiDib:fPDG="<<fPDG<<"(M="<<fMass<<")+sPDG="<<sPDG
            <<"(M="<<sMass<<")="<<sum<<" >? TotM="<<q4M.m()<<q4M<<G4endl;
      //G4cerr<<"***G4QNucl::DecayDiBar: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("***G4QNucleus::DecayDibaryon: DiBaryon DecayIn2 error");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNucleus::DecayAntiDibaryon:(2) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
          <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                 // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                 // Fill "H2" (delete equivalent)
  }
  else if(four)
  {
    q4M=q4M/2.;                                // Divided in 2 !!!
    qM/=2.;                                    // Divide the mass in 2 !
    G4double sum=fMass+sMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cout<<"---Warning---G4QN::DecAntiDib:fPDG="<<fPDG<<"(M="<<fMass<<")+sPDG="<<sPDG
            <<"(M="<<sMass<<")"<<"="<<sum<<">tM="<<q4M.m()<<q4M<<G4endl;
      //G4cerr<<"***G4QN::DecayDibaryon: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("***G4QNucleus::DecDibaryon: Dibaryon DecayIn2 error");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNucleus::DecayAntiDibaryon:(3) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
          <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                      // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                      // Fill "H2" (delete equivalent)
    // Now the second pair mus be decayed
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(!G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      // Should not be here as sum was already compared with qM above for the first delta
      delete qH;
      // G4cerr<<"**G4QNucl::DecAntiDibar:fPDG="<<fPDG<<"(fM="<<fMass<<")+sPDG="<<sPDG<<"(sM="
      //       <<sMass<<")="<<sum<<" >? (DD2,Can't be here) TotM="<<q4M.m()<<q4M<<G4endl;
      // throw G4QException("G4QNucleus::DecayAntiDibaryon: General DecayIn2 error");
      G4ExceptionDescription ed;
      ed << " General DecayIn2 error: fPDG=" << fPDG << "(fM=" << fMass
         << ")+sPDG=" << sPDG << "(sM=" << sMass << ")=" << sum
         << " >? (DD2,Can't be here) TotM=" << q4M.m() << q4M << G4endl;
      G4Exception("G4QNucleus::DecayAntiDibaryon()", "HAD_CHPS_0000",
                  FatalException, ed);
    }
#ifdef debug
    G4cout<<"G4QNucl::DecayAntiDibaryon:(4) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
          <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    G4QHadron* H3 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    evaHV->push_back(H3);                      // Fill "H1" (delete equivalent)
    G4QHadron* H4 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H4);                      // Fill "H2" (delete equivalent)
    delete qH;
  }
  else
  {
    G4double sum=fMass+sMass+tMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
      t4Mom=q4M*(tMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
    {
      G4cout<<"-Warning-G4QN::DecAntiDib:fPDG="<<fPDG<<"(M="<<fMass<<")+sPDG="<<sPDG<<"(M="
            <<sMass<<")+tPDG="<<tPDG<<"(tM="<<tMass<<")="<<sum<<">TotM="<<q4M.m()<<G4endl;
      //G4cerr<<"***G4QNuc::DecayDibaryon:qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("G4QNucleus::DecayDibaryon: diBar DecayIn3 error");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNuc::DecayAbtiDibaryon:(5) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG<<", s4M="
          <<s4Mom<<",sPDG="<<sPDG<<", t4M="<<t4Mom<<",tPDG="<<tPDG<<G4endl;
#endif
    //qH->SetNFragments(2);                    // Fill a#of fragments to decaying Dibaryon
    //evaHV->push_back(qH);               // Fill hadron with nf=2 (delete equivalent)
    // Instead
    delete qH;
    //
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                 // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                 // Fill "H2" (delete equivalent)
    G4QHadron* H3 = new G4QHadron(tPDG,t4Mom); // Create a Hadron for the meson
    evaHV->push_back(H3);                 // Fill "H3" (delete equivalent)
  }
#ifdef qdebug
  if (qH)
  {
    G4cout<<"G4QNucleus::DecayDiBaryon: deleted at end - PDG="<<qH->GetPDGCode()<<G4endl;
    delete qH;
  }
#endif
} // End of DecayAntiDibaryon

//Decay of the nuclear states with antistrangeness (K:/K0)
void G4QNucleus::DecayAntiStrange(G4QHadron* qH, G4QHadronVector* evaHV)
{
  static const G4double mK    = G4QPDGCode(321).GetMass();
  static const G4double mK0   = G4QPDGCode(311).GetMass();
  static const G4double eps   = 0.003;
  G4LorentzVector q4M = qH->Get4Momentum();    // Get 4-mom of the AntiStrangeNuclearState
  G4double         qM = q4M.m();               // Real mass of the AntiStrangeNuclearState
  G4QContent       qQC= qH->GetQC();           // PDGCode of theDecayingAntiStrangeNuclSt.
  G4int            qS = qH->GetStrangeness();  // Strangness of the AntiStrangeNuclearState
  G4int            qB = qH->GetBaryonNumber(); // BaryonNumber of the AntiStrangeNuclearSt.
  G4int            qP = qH->GetCharge();       // Charge of the AntiStranNuclState (a#of p)
#ifdef debug
  G4cout<<"G4QNuc::DecAntS:QC="<<qQC<<",S="<<qS<<",B="<<qB<<",C="<<qP<<",4M="<<q4M<<G4endl;
#endif
  G4int            qN = qB-qP-qS;              // a#of neuterons
  if(qS>=0 || qB<1)
  {
    delete qH;
    // G4cerr<<"G4QNuc::DecayAntiStrange:QC="<<qQC<<",S="<<qS<<",B="<<qB<<",4M="<<q4M<<G4endl;
    // throw G4QException("G4QNucleus::DecayAntiStrange: not an Anti Strange Nucleus");
    G4ExceptionDescription ed;
    ed << "not an Anti Strange Nucleus: QC=" << qQC << ",S=" << qS << ",B="
       << qB << ",4M=" << q4M << G4endl;
    G4Exception("G4QNucleus::DecayAntiStrange()", "HAD_CHPS_0000",
                FatalException, ed);
  }
  G4int n1=1;         // prototype of a#of K0's
  G4double k1M=mK0;
  G4int k1PDG=311;    // K0 (as a prototype)
  G4int n2=0;         // prototype of a#of K+'s
  G4double k2M=mK;
  G4int k2PDG=321;    // K+
  G4int aS=-qS;       // -Strangness = antistrangeness
  G4int sH=aS/2;     // a small half of the antistrangeness
  G4int bH=aS-sH;    // a big half to take into account all the antistrangeness
  if(qP>0 && qP>qN)  // a#of protons > a#of neutrons
  {
    if(qP>=bH)                       // => "Enough protons in nucleus" case
    {
      if(qN>=sH)
      {
        n1=sH;
        n2=bH;
      }
      else
      {
        G4int dPN=qP-qN;
        if(dPN>=aS)
        {
          n1=0;
          n2=aS;
        }
        else
        {
          G4int sS=(aS-dPN)/2;
          G4int bS=aS-dPN-sS;
          sS+=dPN;
          if(qP>=sS && qN>=bS)
          {
            n1=bS;
            n2=sS;
          }
          else if(qP<sS)
          {
            n1=aS-qP;
            n2=qP;
          }
          else
          {
            n1=qN;
            n2=aS-qN;
          }
        }
      }
    }
  }
  else if(qN>=bH)
  {
    if(qP>=sH)
    {
      n2=sH;
      n1=bH;
    }
    else
    {
      G4int dNP=qN-qP;
      if(dNP>=aS)
      {
        n1=aS;
        n2=0;
      }
      else
      {
        G4int sS=(aS-dNP)/2;
        G4int bS=aS-sS;
        if(qN>=bS && qP>=sS)
        {
          n1=bS;
          n2=sS;
        }
        else if(qN<bS)
        {
          n1=qN;
          n2=aS-qN;
        }
        else
        {
          n1=aS-qP;
          n2=qP;
        }
      }
    }
  }
  G4int qPDG=90000000+(qP-n2)*1000+(qN-n1);     // PDG of the Residual Nucleus
  G4double nucM = G4QNucleus(qPDG).GetGSMass(); // Mass of the Residual Nucleus
#ifdef debug
  G4cout<<"G4QNucleus::DecayAnStran:nK0="<<n1<<",nK+="<<n2<<", nucM="<<nucM<<G4endl;
#endif
  G4int m1=0;                        // prototype of a#of K0's
  G4int m2_value=qP;                 // prototype of a#of K+'s
  if(qP>=-qS)   m2_value=-qS;        // Enough charge for K+'s
  else if(qP>0) m1=-qS-qP;           // Anti-Lambdas are partially compensated by neutrons
  G4int sPDG=90000000+(qP-m2_value)*1000+(qN-m1); // PDG of the Residual Nucleus
  G4double mucM = G4QNucleus(sPDG).GetGSMass(); // Mass of the Residual Nucleus
  if(mucM+m1*mK+m2_value*mK0<nucM+n1*mK+n2*mK0) // New is smaller
  {
    qPDG=sPDG;
    nucM=mucM;
    n1=m1;
    n2=m2_value;
  }
#ifdef debug
  G4cout<<"G4QNucleus::DecayAnStran: n1="<<n1<<", n2="<<n2<<", nM="<<nucM<<G4endl;
#endif
  if(!n1||!n2)                            // AntiKaons of only one sort are found
  {
    if(!n1)                               // No K0's only K+'s
    {
      if(n2==1 && mK+nucM>qM+.0001)  // Mass limit: switch to K0
      {
        k1M=mK0;
        n1=1;
        qPDG=90000000+qP*1000+qN-1;       // PDG of the Residual Nucleus
        nucM = G4QNucleus(qPDG).GetGSMass(); // Mass of the Residual Nucleus
      }
      else
      {
        k1M=mK;
        k1PDG=321;                        // Only K+'s (default K0's)
        n1=n2;                            // only n1 is used
      }
    }
    else                                  // No K+'s only K0's
    {
      if(n1==1 && mK0+nucM>qM+.0001) // Mass limit: switch to K+
      {
        k1M=mK;
        k1PDG=321;                        // K+ instead of K0
        qPDG=90000000+(qP-1)*1000+qN;     // PDG of the Residual Nucleus
        nucM = G4QNucleus(qPDG).GetGSMass(); // Mass of the Residual Nucleus
      }
      else k1M=mK0;                      // Only anti-K0's (default k1PDG)
    }
#ifdef debug
    G4int naPDG=90000000+(qP-1)*1000+qN; // Prot PDG of the Alternative Residual Nucleus
    G4double naM=G4QNucleus(naPDG).GetGSMass(); // Prot Mass of the Alt. Residual Nucleus
    G4double kaM=mK;                     // Prot Mass of the Alternative kaon (K+)
    if(k1PDG==321)                       // Calculate alternative to K+
    {
      naPDG=90000000+qP*1000+qN-1;       // PDG of the Alternative Residual Nucleus
      naM=G4QNucleus(naPDG).GetGSMass(); // Mass of the Alt. Residual Nucleus
      kaM=mK0;                           // Prot Mass of the Alternative kaon (K0)
    }
    G4cout<<"G4QNucleus::DecayAnStran:M="<<qM<<",kM="<<k1M<<"+nM="<<nucM<<"="<<k1M+nucM
          <<",m="<<kaM<<"+n="<<naM<<"="<<kaM+naM<<G4endl;
#endif
    G4double n1M=n1*k1M;
    G4LorentzVector f4Mom(0.,0.,0.,n1M);
    G4LorentzVector s4Mom(0.,0.,0.,nucM);
    G4double sum=nucM+n1M;
    if(sum>qM+eps && n1==1)              // Try to use another K if this is the only kaon
    {
      G4int naPDG=90000000+(qP-1)*1000+qN; // Prot PDG of the Alternative Residual Nucleus
      G4double naM=G4QNucleus(naPDG).GetGSMass(); // Prot Mass of the Alt. Residual Nucleus
      G4int akPDG=321;                   // Prototype PDGCode of the AlternativeKaon (K+)
      G4double kaM=mK;                   // Prototype Mass of the AlternativeKaon (K+)
      if(k1PDG==321)                     // Calculate alternative to the K+ meson
      {
        naPDG=90000000+qP*1000+qN-1;     // PDG of the Alternative Residual Nucleus
        naM=G4QNucleus(naPDG).GetGSMass(); // Mass of the Alt. Residual Nucleus
        akPDG=311;                       // PDG Code of the Alternative kaon (K0)
        kaM=mK0;                         // Mass of the Alternative kaon (K0)
      }
      G4double asum=naM+kaM;
      if(asum<sum)                       // Make a KSwap correction
      {
        nucM=naM;
        n1M=kaM;
        k1M=kaM;
        k1PDG=akPDG;
        qPDG=naPDG;
        f4Mom=G4LorentzVector(0.,0.,0.,n1M);
        s4Mom=G4LorentzVector(0.,0.,0.,nucM);
      }
    }
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(n1M/sum);
      s4Mom=q4M*(nucM/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
#ifdef debug
      G4cout<<"--Warning--G4QNuc::DASt:AsItIs, H="<<qQC<<q4M<<qM<<" < sum="<<sum<<"=(F)"
            <<nucM<<"+(kK)"<<n1M<<G4endl;
#endif
      evaHV->push_back(qH);  // @@ Can cause problems with particle conversion in G4
      return;
    }
#ifdef debug
    G4cout<<"G4QNuc::DecAntiS: nK+N "<<n1<<"*K="<<k1PDG<<f4Mom<<",N="<<qPDG<<s4Mom<<G4endl;
#endif
    delete qH;
    //
    f4Mom/=n1;
    for(G4int i1=0; i1<n1; i1++)
    {
      G4QHadron* H1 = new G4QHadron(k1PDG,f4Mom); // Create a Hadron for the Kaon
      evaHV->push_back(H1);                       // Fill "H1" (delete equivalent)
    }
    G4QHadron* H2 = new G4QHadron(qPDG,s4Mom);    // Create a Hadron for the Nucleus
    //evaHV->push_back(H2);                       // Fill "H2" (delete equivalent)
    EvaporateNucleus(H2,evaHV);                   // Fill "H2" (delete equivalent)
#ifdef debug
    G4cout<<"G4QNucleus::DecAntiStr:*** After EvaporateNucleus nH="<<evaHV->size()<<G4endl;
#endif
  }
  else
  {
    G4double n1M=n1*k1M;
    G4double n2M=n2*k2M;
    G4LorentzVector f4Mom(0.,0.,0.,n1M);
    G4LorentzVector s4Mom(0.,0.,0.,n2M);
    G4LorentzVector t4Mom(0.,0.,0.,nucM);
    G4double sum=nucM+n1M+n2M;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(n1M/sum);
      s4Mom=q4M*(n2M/sum);
      t4Mom=q4M*(nucM/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
    {
      G4cout<<"---Warning---G4QN::DASt:nPDG="<<qPDG<<"(M="<<nucM<<")+1="<<k1PDG<<"(M="<<k1M
            <<")+2="<<k2PDG<<"(M="<<k2M<<")="<<nucM+n1*k1M+n2*k2M<<">tM="<<qM<<q4M<<G4endl;
      evaHV->push_back(qH); // @@ Can cause problems with particle conversion in G4
      return;
    }
#ifdef debug
    G4cout<<"G4QNuc::DecAntiS:*DONE* nPDG="<<qPDG<<",1="<<f4Mom<<",2="<<s4Mom<<",n="<<t4Mom
          <<G4endl;
#endif
    delete qH;
    //
    f4Mom/=n1;
    for(G4int i1=0; i1<n1; i1++)
    {
      G4QHadron* H1 = new G4QHadron(k1PDG,f4Mom); // Create a Hadron for the K0
      evaHV->push_back(H1);                       // Fill "H1" (delete equivalent)
    }
    s4Mom/=n2;
    for(G4int i2=0; i2<n2; i2++)
    {
      G4QHadron* H2 = new G4QHadron(k2PDG,s4Mom); // Create a Hadron for the K+
      evaHV->push_back(H2);                       // Fill "H2" (delete equivalent)
    }
    G4QHadron* H3 = new G4QHadron(qPDG,t4Mom);    // Create a Hadron for the nucleus
    //evaHV->push_back(H3);                       // Fill "H3" (delete equivalent)
    EvaporateNucleus(H3,evaHV);                   // Fill "H3" (delete equivalent)
  }
#ifdef qdebug
  if (qH)
  {
    G4cout << "G4QNucleus::DecayAntiStrange: deleted at end - PDG: "
           << qH->GetPDGCode() << G4endl;
    delete qH;
  }
#endif
#ifdef debug
  G4cout<<"G4QNucleus::DecayAntiStrange: ===> End of DecayAntiStrangness"<<G4endl;
#endif
} // End of DecayAntiStrange

//Decay of the excited 3p or 3n systems in three baryons
void G4QNucleus::DecayMultyBaryon(G4QHadron* qH, G4QHadronVector* evaHV)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double eps=0.003;
  G4LorentzVector q4M = qH->Get4Momentum();       // Get 4-momentum of the MultyBaryon
  G4double         qM = q4M.m();                  // Mass of the Multybaryon
  G4int          qPDG = qH->GetPDGCode();         // PDG Code of the decaying multybar
  G4QContent      qQC = qH->GetQC();              // PDG Code of the decaying multibar
#ifdef debug
  G4cout<<"G4QNuc::DecayMultyBaryon: *Called* PDG="<<qPDG<<",4M="<<q4M<<qQC<<G4endl;
#endif
  G4int totS=qQC.GetStrangeness();  //  Total Strangeness       (L)                ^
  G4int totC=qQC.GetCharge();       //  Total Charge            (p)                ^
  G4int totBN=qQC.GetBaryonNumber();// Total Baryon Number      (A)                ^
  G4int totN=totBN-totS-totC;       // Total Number of Neutrons (n)                ^
  G4int          fPDG = 3122;       // Prototype for A lambdas case
  G4double       fMass= mLamb;
  if     (totN==totBN)              // "A-neutron" case
  {
    fPDG = 2112;
    fMass= mNeut;
  }
  else if(totC==totBN)              // "A-protons" case
  {
    fPDG = 2212;
    fMass= mProt;
  }
  else if(totS!=totBN)            // "Bad call" case
  {
    delete qH;
    // G4cerr<<"***G4QNuc::DecayMultyBaryon: PDG="<<qPDG<<G4endl;
    // throw G4QException("***G4QNuc::DecayMultyBaryon: Can not decay this PDG Code");
    G4ExceptionDescription ed;
    ed << "Can not decay this PDG Code: PDG=" << qPDG << G4endl;
    G4Exception("G4QNucleus::DecayMultyBaryon()", "HAD_CHPS_0000",
                FatalException, ed);
  }
#ifdef debug
  else
  {
    delete qH;
    // G4cerr<<"**G4QNucleus::DecayMultyBaryon: PDG="<<qPDG<<G4endl;
    // throw G4QException("***G4QNuc::DecayMultyBaryon: Unknown PDG code of the MultiBaryon");
    G4ExceptionDescription ed;
    ed << "Unknown PDG code of the MultiBaryon: PDG=" << qPDG << G4endl;
    G4Exception("G4QNucleus::DecayMultyBaryon()", "HAD_CHPS_0001",
                FatalException, ed);
  }
#endif
  if(totBN==1) evaHV->push_back(qH);
  else if(totBN==2)
  {
    G4LorentzVector f4Mom(0.,0.,0.,fMass);
    G4LorentzVector s4Mom(0.,0.,0.,fMass);
    G4double sum=fMass+fMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M/2.;
      s4Mom=f4Mom;
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cout<<"---Warning---G4QNucl::DecayMultyBar:fPDG="<<fPDG<<"(fM="<<fMass<<")*2="<<sum
            <<" > TotM="<<q4M.m()<<q4M<<G4endl;
      //G4cerr<<"***G4QNuc::DecayMultyBaryon:qM="<<qM<<"<sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("G4QNuc::DecayMultyBaryon:diBaryon DecayIn2 didn't succeed");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNucleus::DecMulBar:*DONE* fPDG="<<fPDG<<",f="<<f4Mom<<",s="<<s4Mom<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);   // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                   // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(fPDG,s4Mom);   // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                   // Fill "H2" (delete equivalent)
  }
  else if(totBN==3)
  {
    G4LorentzVector f4Mom(0.,0.,0.,fMass);
    G4LorentzVector s4Mom(0.,0.,0.,fMass);
    G4LorentzVector t4Mom(0.,0.,0.,fMass);
    G4double sum=fMass+fMass+fMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M/3.;
      s4Mom=f4Mom;
      t4Mom=f4Mom;
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
    {
      G4cout<<"---Warning---G4QNuc::DecayMultyBaryon: fPDG="<<fPDG<<"(fM="<<fMass<<")*3 = "
            <<3*fMass<<" >? TotM="<<q4M.m()<<q4M<<G4endl;
      //G4cerr<<"***G4QN::DecayMultyBar: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("G4QNucleus::DecayMultyBar:ThreeBaryonDecayIn3 didn't succeed");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNuc::DecMBar:*DONE*, fPDG="<<fPDG<<",f="<<f4Mom<<",s="<<s4Mom<<",t="
          <<t4Mom<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);   // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                   // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(fPDG,s4Mom);   // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                   // Fill "H2" (delete equivalent)
    G4QHadron* H3 = new G4QHadron(fPDG,t4Mom);   // Create a Hadron for the 3-d baryon
    evaHV->push_back(H3);                   // Fill "H3" (delete equivalent)
  }
  else
  {
    // @@It must be checked, that they are not under the mass shell
    // !! OK !! Checked by the warning print that they are mostly in the Ground State !!
    G4LorentzVector f4Mom=q4M/totBN; // @@ Too simple solution (split in two parts!)
#ifdef debug
    // Warning for the future development
    G4cout<<"*G4QNul::DecMulBar:SplitMultiBar inEqParts M="<<totBN<<"*"<<f4Mom.m()<<G4endl;
    G4cout<<"G4QNucleus::DecMultyBaryon: *DONE* fPDG="<<fPDG<<", f="<<f4Mom<<G4endl;
#endif
    delete qH;
    for(G4int h=0; h<totBN; h++)
    {
      G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the baryon
      evaHV->push_back(H1);                 // Fill "H1" (delete equivalent)
    }
  }
#ifdef qdebug
  if (qH)
  {
    G4cout << "G4QNucleus::DecayMultyBaryon: deleted at end - PDG: "
           << qH->GetPDGCode() << G4endl;
    delete qH;
  }
#endif
} // End of DecayMultyBaryon

//Decay of the excited alpha+2p or alpha+2n systems
void G4QNucleus::DecayAlphaDiN(G4QHadron* qH, G4QHadronVector* evaHV)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mHel6= G4QPDGCode(2112).GetNuclMass(2,4,0);
  static const G4double eps=0.003;
  G4LorentzVector q4M = qH->Get4Momentum();      // Get 4-momentum of the AlphaDibaryon
  G4double         qM = q4M.m();                 // Real mass of the AlphaDibaryon
  G4int          qPDG = qH->GetPDGCode();        // PDG Code of the decayin AlphaDybaryon
#ifdef debug
  G4cout<<"G4QNuc::DecayAlphaDiN: *Called* PDG="<<qPDG<<",4M="<<q4M<<G4endl;
#endif
  G4int          fPDG = 2212;                    // Prototype for alpha+pp case
  G4double       fMass= mProt;
  G4int          sPDG = 90002002;
  G4double       sMass= mAlph;
  if     (qPDG==90002004)                        // "alpha+2neutrons" case
  {
    if(fabs(qM-mHel6)<eps)
    {
      evaHV->push_back(qH);                 // Fill as it is (delete equivalent)
      return;
    }
    else if(mNeut+mNeut+mAlph<qM)
    {
      fPDG = 2112;
      fMass= mNeut;
    }
    else
    {
      delete qH;
      // G4cerr<<"***G4QNu::DecAlDiN:M(He6="<<mHel6<<")="<<qM<<"<"<<mNeut+mNeut+mAlph<<G4endl;
      // throw G4QException("G4QNuc::DecayAlphaDiN: Cannot decay excited He6 with this mass");
      G4ExceptionDescription ed;
      ed << "Cannot decay excited He6 with this mass: M(He6=" << mHel6 << ")="
         << qM << "<" << mNeut+mNeut+mAlph << G4endl;
      G4Exception("G4QNucleus::DecayAlphaDiN()", "HAD_CHPS_0000",
                  FatalException, ed);
    }
  }
  else if(qPDG!=90004002)                         // "Bad call" case
  {
    delete qH;
    // G4cerr<<"***G4QNuc::DecayAlphaDiN: PDG="<<qPDG<<G4endl;
    // throw G4QException("G4QNuc::DecayAlphaDiN: Can not decay this PDG Code");
    G4ExceptionDescription ed;
    ed << "Can not decay this PDG Code: PDG=" << qPDG << G4endl;
    G4Exception("G4QNucleus::DecayAlphaDiN()", "HAD_CHPS_0001",
                FatalException, ed);
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,fMass);
  G4LorentzVector t4Mom(0.,0.,0.,sMass);
  G4double sum=fMass+fMass+sMass;
  if(fabs(qM-sum)<eps)
  {
    f4Mom=q4M*(fMass/sum);
    s4Mom=f4Mom;
    t4Mom=q4M*(sMass/sum);
  }
  else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
  {
    G4cout<<"---Warning---G4QNuc::DecayAlphaDiN:fPDG="<<fPDG<<"(M="<<fMass<<")*2+mAlpha = "
          <<sum<<" >? TotM="<<qM<<q4M<<", d="<<sum-qM<<G4endl;
    //G4cerr<<"***G4QNuc::DecayAlphaDiN: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
    //throw G4QException("G4QNuc::DecayAlphaDiN: Alpha+N+N DecayIn3 error");
    evaHV->push_back(qH);
    return;
  }
#ifdef debug
  G4cout<<"G4QNuc::DecAl2N: fPDG="<<fPDG<<",f="<<f4Mom<<",s="<<s4Mom<<",t="<<t4Mom<<G4endl;
#endif
  delete qH;
  G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);    // Create a Hadron for the 1-st baryon
  evaHV->push_back(H1);                    // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(fPDG,s4Mom);    // Create a Hadron for the 2-nd baryon
  evaHV->push_back(H2);                    // Fill "H2" (delete equivalent)
  G4QHadron* H3 = new G4QHadron(sPDG,t4Mom);    // Create a Hadron for the alpha
  evaHV->push_back(H3);                    // Fill "H3" (delete equivalent)
} // End of DecayAlphaDiN

//Decay of the excited alpha+bayon state in alpha and baryons
void G4QNucleus::DecayAlphaBar(G4QHadron* qH, G4QHadronVector* evaHV)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mTrit= G4QPDGCode(2112).GetNuclMass(1,2,0);
  static const G4double mHe3 = G4QPDGCode(2112).GetNuclMass(2,1,0);
  static const G4double eps=0.003;
  G4LorentzVector q4M = qH->Get4Momentum();     // Get 4-momentum of the Alpha-Baryon
  G4double         qM = q4M.m();                // Mass of Alpha-Baryon
  G4int          qPDG = qH->GetPDGCode();       // PDG Code of the decayin Alpha-Baryon
  G4QContent      qQC = qH->GetQC();            // PDG Code of the decaying Alpha-Bar
#ifdef debug
  G4cout<<"G4QNucleus::DecayAlphaBar: *Called* PDG="<<qPDG<<",4M="<<q4M<<qQC<<G4endl;
#endif
  G4int totS=qQC.GetStrangeness();              //  Total Strangeness       (L)
  G4int totC=qQC.GetCharge();                   //  Total Charge            (p)
  G4int totBN=qQC.GetBaryonNumber();            // Total Baryon Number      (A)

  if ( ( (!totS && !totC) || totC == totBN || totS == totBN) 
       && totBN > 1) DecayMultyBaryon(qH,evaHV);
  else if(qPDG==92001002||qPDG==92002001||qPDG==91003001||qPDG==91001003||qPDG==93001001)
    evaHV->push_back(qH);
  else if(qPDG==92000003||qPDG==92003000||qPDG==93000002||qPDG==93002000)
  {
    G4int          fPDG = 3122;                 // 1st Prototype for 2L+3n case
    G4double       fMass= mLamb;
    G4int          sPDG = 2112;
    G4double       sMass= mNeut;
    if     (qPDG==92003000)                     // "2L+3p" case
    {
      sPDG = 2212;
      sMass= mProt;
    }
    else if(qPDG==93000002)                     // "2n+3L" case
    {
      fPDG = 2112;
      fMass= mNeut;
      sPDG = 3122;
      sMass= mLamb;
    }
    else if(qPDG==93002000)                     // "2p+3L" case
    {
      fPDG = 2212;
      fMass= mProt;
      sPDG = 3122;
      sMass= mLamb;
    }
    G4double tfM=fMass+fMass;
    G4double tsM=sMass+sMass+sMass;
    G4LorentzVector f4Mom(0.,0.,0.,tfM);
    G4LorentzVector s4Mom(0.,0.,0.,tsM);
    G4double sum=tfM+tsM;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(tfM/sum);
      s4Mom=q4M*(tsM/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cout<<"--Warning--G4QNuc::DecAlB:fPDG="<<fPDG<<"(M="<<fMass<<")*2="<<2*fMass<<",s="
            <<sPDG<<"(sM="<<sMass<<")*3="<<3*sMass<<"="<<sum<<">M="<<q4M.m()<<q4M<<G4endl;
      //G4cerr<<"***G4QN::DecayAlphaBar: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("G4QNucleus::DecayAlphaBar: DecayIn2 didn't succeed for 3/2");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNucleus::DecAlB:*DONE*, fPDG="<<fPDG<<f4Mom<<",sPDG="<<sPDG<<s4Mom<<G4endl;
#endif
    delete qH;
    G4LorentzVector rf4Mom=f4Mom/2;
    G4QHadron* H1 = new G4QHadron(fPDG,rf4Mom); // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                  // Fill "H1" (delete equivalent)
    evaHV->push_back(H1);                  // Fill "H1" (delete equivalent)
    G4LorentzVector rs4Mom=s4Mom/3;
    G4QHadron* H2 = new G4QHadron(sPDG,rs4Mom); // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                  // Fill "H2" (delete equivalent)
    evaHV->push_back(H2);                  // Fill "H2" (delete equivalent)
    evaHV->push_back(H2);                  // Fill "H2" (delete equivalent)
  }
  else if(qPDG==90004001||qPDG==90001004)
  {
    G4int          fPDG = 90002001;             // Prototype for "He3+2p" case
    G4double       fMass= mHe3;
    G4int          sPDG = 2212;
    G4double       sMass= mProt;
    if     (qPDG==90001004)                     // "t+2n" case
    {
      fPDG = 90001002;
      fMass= mTrit;
      sPDG = 2112;
      sMass= mNeut;
    }
    G4LorentzVector f4Mom(0.,0.,0.,fMass);
    G4LorentzVector s4Mom(0.,0.,0.,sMass);
    G4LorentzVector t4Mom(0.,0.,0.,sMass);
    G4double sum=fMass+sMass+sMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
      t4Mom=s4Mom;
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
    {
      G4cout<<"--Warning--G4QNuc::DecAlB:fPDG="<<fPDG<<",M="<<fMass<<",sPDG="<<sPDG<<",sM="
            <<sMass<<",2sM+fM="<<2*sMass+fMass<<" > TotM="<<q4M.m()<<q4M<<G4endl;
      //G4cerr<<"*G4QNuc::DecayAlphaBar: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("G4QNucleus::DecayAlphaBar: t/nn,He3/pp DecayIn3 didn't");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNucl::DecAlB: *DONE*, f="<<fPDG<<f4Mom<<", s="<<sPDG<<s4Mom<<t4Mom<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);   // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                   // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom);   // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                   // Fill "H2" (delete equivalent)
    G4QHadron* H3 = new G4QHadron(sPDG,t4Mom);   // Create a Hadron for the 3-d baryon
    evaHV->push_back(H3);                   // Fill "H3" (delete equivalent)
  }
  else if(qPDG==94000001||qPDG==94001000||qPDG==91000004||qPDG==91004000)
  {
    G4int          fPDG = 3122;                 // Prototype for "4L+n" case
    G4double       fMass= mLamb+mLamb;
    G4int          sPDG = 2112;
    G4double       sMass= mNeut;
    if     (qPDG==94001000)                     // "4L+p" case
    {
      sPDG = 2212;
      sMass= mProt;
    }
    else if(qPDG==91000004)                     // "4n+L" case
    {
      fPDG = 2112;
      fMass= mNeut+mNeut;
      sPDG = 3122;
      sMass= mLamb;
    }
    else if(qPDG==91004000)                     // "4p+L" case
    {
      fPDG = 2212;
      fMass= mProt+mProt;
      sPDG = 3122;
      sMass= mLamb;
    }
    G4LorentzVector f4Mom(0.,0.,0.,fMass+fMass);
    G4LorentzVector s4Mom(0.,0.,0.,sMass);
    G4double sum=fMass+fMass+sMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*((fMass+fMass)/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cout<<"--Warning--G4QNucl::DecAlphBar:fPDG="<<fPDG<<"(2*fM="<<fMass<<")*2="
            <<2*fMass<<",sPDG="<<sPDG<<"(sM="<<sMass<<" > TotM="<<q4M.m()<<q4M<<G4endl;
      //G4cerr<<"*G4QNuc::DecayAlphaBar: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("G4QNucl::DecayAlphaBar:QuintBaryon DecayIn2 didn't succeed");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNuc::DecAlphaB: *DONE*, fPDG="<<fPDG<<f4Mom<<",sPDG="<<sPDG<<s4Mom<<G4endl;
#endif
    delete qH;
    G4LorentzVector rf4Mom=f4Mom/4;
    G4QHadron* H1 = new G4QHadron(fPDG,rf4Mom); // Create a Hadron for the 1-st baryon
    evaHV->push_back(H1);                  // Fill "H1" (delete equivalent)
    evaHV->push_back(H1);                  // Fill "H1" (delete equivalent)
    evaHV->push_back(H1);                  // Fill "H1" (delete equivalent)
    evaHV->push_back(H1);                  // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom);  // Create a Hadron for the 2-nd baryon
    evaHV->push_back(H2);                  // Fill "H2" (delete equivalent)
  }
  else if(qPDG==90003002||qPDG==90002003||qPDG==91002002)
  {
    G4int          fPDG = 90002002;             // Prototype for "alpha+n" case
    G4int          sPDG = 2112;
    G4double       fMass= mAlph;
    G4double       sMass= mNeut;
    if(qPDG==90003002)                          // "alpha+p" case
    {
      sPDG = 2212;
      sMass= mProt;    
    }
    else if(qPDG==9100202)                      // "alpha+l" case
    {
      sPDG = 3122;
      sMass= mLamb;    
    }
    else if(qPDG!=90002003)
    {
      evaHV->push_back(qH);                     // Fill hadron as it is (delete equivalent)
      //EvaporateNucleus(qH, evaHV);            // Evaporate Nucleus (delete equivivalent)
      return;
    }
    G4double dM=fMass+sMass-qM;
    if(dM>0.&&dM<1.)
    {
#ifdef debug
      G4cout<<"***Corrected***G4QNuc::DecayAlphaBar:fPDG="<<fPDG<<"(fM="<<fMass<<")+ sPDG="
            <<sPDG<<"(sM="<<sMass<<")="<<fMass+sMass<<" > TotM="<<qM<<q4M<<G4endl;
#endif
      G4double hdM=dM/2;
      fMass-=hdM;
      sMass-=hdM;
    }
    G4LorentzVector f4Mom(0.,0.,0.,fMass);
    G4LorentzVector s4Mom(0.,0.,0.,sMass);      // Mass is random since probab. time
    G4double sum=fMass+sMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cout<<"--Warning--G4QNuc::DecAlphaBar:fPDG="<<fPDG<<"(fM="<<fMass<<")+sPDG="<<sPDG
            <<"(sM="<<sMass<<")="<<fMass+sMass<<"="<<sum<<" > TotM="<<q4M.m()<<q4M<<G4endl;
      //G4cout<<"*G4QNuc::DecayAlphaBar: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("***G4QNucl::DecayAlphaBar:Alpha+Baryon DecIn2 didn't succeed");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNucl::DecAlBar:*DONE*a4M="<<f4Mom<<",s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);      // Create a Hadron for the alpha
    evaHV->push_back(H1);                      // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom);      // Create a Hadron for the baryon
    evaHV->push_back(H2);                      // Fill "H2" (delete equivalent)
  }
  else G4cout<<"---Warning---G4QNucleus::DecayAlphaBar: Unknown PDG="<<qPDG<<G4endl;
#ifdef qdebug
  if (qH)
  {
    G4cout << "G4QNucleus::DecayAlphaBar: deleted at end - PDG: "
           << qH->GetPDGCode() << G4endl;
    delete qH;
  }
#endif
} // End of DecayAlphaBar

//Decay of the excited alpha+alpha state in 2 alphas
void G4QNucleus::DecayAlphaAlpha(G4QHadron* qH, G4QHadronVector* evaHV)
{
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double aaGSM= G4QPDGCode(2112).GetNuclMass(4,4,0);
  static const G4double eps=0.003;
  G4int          qPDG = qH->GetPDGCode();         // PDG Code of the decayin dialpha
  if(qPDG!=90004004)
  {
    delete qH;
    // G4cerr<<"***G4QNucleus::DecayAlphaAlpha: qPDG="<<qPDG<<G4endl;
    // throw G4QException("***G4QNucleus::DecayAlphaAlpha: Not Be8 state decais in 2 alphas");
    G4ExceptionDescription ed;
    ed << "Not Be8 state decais in 2 alphas: qPDG=" << qPDG << G4endl;
    G4Exception("G4QNucleus::DecayAlphaAlpha()", "HAD_CHPS_0000",
                FatalException, ed);
  }
  G4LorentzVector q4M = qH->Get4Momentum();       // Get 4-momentum of the Dibaryon
  G4double qM=q4M.m();
#ifdef debug
  G4cout<<"G4QNucleus::DecayAlAl: *Called* PDG="<<qPDG<<",M="<<qM<<q4M<<">"<<aaGSM<<G4endl;
#endif
  //if(qM>aaGSM+.01)  // @@ Be8*->gamma+Be8 (as in evaporation) @@ gamma cooling
  if(2>3)
  {
    G4int          fPDG = 22;
    G4int          sPDG = 90004004;
    G4double       fMass= 0.;
    G4double       sMass= aaGSM;
    G4LorentzVector f4Mom(0.,0.,0.,fMass);
    G4LorentzVector s4Mom(0.,0.,0.,sMass);          // Mass is random since probab. time
    G4double sum=fMass+sMass;
    if(fabs(qM-sum)<eps)
    {
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cout<<"---Warning---G4QNuc::DecayAlphaAlpha:gPDG="<<fPDG<<"(gM="<<fMass<<")+PDG="
            <<sPDG<<"(sM="<<sMass<<")="<<sum<<" > TotM="<<q4M.m()<<q4M<<G4endl;
      //G4cerr<<"***G4QNuc::DecayAlphAlph: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
      //throw G4QException("G4QNucleus::DecayAlphaAlpha:g+diAlph DecayIn2 didn't succeed");
      evaHV->push_back(qH);
      return;
    }
#ifdef debug
    G4cout<<"G4QNucleus::DecayAlphaAlpha: *DONE* gam4M="<<f4Mom<<", aa4M="<<s4Mom<<G4endl;
#endif
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);      // Create a Hadron for the 1-st alpha
    evaHV->push_back(H1);                      // Fill "H1" (delete equivalent)
    qH->Set4Momentum(s4Mom);
    q4M=s4Mom;
  }
  G4int          fPDG = 90002002;
  G4int          sPDG = 90002002;
  G4double       fMass= mAlph;
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,fMass);
  G4double sum=fMass+fMass;
  if(fabs(qM-sum)<eps)
  {
    f4Mom=q4M*(fMass/sum);
    s4Mom=f4Mom;
  }
  else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
  {
    G4cout<<"---Warning---G4QNucl::DecayAlphaAlpha:fPDG="<<fPDG<<"(fM="<<fMass<<")*2="<<sum
          <<" > TotM="<<q4M.m()<<q4M<<G4endl;
    //G4cerr<<"***G4QNuc::DecayAlphAlph: qM="<<qM<<" < sum="<<sum<<",d="<<sum-qM<<G4endl;
    //throw G4QException("G4QNucleus::DecayAlphaAlpha: diAlpha DecayIn2 didn't succeed");
    evaHV->push_back(qH);
    return;
  }
#ifdef debug
  G4cout<<"G4QNucleus::DecayAlphaAlpha: *DONE* fal4M="<<f4Mom<<", sal4M="<<s4Mom<<G4endl;
#endif
  delete qH;
  G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);      // Create a Hadron for the 1-st alpha
  evaHV->push_back(H1);                      // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(sPDG,s4Mom);      // Create a Hadron for the 2-nd alpha
  evaHV->push_back(H2);                      // Fill "H2" (delete equivalent)
} // End of DecayAlphaAlpha
