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
//
// G4 Physics class: G4QDiffractionRatio for N+A Diffraction Interactions
// Created: M.V. Kossov, CERN/ITEP(Moscow), 25-OCT-01
// The last update: M.V. Kossov, CERN/ITEP(Moscow) 10-Nov-09
//
// --------------------------------------------------------------------
// Short description: Difraction excitation is a part of the incoherent
// (inelastic) interaction. This part is calculated in the class.
// --------------------------------------------------------------------

//#define debug
//#define pdebug
//#define fdebug
//#define nandebug

#include "G4QDiffractionRatio.hh"
#include "G4SystemOfUnits.hh"

// Returns Pointer to the G4VQCrossSection class
G4QDiffractionRatio* G4QDiffractionRatio::GetPointer()
{
  static G4QDiffractionRatio theRatios;   // *** Static body of the Diffraction Ratio ***
  return &theRatios;
}

// Calculation of pair(QuasiFree/Inelastic,QuasiElastic/QuasiFree)
G4double G4QDiffractionRatio::GetRatio(G4double pIU, G4int pPDG, G4int tgZ, G4int tgN)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass()/GeV; // in GeV
  static const G4double mProt= G4QPDGCode(2212).GetMass()/GeV; // in GeV
  static const G4double mN=.5*(mNeut+mProt);  // mean nucleon mass in GeV
  static const G4double dmN=mN+mN;            // doubled nuc. mass in GeV
  static const G4double mN2=mN*mN;            // squared nuc. mass in GeV^2
  // Table parameters
  static const G4int    nps=100;        // Number of steps in the R(s) LinTable
  static const G4int    mps=nps+1;      // Number of elements in the R(s) LinTable
  static const G4double sma=6.;         // The first LinTabEl(sqs=0)=1., sqs>sma -> logTab
  static const G4double ds=sma/nps;     // Step of the linear Table
  static const G4int    nls=150;        // Number of steps in the R(lns) logTable
  static const G4int    mls=nls+1;      // Number of elements in the R(lns) logTable
  static const G4double lsi=1.79;       // The min ln(sqs) logTabEl(sqs=5.99 < sma=6.)
  static const G4double lsa=8.;         // The max ln(sqs) logTabEl(sqs=5.99 - 2981 GeV)
  static const G4double mi=std::exp(lsi);// The min s of logTabEl(~ 5.99 GeV)
  static const G4double max_s=std::exp(lsa);// The max s of logTabEl(~ 2981 GeV)
  static const G4double dl=(lsa-lsi)/nls;// Step of the logarithmic Table
  static const G4double edl=std::exp(dl);// Multiplication step of the logarithmic Table
  static const G4double toler=.0001;    // Tolarence (GeV) defining the same sqs
  static G4double lastS=0.;             // Last sqs value for which R was calculated
  static G4double lastR=0.;             // Last ratio R which was calculated
  // Local Associative Data Base:
  static std::vector<G4int>     vA;     // Vector of calculated A
  //static std::vector<G4double>  vH;     // Vector of max sqs initialized in the LinTable
  //static std::vector<G4int>     vN;     // Vector of topBin number initialized in LinTab
  //static std::vector<G4double>  vM;     // Vector of relMax ln(sqs) initialized in LogTab
  //static std::vector<G4int>     vK;     // Vector of topBin number initialized in LogTab
  static std::vector<G4double*> vT;     // Vector of pointers to LinTable in C++ heap
  static std::vector<G4double*> vL;     // Vector of pointers to LogTable in C++ heap
  // Last values of the Associative Data Base:
  //static G4int     lastPDG=0;           // Last PDG for which R was calculated (now fake)
  static G4int     lastA=0;             // theLast of calculated A
  //static G4double  lastH=0.;            // theLast max sqs initialized in the LinTable
  //static G4int     lastN=0;             // theLast of topBin number initialized in LinTab
  //static G4double  lastM=0.;            // theLast relMax ln(sqs) initialized in LogTab
  //static G4int     lastK=0;             // theLast of topBin number initialized in LogTab
  static G4double* lastT=0;             // theLast of pointer to LinTable in the C++ heap
  static G4double* lastL=0;             // theLast of pointer to LogTable in the C++ heap
  // LogTable is created only if necessary. R(sqs>2981GeV) calcul by formula for any nuclei
  G4int A=tgN+tgZ;
  if(pIU<toler || A<1) return 1.;       // Fake use of toler as non zero number
  if(A>238)
  {
    G4cout<<"-*-Warning-*-G4QuasiFreeRatio::GetRatio: A="<<A<<">238, return zero"<<G4endl;
    return 0.;
  }
  //lastPDG=pPDG;                         // @@ at present ratio is PDG independent @@
  // Calculate sqs
  G4double pM=G4QPDGCode(pPDG).GetMass()/GeV; // Projectile mass in GeV
  G4double pM2=pM*pM;
  G4double mom=pIU/GeV;                 // Projectile momentum in GeV
  G4double s_value=std::sqrt(mN2+pM2+dmN*std::sqrt(pM2+mom*mom)); // in GeV
  G4int nDB=vA.size();                  // A number of nuclei already initialized in AMDB
  if(nDB && lastA==A && std::fabs(s_value-lastS)<toler) return lastR;
  if(s_value>max_s)
  {
    lastR=CalcDiff2Prod_Ratio(s_value,A);     // @@ Probably user ought to be notified about bigS
    return lastR;
  }
  G4bool found=false;
  G4int i=-1;
  if(nDB) for (i=0; i<nDB; i++) if(A==vA[i]) // Sirch for this A in AMDB
  {
    found=true;                         // The A value is found
    break;
  }
  if(!nDB || !found)                    // Create new line in the AMDB
  {
    lastA = A;
    lastT = new G4double[mps];          // Create the linear Table
    //lastN = static_cast<int>(s_value/ds)+1;   // MaxBin to be initialized
    //if(lastN>nps)                     // ===> Now initialize all lin table
    //{
    //  lastN=nps;
    //  lastH=sma;
    //}
    //else lastH = lastN*ds;              // Calculate max initialized s for LinTab
    G4double sv=0;
    lastT[0]=1.;
    //for(G4int j=1; j<=lastN; j++)       // Calculate LinTab values
    for(G4int j=1; j<=nps; j++)       // Calculate LinTab values
    {
      sv+=ds;
      lastT[j]=CalcDiff2Prod_Ratio(sv,A);
    }
    lastL = new G4double[mls];          // Create the logarithmic Table
    //G4double ls=std::log(s_value);
    //lastK = static_cast<int>((ls-lsi)/dl)+1; // MaxBin to be initialized in LogTaB
    //if(lastK>nls)                     // ===> Now initialize all lin table
    //{
    //  lastK=nls;
    //  lastM=lsa-lsi;
    //}
    //else lastM = lastK*dl;              // Calculate max initialized ln(s)-lsi for LogTab
    sv=mi;
    //for(G4int j=0; j<=lastK; j++)     // Calculate LogTab values
    for(G4int j=0; j<=nls; j++)     // Calculate LogTab values
    {
      lastL[j]=CalcDiff2Prod_Ratio(sv,A);
      //if(j!=lastK) sv*=edl;
      sv*=edl;
    }
    i++;                                // Make a new record to AMDB and position on it
    vA.push_back(lastA);
    //vH.push_back(lastH);
    //vN.push_back(lastN);
    //vM.push_back(lastM);
    //vK.push_back(lastK);
    vT.push_back(lastT);
    vL.push_back(lastL);
  }
  else                                  // The A value was found in AMDB
  {
    lastA=vA[i];
    //lastH=vH[i];
    //lastN=vN[i];
    //lastM=vM[i];
    //lastK=vK[i];
    lastT=vT[i];
    lastL=vL[i];
    // ==> Now all bins of the tables are initialized immediately for the A
    //if(s_value>lastH)                    // At least LinTab must be updated
    //{
    //  G4int nextN=lastN+1;               // The next bin to be initialized
    //  if(lastN<nps)
    //  {
    //    lastN = static_cast<int>(s_value/ds)+1;// MaxBin to be initialized
    //    G4double sv=lastH;
    //    if(lastN>nps)
    //    {
    //      lastN=nps;
    //      lastH=sma;
    //    }
    //    else lastH = lastN*ds;           // Calculate max initialized s for LinTab
    //    for(G4int j=nextN; j<=lastN; j++)// Calculate LogTab values
    //    {
    //      sv+=ds;
    //      lastT[j]=CalcDiff2Prod_Ratio(sv,A);
    //    }
    //  } // End of LinTab update
    //  if(lastN>=nextN)
    //  {
    //    vH[i]=lastH;
    //    vN[i]=lastN;
    //  }
    //  G4int nextK=lastK+1;
    //  if(s_value>sma && lastK<nls)             // LogTab must be updated
    //  {
    //    G4double sv=std::exp(lastM+lsi); // Define starting poit (lastM will be changed)
    //    G4double ls=std::log(s_value);
    //    lastK = static_cast<int>((ls-lsi)/dl)+1; // MaxBin to be initialized in LogTaB
    //    if(lastK>nls)
    //    {
    //      lastK=nls;
    //      lastM=lsa-lsi;
    //    }
    //    else lastM = lastK*dl;           // Calcul. max initialized ln(s)-lsi for LogTab
    //    for(G4int j=nextK; j<=lastK; j++)// Calculate LogTab values
    //    {
    //      sv*=edl;
    //      lastL[j]=CalcDiff2Prod_Ratio(sv,A);
    //    }
    //  } // End of LogTab update
    //  if(lastK>=nextK)
    //  {
    //    vM[i]=lastM;
    //    vK[i]=lastK;
    //  }
    //}
  }
  // Now one can use tabeles to calculate the value
  if(s_value<sma)                             // Use linear table
  {
    G4int n=static_cast<int>(s_value/ds);     // Low edge number of the bin
    G4double d=s_value-n*ds;                  // Linear shift
    G4double v=lastT[n];                      // Base
    lastR=v+d*(lastT[n+1]-v)/ds;              // Result
  }
  else                                  // Use log table
  {
    G4double ls=std::log(s_value)-lsi;  // ln(s)-l_min
    G4int n=static_cast<int>(ls/dl);    // Low edge number of the bin
    G4double d=ls-n*dl;                 // Log shift
    G4double v=lastL[n];                // Base
    lastR=v+d*(lastL[n+1]-v)/dl;        // Result
  }
  if(lastR<0.) lastR=0.;
  if(lastR>1.) lastR=1.;
  return lastR;
} // End of GetRatio

// Calculate Diffraction/Production Ratio as a function of total sq(s)(hN) (in GeV), A=Z+N
G4double G4QDiffractionRatio::CalcDiff2Prod_Ratio(G4double s_value, G4int A)
{
  static G4int    mA=0;
  static G4double S=.1; // s=SQRT(M_N^2+M_h^2+2*E_h*M_N)
  static G4double R=0.; // Prototype of the result
  static G4double p1=0.;
  static G4double p2=0.;
  static G4double p4=0.;
  static G4double p5=0.;
  static G4double p6=0.;
  static G4double p7=0.;
  if(s_value<=0. || A<=1) return 0.;
  if(A!=mA && A!=1)
  {
    S=s_value;
    mA=A;
    G4double a=mA;
    G4double sa=std::sqrt(a);
    G4double a2=a*a;
    G4double a3=a2*a;
    G4double a4=a3*a;
    G4double a5=a4*a;
    G4double a6=a5*a;
    G4double a7=a6*a;
    G4double a8=a7*a;
    G4double a11=a8*a3;
    G4double a12=a8*a4;
    p1=(.023*std::pow(a,0.37)+3.5/a3+2.1e6/a12+4.e-14*a5)/(1.+7.6e-4*a*sa+2.15e7/a11);
    p2=(1.42*std::pow(a,0.61)+1.6e5/a8+4.5e-8*a4)/(1.+4.e-8*a4+1.2e4/a6);
    G4double q=std::pow(a,0.7);
    p4=(.036/q+.0009*q)/(1.+6./a3+1.e-7*a3);
    p5=1.3*std::pow(a,0.1168)/(1.+1.2e-8*a3);
    p6=.00046*(a+11830./a2);
    p7=1./(1.+6.17/a2+.00406*a);
  }
  else if(A==1 && mA!=1)
  {
    S=s_value;
    p1=.0315;
    p2=.73417;
    p4=.01109;
    p5=1.0972;
    p6=.065787;
    p7=.62976;
  }
  else if(std::fabs(s_value-S)/S<.0001) return R;
  G4double s2=s_value*s_value;
  G4double s4=s2*s2;
  G4double dl=std::log(s_value)-p5;
  R=1./(1.+1./(p1+p2/s4+p4*dl*dl/(1.+p6*std::pow(s_value,p7))));
  return R;
} // End of CalcQF2IN_Ratio


G4QHadronVector* G4QDiffractionRatio::TargFragment(G4int pPDG, G4LorentzVector p4M,
                                                   G4int tgZ, G4int tgN)
{
  static const G4double pFm= 0.; // Fermi momentum in MeV (delta function)
  //static const G4double pFm= 250.; // Fermi momentum in MeV (delta function)
  static const G4double pFm2= pFm*pFm; // Squared Fermi momentum in MeV^2 (delta function)
  static const G4double mPi0= G4QPDGCode(111).GetMass(); // pi0 mass (MeV =min diffraction)
  //static const G4double mPi= G4QPDGCode(211).GetMass();  // pi+- mass (MeV)
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mNeut2=mNeut*mNeut;
  static const G4double dmNeut=mNeut+mNeut;
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mProt2=mProt*mProt;
  static const G4double dmProt=mProt+mProt;
  static const G4double maxDM=mProt*12.;
  //static const G4double mLamb= G4QPDGCode(3122).GetMass();
  //static const G4double mSigZ= G4QPDGCode(3212).GetMass();
  //static const G4double mSigM= G4QPDGCode(3112).GetMass();
  //static const G4double mSigP= G4QPDGCode(3222).GetMass();
  //static const G4double eps=.003;
  static const G4double third=1./3.;
  //
  G4LorentzVector pr4M=p4M/megaelectronvolt;   // Convert 4-momenta in MeV (keep p4M)
  // prepare the DONOTHING answer
  G4QHadronVector* ResHV = new G4QHadronVector;// !! MUST BE DESTROYE/DDELETER by CALLER !!
  G4QHadron* hadron = new G4QHadron(pPDG,p4M); // Hadron for not scattered projectile
  ResHV->push_back(hadron);                // It must be cleaned up for real scattering sec
  // @@ diffraction is simulated as noncoherent (coherent is small)
  G4int tgA=tgZ+tgN;                       // A of the target
  G4int tPDG=90000000+tgZ*1000+tgN;        // PDG code of the targetNucleus/recoilNucleus
  G4double tgM=G4QPDGCode(tPDG).GetMass(); // Mass of the target nucleus
  G4int rPDG=2112;                         // prototype of PDG code of the recoiled nucleon
  if(tgA*G4UniformRand()>tgN)              // Substitute by a proton
  {
    rPDG=2212;                             // PDG code of the recoiled QF nucleon
    tPDG-=1000;                            // PDG code of the recoiled nucleus
  }
  else tPDG-=1;                            // PDG code of the recoiled nucleus
  G4double tM=G4QPDGCode(tPDG).GetMass();  // Mass of the recoiled nucleus
  G4double tE=std::sqrt(tM*tM+pFm2);       // Free energy of the recoil nucleus
  G4ThreeVector tP=pFm*G4RandomDirection();// 3-mom of the recoiled nucleus
  G4LorentzVector t4M(tP,tE);              // 4M of the recoil nucleus
  G4LorentzVector tg4M(0.,0.,0.,tgM);      // Full target 4-momentum
  G4LorentzVector N4M=tg4M-t4M;            // 4-mom of Quasi-free target nucleon
  G4LorentzVector tot4M=N4M+p4M;           // total momentum of quasi-free diffraction
  G4double mT=mNeut;                       // Prototype of mass of QF nucleon
  G4double mT2=mNeut2;                     // Squared mass of a free nucleon to be excited
  G4double dmT=dmNeut;                     // Doubled mass              
  //G4int Z=0;                               // Prototype of the isotope Z
  //G4int N=1;                               // Prototype of the Isotope N
  if(rPDG==2212)                           // Correct it, if this is a proton
  {
    mT=mProt;                              // Prototype of mass of QF nucleon to be excited
    mT2=mProt2;                            // Squared mass of the free nucleon
    dmT=dmProt;                            // Doubled mass              
    //Z=1;                                   // Z of the isotope
    //N=0;                                   // N of the Isotope
  }
  G4double mP2=pr4M.m2();                  // Squared mass of the projectile
  if(mP2<0.) mP2=0.;                       // Can be a problem for photon (m_min = 2*m_pi0)
  G4double s_value=tot4M.m2();             // @@ Check <0 ...
  G4double E=(s_value-mT2-mP2)/dmT;        // Effective interactionEnergy (virtNucl target)
  G4double E2=E*E;
  if(E<0. || E2<mP2)                       // Impossible to fragment: return projectile
  {
#ifdef pdebug
    G4cerr<<"-Warning-G4DifR::TFra:<NegativeEnergy>E="<<E<<",E2="<<E2<<"<M2="<<mP2<<G4endl;
#endif
    return ResHV;                          // *** Do Nothing Action ***
  }
  G4double mP=std::sqrt(mP2);              // Calculate mass of the projectile (to be exc.)
  if(mP<.1) mP=mPi0;                       // For photons minDiffraction is gam+P->P+Pi0
  //G4double dmP=mP+mP;                      // Doubled mass of the projectile
  G4double mMin=mP+mPi0;                   // Minimum diffractive mass
  G4double tA=tgA;                         // Real A of the target
  G4double sA=5./std::pow(tA,third);       // Mass-screaning
  //mMin+=mPi0+G4UniformRand()*(mP*sA+mPi0); // *Experimental*
  mMin+=G4UniformRand()*(mP*sA+mPi0);      // *Experimental*
  G4double ss=std::sqrt(s_value);          // CM compound mass (sqrt(s))
  G4double mMax=ss-mP;                     // Maximum diffraction mass of the projectile
  if(mMax>maxDM) mMax=maxDM;               // Restriction to avoid too big masses
  if(mMin>=mMax)
  {
#ifdef pdebug
    G4cerr<<"-Warning-G4DifR::TFra:ZeroDiffractionMRange, mi="<<mMin<<",ma="<<mMax<<G4endl;
#endif
    return ResHV;                          // Do Nothing Action
  }
  G4double R = G4UniformRand();
  G4double mDif=std::exp(R*std::log(mMax)+(1.-R)*std::log(mMin)); // Low Mass Approximation
  G4double mDif2=mDif*mDif;
  G4double ds=s_value-mP2-mDif2;               
  //G4double e=ds/dmP;
  //G4double P=std::sqrt(e*e-mDif2);      // Momentum in pseudo laboratory system
#ifdef debug
  G4cout<<"G4QDiffR::TargFrag:Before XS, P="<<P<<",Z="<<Z<<",N="<<N<<",PDG="<<pPDG<<G4endl;
#endif
  // @@ Temporary NN t-dependence for all hadrons
  if(pPDG>3400 || pPDG<-3400) G4cout<<"-Warning-G4QDifR::Fragment: pPDG="<<pPDG<<G4endl;
  G4double maxt=(ds*ds-4*mP2*mDif2)/s_value;  // maximum possible -t
  G4double tsl=140000.;                 // slope in MeV^2 
  G4double t=-std::log(G4UniformRand())*tsl;
#ifdef pdebug
  G4cout<<"G4QDifR::TFra:ph="<<pPDG<<",P="<<P<<",t="<<t<<"<"<<maxt<<G4endl;
#endif
#ifdef nandebug
  if(mint>-.0000001);                   // To make the Warning for NAN
  else  G4cout<<"******G4QDiffractionRatio::TargFragment: -t="<<mint<<G4endl;
#endif
  G4double rt=t/maxt;
  G4double cost=1.-rt-rt;               // cos(theta) in CMS
#ifdef ppdebug
  G4cout<<"G4QDiffraRatio::TargFragment: -t="<<t<<", maxt="<<maxt<<", cost="<<cost<<G4endl;
#endif
  if(cost>1. || cost<-1. || !(cost>-1. || cost<=1.))
  {
    if     (cost>1.)  cost=1.;
    else if(cost<-1.) cost=-1.;
    else
    {
      G4cerr<<"G4QDiffRat::TargFragm: *NAN* cost="<<cost<<",t="<<t<<",tmax="<<maxt<<G4endl;
      return ResHV;                     // Do Nothing Action
    }
  }
  G4LorentzVector r4M=G4LorentzVector(0.,0.,0.,mP);      // 4mom of the leading nucleon
  G4LorentzVector d4M=G4LorentzVector(0.,0.,0.,mDif);    // 4mom of the diffract. Quasmon
  G4LorentzVector dir4M=tot4M-G4LorentzVector(0.,0.,0.,(tot4M.e()-mT)*.01);
  if(!G4QHadron(tot4M).RelDecayIn2(r4M, d4M, dir4M, cost, cost))
  {
    G4cerr<<"G4QDifR::TFr:M="<<tot4M.m()<<",T="<<mT<<",D="<<mDif<<",T+D="<<mT+mDif<<G4endl;
    //G4Exception("G4QDifR::Fragm:","009",FatalException,"Decay of ElasticComp");
    return ResHV; // Do Nothing Action
  }
#ifdef debug
  G4cout<<"G4QDifRat::TargFragm:d4M="<<d4M<<"+r4M="<<r4M<<"="<<d4M+r4M<<"="<<tot4M<<G4endl;
#endif
  // Now everything is ready for fragmentation and DoNothing projHadron must be wiped out
  delete hadron;     // Delete the fake (doNothing) projectile hadron
  ResHV->pop_back(); // Clean up pointer to the fake (doNothing) projectile
  hadron = new G4QHadron(pPDG,r4M);     // Hadron for the recoil nucleon
  ResHV->push_back(hadron);             // Fill the recoil nucleon
#ifdef debug
  G4cout<<"G4QDiffractionRatio::TargFragm: *Filled* LeadingNuc="<<r4M<<pPDG<<G4endl;
#endif
  G4QHadronVector* leadhs = 0;   // Prototype of Quasmon Output G4QHadronVector  ---->---*
  G4QContent dQC=G4QPDGCode(rPDG).GetQuarkContent(); // QuarkContent of quasiFreeNucleon | 
  G4Quasmon* quasm = new G4Quasmon(dQC,d4M); // Quasmon=DiffractionExcitationQuasmon-*   |
#ifdef debug
  G4cout<<"G4QDiffRatio::TgFrag:tPDG="<<tPDG<<",rPDG="<<rPDG<<",d4M="<<d4M<<G4endl;//|   |
#endif
  G4QEnvironment* pan= new G4QEnvironment(G4QNucleus(tPDG));// --> DELETED --->---*  |   |
  pan->AddQuasmon(quasm);                    // Add diffractiveQuasmon to Environ.|  |   |
#ifdef debug
  G4cout<<"G4QDiffractionRatio::TargFragment: EnvPDG="<<tPDG<<G4endl; //          |  |   |
#endif
  try                                                           //                |  |   |
  {                                                             //                |  |   |
    leadhs = pan->Fragment();// DESTROYED in the end of the LOOP work space       |  | <-|
  }                                                             //                |  |   |
  catch (G4QException& error)//                                                   |  |   |
  {                                                             //                |  |   |
    //#ifdef pdebug
    G4cerr<<"***G4QDiffractionRatio::TargFrag: G4QException is catched"<<G4endl;//|  |   |
    //#endif
    //  G4Exception("G4QDiffractionRatio::TargFragm:","27",FatalException,"*Nucl");// |  |   |
    G4Exception("G4QDiffractionRatio::TargFragment()","HAD_CHPS_0027",
                FatalException, "Nucl");
  }                                                             //                |  |   |
  delete pan;                              // Delete the Nuclear Environment <-<--*--*   |
  G4int qNH=leadhs->size();                // A#of collected hadrons from diff.frag.     |
  if(qNH) for(G4int iq=0; iq<qNH; iq++)    // Loop over hadrons to fill the result       |
  {                                        //                                            |
    G4QHadron* loh=(*leadhs)[iq];          // Pointer to the output hadron               |
    ResHV->push_back(loh);                 // Fill in the result                         |
  }                                        //                                            |
  leadhs->clear();//                                                                     |
  delete leadhs; // <----<----<----<----<----<----<----<----<----<----<----<----<----<---*
  return ResHV; // Result
} // End of TargFragment


G4QHadronVector* G4QDiffractionRatio::ProjFragment(G4int pPDG, G4LorentzVector p4M,
                                                  G4int tgZ, G4int tgN)
{
  static const G4double pFm= 250.; // Fermi momentum in MeV (delta function)
  static const G4double pFm2= pFm*pFm; // Squared Fermi momentum in MeV^2 (delta function)
  static const G4double mPi0= G4QPDGCode(111).GetMass(); // pi0 mass (MeV =min diffraction)
  static const G4double mPi= G4QPDGCode(211).GetMass();  // pi+- mass (MeV)
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mNeut2=mNeut*mNeut;
  static const G4double dmNeut=mNeut+mNeut;
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mProt2=mProt*mProt;
  static const G4double dmProt=mProt+mProt;
  static const G4double maxDM=mProt*12.;
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mSigZ= G4QPDGCode(3212).GetMass();
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
  static const G4double eps=.003;
  //
  G4LorentzVector pr4M=p4M/megaelectronvolt;   // Convert 4-momenta in MeV (keep p4M)
  // prepare the DONOTHING answer
  G4QHadronVector* ResHV = new G4QHadronVector;// !! MUST BE DESTROYE/DDELETER by CALLER !!
  G4QHadron* hadron = new G4QHadron(pPDG,p4M); // Hadron for not scattered projectile
  ResHV->push_back(hadron);                // It must be cleaned up for real scattering sec
  // @@ diffraction is simulated as noncoherent (coherent is small)
  G4int tgA=tgZ+tgN;                       // A of the target
  G4int tPDG=90000000+tgZ*1000+tgN;        // PDG code of the targetNucleus/recoilNucleus
  G4double tgM=G4QPDGCode(tPDG).GetMass(); // Mass of the target nucleus
  G4int rPDG=2112;                         // prototype of PDG code of the recoiled nucleon
  if(tgA*G4UniformRand()>tgN)              // Substitute by a proton
  {
    rPDG=2212;                             // PDG code of the recoiled QF nucleon
    tPDG-=1000;                            // PDG code of the recoiled nucleus
  }
  else tPDG-=1;                            // PDG code of the recoiled nucleus
  G4double tM=G4QPDGCode(tPDG).GetMass();  // Mass of the recoiled nucleus
  G4double tE=std::sqrt(tM*tM+pFm2);
  G4ThreeVector tP=pFm*G4RandomDirection();
  G4LorentzVector t4M(tP,tE);              // 4M of the recoil nucleus
  G4LorentzVector tg4M(0.,0.,0.,tgM);
  G4LorentzVector N4M=tg4M-t4M;            // Quasi-free target nucleon
  G4LorentzVector tot4M=N4M+p4M;           // total momentum of quasi-free diffraction
  G4double mT=mNeut;
  G4double mT2=mNeut2;                 // Squared mass of the free nucleon spectator
  G4double dmT=dmNeut;
  //G4int Z=0;
  //G4int N=1;
  if(rPDG==2212)
  {
    mT=mProt;
    mT2=mProt2;
    dmT=dmProt;
    //Z=1;
    //N=0;
  }
  G4double mP2=pr4M.m2();               // Squared mass of the projectile
  if(mP2<0.) mP2=0.;                    // A possible problem for photon (m_min = 2*m_pi0)
  G4double s_value=tot4M.m2();          // @@ Check <0 ...
  G4double E=(s_value-mT2-mP2)/dmT;     // Effective interactin energy (virt. nucl. target)
  G4double E2=E*E;
  if(E<0. || E2<mP2)
  {
#ifdef pdebug
    G4cerr<<"-Warning-G4DifR::PFra:<NegativeEnergy>E="<<E<<",E2="<<E2<<"<M2="<<mP2<<G4endl;
#endif
    return ResHV; // Do Nothing Action
  }
  G4double mP=std::sqrt(mP2);
  if(mP<.1)mP=mPi0;                      // For photons min diffraction is gamma+P->Pi0+Pi0
  G4double mMin=mP+mPi0;                 // Minimum diffractive mass
  G4double ss=std::sqrt(s_value);        // CM compound mass (sqrt(s))
  G4double mMax=ss-mT;                   // Maximum diffraction mass
  if(mMax>maxDM) mMax=maxDM;             // Restriction to avoid too big masses
  if(mMin>=mMax)
  {
#ifdef pdebug
    G4cerr<<"-Warning-G4DifR::PFra:ZeroDiffractionMRange, mi="<<mMin<<",ma="<<mMax<<G4endl;
#endif
    return ResHV; // Do Nothing Action
  }
  G4double R = G4UniformRand();
  G4double mDif=std::exp(R*std::log(mMax)+(1.-R)*std::log(mMin)); // LowMassApproximation
  G4double mDif2=mDif*mDif;
  G4double ds=s_value-mT2-mDif2;
  //G4double e=ds/dmT;
  //G4double P=std::sqrt(e*e-mDif2);          // Momentum in pseudo laboratory system
#ifdef debug
  G4cout<<"G4QDiffR::PFra: Before XS, P="<<P<<", Z="<<Z<<", N="<<N<<", PDG="<<pPDG<<G4endl;
#endif
  // @@ Temporary NN t-dependence for all hadrons
  if(pPDG>3400 || pPDG<-3400) G4cout<<"-Warning-G4QDifR::Fragment: pPDG="<<pPDG<<G4endl;
  G4double tsl=140000.;                        // slope in MeV^2 
  G4double t=-std::log(G4UniformRand())*tsl;
  G4double maxt=(ds*ds-4*mT2*mDif2)/s_value;   // maximum possible -t
#ifdef pdebug
  G4cout<<"G4QDifR::PFra:ph="<<pPDG<<",P="<<P<<",t="<<mint<<"<"<<maxt<<G4endl;
#endif
#ifdef nandebug
  if(mint>-.0000001);                          // To make the Warning for NAN
  else  G4cout<<"******G4QDiffractionRatio::ProjFragment: -t="<<mint<<G4endl;
#endif
  G4double rt=t/maxt;
  G4double cost=1.-rt-rt;                      // cos(theta) in CMS
#ifdef ppdebug
  G4cout<<"G4QDiffRatio::ProjFragment: -t="<<t<<", maxt="<<maxt<<", cost="<<cost<<G4endl;
#endif
  if(cost>1. || cost<-1. || !(cost>-1. || cost<=1.))
  {
    if     (cost>1.)  cost=1.;
    else if(cost<-1.) cost=-1.;
    else
    {
      G4cerr<<"G4QDiffRat::ProjFragm: *NAN* cost="<<cost<<",t="<<t<<",tmax="<<maxt<<G4endl;
      return ResHV; // Do Nothing Action
    }
  }
  G4LorentzVector r4M=G4LorentzVector(0.,0.,0.,mT);      // 4mom of the recoil nucleon
  G4LorentzVector d4M=G4LorentzVector(0.,0.,0.,mDif);    // 4mom of the diffract. Quasmon
  G4LorentzVector dir4M=tot4M-G4LorentzVector(0.,0.,0.,(tot4M.e()-mT)*.01);
  if(!G4QHadron(tot4M).RelDecayIn2(d4M, r4M, dir4M, cost, cost))
  {
    G4cerr<<"G4QDifR::PFr:M="<<tot4M.m()<<",T="<<mT<<",D="<<mDif<<",T+D="<<mT+mDif<<G4endl;
    //G4Exception("G4QDifR::Fragm:","009",FatalException,"Decay of ElasticComp");
    return ResHV; // Do Nothing Action
  }
#ifdef debug
  G4cout<<"G4QDiffR::ProjFragm:d4M="<<d4M<<"+r4M="<<r4M<<"="<<d4M+r4M<<"="<<tot4M<<G4endl;
#endif
  // Now everything is ready for fragmentation and DoNothing projHadron must be wiped out
  delete hadron;     // Delete the fake (doNothing) projectile hadron
  ResHV->pop_back(); // Clean up pointer to the fake (doNothing) projectile
  hadron = new G4QHadron(tPDG,t4M);  // Hadron for the recoil neucleus
  ResHV->push_back(hadron);          // Fill the recoil nucleus
#ifdef debug
  G4cout<<"G4QDiffractionRatio::ProjFragment: *Filled* RecNucleus="<<t4M<<tPDG<<G4endl;
#endif
  hadron = new G4QHadron(rPDG,r4M);  // Hadron for the recoil nucleon
  ResHV->push_back(hadron);          // Fill the recoil nucleon
#ifdef debug
  G4cout<<"G4QDiffractionRatio::ProjFragment: *Filled* RecNucleon="<<r4M<<rPDG<<G4endl;
#endif
  G4LorentzVector sum4M(0.,0.,0.,0.);
  // Now the (pPdg,d4M) Quasmon must be fragmented
  G4QHadronVector* leadhs = 0;       // Prototype of QuasmOutput G4QHadronVector
  G4QContent dQC=G4QPDGCode(pPDG).GetQuarkContent(); // Quark Content of the projectile
  G4Quasmon* pan= new G4Quasmon(dQC,d4M); // --->---->---->----->-----> DELETED -->---*
  try                                                           //                    |
  {                                                             //                    |
    G4QNucleus vac(90000000);                                   //                    |
    leadhs=pan->Fragment(vac,1);  // DELETED after it is copied to ResHV vector -->---+-*
  }                                                             //                    | |
  catch (G4QException& error)                                   //                    | |
  {                                                             //                    | |
    G4cerr<<"***G4QDiffractionRatio::ProjFragment: G4Quasmon Exception"<<G4endl;    //| |
    // G4Exception("G4QDiffractionRatio::ProjFragment","72",FatalException,"*Quasmon");//| |
    G4Exception("G4QDiffractionRatio::ProjFragment()", "HAD_CHPS_0072",
                FatalException, "*Quasmon");
  }                                                             //                    | |
  delete pan;                              // Delete the Nuclear Environment <----<---* |
  G4int qNH=leadhs->size();                // A#of collected hadrons from diff.frag.    |
  if(qNH) for(G4int iq=0; iq<qNH; iq++)    // Loop over hadrons to fill the result      |
  {                                        //                                           |
    G4QHadron* loh=(*leadhs)[iq];          // Pointer to the output hadron              |
    G4int nL=loh->GetStrangeness();        // A number of Lambdas in the Hypernucleus   |
    G4int nB=loh->GetBaryonNumber();       // Total Baryon Number of the Hypernucleus   |
    G4int nC = loh->GetCharge();           // Charge of the Hypernucleus                |
    G4int oPDG = loh->GetPDGCode();        // Original CHIPS PDG Code of the hadron     |
    //if((nC>nB || nC<0) && nB>0 && nL>=0 && nL<=nB && oPDG>80000000) // Iso-nucleus    |
    if(2>3) // Closed because "G4QDR::F:90002999,M=-7.768507e-04,B=2,S=0,C=3" is found  |
    {
      G4LorentzVector q4M = loh->Get4Momentum(); // Get 4-momentum of the Isonucleus    |
      G4double qM=q4M.m();                 // Real mass of the Isonucleus
#ifdef fdebug
      G4cout<<"G4QDR::PF:"<<oPDG<<",M="<<qM<<",B="<<nB<<",S="<<nL<<",C="<<nC<<G4endl;// |
#endif
      G4int    qPN=nC-nB;                  // Number of pions in the Isonucleus         |
      G4int    fPDG = 2212;                // Prototype for nP+(Pi+) case               |
      G4int    sPDG = 211;
      tPDG = 3122;                         // @@ Sigma0 (?)                             |
      G4double fMass= mProt;
      G4double sMass= mPi;
      G4double tMass= mLamb;               // @@ Sigma0 (?)                             |
      G4bool   cont=true;                  // Continue flag                             |
      // =--------= Negative state =---------=
      if(nC<0)                             // =----= Only Pi- can help                  |
      {
        if(nL&&nB==nL)                     // --- n*Lamb + k*(Pi-) State ---            |
        {
          sPDG = -211;
          if(-nC==nL && nL==1)             // Only one Sigma- like (nB=1)               |
          {
            if(std::fabs(qM-mSigM)<eps)
            {
              loh->SetQPDG(G4QPDGCode(3112));  // This is Sigma-                        |
              cont=false;                  // Skip decay                                |
            }
            else if(qM>mLamb+mPi)          //(2) Sigma- => Lambda + Pi- decay           |
            {
              fPDG = 3122;
              fMass= mLamb;
            }
            else if(qM>mSigM)              //(2) Sigma+=>Sigma++gamma decay             |
            {
              fPDG = 3112;
              fMass= mSigM;
              sPDG = 22;
              sMass= 0.;
            }
            else                           //(2) Sigma-=>Neutron+Pi- decay              |
            {
              fPDG = 2112;
              fMass= mNeut;
            }
            qPN  = 1;                      // #of (Pi+ or gamma)'s = 1                  |
          }
          else if(-nC==nL)                 //(2) a few Sigma- like                      |
          {
            qPN  = 1;                      // One separated Sigma-                      |
            fPDG = 3112;
            sPDG = 3112;
            sMass= mSigM;
            nB--;
            fMass= mSigM;
          }
          else if(-nC>nL)                  //(2) n*(Sigma-)+m*(Pi-)                     |
          {
            qPN  = -nC-nL;                 // #of Pi-'s                                 |
            fPDG = 3112;
            fMass= mSigM;
          }
          else                             //(2) n*(Sigma-)+m*Lambda(-nC<nL)            |
          {
            nB += nC;                      // #of Lambda's                              |
            fPDG = 3122;
            fMass= mLamb;
            qPN  = -nC;                    // #of Sigma+'s                              |
            sPDG = 3112;
            sMass= mSigM;
          }
          nL   = 0;                        // Only decays in two are above              |
        }
        else if(nL)                        // ->n*Lamb+m*Neut+k*(Pi-) State (nL<nB)     |
        {
          nB -= nL;                        // #of neutrons                              |
          fPDG = 2112;
          fMass= mNeut;
          G4int nPin = -nC;                           // #of Pi-'s                    
          if(nL==nPin)                                //(2) m*Neut+n*Sigma-             |
          {
            qPN  = nL;                                // #of Sigma-                     |
            sPDG = 3112;
            sMass= mSigM;
            nL   = 0;
          }
          else if(nL>nPin)                            //(3) m*P+n*(Sigma+)+k*Lambda     |
          {
            nL-=nPin;                                 // #of Lambdas                    |
            qPN  = nPin;                              // #of Sigma+                     |
            sPDG = 3112;
            sMass= mSigM;
          }
          else                                 //(3) m*N+n*(Sigma-)+k*(Pi-) (nL<nPin)   |
          {
            qPN  = nPin-nL;                           // #of Pi-                        |
            sPDG = -211;
            tPDG = 3112;
            tMass= mSigM;
          }
        }
        else                                          //(2) n*N+m*(Pi-)   (nL=0)        |
        {
          sPDG = -211;
          qPN  = -nC;
          fPDG = 2112;
          fMass= mNeut;
        }
      }
      else if(!nC)                                   // *** Should not be here ***      |
      {
        if(nL && nL<nB)          //(2) n*Lamb+m*N ***Should not be here***              |
        {
          qPN  = nL;
          fPDG = 2112;                               // mN+nL case                      |
          sPDG = 3122;
          sMass= mLamb;
          nB -= nL;
          fMass= mNeut;
          nL   = 0;
        }
        else if(nL>1 && nB==nL)  //(2) m*Lamb(m>1) ***Should not be here***             |
        {
          qPN  = 1;
          fPDG = 3122;
          sPDG = 3122;
          sMass= mLamb;
          nB--;
          fMass= mLamb;
        }
        else if(!nL && nB>1)     //(2) n*Neut(n>1) ***Should not be here***             |
        {
          qPN  = 1;
          fPDG = 2112;
          sPDG = 2112;
          sMass= mNeut;
          nB--;
          fMass= mNeut;
        }
        else G4cout<<"*?*G4QDiffractionRatio::ProjFragment: (1) oPDG="<<oPDG<<G4endl;// |
      }
      else if(nC>0)              // n*Lamb+(m*P)+(k*Pi+)                                |
      {
        if(nL && nL+nC==nB)      //(2) n*Lamb+m*P ***Should not be here***              |
        {
          qPN  = nL;
          nL   = 0;
          fPDG = 2212;
          sPDG = 3122;
          sMass= mLamb;
          nB  = nC;
          fMass= mProt;
        }
        else if(nL  && nC<nB-nL) //(3)n*L+m*P+k*N ***Should not be here***              |
        {
          qPN  = nC;                                  // #of protons                    |
          fPDG = 2112;                                // mP+nL case                     |
          sPDG = 2212;
          sMass= mProt;
          nB -= nL+nC;                                // #of neutrons                   |
          fMass= mNeut;
        }
        else if(nL  && nB==nL)                        // ---> n*L+m*Pi+ State           |
        {
          if(nC==nL && nL==1)                         // Only one Sigma+ like State     |
          {
            if(std::fabs(qM-mSigP)<eps)
            {
              loh->SetQPDG(G4QPDGCode(3222));         // This is GS Sigma+              |
              cont=false;                  // Skip decay                                |
            }
            else if(qM>mLamb+mPi)                     //(2) Sigma+=>Lambda+Pi+ decay    |
            {
              fPDG = 3122;
              fMass= mLamb;
            }
            else if(qM>mNeut+mPi)                     //(2) Sigma+=>Neutron+Pi+ decay   |
            {
              fPDG = 2112;
              fMass= mNeut;
            }
            else if(qM>mSigP)                         //(2) Sigma+=>Sigma++gamma decay  |
            {
              fPDG = 3222;
              fMass= mSigP;
              sPDG = 22;
              sMass= 0.;
            }
            else                                      //(2) Sigma+=>Proton+gamma decay  |
            {
              fPDG = 2212;
              fMass= mProt;
              sPDG = 22;
              sMass= 0.;
            }
            qPN  = 1;                                 // #of (Pi+ or gamma)'s = 1       |
          }
          else if(nC==nL)                             //(2) a few Sigma+ like hyperons  |
          {
            qPN  = 1;
            fPDG = 3222;
            sPDG = 3222;
            sMass= mSigP;
            nB--;
            fMass= mSigP;
          }
          else if(nC>nL)                              //(2) n*(Sigma+)+m*(Pi+)          |
          {
            qPN  = nC-nL;                             // #of Pi+'s                      |
            fPDG = 3222;
            nB  = nL;                                 // #of Sigma+'s                   |
            fMass= mSigP;
          }
          else                                        //(2) n*(Sigma+)+m*Lambda         |
          {
            nB -= nC;                                 // #of Lambda's                   |
            fPDG = 3122;
            fMass= mLamb;
            qPN  = nC;                                // #of Sigma+'s                   |
            sPDG = 3222;
            sMass= mSigP;
          }
          nL   = 0;                                   // All above are decays in 2      |
        }
        else if(nL && nC>nB-nL)                       // n*Lamb+m*P+k*Pi+               |
        {
          nB -= nL;                                   // #of protons                    |
          G4int nPip = nC-nB;                         // #of Pi+'s                      |
          if(nL==nPip)                                //(2) m*P+n*Sigma+                |
          {
            qPN  = nL;                                // #of Sigma+                     |
            sPDG = 3222;
            sMass= mSigP;
            nL   = 0;
          }
          else if(nL>nPip)                            //(3) m*P+n*(Sigma+)+k*Lambda     |
          {
            nL  -= nPip;                              // #of Lambdas                    |
            qPN  = nPip;                              // #of Sigma+                     |
            sPDG = 3222;
            sMass= mSigP;
          }
          else                                        //(3) m*P+n*(Sigma+)+k*(Pi+)      |
          {
            qPN  = nPip-nL;                           // #of Pi+                        |
            tPDG = 3222;
            tMass= mSigP;
          }
        }
        if(nC<nB)                 //(2) n*P+m*N ***Should not be here***                |
        {
          fPDG = 2112;
          fMass= mNeut;
          qPN  = nC;
          sPDG = 2212;
          sMass= mProt;
        }
        else if(nB==nC && nC>1)   //(2) m*Prot(m>1) ***Should not be here***            |
        {
          qPN  = 1;
          fPDG = 2212;
          sPDG = 2212;
          sMass= mProt;
          nB--;
          fMass= mProt;
        }
        else if(nC<=nB||!nB) G4cout<<"*?*G4QDR::ProjFragm: (2) oPDG="<<oPDG<<G4endl; // |
        // !nL && nC>nB                             //(2) Default condition n*P+m*(Pi+) |
      }
      if(cont)                                      // Make a decay                     |
      {
        G4double tfM=nB*fMass;
        G4double tsM=qPN*sMass;
        G4double ttM=0.;
        if(nL) ttM=nL*tMass;
        G4LorentzVector f4Mom(0.,0.,0.,tfM);
        G4LorentzVector s4Mom(0.,0.,0.,tsM);
        G4LorentzVector t4Mom(0.,0.,0.,ttM);
        G4double sum=tfM+tsM+ttM;
        if(std::fabs(qM-sum)<eps)
        {
          f4Mom=q4M*(tfM/sum);
          s4Mom=q4M*(tsM/sum);
          if(nL) t4Mom=q4M*(ttM/sum);
        }
        else if(!nL && (qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))) // Error     |
        {
          //#ifdef fdebug
          G4cout<<"***G4QDR::PrFragm:fPDG="<<fPDG<<"*"<<nB<<"(fM="<<fMass<<")+sPDG="<<sPDG
                <<"*"<<qPN<<"(sM="<<sMass<<")"<<"="<<sum<<" > TM="<<qM<<q4M<<oPDG<<G4endl;
          //#endif
          // throw G4QException("*G4QDiffractionRatio::ProjFragment: Bad decay in 2"); //  |
          G4ExceptionDescription ed;
          ed << "***G4QDR::PrFragm:fPDG=" << fPDG << "*" << nB << "(fM="
             << fMass << ")+sPDG=" << sPDG << "*" << qPN << "(sM=" << sMass
             << ")" << "=" << sum << " > TM=" << qM << q4M << oPDG << G4endl;
          G4Exception("G4QDiffractionRatio::ProjFragment()", "HAD_CHPS_0002",
                      FatalException, ed);
        }
        else if(nL && (qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom)))// Error|
        {
          //#ifdef fdebug
          G4cout<<"***G4DF::PrFrag: "<<fPDG<<"*"<<nB<<"("<<fMass<<")+"<<sPDG<<"*"<<qPN<<"("
                <<sMass<<")+Lamb*"<<nL<<"="<<sum<<" > TotM="<<qM<<q4M<<oPDG<<G4endl;
          //#endif
          // throw G4QException("*G4QDiffractionRatio::ProjFragment: Bad decay in 3"); //  |
          G4ExceptionDescription ed;
          ed << "***G4DF::PrFrag: " << fPDG << "*" << nB << "(" << fMass << ")+"
             << sPDG << "*" << qPN << "(" << sMass << ")+Lamb*" << nL << "="
             << sum << " > TotM=" << qM << q4M << oPDG << G4endl;
          G4Exception("G4QDiffractionRatio::ProjFragment()", "HAD_CHPS_0003",
                      FatalException, ed);
        }
#ifdef fdebug
        G4cout<<"G4QDF::ProjFragm: *DONE* n="<<nB<<f4Mom<<fPDG<<", m="<<qPN<<s4Mom<<sPDG
              <<", l="<<nL<<t4Mom<<G4endl;
#endif
        G4bool notused=true;
        if(nB)                               // There are baryons                       |
        {
          f4Mom/=nB;
          loh->Set4Momentum(f4Mom);          // ! Update the Hadron !                   |
          loh->SetQPDG(G4QPDGCode(fPDG));    // Baryons                                 |
          notused=false;                     // Loh was used                            |
          if(nB>1) for(G4int ih=1; ih<nB; ih++) // Loop over the rest of baryons        |
          {
            G4QHadron* Hi = new G4QHadron(fPDG,f4Mom); // Create a Hadron for Baryon    |
            ResHV->push_back(Hi);            // Fill in the additional nucleon          |
#ifdef fdebug
            sum4M+=r4M;                      // Sum 4-momenta for the EnMom check       |
            G4cout<<"G4QDR::ProjFrag: *additional Nucleon*="<<f4Mom<<fPDG<<G4endl; //   |
#endif
          }
        }
        if(qPN)                              // There are pions                         |
        {
          s4Mom/=qPN;
          G4int min=0;
          if(notused)
          {
            loh->Set4Momentum(s4Mom);        // ! Update the Hadron 4M !                |
            loh->SetQPDG(G4QPDGCode(sPDG));  // Update PDG                              |
            notused=false;                   // loh was used                            |
            min=1;                           // start value                             |
          }
          if(qPN>min) for(G4int ip=min; ip<qPN; ip++) // Loop over pions                |
          {
            G4QHadron* Hj = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the meson |
            ResHV->push_back(Hj);            // Fill in the additional pion             |
#ifdef fdebug
            sum4M+=r4M;                      // Sum 4-momenta for the EnMom check       |
            G4cout<<"G4QDR::ProjFragm: *additional Pion*="<<f4Mom<<fPDG<<G4endl; //     |
#endif
          }
        }
        if(nL)                               // There are Hyperons                      |
        {
          t4Mom/=nL;
          G4int min=0;
          if(notused)
          {
            loh->Set4Momentum(t4Mom);      // ! Update the Hadron 4M !                  |
            loh->SetQPDG(G4QPDGCode(tPDG));// Update PDG                                |
            notused=false;                 // loh was used                              |
            min=1;                         //   
          }
          if(nL>min) for(G4int il=min; il<nL; il++) // Loop over Hyperons               |
          {
            G4QHadron* Hk = new G4QHadron(tPDG,t4Mom); // Create a Hadron for Lambda    |
            ResHV->push_back(Hk);          // Fill in the additional pion               |
#ifdef fdebug
            sum4M+=r4M;                    // Sum 4-momenta for the EnMom check         |
            G4cout<<"G4QDR::ProjFragm: *additional Hyperon*="<<f4Mom<<fPDG<<G4endl; //  |
#endif
          }
        }
      }                                    // --> End of decay                          |
    }                                      // -> End of Iso-nuclear treatment           |
    else if( (nL > 0 && nB > 1) || (nL < 0 && nB < -1) ) 
    {     // Hypernucleus is found                                                      |
      G4bool anti=false;                   // Default=Nucleus (true=antinucleus         |
      if(nB<0)                             // Anti-nucleus                              |
      {
        anti=true;                         // Flag of anti-hypernucleus                 |
        nB=-nB;                            // Reverse the baryon number                 |
        nC=-nC;                            // Reverse the charge                        |
        nL=-nL;                            // Reverse the strangeness                   |
      }
      G4int hPDG = 90000000+nL*999999+nC*999+nB; // CHIPS PDG Code for Hypernucleus     |
      G4int nSM=0;                         // A#0f unavoidable Sigma-                   |
      G4int nSP=0;                         // A#0f unavoidable Sigma+                   |
      if(nC<0)                             // Negative hypernucleus                     |
      {
        if(-nC<=nL)                        // Partial compensation by Sigma-            |
        {
          nSM=-nC;                         // Can be compensated by Sigma-              |
          nL+=nC;                          // Reduce the residual strangeness           |
        }
        else                               // All Charge is compensated by Sigma-       |
        {
          nSM=nL;                          // The maximum number of Sigma-              |
          nL=0;                            // Kill the residual strangeness             |
        }
      }
      else if(nC>nB-nL)                    // Extra positive hypernucleus               |
      {
        if(nC<=nB)                         // Partial compensation by Sigma+            |
        {
          G4int dH=nB-nC;                  // Isotopic shift                            |
          nSP=nL-dH;                       // Can be compensated by Sigma+              |
          nL=dH;                           // Reduce the residual strangeness           |
        }
        else                               // All Charge is compensated by Sigma+       |
        {
          nSP=nL;                          // The maximum number of Sigma+              |
          nL=0;                            // Kill the residual strangeness             |
        }
      }
      r4M=loh->Get4Momentum();                 // Real 4-momentum of the hypernucleus   !
      G4double reM=r4M.m();                    // Real mass of the hypernucleus         |
#ifdef fdebug
      G4cout<<"G4QDiffRatio::PrFrag:oPDG=="<<oPDG<<",hPDG="<<hPDG<<",M="<<reM<<G4endl;//|
#endif
      G4int rlPDG=hPDG-nL*1000000-nSP*1000999-nSM*999001;// Subtract Lamb/Sig from Nucl.|
      G4int    sPDG=3122;                      // Prototype for the Hyperon PDG (Lambda)|
      G4double MLa=mLamb;                      // Prototype for one Hyperon decay       |
#ifdef fdebug
      G4cout<<"G4QDiffRatio::PrFrag:*G4*nS+="<<nSP<<",nS-="<<nSM<<",nL="<<nL<<G4endl;// |
#endif
      if(nSP||nSM)                         // Sigma+/- improvement                      |
      {
        if(nL)                             // By mistake Lambda improvement is found    |
        {
          G4cout<<"***G4QDR::PFr:HypN="<<hPDG<<": bothSigm&Lamb -> ImproveIt"<<G4endl;//|
          //throw G4QException("*G4QDiffractionRatio::Fragment:BothLambda&SigmaInHN");//|
          // @@ Correction, which does not conserv the charge !! (-> add decay in 3)    |
          if(nSP) nL+=nSP;                 // Convert Sigma+ to Lambda                  |
          else    nL+=nSM;                 // Convert Sigma- to Lambda                  |
        }
        if(nSP)                            // Sibma+ should be decayed                  |
        {
          nL=nSP;                          // #of decaying hyperons                     |
          sPDG=3222;                       // PDG code of decaying hyperons             |
          MLa=mSigP;                       // Mass of decaying hyperons                 |
        }
        else                               // Sibma+ should be decayed                  |
        {
          nL=nSM;                          // #of decaying hyperons                     |
          sPDG=3112;                       // PDG code of decaying hyperons             |
          MLa=mSigM;                       // Mass of decaying hyperons                 |
        }
      }
#ifdef fdebug
      G4cout<<"G4QDiffRat::ProjFrag:*G4*mS="<<MLa<<",sPDG="<<sPDG<<",nL="<<nL<<G4endl;//|
#endif
      if(nL>1) MLa*=nL;                    // Total mass of the decaying hyperons       |
      G4double rlM=G4QNucleus(rlPDG).GetMZNS();// Mass of the NonstrangeNucleus         |
      if(!nSP&&!nSM&&nL==1&&reM>rlM+mSigZ&&G4UniformRand()>.5) // Conv Lambda->Sigma0   |
      {
        sPDG=3212;                         // PDG code of a decaying hyperon            |
        MLa=mSigZ;                         // Mass of the decaying hyperon              |
      }
      G4int rnPDG = hPDG-nL*999999;        // Convert Lambdas to neutrons (for convInN) |
      G4QNucleus rnN(rnPDG);               // New nonstrange nucleus                    |
      G4double rnM=rnN.GetMZNS();          // Mass of the new nonstrange nucleus        |
      // @@ In future take into account Iso-Hypernucleus (Add PI+,R & Pi-,R decays)     |
      if(rlPDG==90000000)                  // Multy Hyperon (HyperNuc of only hyperons) |
      {
        if(nL>1) r4M=r4M/nL;               // split the 4-mom for the MultyLambda       |
        for(G4int il=0; il<nL; il++)       // loop over Lambdas                         |
        {
          if(anti) sPDG=-sPDG;             // For anti-nucleus case                     |
          G4QHadron* theLam = new G4QHadron(sPDG,r4M); // Make NewHadr for the Hyperon  |
          ResHV->push_back(theLam);        // Fill in the Lambda                        |
#ifdef fdebug
          sum4M+=r4M;                      // Sum 4-momenta for the EnMom check         |
          G4cout<<"G4QDR::ProjFrag: *additional Lambda*="<<r4M<<sPDG<<G4endl; //        |
#endif
        }
      }
      else if(reM>rlM+MLa-eps)              // Lambda (or Sigma) can be split           |
      {
        G4LorentzVector n4M(0.,0.,0.,rlM);  // 4-mom of the residual nucleus            |
        G4LorentzVector h4M(0.,0.,0.,MLa);  // 4-mom of the Hyperon                     |
        G4double sum=rlM+MLa;               // Safety sum                               |
        if(std::fabs(reM-sum)<eps)          // At rest in CMS                           |
        {
          n4M=r4M*(rlM/sum);                // Split tot 4-mom for resNuc               |
          h4M=r4M*(MLa/sum);                // Split tot 4-mom for Hyperon              |
        }
        else if(reM<sum || !G4QHadron(r4M).DecayIn2(n4M,h4M)) // Error in decay         |
        {
          G4cerr<<"***G4QDF::PF:HypN,M="<<reM<<"<A+n*L="<<sum<<",d="<<sum-reM<<G4endl;//|
          // throw G4QException("***G4QDiffractionRatio::ProjFragment:HypernuclusDecay");//|
          G4Exception("G4QDiffractionRatio::ProjFragment()", "HAD_CHPS_0100",
                      FatalException, "Error in hypernuclus decay");
        }
#ifdef fdebug
        G4cout<<"*G4QDR::PF:HypN="<<r4M<<"->A="<<rlPDG<<n4M<<",n*L="<<nL<<h4M<<G4endl;//|
#endif
        loh->Set4Momentum(n4M);            // ! Update the Hadron !                     |
        if(anti && rlPDG==90000001) rlPDG=-2112; // Convert to anti-neutron             |
        if(anti && rlPDG==90001000) rlPDG=-2212; // Convert to anti-proton              |
        loh->SetQPDG(G4QPDGCode(rlPDG));   // ConvertedHypernucleus to nonstrange(@anti)|
        if(rlPDG==90000002)                // Additional action with loH changed to 2n  |
        {
          G4LorentzVector newLV=n4M/2.;    // Split 4-momentum                          |
          loh->Set4Momentum(newLV);        // Reupdate the hadron                       |
          if(anti) loh->SetQPDG(G4QPDGCode(-2112)); // Make anti-neutron PDG            |
          else loh->SetQPDG(G4QPDGCode(2112)); // Make neutron PDG                      |
          G4QHadron* secHadr = new G4QHadron(loh); // Duplicate the neutron             |
          ResHV->push_back(secHadr);       // Fill in the additional neutron            |
#ifdef fdebug
          sum4M+=r4M;                      // Sum 4-momenta for the EnMom check         |
          G4cout<<"G4QDR::ProgFrag: *additional Neutron*="<<r4M<<sPDG<<G4endl; //       |
#endif
        }
        else if(rlPDG==90002000)           // Additional action with loH change to 2p   |
        {
          G4LorentzVector newLV=n4M/2.;    // Split 4-momentum                          |
          loh->Set4Momentum(newLV);        // Reupdate the hadron                       |
          if(anti) loh->SetQPDG(G4QPDGCode(-2212)); // Make anti-neutron PDG            |
          else loh->SetQPDG(G4QPDGCode(2112)); // Make neutron PDG                      |
          G4QHadron* secHadr = new G4QHadron(loh); // Duplicate the proton              |
          ResHV->push_back(secHadr);       // Fill in the additional neutron            |
#ifdef fdebug
          sum4M+=r4M;                      // Sum 4-momenta for the EnMom check         |
          G4cout<<"G4QDR::ProjFrag: *additional Proton*="<<r4M<<sPDG<<G4endl; //        |
#endif
        }
        // @@(?) Add multybaryon decays if necessary (Now it anyhow is made later)      |
#ifdef fdebug
        G4cout<<"*G4QDiffractionRatio::PrFrag:resNucPDG="<<loh->GetPDGCode()<<G4endl;// |
#endif
        if(nL>1) h4M=h4M/nL;               // split the lambda's 4-mom if necessary     |
        for(G4int il=0; il<nL; il++)       // A loop over excessive hyperons            |
        {
          if(anti) sPDG=-sPDG;             // For anti-nucleus case                     |
          G4QHadron* theLamb = new G4QHadron(sPDG,h4M); // Make NewHadr for the Hyperon |
          ResHV->push_back(theLamb);       // Fill in the additional neutron            |
#ifdef fdebug
          sum4M+=r4M;                      // Sum 4-momenta for the EnMom check         |
          G4cout<<"G4QDR::ProjFrag: *additional Hyperon*="<<r4M<<sPDG<<G4endl; //       |
#endif
        }
      }
      else if(reM>rnM+mPi0-eps&&!nSP&&!nSM)// Lambda->N only if Sigmas are absent       |
      {
        G4int nPi=static_cast<G4int>((reM-rnM)/mPi0); // Calc. pion multiplicity        |
        if (nPi>nL) nPi=nL;                // Cut the pion multiplicity                 |
        G4double npiM=nPi*mPi0;            // Total pion mass                           |
        G4LorentzVector n4M(0.,0.,0.,rnM); // Residual nucleus 4-momentum               |
        G4LorentzVector h4M(0.,0.,0.,npiM);// 4-momentum of pions                       |
        G4double sum=rnM+npiM;             // Safety sum                                |
        if(std::fabs(reM-sum)<eps)         // At rest                                   |
        {
          n4M=r4M*(rnM/sum);               // The residual nucleus part                 |
          h4M=r4M*(npiM/sum);              // The pion part                             |
        }
        else if(reM<sum || !G4QHadron(r4M).DecayIn2(n4M,h4M)) // Error in decay         |
        {
          G4cerr<<"*G4QDR::PF:HypN,M="<<reM<<"<A+n*Pi0="<<sum<<",d="<<sum-reM<<G4endl;//|
          // throw G4QException("***G4QDiffractionRatio::ProjFragment:HypernuclDecay"); // |
          G4Exception("G4QDiffractionRatio::ProjFragment()", "HAD_CHPS_0101",
                       FatalException, "Error in HypernuclDecay"); 
        }
        loh->Set4Momentum(n4M);            // ! Update the Hadron !                     |
        if(anti && rnPDG==90000001) rnPDG=-2112; // Convert to anti-neutron             |
        if(anti && rnPDG==90001000) rnPDG=-2212; // Convert to anti-proton              |
        loh->SetQPDG(G4QPDGCode(rnPDG));   // convert hyperNuc to nonstrangeNuc(@@anti) |
#ifdef fdebug
        G4cout<<"*G4QDR::PF:R="<<r4M<<"->A="<<rnPDG<<n4M<<",n*Pi0="<<nPi<<h4M<<G4endl;//|
#endif
        if(nPi>1) h4M=h4M/nPi;             // Split the 4-mom if necessary              |
        for(G4int ihn=0; ihn<nPi; ihn++)   // A loop over additional pions              |
        {
          G4QHadron* thePion = new G4QHadron(111,h4M); // Make a New Hadr for the pi0   |
          ResHV->push_back(thePion);       // Fill in the Pion                          |
#ifdef fdebug
          sum4M+=r4M;                      // Sum 4-momenta for the EnMom check         |
          G4cout<<"G4QDR::ProjFrag: *additional Pion*="<<r4M<<sPDG<<G4endl; //          |
#endif
        }
        if(rnPDG==90000002)                // Additional action with loH change to 2n   |
        {
          G4LorentzVector newLV=n4M/2.;    // Split 4-momentum                          |
          loh->Set4Momentum(newLV);        // Reupdate the hadron                       |
          if(anti) loh->SetQPDG(G4QPDGCode(-2112)); // Make anti-neutron PDG            |
          else loh->SetQPDG(G4QPDGCode(2112)); // Make neutron PDG                      |
          G4QHadron* secHadr = new G4QHadron(loh); // Duplicate the neutron             |
          ResHV->push_back(secHadr);       // Fill in the additional neutron            |
#ifdef fdebug
          sum4M+=r4M;                      // Sum 4-momenta for the EnMom check         |
          G4cout<<"G4QDR::ProjFrag: *additional Neutron*="<<r4M<<sPDG<<G4endl; //       |
#endif
        }
        else if(rnPDG==90002000)           // Additional action with loH change to 2p   |
        {
          G4LorentzVector newLV=n4M/2.;    // Split 4-momentum                          |
          loh->Set4Momentum(newLV);        // Reupdate the hadron                       |
          if(anti) loh->SetQPDG(G4QPDGCode(-2212)); // Make anti-neutron PDG            |
          else loh->SetQPDG(G4QPDGCode(2112)); // Make neutron PDG                      |
          G4QHadron* secHadr = new G4QHadron(loh); // Duplicate the proton              |
          ResHV->push_back(secHadr);       // Fill in the additional neutron            |
#ifdef fdebug
          sum4M+=r4M;                      // Sum 4-momenta for the EnMom check         |
          G4cout<<"G4QDR::ProjFrag: *additional Proton*="<<r4M<<sPDG<<G4endl; //        |
#endif
        }
        // @@ Add multybaryon decays if necessary                                       |
      }
      else // If this Excepton shows up (lowProbable appearance) => include gamma decay |
      {
        G4double d=rlM+MLa-reM;            // Hyperon Excessive energy                  |
        G4cerr<<"G4QDR::PF:R="<<rlM<<",S+="<<nSP<<",S-="<<nSM<<",L="<<nL<<",d="<<d<<G4endl;
        d=rnM+mPi0-reM;                    // Pion Excessive energy                     |
        G4cerr<<"G4QDR::PF:"<<oPDG<<","<<hPDG<<",M="<<reM<<"<"<<rnM+mPi0<<",d="<<d<<G4endl;
        // throw G4QException("G4QDiffractionRatio::ProjFragment: Hypernuclear conver");// |
        G4Exception("G4QDiffractionRatio::ProjFragment()", "HAD_CHPS_0102",
                    FatalException, "Excessive hypernuclear energy");
      }
    }                                      // => End of G4 Hypernuclear decay           |
    ResHV->push_back(loh);                 // Fill in the result                        |
#ifdef debug
    sum4M+=loh->Get4Momentum();            // Sum 4-momenta for the EnMom check         |
    G4cout<<"G4QDR::PrFra:#"<<iq<<","<<loh->Get4Momentum()<<loh->GetPDGCode()<<G4endl;//|
#endif
  }                                        //                                           |
  leadhs->clear();//                                                                    |
  delete leadhs; // <----<----<----<----<----<----<----<----<----<----<----<----<----<--*
#ifdef debug
  G4cout<<"G4QDiffractionRatio::ProjFragment: *End* Sum="<<sum4M<<" =?= d4M="<<d4M<<G4endl;
#endif
  return ResHV; // Result
} // End of ProjFragment

// Calculates Single Diffraction Taarget Excitation Cross-Section (independent Units)
G4double G4QDiffractionRatio::GetTargSingDiffXS(G4double pIU, G4int pPDG, G4int Z, G4int N)
{
  G4double mom=pIU/gigaelectronvolt;    // Projectile momentum in GeV
  if ( mom < 1. || (pPDG != 2212 && pPDG != 2112) )
    G4cerr<<"G4QDiffractionRatio::GetTargSingDiffXS isn't applicable p="<<mom<<" GeV, PDG="
         <<pPDG<<G4endl;
  G4double A=Z+N;                        // A of the target
  //return 4.5*std::pow(A,.364)*millibarn; // Result
  return 3.7*std::pow(A,.364)*millibarn; // Result after mpi0 correction

} // End of ProjFragment
