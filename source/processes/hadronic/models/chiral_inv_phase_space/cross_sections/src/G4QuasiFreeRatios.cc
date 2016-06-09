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
// G4 Physics class: G4QuasiFreeRatios for N+A elastic cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 28-Oct-11
// 
// ----------------------------------------------------------------------
// Short description: Provides percentage of quasi-free and quasi-elastic
// reactions in the inelastic reactions.
// ----------------------------------------------------------------------

//#define debug
//#define pdebug
//#define ppdebug
//#define nandebug

#include "G4QuasiFreeRatios.hh"
#include "G4SystemOfUnits.hh"

// initialisation of statics
std::vector<G4double*> G4QuasiFreeRatios::vT; // Vector of pointers to LinTable in C++ heap
std::vector<G4double*> G4QuasiFreeRatios::vL; // Vector of pointers to LogTable in C++ heap

G4QuasiFreeRatios::G4QuasiFreeRatios()
{
#ifdef pdebug
  G4cout<<"***^^^*** G4QuasiFreeRatios singletone is created ***^^^***"<<G4endl;
#endif
  QFreeScat=G4QFreeScattering::GetPointer();
}

G4QuasiFreeRatios::~G4QuasiFreeRatios()
{
  std::vector<G4double*>::iterator pos;
  for(pos=vT.begin(); pos<vT.end(); ++pos) delete [] *pos;
  vT.clear();
  for(pos=vL.begin(); pos<vL.end(); ++pos) delete [] *pos;
  vL.clear();
}

// Returns Pointer to the G4VQCrossSection class
G4QuasiFreeRatios* G4QuasiFreeRatios::GetPointer()
{
  static G4QuasiFreeRatios theRatios;   // *** Static body of the QEl Cross Section ***
  return &theRatios;
}

// Calculation of pair(QuasiFree/Inelastic,QuasiElastic/QuasiFree)
std::pair<G4double,G4double> G4QuasiFreeRatios::GetRatios(G4double pIU, G4int pPDG,
                                                           G4int tgZ,    G4int tgN)
{
#ifdef pdebug
  G4cout<<">>>IN>>>G4QFRat::GetQF:P="<<pIU<<",pPDG="<<pPDG<<",Z="<<tgZ<<",N="<<tgN<<G4endl;
#endif
  G4double R=0.;
  G4double QF2In=1.;                        // Prototype of QuasiFree/Inel ratio for hN_tot
  G4int tgA=tgZ+tgN;
  if(tgA<2) return std::make_pair(QF2In,R); // No quasi-elastic on the only nucleon
  std::pair<G4double,G4double> ElTot=GetElTot(pIU, pPDG, tgZ, tgN); // mean (IU)
  //if( ( (pPDG>999 && pIU<227.) || pIU<27.) && tgA>1) R=1.; // @@ TMP to accelerate @lowE
  if(pPDG>999 && pIU<227. && tgZ+tgN>1) R=1.;                // To accelerate @lowE
  else if(ElTot.second>0.)
  {
    R=ElTot.first/ElTot.second;             // El/Total ratio (does not depend on units
    QF2In=GetQF2IN_Ratio(ElTot.second/millibarn, tgZ+tgN);   // QuasiFree/Inelastic ratio
  }
#ifdef pdebug
  G4cout<<">>>OUT>>>G4QuasiFreeRatio::GetQF2IN_Ratio: QF2In="<<QF2In<<", R="<<R<<G4endl;
#endif
  return std::make_pair(QF2In,R);
}

// Calculatio QasiFree/Inelastic Ratio as a function of total hN cross-section (mb) and A
G4double G4QuasiFreeRatios::GetQF2IN_Ratio(G4double s_value, G4int A)
{
  static const G4int    nps=150;        // Number of steps in the R(s) LinTable
  static const G4int    mps=nps+1;      // Number of elements in the R(s) LinTable
  static const G4double sma=150.;       // The first LinTabEl(s=0)=1., s>sma -> LogTable
  static const G4double ds=sma/nps;     // Step of the linear Table
  static const G4int    nls=100;        // Number of steps in the R(lns) LogTable
  static const G4int    mls=nls+1;      // Number of elements in the R(lns) LogTable
  static const G4double lsi=5.;         // The min ln(s) LogTabEl(s=148.4 < sma=150.)
  static const G4double lsa=9.;         // The max ln(s) LogTabEl(s=148.4 - 8103. mb)
  static const G4double mi=std::exp(lsi);// The min s of LogTabEl(~ 148.4 mb)
  static const G4double max_s=std::exp(lsa);// The max s of LogTabEl(~ 8103. mb)
  static const G4double dl=(lsa-lsi)/nls;// A step Log-Length of the Logarithmic Table
  static const G4double edl=std::exp(dl);// Multiplication step of the Logarithmic Table
  static const G4double toler=.01;      // The tolarence mb defining the same cross-section
  static G4double lastS=0.;             // The last sigma value for which R was calculated
  static G4double lastR=0.;             // The last ratio R which was calculated
  // Local Associative Data Base:
  static std::vector<G4int>     vA;     // Vector of calculated A
  static std::vector<G4double>  vH;     // Vector of max s initialized in the LinTable
  static std::vector<G4int>     vN;     // Vector of topBin number initialized in LinTable
  static std::vector<G4double>  vM;     // Vector of rel max ln(s) initialized in LogTable
  static std::vector<G4int>     vK;     // Vector of topBin number initialized in LogTable
  // Last values of the Associative Data Base:
  static G4int     lastA=0;             // theLast of calculated A
  static G4double  lastH=0.;            // theLast of max s initialized in the LinTable
  static G4int     lastN=0;             // theLast of topBin number initialized in LinTable
  static G4double  lastM=0.;            // theLast of rel max ln(s) initialized in LogTable
  static G4int     lastK=0;             // theLast of topBin number initialized in LogTable
  static G4double* lastT=0;             // theLast of pointer to LinTable in the C++ heap
  static G4double* lastL=0;             // theLast of pointer to LogTable in the C++ heap
  // LogTable is created only if necessary. The ratio R(s>8100 mb) = 0 for any nuclei
#ifdef pdebug
  G4cout<<"+++G4QuasiFreeRatio::GetQF2IN_Ratio:A="<<A<<", s="<<s_value<<G4endl;
#endif
  if(s_value<toler || A<2) return 1.;
  if(s_value>max_s) return 0.;
  if(A>238)
  {
    G4cout<<"-Warning-G4QuasiFreeRatio::GetQF2IN_Ratio:A="<<A<<">238, return zero"<<G4endl;
    return 0.;
  }
  G4int nDB=vA.size();                  // A number of nuclei already initialized in AMDB
  //  if(nDB && lastA==A && std::fabs(s_value-lastS)<toler) return lastR;
  if(nDB && lastA==A && s_value==lastS) return lastR;  // VI do not use tolerance
  G4bool found=false;
  G4int i=-1;
  if(nDB) for (i=0; i<nDB; i++) if(A==vA[i]) // Sirch for this A in AMDB
  {
    found=true;                         // The A value is found
    break;
  }
#ifdef pdebug
  G4cout<<"+++G4QuasiFreeRatio::GetQF2IN_Ratio: nDB="<<nDB<<", found="<<found<<G4endl;
#endif
  if(!nDB || !found)                    // Create new line in the AMDB
  {
    lastA = A;
#ifdef pdebug
    G4cout<<"G4QuasiFreeRatios::GetQF2IN_Ratio: NewT, A="<<A<<", nDB="<<nDB<<G4endl;
#endif
    lastT = new G4double[mps];          // Create the Linear Table
    lastN = static_cast<int>(s_value/ds)+1;   // MaxBin to be initialized in Lin Table
    if(lastN>nps)                       // The Lin Table should be completely initialized
    {
      lastN=nps;
      lastH=sma;
    }
    else lastH = lastN*ds;              // Calculate max initialized s for Lin Table
    G4double sv=0;
    lastT[0]=1.;
    for(G4int j=1; j<=lastN; j++)       // Calculate Lin Table values
    {
      sv+=ds;
      lastT[j]=CalcQF2IN_Ratio(sv,A);
    }
    lastL=new G4double[mls];            // Create the Logarithmic Table
    if(s_value>sma)                     // Initialize the Log Table (at least a part of it)
    {
#ifdef pdebug
      G4cout<<"G4QuasiFreeRatios::GetQF2IN_Ratio: NewL, A="<<A<<", nDB="<<nDB<<G4endl;
#endif
      G4double ls=std::log(s_value);
      lastK = static_cast<int>((ls-lsi)/dl)+1; // MaxBin to be initialized in Log Table
      if(lastK>nls)                     // The Log Table should be completely initialized
      {
        lastK=nls;
        lastM=lsa-lsi;                  // lastM is a ln(s)-difference (not a s-difference)
      }
      else lastM = lastK*dl;            // Calculate max initialized ln(s)-lsi for LogTab
      sv=mi;                            // Min s (not ln(s)) for the Log Table
      for(G4int j=0; j<=lastK; j++)     // Calculate Log Table values
      {
        lastL[j]=CalcQF2IN_Ratio(sv,A);
        if(j!=lastK) sv*=edl;
      }
    }
    else                                // Log Tab is not initialized for the firs call
    {
      lastK = 0;
      lastM = 0.;
    }
    i++;                                // Make a new record to AMDB and position on it
    vA.push_back(lastA);
    vH.push_back(lastH);
    vN.push_back(lastN);
    vM.push_back(lastM);
    vK.push_back(lastK);
    vT.push_back(lastT);
    vL.push_back(lastL);
  }
  else                                  // The A value was found in AMDB
  {
    lastA=vA[i];
    lastH=vH[i];
    lastN=vN[i];
    lastM=vM[i];
    lastK=vK[i];
    lastT=vT[i];
    lastL=vL[i];
#ifdef pdebug
    G4cout<<"G4QuasiFreeRatios::GetQF2IN_Ratio: Found, s="<<s_value<<", lastM="<<lastM<<G4endl;
#endif
    if(s_value>lastH)                          // At least Lin Table must be updated
    {
      G4int nextN=lastN+1;               // The next bin to be initialized in Lin Table (p)
#ifdef pdebug
      G4cout<<"G4QuasiFreeRatios::GetQF2IN_Ratio: lastN="<<lastN<<" ?< nps="<<nps<<G4endl;
#endif
      if(lastN < nps)                    // The Lin Table is not completely initialized
      {
        lastN = static_cast<int>(s_value/ds)+1;// MaxBin to be initialized in Lin Table
        G4double sv=lastH;
        if(lastN>nps)
        {
          lastN=nps;
          lastH=sma;
        }
        else lastH = lastN*ds;           // Calculate max initialized s for LinTab
        for(G4int j=nextN; j<=lastN; j++)// Calculate LogTab values
        {
          sv+=ds;
          lastT[j]=CalcQF2IN_Ratio(sv,A);
        }
#ifdef pdebug
        G4cout<<"G4QuasiFreeRatios::GetQF2IN_Ratio: End of LinTab update"<<G4endl;
#endif
      } // End of LinTab update
#ifdef pdebug
      G4cout<<"G4QFRatios::GetQF2IN_Ratio: lN="<<lastN<<", nN="<<nextN<<", i="<<i<<G4endl;
#endif
      if(lastN >= nextN)                 // The Lin Table update really took place
      {
        vH[i]=lastH;
        vN[i]=lastN;
      }
      G4int nextK=lastK+1;               // Possible next element in LogTable to be updated
      if(!lastK) nextK=0;                // The Log Table has not been filled at all
#ifdef pdebug
      G4cout<<"G4QFRat::GetQF2IN_Ratio: sma="<<sma<<", lastK="<<lastK<<" < "<<nls<<G4endl;
#endif
      if(s_value>sma && lastK<nls)             // LogTab must be updated (not filled completely)
      {
        G4double sv=std::exp(lastM+lsi); // Define starting s-poit (lastM will be changed)
        G4double ls=std::log(s_value);
        lastK = static_cast<int>((ls-lsi)/dl)+1; // MaxBin to be initialized in Log Table
        if(lastK>nls)
        {
          lastK=nls;
          lastM=lsa-lsi;
        }
        else lastM = lastK*dl;           // Calculate max initialized ln(s)-lsi for LogTab
#ifdef pdebug
        G4cout<<"G4QFRat::GetQF2IN_Ratio: nK="<<nextK<<", lK="<<lastK<<", sv="<<sv<<G4endl;
#endif
        for(G4int j=nextK; j<=lastK; j++)// Calculate LogTab values
        {
          sv*=edl;
#ifdef pdebug
          G4cout<<"G4QFRat::GetQF2IN_Ratio: j="<<j<<", sv="<<sv<<", A="<<A<<G4endl;
#endif
          lastL[j]=CalcQF2IN_Ratio(sv,A);
        }
#ifdef pdebug
        G4cout<<"G4QuasiFreeRatios::GetQF2IN_Ratio: End of LinTab update"<<G4endl;
#endif
      } // End of LogTab update
#ifdef pdebug
      G4cout<<"G4QFRatios::GetQF2IN_Ratio: lK="<<lastK<<", nK="<<nextK<<", i="<<i<<G4endl;
#endif
      if(lastK >= nextK)                 // The Lin Table update really took place
      {
        vM[i]=lastM;
        vK[i]=lastK;
      }
    }
  }
#ifdef pdebug
  G4cout<<"G4QuasiFreeRatios::GetQF2IN_Ratio: BeforeTab s="<<s_value<<", sma="<<sma<<G4endl;
#endif
  // Now one can use tabeles to calculate the value
  if(s_value<sma)                             // Use linear table
  {
    G4int n=static_cast<int>(s_value/ds);     // Low edge number of the bin
    G4double d=s_value-n*ds;                  // Linear shift
    G4double v=lastT[n];                // Base
    lastR=v+d*(lastT[n+1]-v)/ds;        // Result
  }
  else                                  // Use log table
  {
    G4double ls=std::log(s_value)-lsi;        // ln(s)-l_min
    G4int n=static_cast<int>(ls/dl);    // Low edge number of the bin
    G4double d=ls-n*dl;                 // Log shift
    G4double v=lastL[n];                // Base
    lastR=v+d*(lastL[n+1]-v)/dl;        // Result
  }
  if(lastR<0.) lastR=0.;
  if(lastR>1.) lastR=1.;
#ifdef pdebug
  G4cout<<"G4QuasiFreeRatios::GetQF2IN_Ratio: BeforeRet lastR="<<lastR<<G4endl;
#endif
  return lastR;
} // End of CalcQF2IN_Ratio

// Calculatio QasiFree/Inelastic Ratio as a function of total hN cross-section and A
G4double G4QuasiFreeRatios::CalcQF2IN_Ratio(G4double s_value, G4int A)
{
  static const G4double C=1.246;
  G4double s2=s_value*s_value;
  G4double s4=s2*s2;
  G4double ss=std::sqrt(std::sqrt(s_value));
  G4double P=7.48e-5*s2/(1.+8.77e12/s4/s4/s2);
  G4double E=.2644+.016/(1.+std::exp((29.54-s_value)/2.49));
  G4double F=ss*.1526*std::exp(-s2*ss*.0000859);
  return C*std::exp(-E*std::pow(G4double(A-1.),F))/std::pow(G4double(A),P);
} // End of CalcQF2IN_Ratio

// For hadron PDG with momentum Mom (GeV/c) on N (p/n) calculate <sig_el,sig_tot> pair (mb)
std::pair<G4double,G4double> G4QuasiFreeRatios::GetElTotXS(G4double p, G4int PDG, G4bool F)
{
  G4int ind=0;                                 // Prototype of the reaction index
  G4bool kfl=true;                             // Flag of K0/aK0 oscillation
  G4bool kf=false;
  if(PDG==130||PDG==310)
  {
    kf=true;
    if(G4UniformRand()>.5) kfl=false;
  }
  if      ( (PDG == 2212 && F) || (PDG == 2112 && !F) ) ind=0; // pp/nn
  else if ( (PDG == 2112 && F) || (PDG == 2212 && !F) ) ind=1; // np/pn
  else if ( (PDG == -211 && F) || (PDG == 211 && !F) ) ind=2; // pimp/pipn
  else if ( (PDG == 211 && F) || (PDG == -211 && !F) ) ind=3; // pipp/pimn
  else if ( PDG == -321 || PDG == -311 || (kf && !kfl) ) ind=4; // KmN/K0N
  else if ( PDG == 321 || PDG == 311 || (kf && kfl) ) ind=5; // KpN/aK0N
  else if ( PDG >  3000 && PDG <  3335) ind=6; // @@ for all hyperons - take Lambda
  else if ( PDG > -3335 && PDG < -2000) ind=7; // @@ for all anti-baryons (anti-p/anti-n)
  else {
    G4cout<<"*Error*G4QuasiFreeRatios::CalcElTotXS: PDG="<<PDG
          <<", while it is defined only for p,n,hyperons,anti-baryons,pi,K/antiK"<<G4endl;
    G4Exception("G4QuasiFreeRatio::CalcElTotXS:","22",FatalException,"CHIPScrash");
  }
  return QFreeScat->CalcElTot(p,ind);
}

// (Mean Elastic and Mean Total) Cross-Sections (mb) for PDG+(Z,N) at P=p[GeV/c]
std::pair<G4double,G4double> G4QuasiFreeRatios::GetElTot(G4double pIU, G4int hPDG,
                                                         G4int Z,       G4int N)
{
  G4double pGeV=pIU/gigaelectronvolt;
#ifdef pdebug
  G4cout<<"-->G4QuasiFreeR::GetElTot: P="<<pIU<<",pPDG="<<hPDG<<",Z="<<Z<<",N="<<N<<G4endl;
#endif
  if(Z<1 && N<1)
  {
    G4cout<<"-Warning-G4QuasiFreeRatio::GetElTot:Z="<<Z<<",N="<<N<<", return zero"<<G4endl;
    return std::make_pair(0.,0.);
  }
  std::pair<G4double,G4double> hp=QFreeScat->FetchElTot(pGeV, hPDG, true);
  std::pair<G4double,G4double> hn=QFreeScat->FetchElTot(pGeV, hPDG, false);
#ifdef pdebug
  G4cout<<"-OUT->G4QFRat::GetElTot: hp("<<hp.first<<","<<hp.second<<"), hn("<<hn.first<<","
        <<hn.second<<")"<<G4endl;
#endif
  G4double A=(Z+N)/millibarn;                // To make the result in independent units(IU)
  return std::make_pair((Z*hp.first+N*hn.first)/A,(Z*hp.second+N*hn.second)/A);
} // End of GetElTot

// (Mean Elastic and Mean Total) Cross-Sections (mb) for PDG+(Z,N) at P=p[GeV/c]
std::pair<G4double,G4double> G4QuasiFreeRatios::GetChExFactor(G4double pIU, G4int hPDG,
                                                              G4int Z, G4int N)
{
  G4double pGeV=pIU/gigaelectronvolt;
  G4double resP=0.;
  G4double resN=0.;
  if(Z<1 && N<1)
  {
    G4cout<<"-Warning-G4QuasiFreeRatio::GetChExF:Z="<<Z<<",N="<<N<<", return zero"<<G4endl;
    return std::make_pair(resP,resN);
  }
  G4double A=Z+N;
  G4double pf=0.;                              // Possibility to interact with a proton
  G4double nf=0.;                              // Possibility to interact with a neutron
  if   (hPDG==-211||hPDG==-321||hPDG==3112||hPDG==3212||hPDG==3312||hPDG==3334) pf=Z/(A+N);
  else if(hPDG==211||hPDG==321||hPDG==3222||hPDG==3212||hPDG==3322) nf=N/(A+Z);
  else if(hPDG==-311||hPDG==311||hPDG==130||hPDG==310)
  {
    G4double dA=A+A;
    pf=Z/(dA+N+N);
    nf=N/(dA+Z+Z);
  }
  G4double mult=1.;  // Factor of increasing multiplicity ( ? @@)
  if(pGeV>.5)
  {
    mult=1./(1.+std::log(pGeV+pGeV))/pGeV;
    if(mult>1.) mult=1.;
  }
  if(pf)
  {
    std::pair<G4double,G4double> hp=QFreeScat->FetchElTot(pGeV, hPDG, true);
    resP=pf*(hp.second/hp.first-1.)*mult;
  }
  if(nf)
  {
    std::pair<G4double,G4double> hn=QFreeScat->FetchElTot(pGeV, hPDG, false);
    resN=nf*(hn.second/hn.first-1.)*mult;
  }
  return std::make_pair(resP,resN);
} // End of GetChExFactor

// scatter (pPDG,p4M) on a virtual nucleon (NPDG,N4M), result: final pair(newN4M,newp4M)
// if(newN4M.e()==0.) - below threshold, XS=0, no scattering of the progectile happened
std::pair<G4LorentzVector,G4LorentzVector> G4QuasiFreeRatios::Scatter(G4int NPDG,
                                     G4LorentzVector N4M, G4int pPDG, G4LorentzVector p4M)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);// Mass of deuteron
  static const G4double mTrit= G4QPDGCode(2112).GetNuclMass(1,2,0);// Mass of tritium
  static const G4double mHel3= G4QPDGCode(2112).GetNuclMass(2,1,0);// Mass of Helium3
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);// Mass of alpha
  static const G4double two3d= 2./3.;
  static const G4double two3d2= std::pow(2.,two3d); // The slope reductions for fragments
  static const G4double two3d3= std::pow(3.,two3d);
  static const G4double two3d4= std::pow(4.,two3d);
  G4LorentzVector pr4M=p4M/megaelectronvolt;   // Convert 4-momenta in MeV (keep p4M)
  N4M/=megaelectronvolt;
  G4LorentzVector tot4M=N4M+p4M;
#ifdef ppdebug
  G4cerr<<"->G4QFR::Scat:p4M="<<pr4M<<",N4M="<<N4M<<",t4M="<<tot4M<<",NPDG="<<NPDG<<G4endl;
#endif
  G4double mT=mNeut;
  G4int Z=0;
  G4int N=1;
  if(NPDG==2212||NPDG==90001000)
  {
    mT=mProt;
    Z=1;
    N=0;
  }
  else if(NPDG==90001001)
  {
    mT=mDeut;
    Z=1;
    N=1;
  }
  else if(NPDG==90002001)
  {
    mT=mHel3;
    Z=2;
    N=1;
  }
  else if(NPDG==90001002)
  {
    mT=mTrit;
    Z=1;
    N=2;
  }
  else if(NPDG==90002002)
  {
    mT=mAlph;
    Z=2;
    N=2;
  }
  else if(NPDG!=2112&&NPDG!=90000001)
  {
    G4cout<<"Error:G4QuasiFreeRatios::Scatter:NPDG="<<NPDG<<" is not 2212 or 2112"<<G4endl;
    G4Exception("G4QuasiFreeRatios::Scatter:","21",FatalException,"CHIPScomplain");
    //return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M);// Use this if not exception
  }
  G4double mT2=mT*mT;
  G4double mP2=pr4M.m2();
  G4double E=(tot4M.m2()-mT2-mP2)/(mT+mT);
#ifdef pdebug
  G4cerr<<"G4QFR::Scat:qM="<<mT<<",qM2="<<mT2<<",pM2="<<mP2<<",totM2="<<tot4M.m2()<<G4endl;
#endif
  G4double E2=E*E;
  if(E<0. || E2<mP2)
  {
#ifdef ppdebug
    G4cerr<<"-Warning-G4QFR::Scat:*Negative Energy*E="<<E<<",E2="<<E2<<"<M2="<<mP2<<G4endl;
#endif
    return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
  }
  G4double P=std::sqrt(E2-mP2);                   // Momentum in pseudo laboratory system
  G4VQCrossSection* PCSmanager=G4QProtonElasticCrossSection::GetPointer();
  G4VQCrossSection* NCSmanager=G4QNeutronElasticCrossSection::GetPointer();
#ifdef ppdebug
  G4cout<<"G4QFR::Scatter: Before XS, P="<<P<<", Z="<<Z<<", N="<<N<<", PDG="<<pPDG<<G4endl;
#endif
  // @@ Temporary NN t-dependence for all hadrons
  if(pPDG>3400 || pPDG<-3400) G4cout<<"-Warning-G4QuasiFree::Scatter: pPDG="<<pPDG<<G4endl;
  G4int PDG=2212;                                                // *TMP* instead of pPDG
  if(pPDG==2112||pPDG==-211||pPDG==-321) PDG=2112;               // *TMP* instead of pPDG
  if(!Z && N==1)                 // Change for Quasi-Elastic on neutron
  {
    Z=1;
    N=0;
    if     (PDG==2212) PDG=2112;
    else if(PDG==2112) PDG=2212;
  }
  G4double xSec=0.;                        // Prototype of Recalculated Cross Section *TMP*
  if(PDG==2212) xSec=PCSmanager->GetCrossSection(false, P, Z, N, PDG); // P CrossSect *TMP*
  else          xSec=NCSmanager->GetCrossSection(false, P, Z, N, PDG); // N CrossSect *TMP*
#ifdef ppdebug
  G4cout<<"G4QuasiFreeRat::Scatter: pPDG="<<pPDG<<",P="<<P<<",CS="<<xSec/millibarn<<G4endl;
#endif
#ifdef nandebug
  if(xSec>0. || xSec<0. || xSec==0);
  else  G4cout<<"*Warning*G4QuasiFreeRatios::Scatter: xSec="<<xSec/millibarn<<G4endl;
#endif
  // @@ check a possibility to separate p, n, or alpha (!)
  if(xSec <= 0.)                                    // The cross-section iz 0 -> Do Nothing
  {
#ifdef ppdebug
    G4cerr<<"-Warning-G4QFR::Scat:**Zero XS**PDG="<<pPDG<<",NPDG="<<NPDG<<",P="<<P<<G4endl;
#endif
    return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); //Do Nothing Action
  }
  G4double mint=0.;                        // Prototype of functional rand -t (MeV^2) *TMP*
  if(PDG==2212) mint=PCSmanager->GetExchangeT(Z,N,PDG);// P functional rand -t(MeV^2) *TMP*
  else          mint=NCSmanager->GetExchangeT(Z,N,PDG);// N functional rand -t(MeV^2) *TMP*
  if     (mT==mDeut)              mint/=two3d2;
  else if(mT==mTrit || mT==mHel3) mint/=two3d3; 
  else if(mT==mAlph)              mint/=two3d4; 
  G4double maxt=0.;                                    // Prototype of max possible -t
  if(PDG==2212) maxt=PCSmanager->GetHMaxT();           // max possible -t
  else          maxt=NCSmanager->GetHMaxT();           // max possible -t
#ifdef ppdebug
  G4cout<<"G4QFR::Scat:PDG="<<PDG<<",P="<<P<<",X="<<xSec<<",-t="<<mint<<"<"<<maxt<<", Z="
        <<Z<<",N="<<N<<G4endl;
#endif
#ifdef nandebug
  if(mint>-.0000001);
  else  G4cout<<"*Warning*G4QFR::Scat: -t="<<mint<<G4endl;
#endif
  G4double cost=1.-(mint+mint)/maxt; // cos(theta) in CMS
#ifdef ppdebug
  G4cout<<"G4QFR::Scat:-t="<<mint<<"<"<<maxt<<", cost="<<cost<<", Z="<<Z<<",N="<<N<<G4endl;
#endif
  if(cost>1. || cost<-1. || !(cost>-1. || cost<=1.))
  {
    if     (cost>1.)  cost=1.;
    else if(cost<-1.) cost=-1.;
    else
    {
      G4double tm=0.;
      if(PDG==2212) tm=PCSmanager->GetHMaxT();
      else          tm=NCSmanager->GetHMaxT();
      G4cerr<<"G4QuasiFreeRatio::Scat:*NAN* cost="<<cost<<",-t="<<mint<<",tm="<<tm<<G4endl;
      return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
    }
  }
  G4LorentzVector reco4M=G4LorentzVector(0.,0.,0.,mT);      // 4mom of the recoil nucleon
  G4LorentzVector dir4M=tot4M-G4LorentzVector(0.,0.,0.,(tot4M.e()-mT)*.01);
  if(!G4QHadron(tot4M).RelDecayIn2(pr4M, reco4M, dir4M, cost, cost))
  {
    G4cerr<<"G4QFR::Scat:t="<<tot4M<<tot4M.m()<<",mT="<<mT<<",mP="<<std::sqrt(mP2)<<G4endl;
    //G4Exception("G4QFR::Scat:","009",FatalException,"Decay of ElasticComp");
    return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
  }
#ifdef ppdebug
  G4cout<<"G4QFR::Scat:p4M="<<pr4M<<"+r4M="<<reco4M<<",dr="<<dir4M<<",t4M="<<tot4M<<G4endl;
#endif
  return std::make_pair(reco4M*megaelectronvolt,pr4M*megaelectronvolt); // Result
} // End of Scatter

// scatter (pPDG,p4M) on a virtual nucleon (NPDG,N4M), result: final pair(newN4M,newp4M)
// if(newN4M.e()==0.) - below threshold, XS=0, no scattering of the progectile happened
// User should himself change the charge (PDG) (e.g. pn->np, pi+n->pi0p, pi-p->pi0n etc.)
std::pair<G4LorentzVector,G4LorentzVector> G4QuasiFreeRatios::ChExer(G4int NPDG,
                                     G4LorentzVector N4M, G4int pPDG, G4LorentzVector p4M)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  G4LorentzVector pr4M=p4M/megaelectronvolt;          // Convert 4-momenta in MeV(keep p4M)
  N4M/=megaelectronvolt;
  G4LorentzVector tot4M=N4M+p4M;
  G4int Z=0;
  G4int N=1;
  G4int sPDG=0;                                        // PDG code of the scattered hadron
  G4double mS=0.;                                      // proto of mass of scattered hadron
  G4double mT=mProt;                                   // mass of the recoil nucleon
  G4cout<<"-Warning-G4QFR::ChEx: Impossible for Omega-"<<G4endl; // Omega-
  if(NPDG==2212)
  {
    mT=mNeut;
    Z=1;
    N=0;
    if(pPDG==-211) sPDG=111;                           // pi+    -> pi0
    else if(pPDG==-321)
    {
      sPDG=310;                                        // K+     -> K0S
      if(G4UniformRand()>.5) sPDG=130;                 // K+     -> K0L
    }
    else if(pPDG==-311||pPDG==311||pPDG==130||pPDG==310) sPDG=321;  // K0     -> K+ (?)
    else if(pPDG==3112) sPDG=3212;                     // Sigma- -> Sigma0
    else if(pPDG==3212) sPDG=3222;                     // Sigma0 -> Sigma+
    else if(pPDG==3312) sPDG=3322;                     // Xi-    -> Xi0
  }
  else if(NPDG==2112) // Default
  {
    if(pPDG==211)  sPDG=111;                           // pi+    -> pi0
    else if(pPDG==321)
    {
      sPDG=310;                                        // K+     -> K0S
      if(G4UniformRand()>.5) sPDG=130;                 // K+     -> K0L
    }
    else if(pPDG==-311||pPDG==311||pPDG==130||pPDG==310) sPDG=-321; // K0     -> K- (?)
    else if(pPDG==3222) sPDG=3212;                     // Sigma+ -> Sigma0
    else if(pPDG==3212) sPDG=3112;                     // Sigma0 -> Sigma-
    else if(pPDG==3322) sPDG=3312;                     // Xi0    -> Xi-
  }
  else
  {
    G4cout<<"Error:G4QuasiFreeRatios::ChExer: NPDG="<<NPDG<<" is not 2212 or 2112"<<G4endl;
    G4Exception("G4QuasiFreeRatios::ChExer:","21",FatalException,"CHIPS complain");
    //return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M);// Use this if not exception
  }
  if(sPDG) mS=mNeut;
  else
  {
    G4cout<<"Error:G4QuasiFreeRatios::ChExer: BAD pPDG="<<pPDG<<", NPDG="<<NPDG<<G4endl;
    G4Exception("G4QuasiFreeRatios::ChExer:","21",FatalException,"CHIPS complain");
    //return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M);// Use this if not exception
  }
  G4double mT2=mT*mT;
  G4double mS2=mS*mS;
  G4double E=(tot4M.m2()-mT2-mS2)/(mT+mT);
  G4double E2=E*E;
  if(E<0. || E2<mS2)
  {
#ifdef pdebug
    G4cerr<<"-Warning-G4QFR::ChEx:*Negative Energy*E="<<E<<",E2="<<E2<<"<M2="<<mS2<<G4endl;
#endif
    return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
  }
  G4double P=std::sqrt(E2-mS2);                   // Momentum in pseudo laboratory system
  G4VQCrossSection* PCSmanager=G4QProtonElasticCrossSection::GetPointer();
  G4VQCrossSection* NCSmanager=G4QNeutronElasticCrossSection::GetPointer();
#ifdef debug
  G4cout<<"G4QFR::ChExer: Before XS, P="<<P<<", Z="<<Z<<", N="<<N<<", PDG="<<pPDG<<G4endl;
#endif
  // @@ Temporary NN t-dependence for all hadrons
  G4int PDG=2212;                                                // *TMP* instead of pPDG
  if(pPDG==2112||pPDG==-211||pPDG==-321) PDG=2112;               // *TMP* instead of pPDG
  if(!Z && N==1)                 // Change for Quasi-Elastic on neutron
  {
    Z=1;
    N=0;
    if     (PDG==2212) PDG=2112;
    else if(PDG==2112) PDG=2212;
  }
  G4double xSec=0.;                        // Prototype of Recalculated Cross Section *TMP*
  if(PDG==2212) xSec=PCSmanager->GetCrossSection(false, P, Z, N, PDG); // P CrossSect *TMP*
  else          xSec=NCSmanager->GetCrossSection(false, P, Z, N, PDG); // N CrossSect *TMP*
#ifdef debug
  G4cout<<"G4QuasiFreeRat::ChExer: pPDG="<<pPDG<<",P="<<P<<", CS="<<xSec/millibarn<<G4endl;
#endif
#ifdef nandebug
  if(xSec>0. || xSec<0. || xSec==0);
  else  G4cout<<"*Warning*G4QuasiFreeRatios::ChExer: xSec="<<xSec/millibarn<<G4endl;
#endif
  // @@ check a possibility to separate p, n, or alpha (!)
  if(xSec <= 0.) // The cross-section iz 0 -> Do Nothing
  {
#ifdef pdebug
    G4cerr<<"-Warning-G4QFR::ChEx:**Zero XS**PDG="<<pPDG<<",NPDG="<<NPDG<<",P="<<P<<G4endl;
#endif
    return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); //Do Nothing Action
  }
  G4double mint=0.;                        // Prototype of functional rand -t (MeV^2) *TMP*
  if(PDG==2212) mint=PCSmanager->GetExchangeT(Z,N,PDG);// P functional rand -t(MeV^2) *TMP*
  else          mint=NCSmanager->GetExchangeT(Z,N,PDG);// N functional rand -t(MeV^2) *TMP*
#ifdef pdebug
  G4cout<<"G4QFR::ChEx:PDG="<<pPDG<<", P="<<P<<", CS="<<xSec<<", -t="<<mint<<G4endl;
#endif
#ifdef nandebug
  if(mint>-.0000001);
  else  G4cout<<"*Warning*G4QFR::ChExer: -t="<<mint<<G4endl;
#endif
  G4double maxt=0.;                                    // Prototype of max possible -t
  if(PDG==2212) maxt=PCSmanager->GetHMaxT();           // max possible -t
  else          maxt=NCSmanager->GetHMaxT();           // max possible -t
  G4double cost=1.-mint/maxt;                          // cos(theta) in CMS
#ifdef pdebug
  G4cout<<"G4QuasiFfreeRatio::ChExer: -t="<<mint<<", maxt="<<maxt<<", cost="<<cost<<G4endl;
#endif
  if(cost>1. || cost<-1. || !(cost>-1. || cost<=1.))
  {
    if     (cost>1.)  cost=1.;
    else if(cost<-1.) cost=-1.;
    else
    {
      G4cerr<<"G4QuasiFreeRatio::ChExer:*NAN* c="<<cost<<",t="<<mint<<",tm="<<maxt<<G4endl;
      return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
    }
  }
  G4LorentzVector reco4M=G4LorentzVector(0.,0.,0.,mT);      // 4mom of the recoil nucleon
  pr4M=G4LorentzVector(0.,0.,0.,mS);                        // 4mom of the scattered hadron
  G4LorentzVector dir4M=tot4M-G4LorentzVector(0.,0.,0.,(tot4M.e()-mT)*.01);
  if(!G4QHadron(tot4M).RelDecayIn2(pr4M, reco4M, dir4M, cost, cost))
  {
    G4cerr<<"G4QFR::ChEx:t="<<tot4M<<tot4M.m()<<",mT="<<mT<<",mP="<<mS<<G4endl;
    //G4Exception("G4QFR::ChExer:","009",FatalException,"Decay of ElasticComp");
    return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
  }
#ifdef debug
  G4cout<<"G4QFR::ChEx:p4M="<<p4M<<"+r4M="<<reco4M<<"="<<p4M+reco4M<<"="<<tot4M<<G4endl;
#endif
  return std::make_pair(reco4M*megaelectronvolt,pr4M*megaelectronvolt); // Result
} // End of ChExer

// Calculate ChEx/El ratio (p is in independent units, (Z,N) is target, pPDG is projectile)
G4double G4QuasiFreeRatios::ChExElCoef(G4double p, G4int Z, G4int N, G4int pPDG) 
{
  p/=MeV;                                // Converted from independent units
  G4double A=Z+N;
  if(A<1.5) return 0.;
  G4double C=0.;
  if     (pPDG==2212) C=N/(A+Z);
  else if(pPDG==2112) C=Z/(A+N);
  else G4cout<<"*Warning*G4QCohChrgExchange::ChExElCoef: wrong PDG="<<pPDG<<G4endl;
  C*=C;                         // Coherent processes squares the amplitude
  // @@ This is true only for nucleons: other projectiles must be treated differently
  G4double sp=std::sqrt(p);
  G4double p2=p*p;            
  G4double p4=p2*p2;
  G4double dl1=std::log(p)-5.;
  G4double T=(6.75+.14*dl1*dl1+13./p)/(1.+.14/p4)+.6/(p4+.00013);
  G4double U=(6.25+8.33e-5/p4/p)*(p*sp+.34)/p2/p; 
  G4double R=U/T;
  return C*R*R;
}
