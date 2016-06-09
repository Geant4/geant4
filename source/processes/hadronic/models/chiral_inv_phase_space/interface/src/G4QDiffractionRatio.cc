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
// $Id: G4QDiffractionRatio.cc,v 1.3 2007/06/01 07:22:18 mkossov Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//
// G4 Physics class: G4QDiffractionRatio for N+A elastic cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Oct-06
// 
//================================================================================

//#define debug
//#define pdebug
//#define nandebug

#include "G4QDiffractionRatio.hh"

// Returns Pointer to the G4VQCrossSection class
G4QDiffractionRatio* G4QDiffractionRatio::GetPointer()
{
  static G4QDiffractionRatio theRatios;   // *** Static body of the Diffraction Ratio ***
  return &theRatios;
}

// Calculation of pair(QuasiFree/Inelastic,QuasiElastic/QuasiFree)
G4double G4QDiffractionRatio::GetRatio(G4double pIU, G4int pPDG, G4int tgZ, G4int tgN)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mN=.5*(mNeut+mProt);
  static const G4double dmN=mN+mN;
  static const G4double mN2=mN*mN;
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
  static const G4double ms=std::exp(lsa);// The max s of logTabEl(~ 2981 GeV)
  static const G4double dl=(lsa-lsi)/nls;// Step of the logarithmic Table
  static const G4double edl=std::exp(dl);// Multiplication step of the logarithmic Table
  static const G4double toler=.0001;    // Tolarence (GeV) defining the same sqs
  static G4double lastS=0.;             // Last sqs value for which R was calculated
  static G4double lastR=0.;             // Last ratio R which was calculated
  // Local Associative Data Base:
  static std::vector<G4int>     vA;     // Vector of calculated A
  static std::vector<G4double>  vH;     // Vector of max sqs initialized in the LinTable
  static std::vector<G4int>     vN;     // Vector of topBin number initialized in LinTable
  static std::vector<G4double>  vM;     // Vector of relMax ln(sqs) initialized in LogTable
  static std::vector<G4int>     vK;     // Vector of topBin number initialized in LogTable
  static std::vector<G4double*> vT;     // Vector of pointers to LinTable in C++ heap
  static std::vector<G4double*> vL;     // Vector of pointers to LogTable in C++ heap
  // Last values of the Associative Data Base:
  static G4int     lastPDG=0;           // Last PDG for which R was calculated (now indep)
  static G4int     lastA=0;             // theLast of calculated A
  static G4double  lastH=0.;            // theLast of max sqs initialized in the LinTable
  static G4int     lastN=0;             // theLast of topBin number initialized in LinTable
  static G4double  lastM=0.;            // theLast of relMax ln(sqs) initialized in LogTab.
  static G4int     lastK=0;             // theLast of topBin number initialized in LogTable
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
  lastPDG=pPDG;                         // @@ at present ratio is PDG independent @@
  // Calculate sqs
  G4double pM=G4QPDGCode(pPDG).GetMass()*.001; // Projectile mass in GeV
  G4double pM2=pM*pM;
  G4double mom=pIU/gigaelectronvolt;    // Projectile momentum in GeV
  G4double s=std::sqrt(mN2+pM2+dmN*std::sqrt(pM2+mom*mom));
  G4int nDB=vA.size();                  // A number of nuclei already initialized in AMDB
  if(nDB && lastA==A && std::fabs(s-lastS)<toler) return lastR;
  if(s>ms)
  {
    lastR=CalcDiff2Prod_Ratio(s,A);     // @@ Probably user ought to be notified about bigS
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
    lastN = static_cast<int>(s/ds)+1;   // MaxBin to be initialized
    if(lastN>nps)
    {
      lastN=nps;
      lastH=sma;
    }
    else lastH = lastN*ds;              // Calculate max initialized s for LinTab
    G4double sv=0;
    lastT[0]=1.;
    for(G4int j=1; j<=lastN; j++)       // Calculate LinTab values
    {
      sv+=ds;
      lastT[j]=CalcDiff2Prod_Ratio(sv,A);
    }
    if(s>sma)                           // Initialize the logarithmic Table
    {
      lastL=new G4double[mls];          // Create the logarithmic Table
      G4double ls=std::log(s);
      lastK = static_cast<int>((ls-lsi)/dl)+1; // MaxBin to be initialized in LogTaB
      if(lastK>nls)
      {
        lastK=nls;
        lastM=lsa-lsi;
      }
      else lastM = lastK*dl;            // Calculate max initialized ln(s)-lsi for LogTab
      sv=mi;
      for(G4int j=0; j<=lastK; j++)     // Calculate LogTab values
      {
        lastL[j]=CalcDiff2Prod_Ratio(sv,A);
	       if(j!=lastK) sv*=edl;
      }
    }
    else                                // LogTab is not initialized
    {
      lastL = 0;
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
    if(s>lastM)                          // At least LinTab must be updated
    {
      G4int nextN=lastN+1;               // The next bin to be initialized
      if(lastN<nps)
      {
        lastN = static_cast<int>(s/ds)+1;// MaxBin to be initialized
        if(lastN>nps)
        {
          lastN=nps;
          lastH=sma;
        }
        else lastH = lastN*ds;           // Calculate max initialized s for LinTab
        G4double sv=lastM;
        for(G4int j=nextN; j<=lastN; j++)// Calculate LogTab values
        {
          sv+=ds;
          lastT[j]=CalcDiff2Prod_Ratio(sv,A);
        }
      } // End of LinTab update
      if(lastN>=nextN)
						{
        vH[i]=lastH;
        vN[i]=lastN;
      }
      G4int nextK=lastK+1;
      if(s>sma && lastK<nls)             // LogTab must be updated
						{
        G4double sv=std::exp(lastM+lsi); // Define starting poit (lastM will be changed)
        G4double ls=std::log(s);
        lastK = static_cast<int>((ls-lsi)/dl)+1; // MaxBin to be initialized in LogTaB
        if(lastK>nls)
        {
          lastK=nls;
          lastM=lsa-lsi;
        }
        else lastM = lastK*dl;           // Calculate max initialized ln(s)-lsi for LogTab
        for(G4int j=nextK; j<=lastK; j++)// Calculate LogTab values
        {
	         sv*=edl;
          lastL[j]=CalcDiff2Prod_Ratio(sv,A);
        }
      } // End of LogTab update
      if(lastK>=nextK)
						{
        vM[i]=lastM;
        vK[i]=lastK;
      }
    }
  }
  // Now one can use tabeles to calculate the value
  if(s<sma)                             // Use linear table
		{
    G4int n=static_cast<int>(s/ds);     // Low edge number of the bin
    G4double d=s-n*ds;                  // Linear shift
    G4double v=lastT[n];                // Base
    lastR=v+d*(lastT[n+1]-v)/ds;        // Result
  }
		else                                  // Use log table
		{
    G4double ls=std::log(s)-lsi;        // ln(s)-l_min
    G4int n=static_cast<int>(ls/dl);    // Low edge number of the bin
    G4double d=ls-n*dl;                 // Log shift
    G4double v=lastL[n];                // Base
    lastR=v+d*(lastL[n+1]-v)/dl;        // Result
  }
  if(lastR<0.) lastR=0.;
  if(lastR>1.) lastR=1.;
  return lastR;
} // End of CalcDiff2Prod_Ratio

// Calculate Diffraction/Production Ratio as a function of total sq(s)(hN) (in GeV), A=Z+N
G4double G4QDiffractionRatio::CalcDiff2Prod_Ratio(G4double s, G4int A)
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
  if(s<=0. || A<=1) return 0.;
  if(A!=mA && A!=1)
		{
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
    G4double p=std::pow(a,0.37);
    p1=(.023*p+3.5/a3+2.1e6/a12+4.e-14*a5)/(1.+7.6e-4*a*sa+2.15e7/a11);
    p2=(1.42*std::pow(a,0.61)+1.6e5/a8+4.5e-8*a4)/(1.+4.e-8*a4+1.2e4/a6);
    G4double q=std::pow(a,0.7);
    p4=(.036/q+.0009*q)/(1.+6./a3+1.e-7*a3);
				p5=1.3*std::pow(a,0.1168)/(1.+1.2e-8*a3);
    p6=.00046*(a+11830./a2);
    p7=1./(1.+6.17/a2+.00406*a);
  }
  else if(A==1 && mA!=1)
  {
    p1=.0315;
    p2=.73417;
    p4=.01109;
    p5=1.0972;
    p6=.065787;
    p7=.62976;
  }
  else if(std::fabs(s-S)/S<.0001) return R;
  G4double s2=s*s;
  G4double s4=s2*s2;
		G4double dl=std::log(s)-p5;
  R=1./(1.+1./(p1+p2/s4+p4*dl*dl/(1.+p6*std::pow(s,p7))));
	 return R;
} // End of CalcQF2IN_Ratio

G4QHadronVector* G4QDiffractionRatio::Fragment(G4int pPDG, G4LorentzVector p4M,
                                               G4int tgZ, G4int tgN)
{
  static const G4double pFm= 250.; // Fermi momentum in MeV (delta function)
  static const G4double pFm2= pFm*pFm; // Squared Fermi momentum in MeV^2 (delta function)
  static const G4double mPi0= G4QPDGCode(111).GetMass(); // pi0 mass (MeV =min diffraction)
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mNeut2=mNeut*mNeut;
  static const G4double dmNeut=mNeut+mNeut;
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mProt2=mProt*mProt;
  static const G4double dmProt=mProt+mProt;
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
  G4int Z=0;
  G4int N=1;
  if(rPDG==2212)
  {
    mT=mProt;
    mT2=mProt2;
    dmT=dmProt;
    Z=1;
    N=0;
  }
  G4double mP2=pr4M.m2();               // Squared mass of the projectile
  if(mP2<0.) mP2=0.;                    // A possible problem for photon (m_min = 2*m_pi0)
  G4double s=tot4M.m2();                 // @@ Check <0 ...
  G4double E=(s-mT2-mP2)/dmT;  // Effective interactin energy (virt. nucl. target)
  G4double E2=E*E;
  if(E<0. || E2<mP2)
  {
#ifdef pdebug
    G4cerr<<"Warning-G4DifR::Frag:*Negative Energy*E="<<E<<",E2="<<E2<<"<M2="<<mP2<<G4endl;
#endif
    return ResHV; // Do Nothing Action
  }
  G4double mP=std::sqrt(mP2);
  if(mP<.1)mP=mPi0;                      // For photons min diffraction is gamma+P->Pi0+Pi0
  G4double mMin=mP+mPi0;                 // Minimum diffractive mass
  G4double ss=std::sqrt(s);              // CM compound mass (sqrt(s))
  G4double mMax=ss-mT;                   // Maximum diffraction mass
  if(mMin>=mMax)
  {
#ifdef pdebug
    G4cerr<<"Warning-G4DifR::Frag:ZeroDiffractionMRange, mi="<<mMin<<", ma="<<mMax<<G4endl;
#endif
    return ResHV; // Do Nothing Action
  }
  G4double mDif=mMin/(1.-G4UniformRand()*(1.-mMin/mMax)); // Low Mass Approximation
  G4double mDif2=mDif*mDif;
  G4double ds=s-mT2-mDif2;
  G4double e=ds/dmT;
		G4double P=std::sqrt(e*e-mDif2);          // Momentum in pseudo laboratory system
  G4VQCrossSection* CSmanager=G4QElasticCrossSection::GetPointer();
#ifdef debug
  G4cout<<"G4QDiffR::Frag: Before XS, P="<<P<<", Z="<<Z<<", N="<<N<<", PDG="<<pPDG<<G4endl;
#endif
  // @@ Temporary NN t-dependence for all hadrons
  if(pPDG>3400 || pPDG<-3400) G4cout<<"-Warning-G4QDifR::Fragment: pPDG="<<pPDG<<G4endl;
  G4int PDG=2212;                                                  // *TMP* instead of pPDG
  G4double xSec=CSmanager->GetCrossSection(false, P, tgZ, tgN, PDG);// Rec.CrossSect *TMP*
  //G4double xSec=CSmanager->GetCrossSection(false, P, tgZ, tgN, pPDG); // Rec.CrossSect
#ifdef debug
  G4cout<<"G4QElast::Scatter:pPDG="<<pPDG<<",P="<<P<<",CS="<<xSec/millibarn<<G4endl;
#endif
#ifdef nandebug
  if(xSec>0. || xSec<0. || xSec==0);
  else  G4cout<<"***NAN***G4QDiffR::Fragment: xSec="<<xSec/millibarn<<G4endl;
#endif
  // @@ check a possibility to separate p, n, or alpha (!)
  if(xSec <= 0.) // The cross-section iz 0 -> Do Nothing
  {
#ifdef pdebug
    G4cerr<<"-Warning-G4QDiffR::Fragment:**Zero XS**PDG="<<pPDG<<",P="<<P<<G4endl;
#endif
    return ResHV; //Do Nothing Action
  }
  G4double t=CSmanager->GetExchangeT(tgZ,tgN,pPDG); // functional randomized -t (MeV^2)
  G4double maxt=(ds*ds-4*mT2*mDif2)/s;                 // maximum possible -t
#ifdef pdebug
  G4cout<<"G4QDifR::Frag:ph="<<pPDG<<",P="<<P<<",X="<<xSec<<",t="<<mint<<"<"<<maxt<<G4endl;
#endif
#ifdef nandebug
  if(mint>-.0000001);                          // To make the Warning for NAN
  else  G4cout<<"******G4QDiffR::Frag:-t="<<mint<<G4endl;
#endif
  G4double rt=t/maxt;
  G4double cost=1.-rt-rt;                          // cos(theta) in CMS
#ifdef ppdebug
  G4cout<<"G4QDiffR::Fragment: -t="<<t<<", maxt="<<maxt<<", cost="<<cost<<G4endl;
#endif
  if(cost>1. || cost<-1. || !(cost>-1. || cost<=1.))
  {
    if     (cost>1.)  cost=1.;
    else if(cost<-1.) cost=-1.;
    else
				{
      G4cerr<<"G4QDiffR::Fragm: *NAN* cost="<<cost<<",t="<<t<<",tmax="<<maxt<<G4endl;
      return ResHV; // Do Nothing Action
    }
  }
  G4LorentzVector r4M=G4LorentzVector(0.,0.,0.,mT);      // 4mom of the recoil nucleon
  G4LorentzVector d4M=G4LorentzVector(0.,0.,0.,mDif);    // 4mom of the diffract. Quasmon
  G4LorentzVector dir4M=tot4M-G4LorentzVector(0.,0.,0.,(tot4M.e()-mT)*.01);
  if(!G4QHadron(tot4M).RelDecayIn2(pr4M, r4M, dir4M, cost, cost))
  {
    G4cerr<<"G4QDifR::Fra:M="<<tot4M.m()<<",T="<<mT<<",D="<<mDif<<",T+D="<<mT+mDif<<G4endl;
    //G4Exception("G4QDifR::Fragm:","009",FatalException,"Decay of ElasticComp");
    return ResHV; // Do Nothing Action
  }
#ifdef debug
		G4cout<<"G4QFR::Scat:p4M="<<p4M<<"+r4M="<<reco4M<<"="<<scat4M+reco4M<<"="<<tot4M<<G4endl;
#endif
  // Now everything is ready for fragmentation and DoNothing projHadron must be wiped out
  ResHV->pop_back(); // Clean up pointer to the fake (doNothing) projectile
  delete hadron;     // Deklete the fake (doNothing) projectile hadron
  hadron = new G4QHadron(tPDG,t4M);  // Hadron for the recoil neucleus
  ResHV->push_back(hadron);          // Fill the recoil nucleus
  hadron = new G4QHadron(rPDG,r4M);  // Hadron for the recoil neucleus
  // Now the (pPdg,d4M) Quasmon must be fragmented
  G4QHadronVector* leadhs=new G4QHadronVector;// Prototype of QuasmOutput G4QHadronVector
  G4QContent dQC=G4QPDGCode(pPDG).GetQuarkContent(); // Quark Content of the projectile
  G4Quasmon* pan= new G4Quasmon(dQC,d4M); // --->---->---->----->-----> DELETED -->---*
  try                                                           //                    |
	 {                                                             //                    |
    G4QNucleus vac(90000000);                                   //                    |
    leadhs=pan->Fragment(vac,1);  // DELETED after it is copied to ResHV vector -->---+-*
  }                                                             //                    | |
  catch (G4QException& error)                                   //                    | |
	 {                                                             //                    | |
    G4cerr<<"***G4QDiffractionRatio::Fragment: G4Quasmon Exception"<<G4endl;        //| |
    G4Exception("G4QDiffractionRatio::Fragment","72",FatalException,"QuasmonCrash");//| |
  }                                                             //                    | |
  delete pan;                              // Delete the Nuclear Environment <----<---* |
  G4int qNH=leadhs->size();                // A#of collected hadrons from diff.frag.    |
  if(qNH) for(G4int iq=0; iq<qNH; iq++)    // Loop over hadrons to fill the result      |
  {                                        //                                           |
    G4QHadron* loh=(*leadhs)[iq];          // Pointer to the output hadron              |
    ResHV->push_back(loh);                 // Fill in the result                        |
  }                                        //                                           |
  delete leadhs; // <----<----<----<----<----<----<----<----<----<----<----<----<----<--*

		return ResHV; // Result
} // End of Scatter
