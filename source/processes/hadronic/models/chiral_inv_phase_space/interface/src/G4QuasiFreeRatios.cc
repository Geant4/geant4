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
// $Id: G4QuasiFreeRatios.cc,v 1.11.2.1 2007/08/14 09:28:06 mkossov Exp $
// GEANT4 tag $Name: geant4-08-03-patch-01 $
//
//
// G4 Physics class: G4QuasiFreeRatios for N+A elastic cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Oct-06
// 
//================================================================================

//#define debug
//#define pdebug
//#define nandebug

#include "G4QuasiFreeRatios.hh"

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
		std::pair<G4double,G4double> ElTot=GetElTot(pIU, pPDG, tgZ, tgN); // mean hN El&Tot(IU)
  G4double R=0.;
  G4double QF2In=1.;               // Prototype of QuasiFree/Inelastic ratio for hN_tot
  if(ElTot.second>0.)
  {
    R=ElTot.first/ElTot.second;    // Elastic/Total ratio (does not depend on units
    QF2In=GetQF2IN_Ratio(ElTot.second/millibarn, tgZ+tgN); // QuasiFree/Inelastic ratio
  }
  return std::make_pair(QF2In,R);
}

// Calculatio QasiFree/Inelastic Ratio as a function of total hN cross-section (mb) and A
G4double G4QuasiFreeRatios::GetQF2IN_Ratio(G4double s, G4int A)
{
  static const G4int    nps=150;        // Number of steps in the R(s) LinTable
  static const G4int    mps=nps+1;      // Number of elements in the R(s) LinTable
  static const G4double sma=150.;       // The first LinTabEl(s=0)=1., s>sma -> logTab
  static const G4double ds=sma/nps;     // Step of the linear Table
  static const G4int    nls=100;        // Number of steps in the R(lns) logTable
  static const G4int    mls=nls+1;      // Number of elements in the R(lns) logTable
  static const G4double lsi=5.;         // The min ln(s) logTabEl(s=148.4 < sma=150.)
  static const G4double lsa=9.;         // The max ln(s) logTabEl(s=148.4 - 8103. mb)
  static const G4double mi=std::exp(lsi);// The min s of logTabEl(~ 148.4 mb)
  static const G4double ms=std::exp(lsa);// The max s of logTabEl(~ 8103. mb)
  static const G4double dl=(lsa-lsi)/nls;// Step of the logarithmic Table
  static const G4double edl=std::exp(dl);// Multiplication step of the logarithmic Table
  static const G4double toler=.01;      // The tolarence mb defining the same cross-section
  static G4double lastS=0.;             // The last sigma value for which R was calculated
  static G4double lastR=0.;             // The last ratio R which was calculated
  // Local Associative Data Base:
  static std::vector<G4int>     vA;     // Vector of calculated A
  static std::vector<G4double>  vH;     // Vector of max s initialized in the LinTable
  static std::vector<G4int>     vN;     // Vector of topBin number initialized in LinTable
  static std::vector<G4double>  vM;     // Vector of rel max ln(s) initialized in LogTable
  static std::vector<G4int>     vK;     // Vector of topBin number initialized in LogTable
  static std::vector<G4double*> vT;     // Vector of pointers to LinTable in C++ heap
  static std::vector<G4double*> vL;     // Vector of pointers to LogTable in C++ heap
  // Last values of the Associative Data Base:
  static G4int     lastA=0;             // theLast of calculated A
  static G4double  lastH=0.;            // theLast of max s initialized in the LinTable
  static G4int     lastN=0;             // theLast of topBin number initialized in LinTable
  static G4double  lastM=0.;            // theLast of rel max ln(s) initialized in LogTable
  static G4int     lastK=0;             // theLast of topBin number initialized in LogTable
  static G4double* lastT=0;             // theLast of pointer to LinTable in the C++ heap
  static G4double* lastL=0;             // theLast of pointer to LogTable in the C++ heap
  // LogTable is created only if necessary. The ratio R(s>8100 mb) = 0 for any nuclei
  if(s<toler || A<2) return 1.;
  if(s>ms) return 0.;
  if(A>238)
  {
    G4cout<<"-Warning-G4QuasiFreeRatio::GetQF2IN_Ratio:A="<<A<<">238, return zero"<<G4endl;
    return 0.;
  }
  G4int nDB=vA.size();                  // A number of nuclei already initialized in AMDB
  if(nDB && lastA==A && std::fabs(s-lastS)<toler) return lastR;
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
    for(G4int j=1; j<=lastN; j++)       // Calculate LogTab values
    {
      sv+=ds;
      lastT[j]=CalcQF2IN_Ratio(sv,A);
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
        lastL[j]=CalcQF2IN_Ratio(sv,A);
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
          lastT[j]=CalcQF2IN_Ratio(sv,A);
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
          lastL[j]=CalcQF2IN_Ratio(sv,A);
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
} // End of CalcQF2IN_Ratio

// Calculatio QasiFree/Inelastic Ratio as a function of total hN cross-section and A
G4double G4QuasiFreeRatios::CalcQF2IN_Ratio(G4double s, G4int A)
{
  static const G4double C=1.246;
		G4double s2=s*s;
  G4double s4=s2*s2;
		G4double ss=std::sqrt(std::sqrt(s));
  G4double P=7.48e-5*s2/(1.+8.77e12/s4/s4/s2);
  G4double E=.2644+.016/(1.+std::exp((29.54-s)/2.49));
  G4double F=ss*.1526*std::exp(-s2*ss*.0000859);
	 return C*std::exp(-E*std::pow(G4double(A-1.),F))/std::pow(G4double(A),P);
} // End of CalcQF2IN_Ratio

// Calculatio pair(hN_el,hN_tot) (mb): p in GeV/c, index(PDG,F) (see FetchElTot)
std::pair<G4double,G4double> G4QuasiFreeRatios::CalcElTot(G4double p, G4int I)
{
  // ---------> Each parameter set can have not more than nPoints=128 parameters
  static const G4double lmi=3.5;       // min of (lnP-lmi)^2 parabola
  static const G4double pbe=.0557;     // elastic (lnP-lmi)^2 parabola coefficient
  static const G4double pbt=.3;        // total (lnP-lmi)^2 parabola coefficient
  static const G4double pmi=.1;        // Below that fast LE calculation is made
  static const G4double pma=1000.;     // Above that fast HE calculation is made
  // ======================================================================================
  G4double El=0.;                      // prototype of the elastic hN cross-section
  G4double To=0.;                      // prototype of the total hN cross-section
  if(p<=0.)
  {
    G4cout<<"-Warning-G4QuasiFreeRatios::CalcElTot: p="<<p<<" is zero or negative"<<G4endl;
    return std::make_pair(El,To);
  }
  if     (!I)                          // pp/nn
		{
#ifdef debug
				G4cout<<"G4QuasiFreeR::CalcElTot:I=0, p="<<p<<", pmi="<<pmi<<", pma="<<pma<<G4endl;
#endif
    if(p<pmi)
    {
      G4double p2=p*p;
      El=1./(.00012+p2*.2);
      To=El;
#ifdef debug
				  G4cout<<"G4QuasiFreeR::CalcElTot:I=0i, El="<<El<<", To="<<To<<", p2="<<p2<<G4endl;
#endif
    }
    else if(p>pma)
    {
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      El=pbe*lp2+6.72;
      To=pbt*lp2+38.2;
#ifdef debug
				  G4cout<<"G4QuasiFreeR::CalcElTot:I=0a, El="<<El<<", To="<<To<<", lp2="<<lp2<<G4endl;
#endif
    }
    else
    {
      G4double p2=p*p;
      G4double LE=1./(.00012+p2*.2);
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      G4double rp2=1./p2;
      El=LE+(pbe*lp2+6.72+32.6/p)/(1.+rp2/p);
      To=LE+(pbt*lp2+38.2+52.7*rp2)/(1.+2.72*rp2*rp2);
#ifdef debug
				  G4cout<<"G4QuasiFreeR::CalcElTot:0,E="<<El<<",T="<<To<<",s="<<p2<<",l="<<lp2<<G4endl;
#endif
    }
  }
  else if(I==1)                        // np/pn
		{
    if(p<pmi)
    {
      G4double p2=p*p;
      El=1./(.00012+p2*(.051+.1*p2));
      To=El;
    }
    else if(p>pma)
    {
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      El=pbe*lp2+6.72;
      To=pbt*lp2+38.2;
    }
    else
    {
      G4double p2=p*p;
      G4double LE=1./(.00012+p2*(.051+.1*p2));
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      G4double rp2=1./p2;
      El=LE+(pbe*lp2+6.72+30./p)/(1.+.49*rp2/p);
      To=LE+(pbt*lp2+38.2)/(1.+.54*rp2*rp2);
    }
  }
  else if(I==2)                        // pimp/pipn
		{
    G4double lp=std::log(p);
    if(p<pmi)
    {
      G4double lr=lp+1.27;
      El=1.53/(lr*lr+.0676);
      To=El*3;
    }
    else if(p>pma)
    {
      G4double ld=lp-lmi;
      G4double ld2=ld*ld;
      G4double sp=std::sqrt(p);
      El=pbe*ld2+2.4+7./sp;
      To=pbt*ld2+22.3+12./sp;
    }
    else
    {
      G4double lr=lp+1.27;
      G4double LE=1.53/(lr*lr+.0676);
      G4double ld=lp-lmi;
      G4double ld2=ld*ld;
      G4double p2=p*p;
      G4double p4=p2*p2;
      G4double sp=std::sqrt(p);
      G4double lm=lp+.36;
      G4double md=lm*lm+.04;
      G4double lh=lp-.017;
      G4double hd=lh*lh+.0025;
      El=LE+(pbe*ld2+2.4+7./sp)/(1.+.7/p4)+.6/md+.05/hd;
      To=LE*3+(pbt*ld2+22.3+12./sp)/(1.+.4/p4)+1./md+.06/hd;
    }
  }
  else if(I==3)                        // pipp/pimn
		{
    G4double lp=std::log(p);
    if(p<pmi)
    {
      G4double lr=lp+1.27;
      G4double lr2=lr*lr;
      El=13./(lr2+lr2*lr2+.0676);
      To=El;
    }
    else if(p>pma)
    {
      G4double ld=lp-lmi;
      G4double ld2=ld*ld;
      G4double sp=std::sqrt(p);
      El=pbe*ld2+2.4+6./sp;
      To=pbt*ld2+22.3+5./sp;
    }
    else
    {
      G4double lr=lp+1.27;
      G4double lr2=lr*lr;
      G4double LE=13./(lr2+lr2*lr2+.0676);
      G4double ld=lp-lmi;
      G4double ld2=ld*ld;
      G4double p2=p*p;
      G4double p4=p2*p2;
      G4double sp=std::sqrt(p);
      G4double lm=lp-.32;
      G4double md=lm*lm+.0576;
      El=LE+(pbe*ld2+2.4+6./sp)/(1.+3./p4)+.7/md;
      To=LE+(pbt*ld2+22.3+5./sp)/(1.+1./p4)+.8/md;
    }
  }
		else if(I==4)                        // Kmp/Kmn/K0p/K0n
		{

    if(p<pmi)
    {
      G4double psp=p*std::sqrt(p);
      El=5.2/psp;
      To=14./psp;
    }
    else if(p>pma)
    {
      G4double ld=std::log(p)-lmi;
      G4double ld2=ld*ld;
      El=pbe*ld2+2.23;
      To=pbt*ld2+19.5;
    }
    else
    {
      G4double ld=std::log(p)-lmi;
      G4double ld2=ld*ld;
      G4double sp=std::sqrt(p);
      G4double psp=p*sp;
      G4double p2=p*p;
      G4double p4=p2*p2;
      G4double lm=p-.39;
      G4double md=lm*lm+.000156;
      G4double lh=p-1.;
      G4double hd=lh*lh+.0156;
      El=5.2/psp+(pbe*ld2+2.23)/(1.-.7/sp+.075/p4)+.004/md+.15/hd;
      To=14./psp+(pbt*ld2+19.5)/(1.-.21/sp+.52/p4)+.006/md+.30/hd;
    }
  }
  else if(I==5)                        // Kpp/Kpn/aKp/aKn
		{
    if(p<pmi)
    {
      G4double lr=p-.38;
      G4double lm=p-1.;
      G4double md=lm*lm+.372;   
      El=.7/(lr*lr+.0676)+2./md;
      To=El+.6/md;
    }
    else if(p>pma)
    {
      G4double ld=std::log(p)-lmi;
      G4double ld2=ld*ld;
      El=pbe*ld2+2.23;
      To=pbt*ld2+19.5;
    }
    else
    {
      G4double ld=std::log(p)-lmi;
      G4double ld2=ld*ld;
      G4double lr=p-.38;
      G4double LE=.7/(lr*lr+.0676);
      G4double sp=std::sqrt(p);
      G4double p2=p*p;
      G4double p4=p2*p2;
      G4double lm=p-1.;
      G4double md=lm*lm+.372;
      El=LE+(pbe*ld2+2.23)/(1.-.7/sp+.1/p4)+2./md;
      To=LE+(pbt*ld2+19.5)/(1.+.46/sp+1.6/p4)+2.6/md;
    }
  }
  else if(I==6)                        // hyperon-N
		{
    if(p<pmi)
    {
      G4double p2=p*p;
      El=1./(.002+p2*(.12+p2));
      To=El;
    }
    else if(p>pma)
    {
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      G4double sp=std::sqrt(p);
      El=(pbe*lp2+6.72)/(1.+2./sp);
      To=(pbt*lp2+38.2+900./sp)/(1.+27./sp);
    }
    else
    {
      G4double p2=p*p;
      G4double LE=1./(.002+p2*(.12+p2));
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      G4double p4=p2*p2;
      G4double sp=std::sqrt(p);
      El=LE+(pbe*lp2+6.72+99./p2)/(1.+2./sp+2./p4);
      To=LE+(pbt*lp2+38.2+900./sp)/(1.+27./sp+3./p4);
    }
  }
  else if(I==7)                        // antibaryon-N
		{
    if(p>pma)
    {
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      El=pbe*lp2+6.72;
      To=pbt*lp2+38.2;
    }
    else
    {
      G4double ye=std::pow(p,1.25);
      G4double yt=std::pow(p,.35);
      G4double lp=std::log(p)-lmi;
      G4double lp2=lp*lp;
      El=80./(ye+1.)+pbe*lp2+6.72;
      To=(80./yt+.3)/yt+pbt*lp2+38.2;
    }
  }
  else
  {
    G4cout<<"*Error*G4QuasiFreeRatios::CalcElTot:ind="<<I<<" is not defined (0-7)"<<G4endl;
    G4Exception("G4QuasiFreeRatios::CalcElTot:","23",FatalException,"CHIPScrash");
  }
  if(El>To) El=To;
  return std::make_pair(El,To);
} // End of CalcElTot

// Calculatio pair(hN_el,hN_tot)(mb): p in GeV/c, F=true -> N=proton, F=false -> N=neutron
std::pair<G4double,G4double> G4QuasiFreeRatios::FetchElTot(G4double p, G4int PDG, G4bool F)
{
  static const G4int    nlp=300;         // Number of steps in the S(lnp) logTable(5% step)
  static const G4int    mlp=nlp+1;       // Number of elements in the S(lnp) logTable
  static const G4double lpi=-5.;         // The min ln(p) logTabEl(p=6.7 MeV/c - 22. TeV/c)
  static const G4double lpa=10.;         // The max ln(p) logTabEl(p=6.7 MeV/c - 22. TeV/c)
  static const G4double mi=std::exp(lpi);// The min p of logTabEl(~ 6.7 MeV/c)
  static const G4double ma=std::exp(lpa);// The max p of logTabEl(~ 22. TeV)
  static const G4double dl=(lpa-lpi)/nlp;// Step of the logarithmic Table
  static const G4double edl=std::exp(dl);// Multiplication step of the logarithmic Table
  static const G4double toler=.001;      // Relative tolarence defining "the same momentum"
  static G4double lastP=0.;              // The last momentum for which XS was calculated
  static G4int    lastH=0;               // The last projPDG for which XS was calculated
  static G4bool   lastF=true;            // The last nucleon for which XS was calculated
  static std::pair<G4double,G4double> lastR=std::make_pair(0.,0.); // The last result
  // Local Associative Data Base:
  static std::vector<G4int>     vI;      // Vector of index for which XS was calculated
  static std::vector<G4double>  vM;      // Vector of rel max ln(p) initialized in LogTable
  static std::vector<G4int>     vK;      // Vector of topBin number initialized in LogTable
  static std::vector<G4double*> vE;      // Vector of ElastPointers to LogTable in C++ heap
  static std::vector<std::pair<G4double,G4double>*> vX; // Vector of ETPointers to LogTable
  // Last values of the Associative Data Base:
  static G4int     lastI=0;              // The Last index for which XS was calculated
  static G4double  lastM=0.;             // The Last rel max ln(p) initialized in LogTable
  static G4int     lastK=0;             // The Last topBin number initialized in LogTable
  static std::pair<G4double,G4double>* lastX=0; // The Last ETPointers to LogTable in heap
  // LogTable is created only if necessary. The ratio R(s>8100 mb) = 0 for any nuclei
  G4int nDB=vI.size();                   // A number of hadrons already initialized in AMDB
#ifdef pdebug
		G4cout<<"G4QuasiFreeR::FetchElTot:p="<<p<<",PDG="<<PDG<<",F="<<F<<",nDB="<<nDB<<G4endl;
#endif
  if(nDB && lastH==PDG && lastF==F && p>0. && std::fabs(p-lastP)/p<toler) return lastR;
  lastH=PDG;
  lastF=F;
  G4int ind=-1;                          // Prototipe of the index of the PDG/F combination
  // i=0: pp(nn), i=1: np(pn), i=2: pimp(pipn), i=3: pipp(pimn), i=4: Kmp(Kmn,K0n,K0p),
		// i=5: Kpp(Kpn,aK0n,aK0p), i=6: Hp(Hn), i=7: app(apn,ann,anp) 
  G4bool kfl=true;                             // Flag of K0/aK0 oscillation
  G4bool kf=false;
  if(PDG==130||PDG==310)
  {
    kf=true;
    if(G4UniformRand()>.5) kfl=false;
  }
  if     (PDG==2212&&F || PDG==2112&&!F) ind=0; // pp/nn
  else if(PDG==2112&&F || PDG==2212&&!F) ind=1; // np/pn
  else if(PDG==-211&&F || PDG== 211&&!F) ind=2; // pimp/pipn
  else if(PDG== 211&&F || PDG==-211&&!F) ind=3; // pipp/pimn
  else if(PDG==-321 || PDG==-311 || (kf &&!kfl)) ind=4; // KmN/K0N
		else if(PDG== 321 || PDG== 311 || (kf && kfl)) ind=5; // KpN/aK0N
  else if(PDG> 3000 && PDG< 3335) ind=6;        // @@ for all hyperons - take Lambda
  else if(PDG<-2000 && PDG>-3335) ind=7;        // @@ for all anti-baryons - anti-p/anti-n
  else
  {
    G4cout<<"*Error*G4QuasiFreeRatios::FetchElTot: PDG="<<PDG
          <<", while it is defined only for p,n,hyperons,anti-baryons,pi,K/antiK"<<G4endl;
    G4Exception("G4QuasiFreeRatio::FetchElTot:","22",FatalException,"CHIPScrash");
  }
  if(nDB && lastI==ind && p>0. && std::fabs(p-lastP)/p<toler) return lastR;
  if(p<=mi || p>=ma) return CalcElTot(p,ind);   // @@ Slow calculation ! (Warning?)
  G4bool found=false;
  G4int i=-1;
		if(nDB) for (i=0; i<nDB; i++) if(ind==vI[i])  // Sirch for this index in AMDB
  {
    found=true;                                 // The index is found
    break;
  }
  G4double lp=std::log(p);
#ifdef pdebug
		G4cout<<"G4QuasiFreeR::FetchElTot:I="<<ind<<",i="<<i<<",fd="<<found<<",lp="<<lp<<G4endl;
#endif
  if(!nDB || !found)                            // Create new line in the AMDB
	 {
    lastX = new std::pair<G4double,G4double>[mlp]; // Create logarithmic Table for ElTot
    lastI = ind;                                // Remember the initialized inex
    lastK = static_cast<int>((lp-lpi)/dl)+1;    // MaxBin to be initialized in LogTaB
    if(lastK>nlp)
    {
      lastK=nlp;
      lastM=lpa-lpi;
    }
    else lastM = lastK*dl;               // Calculate max initialized ln(p)-lpi for LogTab
    G4double pv=mi;
    for(G4int j=0; j<=lastK; j++)        // Calculate LogTab values
    {
      lastX[j]=CalcElTot(pv,ind);
#ifdef pdebug
		    G4cout<<"G4QuasiFreeR::FetchElTot:I,j="<<j<<",pv="<<pv<<",E="<<lastX[j].first<<",T="
            <<lastX[j].second<<G4endl;
#endif
      if(j!=lastK) pv*=edl;
    }
    i++;                                 // Make a new record to AMDB and position on it
    vI.push_back(lastI);
    vM.push_back(lastM);
    vK.push_back(lastK);
    vX.push_back(lastX);
  }
  else                                   // The A value was found in AMDB
	 {
    lastI=vI[i];
    lastM=vM[i];
    lastK=vK[i];
    lastX=vX[i];
    G4int nextK=lastK+1;
    G4double lpM=lastM+lpi;
#ifdef pdebug
		  G4cout<<"G4QuasiFreeR::FetchElTo:M="<<lpM<<",l="<<lp<<",K="<<lastK<<",n="<<nlp<<G4endl;
#endif
    if(lp>lpM && lastK<nlp)              // LogTab must be updated
				{
      lastK = static_cast<int>((lp-lpi)/dl)+1; // MaxBin to be initialized in LogTab
#ifdef pdebug
		    G4cout<<"G4QuasiFreeR::FetET:K="<<lastK<<",lp="<<lp<<",li="<<lpi<<",dl="<<dl<<G4endl;
#endif
      if(lastK>nlp)
      {
        lastK=nlp;
        lastM=lpa-lpi;
      }
      else lastM = lastK*dl;           // Calculate max initialized ln(p)-lpi for LogTab
      G4double pv=std::exp(lpM);       // momentum of the last calculated beam
      for(G4int j=nextK; j<=lastK; j++)// Calculate LogTab values
      {
        pv*=edl;
        lastX[j]=CalcElTot(pv,ind);
#ifdef pdebug
		      G4cout<<"G4QuasiFreeR::FetchElTot:U:j="<<j<<",p="<<pv<<",E="<<lastX[j].first<<",T="
              <<lastX[j].second<<G4endl;
#endif
      }
    } // End of LogTab update
    if(lastK>=nextK)                   // The AMDB was apdated
				{
      vM[i]=lastM;
      vK[i]=lastK;
    }
  }
  // Now one can use tabeles to calculate the value
  G4double dlp=lp-lpi;                       // Shifted log(p) value
  G4int n=static_cast<int>(dlp/dl);          // Low edge number of the bin
  G4double d=dlp-n*dl;                       // Log shift
  G4double e=lastX[n].first;                 // E-Base
  lastR.first=e+d*(lastX[n+1].first-e)/dl;   // E-Result
  if(lastR.first<0.)  lastR.first = 0.;
  G4double t=lastX[n].second;                // T-Base
  lastR.second=t+d*(lastX[n+1].second-t)/dl; // T-Result
  if(lastR.second<0.) lastR.second= 0.;
  if(lastR.first>lastR.second) lastR.first = lastR.second;
  return lastR;
} // End of FetchElTot

// (Mean Elastic and Mean Total) Cross-Sections (mb) for PDG+(Z,N) at P=p[GeV/c]
std::pair<G4double,G4double> G4QuasiFreeRatios::GetElTot(G4double pIU, G4int hPDG,
                                                         G4int Z,       G4int N)
{
  G4double pGeV=pIU/gigaelectronvolt;
  if(Z<1 && N<1)
  {
    G4cout<<"-Warning-G4QuasiFreeRatio::GetElTot:Z="<<Z<<",N="<<N<<", return zero"<<G4endl;
    return std::make_pair(0.,0.);
  }
  std::pair<G4double,G4double> pA=FetchElTot(pGeV, hPDG, true);
  std::pair<G4double,G4double> nA=FetchElTot(pGeV, hPDG, false);
  G4double A=(Z+N)/millibarn;                // To make the result in independent units(IU)
  return std::make_pair((Z*pA.first+N*nA.first)/A,(Z*pA.second+N*nA.second)/A);
} // End of GetElTot

// scatter (pPDG,p4M) on a virtual nucleon (NPDG,N4M), result: final pair(newN4M,newp4M)
// if(newN4M.e()==0.) - below threshold, XS=0, no scattering of the progectile happened
std::pair<G4LorentzVector,G4LorentzVector> G4QuasiFreeRatios::Scatter(G4int NPDG,
                                     G4LorentzVector N4M, G4int pPDG, G4LorentzVector p4M)
{
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  G4LorentzVector pr4M=p4M/megaelectronvolt;   // Convert 4-momenta in MeV (keep p4M)
  N4M/=megaelectronvolt;
  G4LorentzVector tot4M=N4M+p4M;
  G4double mT=mNeut;
  G4int Z=0;
  G4int N=1;
  if(NPDG==2212)
  {
    mT=mProt;
    Z=1;
    N=0;
  }
  else if(NPDG!=2112)
  {
    G4cout<<"Error:G4QuasiFreeRatios::Scatter:NPDG="<<NPDG<<" is not 2212 or 2112"<<G4endl;
    G4Exception("G4QuasiFreeRatios::Scatter:","21",FatalException,"CHIPScomplain");
    //return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M);// Use this if not exception
  }
  G4double mT2=mT*mT;
  G4double mP2=pr4M.m2();
  G4double E=(tot4M.m2()-mT2-mP2)/(mT+mT);
  G4double E2=E*E;
  if(E<0. || E2<mP2)
  {
#ifdef pdebug
    G4cerr<<"-Warning-G4QFR::Scat:*Negative Energy*E="<<E<<",E2="<<E2<<"<M2="<<mP2<<G4endl;
#endif
    return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); // Do Nothing Action
  }
		G4double P=std::sqrt(E2-mP2);                   // Momentum in pseudo laboratory system
  G4VQCrossSection* CSmanager=G4QElasticCrossSection::GetPointer();
#ifdef debug
  G4cout<<"G4QFR::Scatter: Before XS, P="<<P<<", Z="<<Z<<", N="<<N<<", PDG="<<pPDG<<G4endl;
#endif
  // @@ Temporary NN t-dependence for all hadrons
  if(pPDG>3400 || pPDG<-3400) G4cout<<"-Warning-G4QElast::Scatter: pPDG="<<pPDG<<G4endl;
  G4int PDG=2212;                                                // *TMP* instead of pPDG
  G4double xSec=CSmanager->GetCrossSection(false, P, Z, N, PDG); // Rec.CrossSect *TMP*
  //G4double xSec=CSmanager->GetCrossSection(false, P, Z, N, pPDG); // Rec.CrossSect
#ifdef debug
  G4cout<<"G4QElast::Scatter:pPDG="<<pPDG<<",P="<<P<<",CS="<<xSec/millibarn<<G4endl;
#endif
#ifdef nandebug
  if(xSec>0. || xSec<0. || xSec==0);
  else  G4cout<<"******G4QElast::Scatter:xSec="<<xSec/millibarn<<G4endl;
#endif
  // @@ check a possibility to separate p, n, or alpha (!)
  if(xSec <= 0.) // The cross-section iz 0 -> Do Nothing
  {
#ifdef pdebug
    G4cerr<<"-Warning-G4QFR::Scat:**Zero XS**PDG="<<pPDG<<",NPDG="<<NPDG<<",P="<<P<<G4endl;
#endif
    return std::make_pair(G4LorentzVector(0.,0.,0.,0.),p4M); //Do Nothing Action
  }
  G4double mint=CSmanager->GetExchangeT(Z,N,PDG); // functional randomized -t (MeV^2) *TMP*
  //G4double mint=CSmanager->GetExchangeT(Z,N,pPDG); // functional randomized -t in MeV^2
#ifdef pdebug
  G4cout<<"G4QFR::SCAT:PDG="<<pPDG<<", P="<<Momentum<<", CS="<<xSec<<", -t="<<mint<<G4endl;
#endif
#ifdef nandebug
  if(mint>-.0000001);
  else  G4cout<<"******G4QFR::Scat:-t="<<mint<<G4endl;
#endif
  G4double cost=1.-(mint+mint)/CSmanager->GetHMaxT(); // cos(theta) in CMS
#ifdef ppdebug
  G4cout<<"G4QFR::Scat:-t="<<mint<<",dpc2="<<CSmanager->GetHMaxT()<<",cost="<<cost<<G4endl;
#endif
  if(cost>1. || cost<-1. || !(cost>-1. || cost<=1.))
  {
    if     (cost>1.)  cost=1.;
    else if(cost<-1.) cost=-1.;
    else
				{
      G4cerr<<"G4QFR::S:*NAN*c="<<cost<<",t="<<mint<<",tm="<<CSmanager->GetHMaxT()<<G4endl;
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
#ifdef debug
		G4cout<<"G4QFR::Scat:p4M="<<p4M<<"+r4M="<<reco4M<<"="<<scat4M+reco4M<<"="<<tot4M<<G4endl;
#endif
		return std::make_pair(reco4M*megaelectronvolt,pr4M*megaelectronvolt); // Result
} // End of Scatter
