//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4QElasticCrossSection.cc,v 1.1 2006-02-06 09:35:57 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: G4QElasticCrossSection for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-Oct-03
// 
//================================================================================

//#define debug
#define edebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4QElasticCrossSection.hh"

// Initialization of the static parameters
const G4int G4QElasticCrossSection::nPoints=50;// #of points in the AMDB tables(>anyPar)(D)
const G4int G4QElasticCrossSection::nLast=nPoints-1; // the Last element in the table   (D)
G4double  G4QElasticCrossSection::lPMin=-2.;  // Min tabulated logarithmic Momentum     (D)
G4double  G4QElasticCrossSection::lPMax= 9.;  // Max tabulated logarithmic Momentum     (D)
G4double  G4QElasticCrossSection::dlnP=(lPMax-lPMin)/nLast;// Log step in the table     (D)
G4bool    G4QElasticCrossSection::onlyCS=true;// Flag to calculate only CS (not Si/Bi)		(L)
G4double  G4QElasticCrossSection::lastSIG=0.; // Last calculated cross section								  (L)
G4double  G4QElasticCrossSection::lastLP=-10.;// Last log(mom_of_the_incident_hadron)   (L)
G4double  G4QElasticCrossSection::lastTM=0.;  // Last t_maximum                         (L)
G4double  G4QElasticCrossSection::theS1=0.;   // The Last mantissa of first difruction  (L)
G4double  G4QElasticCrossSection::theB1=0.;   // The Last slope of first difruction     (L)
G4double  G4QElasticCrossSection::theS2=0.;   // The Last mantissa of second difruction (L)
G4double  G4QElasticCrossSection::theB2=0.;   // The Last slope of second difruction    (L)
G4double  G4QElasticCrossSection::theS3=0.;   // The Last mantissa of third difruction  (L)
G4double  G4QElasticCrossSection::theB3=0.;   // The Last slope of third difruction     (L)
G4int     G4QElasticCrossSection::lastPDG=0;  // Last PDG code of the projectile					
G4int     G4QElasticCrossSection::lastTZ=0;   // Last atomic number of the target				
G4int     G4QElasticCrossSection::lastTN=0;   // Last number of neutrons of the target
G4double  G4QElasticCrossSection::lastPIN=0.; // Last initialized max momentum
G4double* G4QElasticCrossSection::lastCST=0;  // Elastic cross-section table															
G4double* G4QElasticCrossSection::lastPAR=0;  // Parameters for functional calculation					
G4double* G4QElasticCrossSection::lastS1T=0;  // E-dep of mantissa of the first difruction	
G4double* G4QElasticCrossSection::lastB1T=0;  // E-dep of the slope of the first difruction
G4double* G4QElasticCrossSection::lastS2T=0;  // E-dep of mantissa of the second difruction
G4double* G4QElasticCrossSection::lastB2T=0;  // E-dep of the slope of theSecond difruction
G4double* G4QElasticCrossSection::lastS3T=0;  // E-dep of mantissa of the third difruction	
G4double* G4QElasticCrossSection::lastB3T=0;  // E-dep of the slope of the third difruction

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QElasticCrossSection::GetPointer()
{
  static G4QElasticCrossSection theCrossSection; //**Static body of the QEl Cross Section**
  return &theCrossSection;
}


// Calculation of total elastic cross section (p in IU, CS in mb) @@ Units (?)
// F=0 - create AMDB, F=-1 - read&update AMDB, F=1 - update AMDB (sinchro with higher AMDB)
G4double G4QElasticCrossSection::CalculateCrossSection(G4bool CS,G4int F,G4int I,G4int PDG,
                                                       G4int tgZ, G4int tgN, G4double pIU)
{
  // *** Begin of Associative Memory DB for acceleration of the cross section calculations
  static std::vector <G4int>    pPDG;   // Vector of PDG code of the projectile					
  static std::vector <G4double>  PIN;   // Vector of max initialized log(P) in the table
  static std::vector <G4double*> PAR;   // Vector of parameters for functional calculations
  static std::vector <G4double*> CST;   // Vector of cross-section table
  static std::vector <G4double*> S1T;   // Vector of the first mantissa
  static std::vector <G4double*> B1T;   // Vector of the first slope
  static std::vector <G4double*> S2T;   // Vector of the first mantissa
  static std::vector <G4double*> B2T;   // Vector of the first slope
  static std::vector <G4double*> S3T;   // Vector of the first mantissa
  static std::vector <G4double*> B3T;   // Vector of the first slope
  // *** End of Static Definitions (Associative Memory Data Base) ***
  G4double pMom=pIU/GeV;                // All calculations are in GeV
  onlyCS=CS;                            // Flag to calculate only CS (not Si/Bi)
  lastLP=std::log(pMom);                // Make a logarithm of the momentum for calculation
  if(F)                                 // This isotope was found in AMDB =>RETRIEVE/UPDATE
		{
    if(F<0)                             // the AMDB must be loded
    {
      lastPIN = PIN[I];                 // Max log(P) initialised for this table set
      lastPAR = PAR[I];                 // Pointer to the parameter set
      lastCST = CST[I];                 // Pointer to the total sross-section table
      lastS1T = S1T[I];                 // Pointer to the first mantissa
      lastB1T = B1T[I];                 // Pointer to the first slope
      lastS2T = S2T[I];                 // Pointer to the second mantissa
      lastB2T = B2T[I];                 // Pointer to the second slope
      lastS3T = S3T[I];                 // Pointer to the third mantissa
      lastB3T = B3T[I];                 // Pointer to the rhird slope
    }
#ifdef pdebug
    G4cout<<"G4QElasticCrossSection::CalcCS:*read*, LP="<<lastLP<<",PIN="<<lastPIN<<G4endl;
#endif
    lastPIN = GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);// Can update upper logP-Limit in tabs
    PIN[I]=lastPIN;                     // Remember the new P-Limit of the tables
	 }
	 else                                  // This isotope wasn't initialized => CREATE
	 {
    lastPAR = new G4double[nPoints];    // Allocate memory for parameters of CS function
    lastCST = new G4double[nPoints];    // Allocate memory for Tabulated CS function				
    lastS1T = new G4double[nPoints];    // Allocate memory for Tabulated first mantissa	
    lastB1T = new G4double[nPoints];    // Allocate memory for Tabulated first slope				
    lastS2T = new G4double[nPoints];    // Allocate memory for Tabulated second mantissa
    lastB2T = new G4double[nPoints];    // Allocate memory for Tabulated second slope			
    lastS3T = new G4double[nPoints];    // Allocate memory for Tabulated third mantissa	
    lastB3T = new G4double[nPoints];    // Allocate memory for Tabulated third slope    
#ifdef pdebug
    G4cout<<"G4QElasticCrossSection::CalcCS:*ini*, lastLP="<<lastLP<<", -27."<<G4endl;
#endif
    lastPIN = GetPTables(lastLP,-27.,PDG,tgZ,tgN); // Returns the new P-limit for tables
#ifdef pdebug
    G4cout<<"G4QElasticCrossSection::CalcCS: Z="<<tgZ<<", N="<<tgN<<", PDG="<<PDG<<G4endl;
#endif
    pPDG.push_back(lastPDG);            // Fill the associative memmory for the projPDG
    PIN.push_back(lastPIN);             // Fill parameters of CS function to AMDB
    PAR.push_back(lastPAR);             // Fill parameters of CS function to AMDB
    CST.push_back(lastCST);             // Fill Tabulated CS function to AMDB				
    S1T.push_back(lastS1T);             // Fill Tabulated first mantissa to AMDB	
    B1T.push_back(lastB1T);             // Fill Tabulated first slope to AMDB    
    S2T.push_back(lastS2T);             // Fill Tabulated second mantissa to AMDB	
    B2T.push_back(lastB2T);             // Fill Tabulated second slope to AMDB    
    S3T.push_back(lastS3T);             // Fill Tabulated third mantissa to AMDB	
    B3T.push_back(lastB3T);             // Fill Tabulated third slope to AMDB    
	 } // End of creation/update of the new set of parameters and tables
  // ============= NOW Update (if necessary) and Calculate the Cross Section ===========
  if(lastLP>lastPIN && lastLP<=lPMax)
  {
#ifdef pdebug
    G4cout<<"G4QElasticCrossSection::CalcCS:*update*,LP="<<lastLP<<",IN="<<lastPIN<<G4endl;
#endif
    lastPIN = GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);
  }
  if(lastLP>lPMin && lastLP<=lastPIN)   // Linear fit is made using precalculated tables
  {
    if(!onlyCS) lastTM=GetQ2max(PDG, tgZ, tgN, pMom); // Calculate (-t)_max=Q2_max (GeV2)
#ifdef pdebug
    G4cout<<"G4QElasticCrosSec::CalcCS:F="<<onlyCS<<",-t="<<lastTM<<", p="<<lastLP<<G4endl;
#endif
    if(std::fabs(lastLP-lastPIN)<.0001) // Just take the highest tabulated value
    {
      G4double shift=(lastLP-lPMin)/dlnP+.000001; // Log distance from lPMin
      G4int    blast=static_cast<int>(shift); // this is a bin number of the lower edge (0)
      if(blast<0 || blast>=nLast) G4cout<<"G4QEleastCS::CCS:b="<<blast<<","<<nLast<<G4endl;
      lastSIG = lastCST[blast];
      if(!onlyCS)                       // Skip the differential cross-section parameters
      {
        theS1  = lastS1T[blast];
        theB1  = lastB1T[blast];
        theS2  = lastS2T[blast];
        theB2  = lastB2T[blast];
        theS3  = lastS3T[blast];
        theB3  = lastB3T[blast];
      }
#ifdef pdebug
      G4cout<<"G4QElasticCrossSection::CalculateCS:(E) S1="<<theS1<<", B1="<<theB1<<G4endl;
#endif
				}
    else
    {
      G4double shift=(lastLP-lPMin)/dlnP;        // a shift from the beginning of the table
      G4int    blast=static_cast<int>(shift);    // the lower bin number
      if(blast<0)   blast=0;
      if(blast>=nLast) blast=nLast-1;            // low edge of the last bin
      shift-=blast;                              // step inside the unit bin
      G4int lastL=blast+1;                       // the upper bin number
      G4double SIGL=lastCST[blast];              // the basic value of the cross-section
      lastSIG= SIGL+shift*(lastCST[lastL]-SIGL); // calculated total elastic cross-section
#ifdef pdebug
      G4cout<<"G4QElCS::CalcCrossSection: Sig="<<lastSIG<<", P="<<pMom<<", Z="<<tgZ<<", N="
            <<tgN<<", PDG="<<PDG<<", onlyCS="<<onlyCS<<G4endl;
#endif
      if(!onlyCS)                       // Skip the differential cross-section parameters
      {
        G4double S1TL=lastS1T[blast];           // the low bin of the first mantissa
        theS1=S1TL+shift*(lastS1T[lastL]-S1TL); // the basic value of the first mantissa
        G4double B1TL=lastB1T[blast];           // the low bin of the first slope
#ifdef pdebug
        G4cout<<"G4QElCS::CalcCrossSection:bl="<<blast<<",ls="<<lastL<<",SL="<<S1TL<<",SU="
              <<lastB1T[lastL]<<",BL="<<B1TL<<",BU="<<lastB1T[lastL]<<G4endl;
#endif
        theB1=B1TL+shift*(lastB1T[lastL]-B1TL); // the basic value of the first slope
        G4double S2TL=lastS2T[blast];           // the low bin of the second mantissa
        theS2=S2TL+shift*(lastS2T[lastL]-S2TL); // the basic value of the second mantissa
        G4double B2TL=lastB2T[blast];           // the low bin of the second slope
        theB2=B2TL+shift*(lastB2T[lastL]-B2TL); // the basic value of the second slope
        G4double S3TL=lastS3T[blast];           // the low bin of the third mantissa
        theS3=S3TL+shift*(lastS3T[lastL]-S3TL); // the basic value of the third mantissa
        G4double B3TL=lastB3T[blast];           // the low bin of the third slope
        theB3=B3TL+shift*(lastB3T[lastL]-B3TL); // the basic value of the third slope
      }
#ifdef pdebug
      G4cout<<"G4QElasticCrossSection::CalculateCS:(I) S1="<<theS1<<", B1="<<theB1<<G4endl;
#endif
    }
  }
  else lastSIG=GetTabValues(lastLP, PDG, tgZ, tgN); // Direct calculation beyond the table
  if(lastSIG<0.) lastSIG = 0.;                   // @@ a Warning print can be added
  return lastSIG;
}

// It has parameter sets for all tZ/tN/PDG, using them the tables can be created/updated
G4double G4QElasticCrossSection::GetPTables(G4double LP,G4double ILP, G4int PDG, G4int tgZ,
                                                                                 G4int tgN)
{
  static const G4double pwd=2727;
  const G4int n_ppel=22;
  //                       -0-   -1-   -2- -3- -4- -5-   -6-  -7-   -8-    -9-   -10-
  G4double pp_el[n_ppel]={2.648,18.73,.6351,3.,9.,.4186,.3953,4.517,75.07,2.585,.8392,
                          6.132,.8539,.7928,20.,1.,1.724,1.3165,80.,.1304,.32321,81.9};
  //                       -11- -12-  -13- -14--15- -16-  -17- -18-  -19-  -20-  -21-
  if(PDG==2212 && tgZ==1 && tgN==0)
  {
    // -- Total pp elastic cross section cs & s1/b1 (main), s2/b2 (tail1), s3/b3 (tail2) --
    //p2=p*p; p3=p2*p; sp=sqrt(p); lp=log(p); dl1=lp-(3.=par(3)); p4=p2*p2; p=momentum
				//CS=2.648/p2/sp+(18.73+.6351*dl1*dl1+9./p)/(1.+.4186*lp)/(1.+.3953/p4);
    //   par(0)      par(1) par(2)      par(4)      par(5)         par(6)
    //dl2=lp-4.517, s1=(75.07+2.585*dl2*dl2)/(1+.8389/p4), b1=(6.132+.8539*lp)/(1+.7928/p3)
    //       par(7)    par(8) par(9)            par(10)       par(11)par(12)        par(13)
    //s2=(20.+1./p3)/p2; b2=1.724
    //par(14) par(15)       par(16)
    //s3=exp(p/(1.+p^1.3165/80.)); b3=.1304*p^.32321/(1+81.9/p4)
    //              par(17) par(18)   par(19) par(20)  par(21)
    //
    if(lastPAR[nLast]!=pwd)
    {
      for(G4int ip=0; ip<n_ppel; ip++) lastPAR[ip]=pp_el[ip]; // fill the parameters
      lastPAR[nLast]=pwd;
    }
    if(LP>ILP)
				{
      G4int ini = static_cast<int>((ILP-lPMin+.000001)/dlnP)+1; // already initialized till
      if(ini<0)ini=0;
      if(ini<nPoints)
      {
        G4int fin = static_cast<int>((LP-lPMin)/dlnP)+1; // final bin of initialization
        if(fin>nPoints) fin=nLast;                // Limit of the tabular initialization
        if(fin>=ini)
        {
          G4double lp=0.;
          for(G4int ip=ini; ip<=fin; ip++)        // Calculate tabular CS,S1,B1,S2,B2,S3,B3
										{
            lp=lPMin+ip*dlnP;                     // ln(momentum)
            G4bool memCS=onlyCS;
            onlyCS=true;
            lastCST[ip]=GetTabValues(lp, PDG, tgZ, tgN); // Calculate AMDB tables
            onlyCS=memCS;
            lastS1T[ip]=theS1;
            lastB1T[ip]=theB1;
            lastS2T[ip]=theS2;
            lastB2T[ip]=theB2;
            lastS3T[ip]=theS3;
            lastB3T[ip]=theB3;
#ifdef pdebug
            G4cout<<"G4QElasticCrossSection::GetPTables: ip="<<ip<<",lp="<<lp
																		<<",S1="<<theS1<<",B1="<<theB1<<",S2="<<theS2<<",S3="<<theS3<<G4endl;
#endif
          }
          return lp;
        }
        else G4cout<<"*Warning*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
                   <<", N="<<tgN<<", i="<<ini<<" > fin="<<fin<<" nothing is done!"<<G4endl;
      }
      else G4cout<<"*Warning*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
                 <<", N="<<tgN<<", i="<<ini<<">= max="<<nLast<<" nothing is done!"<<G4endl;
    }
#ifdef pdebug
    else G4cout<<"*Warning*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
               <<", N="<<tgN<<", Lp="<<LP<<" <= ILP="<<ILP<<" nothing is done!"<<G4endl;
#endif
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212, Z=1, N=0"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetPTables: only pp is implemented");
  }
  return ILP;
}

// Returns Q2=-t in independent units (all internal calculations are in GeV)
G4double G4QElasticCrossSection::GetExchangeT(G4int tgZ, G4int tgN, G4int PDG)
{
  static const G4double GeVSQ=gigaelectronvolt*gigaelectronvolt;
#ifdef pdebug
  G4cout<<"G4QElasticCS::GetExchangeT:F="<<onlyCS<<",Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<G4endl;
#endif
  if(onlyCS) G4cout<<"*Warning*G4QElasticCrossSection::GetExchangeQ2: onlyCS=true"<<G4endl;
  G4double q2=0.;
  if(PDG==2212 && tgZ==1 && tgN==0)
  {
#ifdef pdebug
    G4cout<<"G4QElasticCS::GetExchangeT: TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1<<",S2="
          <<theS2<<",B2="<<theB2<<",S3="<<theS3<<",B3="<<theB3<<",GeV2="<<GeVSQ<<G4endl;
#endif
    G4double E1=lastTM*theB1;
  		G4double R1=(1.-std::exp(-E1));
#ifdef ppdebug
    G4double ts1=-std::log(1.-R1)/theB1;
    G4double ds1=std::fabs(ts1-lastTM)/lastTM;
    if(ds1>.0001)G4cout<<"*Warn*G4QElCS::GetExT:1 "<<ts1<<"#"<<lastTM<<", d="<<ds1<<G4endl;
#endif
    G4double E2=lastTM*theB2;
  		G4double R2=(1.-std::exp(-E2));
#ifdef ppdebug
    G4double ts2=-std::log(1.-R2)/theB2;
    G4double ds2=std::fabs(ts2-lastTM)/lastTM;
    if(ds2>.0001)G4cout<<"*Warn*G4QElCS::GetExT:2 "<<ts2<<"#"<<lastTM<<", d="<<ds2<<G4endl;
#endif
    G4double E3=lastTM*theB3;
  		G4double R3=(1.-std::exp(-E3));
#ifdef ppdebug
    G4double ts3=-std::log(1.-R3)/theB3;
    G4double ds3=std::fabs(ts3-lastTM)/lastTM;
    if(ds3>.0001)G4cout<<"*Warn*G4QElCS::GetExT:3 "<<ts3<<"#"<<lastTM<<", d="<<ds3<<G4endl;
#endif
  		G4double I1=R1*theS1/theB1;
  		G4double I2=R2*theS2/theB2;
				G4double I3=R3*theS3/theB3;
    G4double I12=I1+I2;
    G4double rand=(I12+I3)*G4UniformRand();
    if     (rand<I1 ) q2=-std::log(1.-R1*G4UniformRand())/theB1;
				else if(rand<I12) q2=-std::log(1.-R2*G4UniformRand())/theB2;
    else              q2=-std::log(1.-R3*G4UniformRand())/theB3;
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetExchangeT: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212, Z=1, N=0"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetExchangeT: only pp is implemented");
  }
  if(q2>lastTM)G4cout<<"*Warning*G4QElasticCrossSect::GetExT:-t="<<q2<<">"<<lastTM<<G4endl;
  return q2*GeVSQ;
}

// lastLP is used, so calculating tables, one need to remember and then recover lastLP
G4double G4QElasticCrossSection::GetTabValues(G4double lp, G4int PDG, G4int tgZ,G4int tgN)
{
#ifdef pdebug
  G4cout<<"G4QElasticCS::GetTabVal: lp="<<lp<<",Z="<<tgZ<<",N="<<tgN<<",PDG="<<PDG<<G4endl;
#endif
  if(PDG==2212 && tgZ==1 && tgN==0)
  {
    G4double p=std::exp(lp);              // momentum
    G4double sp=std::sqrt(p);             // sqrt(p)
    G4double p2=p*p;            
    G4double p4=p2*p2;            
				G4double dl1=lp-lastPAR[3];
    G4double p3=p2*p;
				G4double dl2=lp-lastPAR[7];
    theS1=(lastPAR[8]+lastPAR[9]*dl2*dl2)/(1.+lastPAR[10]/p4);
    theB1=(lastPAR[11]+lastPAR[12]*lp)/(1.+lastPAR[13]/p3);
    theS2=(lastPAR[14]+lastPAR[15]/p3)/p2;
    theB2=lastPAR[16]; 
    theS3=std::exp(-p/(1.+std::pow(p,lastPAR[17])/lastPAR[18]));
    theB3=lastPAR[19]*std::pow(p,lastPAR[20])/(1.+lastPAR[21]/p4); 
#ifdef pdebug
    G4cout<<"G4QElasticCS::GetTableValues: TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1
          <<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS1<<",B3="<<theB1<<G4endl;
#endif
    // Returns the total elastic pp cross-section (to avoid spoiling lastSIG)
    return lastPAR[0]/p2/sp+(lastPAR[1]+lastPAR[2]*dl1*dl1+lastPAR[4]/p)/(1.+lastPAR[5]*lp)
                                                                       /(1.+lastPAR[6]/p4);
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetTabValues: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212, Z=1, N=0"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetTabValues: only pp is implemented");
  }
  return 0.;
}

// Returns max -t=Q2 (GeV^2) for the momentum pP(GeV)PDG and the target nucleus (tgN,tgZ)
G4double G4QElasticCrossSection::GetQ2max(G4int PDG, G4int tgZ, G4int tgN, G4double pP)
{
  //static const G4double mNeut= G4QPDGCode(2112).GetMass()*.001; // MeV to GeV
  static const G4double mProt= G4QPDGCode(2212).GetMass()*.001; // MeV to GeV
  //static const G4double mLamb= G4QPDGCode(3122).GetMass()*.001; // MeV to GeV
  static const G4double mProt2= mProt*mProt;
  if(PDG==2212 && tgZ==1 && tgN==0)
  {
    G4double tMid=std::sqrt(pP*pP+mProt2)*mProt-mProt2; // CMS 90deg value of -t=Q2 (GeV^2)
    return tMid+tMid;
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetQ2max: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212, Z=1, N=0"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetQ2max: only pp is implemented");
  }
}
