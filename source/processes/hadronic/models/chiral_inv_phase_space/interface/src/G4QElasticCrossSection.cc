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
// $Id: G4QElasticCrossSection.cc,v 1.9 2006/06/29 20:08:34 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
G4double  G4QElasticCrossSection::lPMin=-3.;  // Min tabulated logarithmic Momentum     (D)
G4double  G4QElasticCrossSection::lPMax= 9.;  // Max tabulated logarithmic Momentum     (D)
G4double  G4QElasticCrossSection::dlnP=(lPMax-lPMin)/nLast;// Log step in the table     (D)
G4bool    G4QElasticCrossSection::onlyCS=true;// Flag to calculate only CS (not Si/Bi)		(L)
G4double  G4QElasticCrossSection::lastSIG=0.; // Last calculated cross section								  (L)
G4double  G4QElasticCrossSection::lastLP=-10.;// Last log(mom_of_the_incident_hadron)   (L)
G4double  G4QElasticCrossSection::lastTM=0.;  // Last t_maximum                         (L)
G4double  G4QElasticCrossSection::theSS=0.;   // The Last sq.slope of first difruction  (L)
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
G4double* G4QElasticCrossSection::lastSST=0;  // E-dep of sq.slope of the first difruction	
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
  static std::vector <G4double*> SST;   // Vector of the first squared slope
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
      lastSST = SST[I];                 // Pointer to the first squared slope
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
    if(lastLP>lastPIN && lastLP<lPMax)
    {
      lastPIN = GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);// Can update upper logP-Limit in tabs
      PIN[I]=lastPIN;                     // Remember the new P-Limit of the tables
    }
	 }
	 else                                  // This isotope wasn't initialized => CREATE
	 {
    lastPAR = new G4double[nPoints];    // Allocate memory for parameters of CS function
    lastPAR[nLast]=0;                   // Initialization for VALGRIND
    lastCST = new G4double[nPoints];    // Allocate memory for Tabulated CS function				
    lastSST = new G4double[nPoints];    // Allocate memory for Tabulated first sqaredSlope	
    lastS1T = new G4double[nPoints];    // Allocate memory for Tabulated first mantissa	
    lastB1T = new G4double[nPoints];    // Allocate memory for Tabulated first slope				
    lastS2T = new G4double[nPoints];    // Allocate memory for Tabulated second mantissa
    lastB2T = new G4double[nPoints];    // Allocate memory for Tabulated second slope			
    lastS3T = new G4double[nPoints];    // Allocate memory for Tabulated third mantissa	
    lastB3T = new G4double[nPoints];    // Allocate memory for Tabulated third slope    
#ifdef pdebug
    G4cout<<"G4QElasticCrossSection::CalcCS:*ini*,lastLP="<<lastLP<<",min="<<lPMin<<G4endl;
#endif
    lastPIN = GetPTables(lastLP,lPMin,PDG,tgZ,tgN); // Returns the new P-limit for tables
#ifdef pdebug
    G4cout<<"G4QElasticCrossSection::CalcCS: Z="<<tgZ<<", N="<<tgN<<", PDG="<<PDG<<G4endl;
#endif
    pPDG.push_back(lastPDG);            // Fill the associative memmory for the projPDG
    PIN.push_back(lastPIN);             // Fill parameters of CS function to AMDB
    PAR.push_back(lastPAR);             // Fill parameters of CS function to AMDB
    CST.push_back(lastCST);             // Fill Tabulated CS function to AMDB				
    SST.push_back(lastSST);             // Fill Tabulated first sq.slope to AMDB	
    S1T.push_back(lastS1T);             // Fill Tabulated first mantissa to AMDB	
    B1T.push_back(lastB1T);             // Fill Tabulated first slope to AMDB    
    S2T.push_back(lastS2T);             // Fill Tabulated second mantissa to AMDB	
    B2T.push_back(lastB2T);             // Fill Tabulated second slope to AMDB    
    S3T.push_back(lastS3T);             // Fill Tabulated third mantissa to AMDB	
    B3T.push_back(lastB3T);             // Fill Tabulated third slope to AMDB    
	 } // End of creation/update of the new set of parameters and tables
  // ============= NOW Update (if necessary) and Calculate the Cross Section ===========
  if(lastLP>lastPIN && lastLP<lPMax)
  {
#ifdef pdebug
    G4cout<<"G4QElasticCrossSection::CalcCS:*update*,LP="<<lastLP<<",IN="<<lastPIN<<G4endl;
#endif
    lastPIN = GetPTables(lastLP,lastPIN,PDG,tgZ,tgN);
  }
#ifdef pdebug
  G4cout<<"G4QElastCS::CalcCS: lastLP="<<lastLP<<",lPM="<<lPMin<<",lPIN="<<lastPIN<<G4endl;
#endif
  if(!onlyCS) lastTM=GetQ2max(PDG, tgZ, tgN, pMom); // Calculate (-t)_max=Q2_max (GeV2)
#ifdef pdebug
  G4cout<<"G4QElasticCrosSec::CalcCS:F="<<onlyCS<<",-t="<<lastTM<<", p="<<lastLP<<G4endl;
#endif
  if(lastLP>lPMin && lastLP<=lastPIN)   // Linear fit is made using precalculated tables
  {
    if(std::fabs(lastLP-lastPIN)<.0001) // Just take the highest tabulated value
    {
      G4double shift=(lastLP-lPMin)/dlnP+.000001; // Log distance from lPMin
      G4int    blast=static_cast<int>(shift); // this is a bin number of the lower edge (0)
      if(blast<0 || blast>=nLast) G4cout<<"G4QEleastCS::CCS:b="<<blast<<","<<nLast<<G4endl;
      lastSIG = lastCST[blast];
      if(!onlyCS)                       // Skip the differential cross-section parameters
      {
        theSS  = lastSST[blast];
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
        G4double SSTL=lastSST[blast];           // the low bin of the first squared slope
        theSS=SSTL+shift*(lastSST[lastL]-SSTL); // the basic value of the first sq.slope
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
#ifdef pdebug
        G4cout<<"G4QElCS::CCS: s3l="<<S3TL<<",sh3="<<shift<<",s3h="<<lastS3T[lastL]<<",b="
              <<blast<<",l="<<lastL<<G4endl;
#endif
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
  const G4int n_npel=24;                // #of parameters for np-elastic (<nPoints=50)
  const G4int n_ppel=32;                // #of parameters for pp-elastic (<nPoints=50)
  const G4int n_pdel=30;                // #of parameters for pp-elastic (<nPoints=50)
  const G4int n_phe4=32;                // #of parameters for phe4-elastic (<nPoints=50)
  //                      -0- -1-  -2- -3- -4-  -5- -6- -7- -8- -9--10--11--12--13- -14-
  G4double np_el[n_npel]={12.,.05,.0001,5.,.35,6.75,.14,19.,.6,6.75,.14,13.,.14,.6,.00013,
                          75.,.001,7.2,4.32,.012,2.5,0.0,12.,.34};
  //                      -15--16--17- -18- -19--20--21--22--23-
  //                       -0-   -1-  -2- -3- -4- -5-  -6-  -7-  -8--9--10--11--12--13-
  G4double pp_el[n_ppel]={2.865,18.9,.6461,3.,9.,.425,.4276,.0022,5.,74.,3.,3.4,.2,.17,
                          .001,8.,.055,3.64,5.e-5,4000.,1500.,.46,1.2e6,3.5e6,5.e-5,1.e10,
                          8.5e8,1.e10,1.1,3.4e6,6.8e6,0.};
  //                      -14--15- -16- -17- -18-  -19- -20- -21- -22-  -23-   -24-  -25-
  //                       -26- -27-  -28- -29- -30- -31-
  //                      -0- -1--2- -3-  -4-  -5-  -6-    -7-  -8- -9- -10--11--12--13-
  G4double pd_el[n_pdel]={7.5,.09,5.,.06,.0013,.1,.000006,1.287,.001,34.,.4,4.726,7.5,.1,
                          .05,.000017,.0004,1.15,5.5,.13,.02,.1911,4.857,46.,40.,2.,.01,
                          .0000007,13.,.1};
  //                       -14- -15-  -16-  -17--18- -19--20- -21- -22- -23--24--25--26-
  //                        -27-  -28- -29-
  //                      -0- -1- -2- -3-  -4-   -5-   -6-  -7-  -8- -9- -10--11--12-
  G4double p_he4[n_phe4]={22.5,.3,5.,.037,.004,.00003,12.9,2500.,.05,22.5,.3,.04,.0065,
																										12.5,2500.,.053,.035,2.5e-8,29.,.6,.7,.1,1.5,.03,.00015,4.24,
                          .004,.005,1.5e-11,2.3,.002,23.};
  //                      -13- -14-  -15- -16- -17- -18--19--20-21--22--23-  -24- -25-
  //                      -26- -27-   -28- -29- -30- -31-
  if((PDG==2212 || PDG==2112) && (tgZ==1&&tgN==1 || tgZ==2&&tgN==2 || tgZ==1&&!tgN))
  {
    // --- Total np elastic cross section cs & s1/b1 (t), s2/b2 (u) --- NotTuned for highE
    //p2=p*p;p3=p2*p;sp=sqrt(p);p2s=p2*sp;lp=log(p);dl1=lp-(5.=par(3));p4=p2*p2; p=|3-mom|
				//CS=12./(p2s+.05*p+.0001/sqrt(sp))+.35/p+(6.75+.14*dl1*dl1+19./p)/(1.+.6/p4);
    //  par(0)   par(1) par(2)        par(4) par(5) par(6)     par(7)     par(8)
    //s1=(6.75+.14*dl2*dl2+13./p)/(1.+.14/p4)+.6/(p4+.00013), s2=(75.+.001/p4/p)/p3
    //  par(9) par(10)   par(11)     par(12) par(13) par(14)  par(15) par(16)
    //b1=(7.2+4.32/(p4*p4+.012*p3))/(1.+2.5/p4), ss=0., b2=12./(p*sp+.34)
    //par(17) par(18)    par(19)       par(20)  par(21)   par(22)   par(23)
    //
    // -- Total pp elastic cross section cs & s1/b1 (main), s2/b2 (tail1), s3/b3 (tail2) --
    //p2=p*p;p3=p2*p;sp=sqrt(p);p2s=p2*sp;lp=log(p);dl1=lp-(3.=par(3));p4=p2*p2; p=|3-mom|
				//CS=2.865/p2s/(1+.0022/p2s)+(18.9+.6461*dl1*dl1+9./p)/(1.+.425*lp)/(1.+.4276/p4);
    //   par(0)       par(7)     par(1) par(2)      par(4)      par(5)         par(6)
    //dl2=lp-5., s1=(74.+3.*dl2*dl2)/(1+3.4/p4/p)+(.2/p2+17.*p)/(p4+.001*sp),
    //     par(8) par(9) par(10)        par(11)   par(12)par(13)    par(14)
    // b1=8.*p**.055/(1.+3.64/p3); s2=5.e-5+4000./(p4+1500.*p); b2=.46+1.2e6/(p4+3.5e6/sp);
    // par(15) par(16)  par(17)     par(18) par(19)  par(20)   par(21) par(22)  par(23)
    // s3=5.e-5+1.e10/(p4*p4+8.5e8*p2+1.e10); b3=1.1+3.4e6/(p4+6.8e6); ss=0.
    //  par(24) par(25)     par(26)  par(27) par(28) par(29)  par(30)   par(31)
    //
    //-- Total pd elastic cross section cs & s1/b1/ss(main), s2/b2(tail1), s3/b3(u-chan) --
    //p2=p*p;p3=p2*p;sp=sqrt(p);p2s=p2*sp;lp=log(p);dl1=lp-(5.=par(2));p4=p2*p2; p=|3-mom|
    //CS=(7.5+.09*dl1*dl1)*(1+.06/(p4+.0013*sp)+.1/p)/(1+.000006/p3),
    // par(0) par(1)         par(3)   par(4)    par(5)     par(6)
    // b1=1.287/p2+.001/p4+34./(1.+.4/p)
    //    par(7)  par(8)  par(9)  par(10)
    //dl2=lp-4.726, s1=(7.5+.1*dl2*dl2)*(1.+.05/(p4+.000017/p+.0004)), ss=1.15/p4-5.5/p
    //      par(11) par(12) par(13)        par(14)  par(15)   par(16)    par(17)  par(18)
    //psp=p*sp, s2=.13+.02*lp+(1.+.1911/p2)/p2, b2=(4.857+46./psp)/(1.+40./psp)
    //         par(19) par(20)   par(21)          par(22) par(23)      par(24)
    // s3=2./((p3+.01)*p3+.0000007), b3=13./(p3+.1*sp)
    //  par(25)   par(26) par(27)      par(28)  par(29)
    //
    //-- Total pHe4 elastic cross section cs & s1/b1/ss(main), s2/b2(tail), s3/b3(u-chan)--
    //p2=p*p;p3=p2*p;sp=sqrt(p);p2s=p2*sp;lp=log(p);dl1=lp-(5.=par(2));p4=p2*p2; p=|3-mom|
    //CS=(22.5+.3*dl1*dl1)*(1+.037/(p4+.004*p+.00003))-12.9/(1+2500*dl2*dl2)/p2, dl2=p-.05
    // par(0) par(1)         par(3)    par(4) par(5)   par(6)  par(7)                par(8)
    // s1=(22.5+.3*dl1*dl1)*(1.+.04/(p4+.0065*p))-12.5/(1+2500*dl3*dl3)/p2, dl3=p-.053
    //   par(9) par(10)      par(11)    par(12)  par(13) par(14)                 par(15)
    // p3sp=p3*sp, b1=.035/(p3sp+2.5e-8/p2)+29.+.6*lp, ss=.7/p4
    //                par(16)   par(17) par(18) par(19)  par(20)  
    // s2=.1+1.5/p+.03/(p4+.00015)/p2, b2=4.24,
    // par(21)par(22)par(23)par(24)     par(25)
    // p5=p4*p; p6=p4*p2, s3=.004/((p6+.005)*p5+1.5e-11/p2), b3=2.3/(p2+.002)+23.
    //                    par(26)    par(27)    par(28)       par(29) par(30) par(31)
    if(lastPAR[nLast]!=pwd) // A unique flag to avoid the repeatable definition
    {
      if     (PDG==2112&&tgZ==1&&tgN==0)
                                for(G4int ip=0;ip<n_npel;ip++)lastPAR[ip]=np_el[ip];// np
						else if(PDG==2212&&tgZ==1&&tgN==0)
                                for(G4int ip=0;ip<n_ppel;ip++)lastPAR[ip]=pp_el[ip];// pp
      else if((PDG==2212||PDG==2112)&&tgZ==1&&tgN==1)
                                for(G4int ip=0;ip<n_pdel;ip++)lastPAR[ip]=pd_el[ip];// pd
      else if((PDG==2212||PDG==2112)&&tgZ==2&&tgN==2)
                                for(G4int ip=0;ip<n_phe4;ip++)lastPAR[ip]=p_he4[ip];// phe4
      else G4cout<<"*Warning*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
               <<", N="<<tgN<<" isn't supported. No initialization of parameters!"<<G4endl;
      lastPAR[nLast]=pwd;
      // and initialize the zero element of the table
      G4double lp=lPMin;                                      // ln(momentum)
      G4bool memCS=onlyCS;                                    // ??
      onlyCS=true;
      lastCST[0]=GetTabValues(lp, PDG, tgZ, tgN); // Calculate AMDB tables
      onlyCS=memCS;
      lastSST[0]=theSS;
      lastS1T[0]=theS1;
      lastB1T[0]=theB1;
      lastS2T[0]=theS2;
      lastB2T[0]=theB2;
      lastS3T[0]=theS3;
      lastB3T[0]=theB3;
#ifdef pdebug
      G4cout<<"G4QElasticCrossSection::GetPTables:ip=0(init), lp="<<lp<<",S1="<<theS1
												<<",B1="<<theB1<<",S2="<<theS2<<",S3="<<theS3<<",B3="<<theB3<<G4endl;
#endif
    }
    if(LP>ILP)
				{
      G4int ini = static_cast<int>((ILP-lPMin+.000001)/dlnP)+1; // already inited till this
      if(ini<0) ini=0;
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
            lastSST[ip]=theSS;
            lastS1T[ip]=theS1;
            lastB1T[ip]=theB1;
            lastS2T[ip]=theS2;
            lastB2T[ip]=theB2;
            lastS3T[ip]=theS3;
            lastB3T[ip]=theB3;
#ifdef pdebug
            G4cout<<"G4QElasticCrossSection::GetPTables:ip="<<ip<<",lp="<<lp<<",S1="<<theS1
                  <<",B1="<<theB1<<",S2="<<theS2<<",S3="<<theS3<<",B3="<<theB3<<G4endl;
#endif
          }
          return lp;
        }
        else G4cout<<"*Warning*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
                   <<", N="<<tgN<<", i="<<ini<<" > fin="<<fin<<", LP="<<LP<<" > ILP="<<ILP
                   <<" nothing is done!"<<G4endl;
      }
      else G4cout<<"*Warning*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
                 <<", N="<<tgN<<", i="<<ini<<">= max="<<nPoints<<", LP="<<LP<<" > ILP="
                 <<ILP<<", lPMax="<<lPMax<<" nothing is done!"<<G4endl;
    }
#ifdef pdebug
    else G4cout<<"*Warning*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ
               <<", N="<<tgN<<", LP="<<LP<<" <= ILP="<<ILP<<" nothing is done!"<<G4endl;
#endif
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetPTables: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212/2112, Z=1, N=0,1"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetPTables: np,pp,pd,pHe are implemented");
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
  if(PDG==2112 && tgZ==1 && tgN==0)                // ===> n+p=n+p
  {
#ifdef pdebug
    G4cout<<"G4QElasticCS::GetExchangeT: TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1<<",S2="
          <<theS2<<",B2="<<theB2<<",GeV2="<<GeVSQ<<G4endl;
#endif
    G4double E1=lastTM*theB1;
  		G4double R1=(1.-std::exp(-E1));
#ifdef ppdebug
    G4double ts1=-std::log(1.-R1)/theB1;
    G4double ds1=std::fabs(ts1-lastTM)/lastTM;
    if(ds1>.0001)G4cout<<"*Warn*G4QElCS::GetExT:1n "<<ts1<<"#"<<lastTM<<",d="<<ds1<<G4endl;
#endif
    G4double E2=lastTM*theB2;
  		G4double R2=(1.-std::exp(-E2));
#ifdef ppdebug
    G4double ts2=-std::log(1.-R2)/theB2;
    G4double ds2=std::fabs(ts2-lastTM)/lastTM;
    if(ds2>.0001)G4cout<<"*Warn*G4QElCS::GetExT:2n "<<ts2<<"#"<<lastTM<<",d="<<ds2<<G4endl;
#endif
    //G4double E3=lastTM*theB3;
  		//G4double R3=(1.-std::exp(-E3));
#ifdef ppdebug
    //G4double ts3=-std::log(1.-R3)/theB3;
    //G4double ds3=std::fabs(ts3-lastTM)/lastTM;
    //if(ds3>.01)G4cout<<"*Warn*G4QElCS::GetExT:3n "<<ts3<<"#"<<lastTM<<",d="<<ds3<<G4endl;
#endif
  		G4double I1=R1*theS1;
  		G4double I2=R2*theS2/theB2;
				//G4double I3=R3*theS3/theB3;
    G4double I12=I1+I2;
    //G4double rand=(I12+I3)*G4UniformRand();
    G4double rand=I12*G4UniformRand();
    if     (rand<I1 ) q2=-std::log(1.-R1*G4UniformRand())/theB1;       // t-chan
				else              q2=lastTM+std::log(1.-R2*G4UniformRand())/theB2; // u-chan (ChEx)
  }
  else if(PDG==2212 && tgZ==1 && tgN==0)                // ===> p+p=p+p
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
    if(ds1>.0001)G4cout<<"*Warn*G4QElCS::GetExT:1p "<<ts1<<"#"<<lastTM<<",d="<<ds1<<G4endl;
#endif
    G4double E2=lastTM*theB2;
  		G4double R2=(1.-std::exp(-E2*E2*E2));
#ifdef ppdebug
    G4double ts2=std::pow(-std::log(1.-R2),.333333333)/theB2;
    G4double ds2=std::fabs(ts2-lastTM)/lastTM;
    if(ds2>.0001)G4cout<<"*Warn*G4QElCS::GetExT:2p "<<ts2<<"#"<<lastTM<<",d="<<ds2<<G4endl;
#endif
    G4double E3=lastTM*theB3;
  		G4double R3=(1.-std::exp(-E3));
#ifdef ppdebug
    G4double ts3=-std::log(1.-R3)/theB3;
    G4double ds3=std::fabs(ts3-lastTM)/lastTM;
    if(ds3>.0001)G4cout<<"*Warn*G4QElCS::GetExT:3p "<<ts3<<"#"<<lastTM<<",d="<<ds3<<G4endl;
#endif
  		G4double I1=R1*theS1/theB1;
  		G4double I2=R2*theS2;
				G4double I3=R3*theS3;
    G4double I12=I1+I2;
    G4double rand=(I12+I3)*G4UniformRand();
    if     (rand<I1 ) q2=-std::log(1.-R1*G4UniformRand())/theB1;
				else if(rand<I12) q2=std::pow(-std::log(1.-R2*G4UniformRand()),.33333333333)/theB2;
    else              q2=-std::log(1.-R3*G4UniformRand())/theB3;
  }
  else if((PDG==2212 || PDG==2112) && (tgZ==1 && tgN==1 || tgZ==2 && tgN==2))
  {
#ifdef tdebug
    G4cout<<"G4QElasticCS::GetExchangeT: TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1<<",S2="
          <<theS2<<",B2="<<theB2<<",S3="<<theS3<<",B3="<<theB3<<",GeV2="<<GeVSQ<<G4endl;
#endif
    G4double E1=lastTM*(theB1+lastTM*theSS);
  		G4double R1=(1.-std::exp(-E1));
    G4double tss=theSS+theSS;
#ifdef tdebug
    G4double ts1=-std::log(1.-R1)/theB1;
    if(std::fabs(tss)>1.e-7)ts1=(std::sqrt(theB1*(theB1+(tss+tss)*ts1))-theB1)/tss;
    G4double ds1=(ts1-lastTM)/lastTM;
    if(ds1>.0001)G4cout<<"*Warn*G4QElCS::GetExT:1a "<<ts1<<"#"<<lastTM<<",d="<<ds1<<G4endl;
#endif
    G4double E2=lastTM*theB2;
  		G4double R2=(1.-std::exp(-E2));
#ifdef tdebug
    G4double ts2=-std::log(1.-R2)/theB2;
    G4double ds2=std::fabs(ts2-lastTM)/lastTM;
    if(ds2>.0001)G4cout<<"*Warn*G4QElCS::GetExT:2a "<<ts2<<"#"<<lastTM<<",d="<<ds2<<G4endl;
#endif
    G4double E3=lastTM*theB3;
  		G4double R3=(1.-std::exp(-E3));
#ifdef tdebug
    G4double ts3=-std::log(1.-R3)/theB3;
    G4double ds3=std::fabs(ts3-lastTM)/lastTM;
    if(ds3>.0001)G4cout<<"*Warn*G4QElCS::GetExT:3a "<<ts3<<"#"<<lastTM<<",d="<<ds3<<G4endl;
#endif
  		G4double I1=R1*theS1;
  		G4double I2=R2*theS2/theB2;
				G4double I3=R3*theS3/theB3;
    G4double I12=I1+I2;
    G4double rand=(I12+I3)*G4UniformRand();
#ifdef tdebug
    G4cout<<"G4QElCS::GtExT:1="<<I1<<",2="<<I2<<",3="<<I3<<",R3="<<R3<<",r="<<rand<<G4endl;
#endif
    if(rand<I1)
    {
      q2=(std::sqrt(theB1*theB1-(tss+tss)*std::log(1.-R1*G4UniformRand()))-theB1)/tss;
#ifdef tdebug
      G4cout<<"G4QElCS::GetExT:Q2="<<q2<<",ss="<<tss/2<<",b1="<<theB1<<",t1="<<ts1<<G4endl;
#endif
    }
				else if(rand<I12)
    {
      q2=-std::log(1.-R2*G4UniformRand())/theB2;
#ifdef tdebug
      G4cout<<"G4QElCS::GetExT: Q2="<<q2<<", r2="<<R2<<", b2="<<theB2<<",t2="<<ts2<<G4endl;
#endif
    }
    else
    {
      q2=lastTM+std::log(1.-R3*G4UniformRand())/theB3;
#ifdef tdebug
      G4cout<<"G4QElCS::GetExT:Q2="<<q2<<",m="<<lastTM<<",b3="<<theB3<<",t3="<<ts3<<G4endl;
#endif
    }
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetExchangeT: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212, Z=1,N=(0,1) & Z=2,N=2"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetExchangeT: n-p,p-H/He are implemented");
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
  G4double p=std::exp(lp);              // momentum
  G4double sp=std::sqrt(p);             // sqrt(p)
  G4double p2=p*p;            
  G4double p4=p2*p2;
  G4double p3=p2*p;

  if(PDG==2112 && tgZ==1 && tgN==0)
  {
    G4double ssp=std::sqrt(sp);           // sqrt(sqrt(p))=p^.25
    G4double p2s=p2*sp;
		  G4double dl1=lp-lastPAR[3];
    theSS=lastPAR[21];
    theS1=(lastPAR[9]+lastPAR[10]*dl1*dl1+lastPAR[11]/p)/(1.+lastPAR[12]/p4)
          +lastPAR[13]/(p4+lastPAR[14]);
    theB1=(lastPAR[17]+lastPAR[18]/(p4*p4+lastPAR[19]*p3))/(1.+lastPAR[20]/p4);
    theS2=(lastPAR[15]+lastPAR[16]/p4/p)/p3;
    theB2=lastPAR[22]/(p*sp+lastPAR[23]); 
    theS3=0.;
    theB3=0.; 
#ifdef tdebug
    G4cout<<"G4QElasticCS::GetTableValues:(pp) TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1
          <<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS1<<",B3="<<theB1<<G4endl;
#endif
    // Returns the total elastic pp cross-section (to avoid spoiling lastSIG)
    return lastPAR[0]/(p2s+lastPAR[1]*p+lastPAR[2]/ssp)+lastPAR[4]/p
           +(lastPAR[5]+lastPAR[6]*dl1*dl1+lastPAR[7]/p)/(1.+lastPAR[8]/p4);
  }
  else if(PDG==2212 && tgZ==1 && tgN==0)
  {
    G4double p2s=p2*sp;
		  G4double dl1=lp-lastPAR[3];
		  G4double dl2=lp-lastPAR[8];
    theSS=lastPAR[31];
    theS1=(lastPAR[9]+lastPAR[10]*dl2*dl2)/(1.+lastPAR[11]/p4/p)+
          (lastPAR[12]/p2+lastPAR[13]*p)/(p4+lastPAR[14]*sp);
    theB1=lastPAR[15]*std::pow(p,lastPAR[16])/(1.+lastPAR[17]/p3);
    theS2=lastPAR[18]+lastPAR[19]/(p4+lastPAR[20]*p);
    theB2=lastPAR[21]+lastPAR[22]/(p4+lastPAR[23]/sp); 
    theS3=lastPAR[24]+lastPAR[25]/(p4*p4+lastPAR[26]*p2+lastPAR[27]);
    theB3=lastPAR[28]+lastPAR[29]/(p4+lastPAR[30]); 
#ifdef tdebug
    G4cout<<"G4QElasticCS::GetTableValues:(pp) TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1
          <<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS1<<",B3="<<theB1<<G4endl;
#endif
    // Returns the total elastic pp cross-section (to avoid spoiling lastSIG)
    return lastPAR[0]/p2s/(1.+lastPAR[7]/p2s)+(lastPAR[1]+lastPAR[2]*dl1*dl1+lastPAR[4]/p)
                                                   /(1.+lastPAR[5]*lp)/(1.+lastPAR[6]/p4);
  }
  else if((PDG==2212 || PDG==2112) && tgZ==1 && tgN==1) // n/p+d (isotope invariant)
  {
    G4double psp=p*sp;                      // p*sqrt(p)
		  G4double dl1=lp-lastPAR[2];
		  G4double dl2=lp-lastPAR[11];
    theSS=lastPAR[17]/p4-lastPAR[18]/p;
    theS1=(lastPAR[12]+lastPAR[13]*dl2*dl2)*(1+lastPAR[14]/(p4+lastPAR[15]/p+lastPAR[16]));
    theB1=lastPAR[7]/p2+lastPAR[8]/p4+lastPAR[9]/(1.+lastPAR[10]/p);
    theS2=lastPAR[19]+lastPAR[20]*lp+(1.+lastPAR[21]/p2)/p2;
    theB2=(lastPAR[22]+lastPAR[23]/psp)/(1.+lastPAR[24]/psp); 
    theS3=lastPAR[25]/((p3+lastPAR[26])*p3+lastPAR[27]);
    theB3=lastPAR[28]/(p3+lastPAR[29]*sp); 
#ifdef tdebug
    G4cout<<"G4QElasticCS::GetTableValues:(pd) p0="<<lastPAR[0]<<",p1="<<lastPAR[1]<<",p2="
          <<lastPAR[2]<<",p3="<<lastPAR[3]<<",p4="<<lastPAR[4]<<",p5="<<lastPAR[5]<<G4endl;
    G4cout<<"G4QElasticCS::GetTableValues:(pd) TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1<<
         ",SS="<<theSS<<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS3<<",B3="<<theB3<<G4endl;
#endif
    // Returns the total elastic pd cross-section (to avoid spoiling lastSIG)
    return (lastPAR[0]+lastPAR[1]*dl1*dl1)*(1.+lastPAR[3]/(p4+lastPAR[4]*sp)+lastPAR[5]/p)
           /(1.+lastPAR[6]/p3);
  }
  else if((PDG==2212 || PDG==2112) && tgZ==2 && tgN==2) // n/p+He4 (isotope invariant)
  {
    G4double p3sp=p3*sp;                      // p^3*sqrt(p)
    G4double p5=p4*p;
    G4double p6=p5*p;
		  G4double dl1=lp-lastPAR[2];
		  G4double dl2=p-lastPAR[8];
		  G4double dl3=p-lastPAR[15];
    theSS=lastPAR[20]/p4;
    theS1=(lastPAR[9]+lastPAR[10]*dl1*dl1)*(1.+lastPAR[11]/(p4+lastPAR[12]*p))-
                                                   lastPAR[13]/(1.+lastPAR[14]*dl3*dl3)/p2;
    theB1=lastPAR[16]/(p3sp+lastPAR[17]/p2)+lastPAR[18]+lastPAR[19]*lp;
    theS2=lastPAR[21]+lastPAR[22]/p+lastPAR[23]/(p4+lastPAR[24])/p2;
    theB2=lastPAR[25]; 
    theS3=lastPAR[26]/((p6+lastPAR[27])*p5+lastPAR[28]/p2);
    theB3=lastPAR[29]/(p2+lastPAR[30])+lastPAR[31]; 
#ifdef tdebug
    G4cout<<"G4QElasticCS::GetTableValues:(pa) p0="<<lastPAR[0]<<",p1="<<lastPAR[1]<<",p2="
          <<lastPAR[2]<<",p3="<<lastPAR[3]<<",p4="<<lastPAR[4]<<",p5="<<lastPAR[5]<<G4endl;
    G4cout<<"G4QElasticCS::GetTableValues:(pd) TM="<<lastTM<<",S1="<<theS1<<",B1="<<theB1<<
         ",SS="<<theSS<<",S2="<<theS2<<",B2="<<theB2<<",S3="<<theS3<<",B3="<<theB3<<G4endl;
#endif
    // Returns the total elastic pd cross-section (to avoid spoiling lastSIG)
    return (lastPAR[0]+lastPAR[1]*dl1*dl1)*(1.+lastPAR[3]/(p4+lastPAR[4]*p+lastPAR[5]))-
                                                     lastPAR[6]/(1.+lastPAR[7]*dl2*dl2)/p2;
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
  static const G4double mNeut= G4QPDGCode(2112).GetMass()*.001; // MeV to GeV
  static const G4double mProt= G4QPDGCode(2212).GetMass()*.001; // MeV to GeV
  //static const G4double mLamb= G4QPDGCode(3122).GetMass()*.001; // MeV to GeV
  //static const G4double mHe3 = G4QPDGCode(2112).GetNuclMass(2,1,0)*.001; // MeV to GeV
  //static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0)*.001; // MeV to GeV
  static const G4double mDeut = G4QPDGCode(2112).GetNuclMass(1,1,0)*.001; // MeV to GeV
  static const G4double mProt2= mProt*mProt;
  static const G4double mNeut2= mNeut*mNeut;
  static const G4double mDeut2= mDeut*mDeut;
  if(PDG==2212 && tgZ==1 && tgN==0)
  {
    G4double tMid=std::sqrt(pP*pP+mProt2)*mProt-mProt2; // CMS 90deg value of -t=Q2 (GeV^2)
    return tMid+tMid;
  }
  else if(PDG==2112 && tgZ==1&&tgN==0)
  {
    G4double pP2=pP*pP;
    G4double sM=(mProt+mProt)*std::sqrt(pP*pP+mNeut2)+mNeut2+mProt2; // Mondelstam s
    return 4*mProt2*pP2/sM;
  }
  else if((PDG==2212 || PDG==2112) && tgZ==1&&tgN==1)                // n/p+He4
  {
    G4double pP2=pP*pP;
    G4double sM=(mDeut+mDeut)*std::sqrt(pP*pP+mProt2)+mProt2+mDeut2; // Mondelstam s
    return 4*mDeut2*pP2/sM;
  }
  else if((PDG==2212 || PDG==2112) && tgZ==2&&tgN==2)                // n/p+He4
  {
    G4double pP2=pP*pP;
    G4double mNuc=G4QPDGCode(2112).GetNuclMass(tgZ,tgN,0)*.001;      // MeV to GeV
    G4double mNuc2=mNuc*mNuc;
    G4double sM=(mNuc+mNuc)*std::sqrt(pP*pP+mProt2)+mProt2+mNuc2;    // Mondelstam s
    return 4*mNuc2*pP2/sM;
  }
  else
  {
    G4cout<<"*Error*G4QElasticCrossSection::GetQ2max: PDG="<<PDG<<", Z="<<tgZ<<", N="
          <<tgN<<", while it is defined only for PDG=2212/2112, Z=1, N=0"<<G4endl;
    throw G4QException("G4QElasticCrossSection::GetQ2max: np,pp,pd,pHe are implemented");
  }
}
