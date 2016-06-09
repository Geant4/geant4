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
// G4 Physics class: G4QNeutronCaptureRatio for N+A Diffraction Interactions
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Oct-06
//
// ----------------------------------------------------------------------
// Short description: (n,gamma) capture is a part of the incoherent
// (inelastic) interaction. This part is calculated in the class.
// ----------------------------------------------------------------------

//#define debug
//#define pdebug
//#define fdebug
//#define nandebug

#include "G4QNeutronCaptureRatio.hh"
#include "G4SystemOfUnits.hh"

// Returns Pointer to the  class G4QNeutronCaptureRatio
G4QNeutronCaptureRatio* G4QNeutronCaptureRatio::GetPointer()
{
  static G4QNeutronCaptureRatio theRatios; // *** Static body of the NeutCaptureRatio ***
  return &theRatios;
}

// Calculation of netronCapture/inelastic ratio
G4double G4QNeutronCaptureRatio::GetRatio(G4double pIU, G4int tgZ, G4int tgN)
{
  // Table parameters
  static const G4int    npp=100;         // Number of steps in the R(s) LinTable
  static const G4int    mpp=npp+1;       // Number of elements in the R(s) LinTable
  static const G4double pma=6.;          // The first LinTabEl(mom=0)=1., mom>pma -> logTab
  static const G4double dp=pma/npp;      // Step of the linear Table
  static const G4int    nls=150;         // Number of steps in the R(lns) logTable
  static const G4int    mls=nls+1;       // Number of elements in the R(lns) logTable
  static const G4double lpi=1.79;        // The min ln(p) logTabEl(p=5.99 < pma=6.)
  static const G4double lpa=8.;          // The max ln(p) logTabEl(p=5.99 - 2981 GeV)
  static const G4double mi=std::exp(lpi);// The min mom of logTabEl(~ 5.99 GeV)
  static const G4double max_s=std::exp(lpa);// The max mom of logTabEl(~ 2981 GeV)
  static const G4double dl=(lpa-lpi)/nls;// Step of the logarithmic Table
  static const G4double edl=std::exp(dl);// Multiplication step of the logarithmic Table
  static const G4double toler=.0001;     // Tolarence (GeV) defining the same momentum
  static G4double lastP=0.;              // Last mometum value for which R was calculated
  static G4double lastR=0.;              // Last ratio R which was calculated
  // Local Associative Data Base:
  static std::vector<G4int>     vZ;      // Vector of calculated Z (target)
  static std::vector<G4int>     vN;      // Vector of calculated N (target)
  static std::vector<G4double>  vH;      // Vector of max mom initialized in the LinTable
  static std::vector<G4int>     vJ;      // Vector of topBin number initialized in LinTable
  static std::vector<G4double>  vM;      // Vector of relMax ln(p) initialized in LogTable
  static std::vector<G4int>     vK;      // Vector of topBin number initialized in LogTable
  static std::vector<G4double*> vT;      // Vector of pointers to LinTable in C++ heap
  static std::vector<G4double*> vL;      // Vector of pointers to LogTable in C++ heap
  // Last values of the Associative Data Base:
  static G4int     lastZ=0;              // theLast of calculated A
  static G4int     lastN=0;              // theLast of calculated A
  static G4double  lastH=0.;             // theLast of max mom initialized in the LinTable
  static G4int     lastJ=0;              // theLast of topBinNumber initialized in LinTable
  static G4double  lastM=0.;             // theLast of relMax ln(p) initialized in LogTab.
  static G4int     lastK=0;              // theLast of topBinNumber initialized in LogTable
  static G4double* lastT=0;              // theLast of pointer to LinTable in the C++ heap
  static G4double* lastL=0;              // theLast of pointer to LogTable in the C++ heap
  // LogTable is created only if necessary. R(p>2981GeV) calcul by formula for any nuclei
  G4int A=tgN+tgZ;
  if(pIU > 50) return 0.;
  if(pIU > 30 && ((tgN==1 && tgZ==1) || (tgN==8 && tgZ==7))) return 0.;
  if(pIU > 20 && tgN==2 && tgZ==1) return 0.;
  if(pIU > 15 && ((tgN==1 && tgZ==2) || (tgN==8 && tgZ==8))) return 0.;
  if(pIU<toler || A<1) return 1.;        // Fake use of toler as non zero number
  if(A>247)
  {
    G4cout<<"-*-Warning-*-G4NeutronCaptureRatio::GetRatio:A="<<A<<">247, return 0"<<G4endl;
    return 0.;
  }
  G4int nDB=vZ.size();                   // A number of nuclei already initialized in AMDB
  if(nDB && lastZ==tgZ && lastN==tgN && std::fabs(pIU-lastP)<toler) return lastR;
  if(pIU>max_s)
  {
    lastR=CalcCap2In_Ratio(s,tgZ,tgN);   //@@ Probably user ought to be notified about bigP
    return lastR;
  }
  G4bool found=false;
  G4int i=-1;
  if(nDB) for (i=0; i<nDB; i++) if(tgZ==vZ[i] && tgN==vN[i]) // Sirch for this Z,N in AMDB
  {
    found=true;                          // The (Z,N) is found
    break;
  }
  if(!nDB || !found)                     // Create new line in the AMDB
  {
    lastZ = tgZ;
    lastN = tgN;
    lastT = new G4double[mpp];           // Create the linear Table
    lastJ = static_cast<int>(pIU/dp)+1;  // MaxBin to be initialized
    if(lastJ>npp)
    {
      lastJ=npp;
      lastH=pma;
    }
    else lastH = lastJ*dp;               // Calculate max initialized s for LinTab
    G4double pv=0;
    lastT[0]=1.;
    for(G4int j=1; j<=lastJ; j++)        // Calculate LinTab values
    {
      pv+=dp;
      lastT[j]=CalcCap2In_Ratio(pv,tgZ,tgN); // ??
    }
    lastL=new G4double[mls];           // Create the logarithmic Table
    G4double ls=std::log(s);
    lastK = static_cast<int>((ls-lpi)/dl)+1; // MaxBin to be initialized in LogTaB
    if(lastK>nls)
    {
      lastK=nls;
      lastM=lpa-lpi;
    }
    else lastM = lastK*dl;             // Calculate max initialized ln(s)-lpi for LogTab
    pv=mi;
    for(G4int j=0; j<=lastK; j++)      // Calculate LogTab values
    {
      lastL[j]=CalcCap2In_Ratio(pv,tgZ,tgN);
      if(j!=lastK) pv*=edl;
    }
    i++;                                 // Make a new record to AMDB and position on it
    vZ.push_back(lastZ);
    vN.push_back(lastN);
    vH.push_back(lastH);
    vJ.push_back(lastJ);
    vM.push_back(lastM);
    vK.push_back(lastK);
    vT.push_back(lastT);
    vL.push_back(lastL);
  }
  else                                   // The A value was found in AMDB
  {
    lastZ=vZ[i];
    lastN=vN[i];
    lastH=vH[i];
    lastJ=vJ[i];
    lastM=vM[i];
    lastK=vK[i];
    lastT=vT[i];
    lastL=vL[i];
    if(s>lastH)                          // At least LinTab must be updated
    {
      G4int nextN=lastJ+1;               // The next bin to be initialized
      if(lastJ<npp)
      {
        lastJ = static_cast<int>(pIU/dp)+1;// MaxBin to be initialized
        G4double pv=lastH;
        if(lastJ>npp)
        {
          lastJ=npp;
          lastH=pma;
        }
        else lastH = lastJ*dp;           // Calculate max initialized s for LinTab
        for(G4int j=nextN; j<=lastJ; j++)// Calculate LogTab values
        {
          pv+=dp;
          lastT[j]=CalcCap2In_Ratio(pv,tgZ,tgN);
        }
      } // End of LinTab update
      if(lastJ>=nextN)
      {
        vH[i]=lastH;
        vJ[i]=lastJ;
      }
      G4int nextK=lastK+1;
      if(pIU>pma && lastK<nls)           // LogTab must be updated
      {
        G4double pv=std::exp(lastM+lpi); // Define starting poit (lastM will be changed)
        G4double ls=std::log(s);
        lastK = static_cast<int>((ls-lpi)/dl)+1; // MaxBin to be initialized in LogTaB
        if(lastK>nls)
        {
          lastK=nls;
          lastM=lpa-lpi;
        }
        else lastM = lastK*dl;           // Calculate max initialized ln(p)-lpi for LogTab
        for(G4int j=nextK; j<=lastK; j++)// Calculate LogTab values
        {
          pv*=edl;
          lastL[j]=CalcCap2In_Ratio(pv,tgZ,tgN);
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
  if(pIU<pma)                           // Use linear table
  {
    G4int n=static_cast<int>(pIU/dp);   // Low edge number of the bin
    G4double d=s-n*dp;                  // Linear shift
    G4double v=lastT[n];                // Base
    lastR=v+d*(lastT[n+1]-v)/dp;        // Result
  }
  else                                  // Use log table
  {
    G4double ls=std::log(pIU)-lpi;      // ln(p)-l_min
    G4int n=static_cast<int>(ls/dl);    // Low edge number of the bin
    G4double d=ls-n*dl;                 // Log shift
    G4double v=lastL[n];                // Base
    lastR=v+d*(lastL[n+1]-v)/dl;        // Result
  }
  if(lastR<0.) lastR=0.;
  if(lastR>1.) lastR=1.;
  return lastR;
} // End of GetRatio

// Calculate Capture/Inelastic Ratio as a function of total momentum (in GeV/c)
G4double G4QNeutronCaptureRatio::CalcCap2In_Ratio(G4double p, G4int Z, G4int N)
{
  //==> n (Z=0)
  static const G4int N0=1;
  static const G4double pZ0N1[5]={.0001, 40., 0., 0., 1.}; // *** No Capture ?
  static const std::pair<G4int, const G4double*> Z0N1(1,pZ0N1);
  static const std::pair<G4int, const G4double*> Z0[N0]={Z0N1};
  //==> H (Z=1) *** no protons, which are treated separately ***
  static const G4int N1=2;
  static const G4double pZ1N1[5]={.07, 40., .0006, .0001, .01};
  static const std::pair<G4int, const G4double*> Z1N1(1,pZ1N1);
  static const G4double pZ1N2[5]={.0001, 40., 0., 0., 1.}; // *** No Capture ?
  static const std::pair<G4int, const G4double*> Z1N2(2,pZ1N2);
  static const std::pair<G4int, const G4double*> Z1[N1]={Z1N1, Z1N2};
  //==> He(Z=2)
  static const G4int N2=2;
  static const G4double pZ2N1[5]={.1, 40., 0., 0., 1.};   // *** Unknown threshold
  static const std::pair<G4int, const G4double*> Z2N1(1,pZ2N1);
  static const G4double pZ2N2[5]={.0001, 40., 0., 0., 1.}; // *** No Capture ?
  static const std::pair<G4int, const G4double*> Z2N2(2,pZ2N2);
  static const std::pair<G4int, const G4double*> Z2[N2]={Z2N1, Z2N2};
  //==> Li(Z=3)
  static const G4int N3=2;
  static const G4double pZ3N3[5]={.001, 40., 3.E-5, .0001, .1};
  static const std::pair<G4int, const G4double*> Z3N1(3,pZ3N3);
  static const G4double pZ3N4[5]={.022, 19., 3.E-5, .0001, .04};
  static const std::pair<G4int, const G4double*> Z3N2(4,pZ3N4);
  static const std::pair<G4int, const G4double*> Z3[N3]={Z3N1, Z3N2};
  //==> Be(Z=4)
  static const G4int N4=1;
  static const G4double pZ4N5[5]={.0004, 40., 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z4N5(5,pZ4N5);
  static const std::pair<G4int, const G4double*> Z4[N4]={Z4N5};
  //==> B (Z=5)
  static const G4int N5=2;
  static const G4double pZ5N5[5]={.011, 9., .0002, .0001, .002};
  static const std::pair<G4int, const G4double*> Z5N5(5,pZ5N5);
  static const G4double pZ5N6[5]={.027, 9., 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z5N6(6,pZ5N6);
  static const std::pair<G4int, const G4double*> Z5[N5]={Z5N5, Z5N6};
  //==> C (Z=6)
  static const G4int N6=2;
  static const G4double pZ6N6[5]={.08, 40., .0003, .0001, .07}; // *** Only Nat Mix ***
  static const std::pair<G4int, const G4double*> Z6N6(6,pZ6N6);
  static const G4double pZ6N7[5]={.08, 40., .0003, .0001, .07}; // *** Only Nat Mix ***
  static const std::pair<G4int, const G4double*> Z6N7(7,pZ6N7);
  static const std::pair<G4int, const G4double*> Z6[N6]={Z6N6, Z6N7};
  //==> N (Z=7)
  static const G4int N7=2;
  static const G4double pZ7N7[5]={.005, 3., 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z7N7(7,pZ7N7);
  static const G4double pZ7N8[5]={.084, 40., .0001, .0001, .015};
  static const std::pair<G4int, const G4double*> Z7N8(8,pZ7N8);
  static const std::pair<G4int, const G4double*> Z7[N7]={Z7N7, Z7N8};
  //==> O (Z=8)
  static const G4int N8=3;
  static const G4double pZ8N8[5]={.08, 40., 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z8N8(8,pZ8N8);
  static const G4double pZ8N9[5]={.0065, 5., .0013, .0001, .02};
  static const std::pair<G4int, const G4double*> Z8N9(9,pZ8N9);
  static const G4double pZ8N10[5]={.01, 27., 0., 0., 1.}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z8N10(10,pZ8N10);
  static const std::pair<G4int, const G4double*> Z8[N8]={Z8N8, Z8N9, Z8N10};
  //==> F (Z=9)
  static const G4int N9=1;
  static const G4double pZ9N10[5]={.013, 27., .0001, .0001, .02};
  static const std::pair<G4int, const G4double*> Z9N10(10,pZ9N10);
  static const std::pair<G4int, const G4double*> Z9[N9]={Z9N10};
  //==> Ne(Z=10)
  static const G4int N10=3;
  static const G4double pZ10N10[5]={.01, 27., 0., 0., 1.}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z10N10(10,pZ10N10);
  static const G4double pZ10N11[5]={.01, 27., 0., 0., 1.}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z10N11(11,pZ10N11);
  static const G4double pZ10N12[5]={.01, 27., 0., 0., 1.}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z10N12(12,pZ10N12);
  static const std::pair<G4int, const G4double*> Z10[N10]={Z10N10, Z10N11, Z10N12};
  //==> Na(Z=11)
  static const G4int N11=1;
  static const G4double pZ11N12[5]={.024, 17., .0005, .0001, .03};
  static const std::pair<G4int, const G4double*> Z11N12(12,pZ11N12);
  static const std::pair<G4int, const G4double*> Z11[N11]={Z11N12};
  //==> Mg(Z=12)
  static const G4int N12=3;
  static const G4double pZ12N12[5]={.045, 40., .0003, .0001, .02};
  static const std::pair<G4int, const G4double*> Z12N12(12,pZ12N12);
  static const G4double pZ12N13[5]={.019, 7., .0002, .0001, .01};
  static const std::pair<G4int, const G4double*> Z12N13(13,pZ12N13);
  static const G4double pZ12N14[5]={.053, 40., .0006, .0001, .007};
  static const std::pair<G4int, const G4double*> Z12N14(14,pZ12N14);
  static const std::pair<G4int, const G4double*> Z12[N12]={Z12N12, Z12N13, Z12N14};
  //==> Al(Z=13)
  static const G4int N13=1;
  static const G4double pZ13N14[5]={.035, 17., .001, .03, 1.};
  static const std::pair<G4int, const G4double*> Z13N14(14,pZ13N14);
  static const std::pair<G4int, const G4double*> Z13[N13]={Z13N14};
  //==> Si(Z=14)
  static const G4int N14=3;
  static const G4double pZ14N14[5]={.052, 40., .002, .0001, .008};
  static const std::pair<G4int, const G4double*> Z14N14(14,pZ14N14);
  static const G4double pZ14N15[5]={.048, 40., .0001, .0004, .02};
  static const std::pair<G4int, const G4double*> Z14N15(15,pZ14N15);
  static const G4double pZ14N16[5]={.06, 40., .0015, .0001, .01};
  static const std::pair<G4int, const G4double*> Z14N16(16,pZ14N16);
  static const std::pair<G4int, const G4double*> Z14[N14]={Z14N14, Z14N15, Z14N16};
  //==> P (Z=15)
  static const G4int N15=1;
  static const G4double pZ15N16[5]={.024, 7., .0008, .0001, .03};
  static const std::pair<G4int, const G4double*> Z15N16(16,pZ15N16);
  static const std::pair<G4int, const G4double*> Z15[N15]={Z15N16};
  //==> S (Z=16)
  static const G4int N16=4;
  static const G4double pZ16N16[5]={.036, 12., .0003, .03, .004};
  static const std::pair<G4int, const G4double*> Z16N16(16,pZ16N16);
  static const G4double pZ16N17[5]={.018, 40., .0033, .0001, .002};
  static const std::pair<G4int, const G4double*> Z16N17(17,pZ16N17);
  static const G4double pZ16N18[5]={.053, 25., .002, .0001, .0043};
  static const std::pair<G4int, const G4double*> Z16N18(18,pZ16N18);
  static const G4double pZ16N20[5]={.065, 25., .002, .0001, .0043};
  static const std::pair<G4int, const G4double*> Z16N20(20,pZ16N20);
  static const std::pair<G4int, const G4double*> Z16[N16]={Z16N16, Z16N17, Z16N18, Z16N20};
  //==> Cl(Z=17)
  static const G4int N17=2;
  static const G4double pZ17N18[5]={.014, 4., .0004, .175, .002};
  static const std::pair<G4int, const G4double*> Z17N18(18,pZ17N18);
  static const G4double pZ17N20[5]={.035, 8., .008, 18., .0005};
  static const std::pair<G4int, const G4double*> Z17N20(20,pZ17N20);
  static const std::pair<G4int, const G4double*> Z17[N17]={Z17N18, Z17N20};
  //==> Ar(Z=18)
  static const G4int N18=3;
  static const G4double pZ18N18[5]={.036, 8., .0005, .1, .01};
  static const std::pair<G4int, const G4double*> Z18N18(18,pZ18N18);
  static const G4double pZ18N20[5]={.025, 6., .0027, .19, .0003};
  static const std::pair<G4int, const G4double*> Z18N20(20,pZ18N20);
  static const G4double pZ18N22[5]={.028, 6., .001, .19, .0003};
  static const std::pair<G4int, const G4double*> Z18N22(22,pZ18N22);
  static const std::pair<G4int, const G4double*> Z18[N18]={Z18N18, Z18N20, Z18N22};
  //==> K (Z=19)
  static const G4int N19=3;
  static const G4double pZ19N20[5]={.04, 8., 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z19N20(20,pZ19N20);
  static const G4double pZ19N21[5]={.049, 5., 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z19N21(21,pZ19N21);
  static const G4double pZ19N22[5]={.04, 11., .005, .0001, .005};
  static const std::pair<G4int, const G4double*> Z19N22(22,pZ19N22);
  static const std::pair<G4int, const G4double*> Z19[N19]={Z19N20, Z19N21, Z19N22};
  //==> Ca(Z=20)
  static const G4int N20=6;
  static const G4double pZ20N20[5]={.05, 14., .0006, .09, .009};
  static const std::pair<G4int, const G4double*> Z20N20(20,pZ20N20);
  static const G4double pZ20N22[5]={.047, 30., .003, .0001, .014};
  static const std::pair<G4int, const G4double*> Z20N22(22,pZ20N22);
  static const G4double pZ20N23[5]={.01, 3.5, .0015, .16, .002};
  static const std::pair<G4int, const G4double*> Z20N23(23,pZ20N23);
  static const G4double pZ20N24[5]={.04, 30., .002, .0001, .008};
  static const std::pair<G4int, const G4double*> Z20N24(24,pZ20N24);
  static const G4double pZ20N26[5]={.044, 40., .001, .0001, .01};
  static const std::pair<G4int, const G4double*> Z20N26(26,pZ20N26);
  static const G4double pZ20N28[5]={.055, 14., .001, .18, .001};
  static const std::pair<G4int, const G4double*> Z20N28(28,pZ20N28);
  static const std::pair<G4int, const G4double*> Z20[N20]={Z20N20, Z20N22, Z20N23,
                                                           Z20N24, Z20N26, Z20N28};
  //==> Sc(Z=21)
  static const G4int N21=1;
  static const G4double pZ21N24[5]={.014, 4., 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z21N24(24,pZ21N24);
  static const std::pair<G4int, const G4double*> Z21[N21]={Z21N24};
  //==> Ti(Z=22)
  static const G4int N22=5;
  static const G4double pZ22N24[5]={.036, 27., .007, .0001, .005};
  static const std::pair<G4int, const G4double*> Z22N24(24,pZ22N24);
  static const G4double pZ22N25[5]={.013, 9., .017, .0001, .005};
  static const std::pair<G4int, const G4double*> Z22N25(25,pZ22N25);
  static const G4double pZ22N26[5]={.043, 40., .002, .0001, .01};
  static const std::pair<G4int, const G4double*> Z22N26(26,pZ22N26);
  static const G4double pZ22N27[5]={.047, 30., .007, .0001, .01};
  static const std::pair<G4int, const G4double*> Z22N27(27,pZ22N27);
  static const G4double pZ22N28[5]={.052, 40., .0005, .0001, .01};
  static const std::pair<G4int, const G4double*> Z22N28(28,pZ22N28);
  static const std::pair<G4int, const G4double*> Z22[N22]={Z22N24, Z22N25, Z22N26,
                                                         Z22N27, Z22N28};
  //==> V (Z=23)
  static const G4int N23=2;
  static const G4double pZ23N27[5]={.023, 30., .01, .0001, .003}; // *** Only Nat mix ***
  static const std::pair<G4int, const G4double*> Z23N27(27,pZ23N27);
  static const G4double pZ23N28[5]={.023, 30., .01, .0001, .003}; // *** Only Nat mix ***
  static const std::pair<G4int, const G4double*> Z23N28(28,pZ23N28);
  static const std::pair<G4int, const G4double*> Z23[N23]={Z23N27, Z23N28};
  //==> Cr(Z=24)
  static const G4int N24=4;
  static const G4double pZ24N26[5]={.035, 27., .004, .0001, .01};
  static const std::pair<G4int, const G4double*> Z24N26(26,pZ24N26);
  static const G4double pZ24N28[5]={.049, 40., .001, .0001, .016};
  static const std::pair<G4int, const G4double*> Z24N28(28,pZ24N28);
  static const G4double pZ24N29[5]={.032, 30., .005, .0001, .005};
  static const std::pair<G4int, const G4double*> Z24N29(29,pZ24N29);
  static const G4double pZ24N30[5]={.034, 30., .002, .0001, .008};
  static const std::pair<G4int, const G4double*> Z24N30(30,pZ24N30);
  static const std::pair<G4int, const G4double*> Z24[N24]={Z24N26, Z24N28, Z24N29, Z24N30};
  //==> Mn(Z=25)
  static const G4int N25=1;
  static const G4double pZ25N30[5]={.0145, 10., .01, .0001, .003};
  static const std::pair<G4int, const G4double*> Z25N30(30,pZ25N30);
  static const std::pair<G4int, const G4double*> Z25[N25]={Z25N30};
  //==> Fe(Z=26)
  static const G4int N26=4;
  static const G4double pZ26N28[5]={.048, 27., .0016, .0001, .016};
  static const std::pair<G4int, const G4double*> Z26N28(28,pZ26N28);
  static const G4double pZ26N30[5]={.036, 27., .004, .0001, .006};
  static const std::pair<G4int, const G4double*> Z26N30(30,pZ26N30);
  static const G4double pZ26N31[5]={.036, 27., .005, .0001, .005};
  static const std::pair<G4int, const G4double*> Z26N31(31,pZ26N31);
  static const G4double pZ26N32[5]={.036, 27., .005, .0001, .007};
  static const std::pair<G4int, const G4double*> Z26N32(32,pZ26N32);
  static const std::pair<G4int, const G4double*> Z26[N26]={Z26N28, Z26N30, Z26N31, Z26N32};
  //==> Co(Z=27)
  static const G4int N27=1;
  static const G4double pZ27N32[5]={.044, 22., .002, .0001, .016};
  static const std::pair<G4int, const G4double*> Z27N32(32,pZ27N32);
  static const std::pair<G4int, const G4double*> Z27[N27]={Z27N32};
  //==> Ni(Z=28)
  static const G4int N28=5;
  static const G4double pZ28N30[5]={.045, 20., .003, .0001, .01};
  static const std::pair<G4int, const G4double*> Z28N30(30,pZ28N30);
  static const G4double pZ28N32[5]={.046, 20., .016, .0001, .01};
  static const std::pair<G4int, const G4double*> Z28N32(32,pZ28N32);
  static const G4double pZ28N33[5]={.046, 20., .016, .0001, .01};
  static const std::pair<G4int, const G4double*> Z28N33(33,pZ28N33);
  static const G4double pZ28N34[5]={.045, 20., .005, .0001, .007};
  static const std::pair<G4int, const G4double*> Z28N34(34,pZ28N34);
  static const G4double pZ28N36[5]={.045, 20., .005, .0001, .007};
  static const std::pair<G4int, const G4double*> Z28N36(36,pZ28N36);
  static const std::pair<G4int, const G4double*> Z28[N28]={Z28N30, Z28N32, Z28N33,
                                                         Z28N34, Z28N36};
  //==> Cu(Z=29)
  static const G4int N29=2;
  static const G4double pZ29N34[5]={.035, 15., .008, .0001, .015};
  static const std::pair<G4int, const G4double*> Z29N34(34,pZ29N34);
  static const G4double pZ29N36[5]={.036, 15., .003, .0001, .013};
  static const std::pair<G4int, const G4double*> Z29N36(36,pZ29N36);
  static const std::pair<G4int, const G4double*> Z29[N29]={Z29N34, Z29N36};
  //==> Zn(Z=30)
  static const G4int N30=5;
  static const G4double pZ30N34[5]={.041, 20., .008, .0001, .02}; // *** only NAT mix ***
  static const std::pair<G4int, const G4double*> Z30N34(34,pZ30N34);
  static const G4double pZ30N36[5]={.041, 20., .008, .0001, .02}; // *** only NAT mix ***
  static const std::pair<G4int, const G4double*> Z30N36(36,pZ30N36);
  static const G4double pZ30N37[5]={.041, 20., .008, .0001, .02}; // *** only NAT mix ***
  static const std::pair<G4int, const G4double*> Z30N37(37,pZ30N37);
  static const G4double pZ30N38[5]={.041, 20., .008, .0001, .02}; // *** only NAT mix ***
  static const std::pair<G4int, const G4double*> Z30N38(38,pZ30N38);
  static const G4double pZ30N40[5]={.041, 20., .008, .0001, .02}; // *** only NAT mix ***
  static const std::pair<G4int, const G4double*> Z30N40(40,pZ30N40);
  static const std::pair<G4int, const G4double*> Z30[N30]={Z30N34, Z30N36, Z30N37,
                                                           Z30N38, Z30N40};
  //==> Ga(Z=31)
  static const G4int N31=2;
  static const G4double pZ31N38[5]={.024, 7., .03, .01, .003};
  static const std::pair<G4int, const G4double*> Z31N38(38,pZ31N38);
  static const G4double pZ31N40[5]={.026, 9., .015, .01, .003};
  static const std::pair<G4int, const G4double*> Z31N40(40,pZ31N40);
  static const std::pair<G4int, const G4double*> Z31[N31]={Z31N38, Z31N40};
  //==> Ge(Z=32)
  static const G4int N32=5;
  static const G4double pZ32N38[5]={.037, 12., .15, .025, .003};
  static const std::pair<G4int, const G4double*> Z32N38(38,pZ32N38);
  static const G4double pZ32N40[5]={.035, 20., .015, .01, .0035};
  static const std::pair<G4int, const G4double*> Z32N40(40,pZ32N40);
  static const G4double pZ32N41[5]={.009, 3., .02, .03, .0001};
  static const std::pair<G4int, const G4double*> Z32N41(41,pZ32N41);
  static const G4double pZ32N42[5]={.027, 12., .003, .0001, .01};
  static const std::pair<G4int, const G4double*> Z32N42(42,pZ32N42);
  static const G4double pZ32N44[5]={.031, 20., .025, .0005, .0045};
  static const std::pair<G4int, const G4double*> Z32N44(44,pZ32N44);
  static const std::pair<G4int, const G4double*> Z32[N32]={Z32N38, Z32N40, Z32N41,
                                                           Z32N42, Z32N44};
  //==> As(Z=33)
  static const G4int N33=1;
  static const G4double pZ33N42[5]={.017, 5., .004, .05, .006};
  static const std::pair<G4int, const G4double*> Z33N42(42,pZ33N42);
  static const std::pair<G4int, const G4double*> Z33[N33]={Z33N42};
  //==> Se(Z=34)
  static const G4int N34=6;
  static const G4double pZ34N40[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N40(40,pZ34N40);
  static const G4double pZ34N42[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N42(42,pZ34N42);
  static const G4double pZ34N43[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N43(43,pZ34N43);
  static const G4double pZ34N44[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N44(44,pZ34N44);
  static const G4double pZ34N46[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N46(46,pZ34N46);
  static const G4double pZ34N48[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z34N48(48,pZ34N48);
  static const std::pair<G4int, const G4double*> Z34[N34]={Z34N40, Z34N42, Z34N43,
                                                           Z34N44, Z34N46, Z34N48};
  //==> Br(Z=35)
  static const G4int N35=2;
  static const G4double pZ35N44[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z35N44(44,pZ35N44);
  static const G4double pZ35N46[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z35N46(46,pZ35N46);
  static const std::pair<G4int, const G4double*> Z35[N35]={Z35N44, Z35N46};
  //==> Kr(Z=36)
  static const G4int N36=6;
  static const G4double pZ36N42[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N42(42,pZ36N42);
  static const G4double pZ36N44[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N44(44,pZ36N44);
  static const G4double pZ36N46[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N46(46,pZ36N46);
  static const G4double pZ36N47[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N47(47,pZ36N47);
  static const G4double pZ36N48[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N48(48,pZ36N48);
  static const G4double pZ36N50[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z36N50(50,pZ36N50);
  static const std::pair<G4int, const G4double*> Z36[N36]={Z36N42, Z36N44, Z36N46,
                                                           Z36N47, Z36N48, Z36N50};
  //==> Rb(Z=37)
  static const G4int N37=2;
  static const G4double pZ37N48[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z37N48(48,pZ37N48);
  static const G4double pZ37N50[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z37N50(50,pZ37N50);
  static const std::pair<G4int, const G4double*> Z37[N37]={Z37N48, Z37N50};
  //==> Sr(Z=38)
  static const G4int N38=4;
  static const G4double pZ38N46[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N46(46,pZ38N46);
  static const G4double pZ38N48[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N48(48,pZ38N48);
  static const G4double pZ38N49[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N49(49,pZ38N49);
  static const G4double pZ38N50[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z38N50(50,pZ38N50);
  static const std::pair<G4int, const G4double*> Z38[N38]={Z38N46, Z38N48, Z38N49, Z38N50};
  //==> Y (Z=39)
  static const G4int N39=1;
  static const G4double pZ39N50[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z39N50(50,pZ39N50);
  static const std::pair<G4int, const G4double*> Z39[N39]={Z39N50};
  //==> Zr(Z=40)
  static const G4int N40=5;
  static const G4double pZ40N50[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N50(50,pZ40N50);
  static const G4double pZ40N51[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N51(51,pZ40N51);
  static const G4double pZ40N52[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N52(52,pZ40N52);
  static const G4double pZ40N54[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N54(54,pZ40N54);
  static const G4double pZ40N56[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z40N56(56,pZ40N56);
  static const std::pair<G4int, const G4double*> Z40[N40]={Z40N50, Z40N51, Z40N52,
                                                           Z40N54, Z40N56};
  //==> Nb(Z=41)
  static const G4int N41=1;
  static const G4double pZ41N52[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z41N52(52,pZ41N52);
  static const std::pair<G4int, const G4double*> Z41[N41]={Z41N52};
  //==> Mo(Z=42)
  static const G4int N42=7;
  static const G4double pZ42N50[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N50(50,pZ42N50);
  static const G4double pZ42N52[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N52(52,pZ42N52);
  static const G4double pZ42N53[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N53(53,pZ42N53);
  static const G4double pZ42N54[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N54(54,pZ42N54);
  static const G4double pZ42N55[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N55(55,pZ42N55);
  static const G4double pZ42N56[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N56(56,pZ42N56);
  static const G4double pZ42N58[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z42N58(58,pZ42N58);
  static const std::pair<G4int, const G4double*> Z42[N42]={Z42N50, Z42N52, Z42N53, Z42N54,
                                                           Z42N55, Z42N56, Z42N58};
  //==> Mo(Z=43)
  static const G4int N43=1;
  static const G4double pZ43N0[5]={.018, 5., .1, .03, .002}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z43N0(0,pZ43N0);
  static const std::pair<G4int, const G4double*> Z43[N43]={Z43N0};
  //==> Ru(Z=44)
  static const G4int N44=7;
  static const G4double pZ44N52[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N52(52,pZ44N52);
  static const G4double pZ44N54[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N54(54,pZ44N54);
  static const G4double pZ44N55[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N55(55,pZ44N55);
  static const G4double pZ44N56[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N56(56,pZ44N56);
  static const G4double pZ44N57[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N57(57,pZ44N57);
  static const G4double pZ44N58[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N58(58,pZ44N58);
  static const G4double pZ44N60[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z44N60(60,pZ44N60);
  static const std::pair<G4int, const G4double*> Z44[N44]={Z44N52, Z44N54, Z44N55, Z44N56,
                                                           Z44N57, Z44N58, Z44N60};
  //==> Rh(Z=45)
  static const G4int N45=1;
  static const G4double pZ45N58[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z45N58(58,pZ45N58);
  static const std::pair<G4int, const G4double*> Z45[N45]={Z45N58};
  //==> Pd(Z=46)
  static const G4int N46=6;
  static const G4double pZ46N56[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N56(56,pZ46N56);
  static const G4double pZ46N58[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N58(58,pZ46N58);
  static const G4double pZ46N59[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N59(59,pZ46N59);
  static const G4double pZ46N60[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N60(60,pZ46N60);
  static const G4double pZ46N62[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N62(62,pZ46N62);
  static const G4double pZ46N64[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z46N64(64,pZ46N64);
  static const std::pair<G4int, const G4double*> Z46[N46]={Z46N56, Z46N58, Z46N59,
                                                           Z46N60, Z46N62, Z46N64};
  //==> Ag(Z=47)
  static const G4int N47=2;
  static const G4double pZ47N60[5]={.018, 5., .1, .004, .003};
  static const std::pair<G4int, const G4double*> Z47N60(60,pZ47N60);
  static const G4double pZ47N62[5]={.018, 4., .015, .06, .0008};
  static const std::pair<G4int, const G4double*> Z47N62(62,pZ47N62);
  static const std::pair<G4int, const G4double*> Z47[N47]={Z47N60, Z47N62};
  //==> Cd(Z=48)
  static const G4int N48=8;
  static const G4double pZ48N58[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N58(58,pZ48N58);
  static const G4double pZ48N60[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N60(60,pZ48N60);
  static const G4double pZ48N62[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N62(62,pZ48N62);
  static const G4double pZ48N63[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N63(63,pZ48N63);
  static const G4double pZ48N64[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N64(64,pZ48N64);
  static const G4double pZ48N65[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N65(65,pZ48N65);
  static const G4double pZ48N66[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N66(66,pZ48N66);
  static const G4double pZ48N68[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z48N68(68,pZ48N68);
  static const std::pair<G4int, const G4double*> Z48[N48]={Z48N58, Z48N60, Z48N62, Z48N63,
                                                           Z48N64, Z48N65, Z48N66, Z48N68};
  //==> In(Z=49)
  static const G4int N49=2;
  static const G4double pZ49N64[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z49N64(64,pZ49N64);
  static const G4double pZ49N66[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z49N66(66,pZ49N66);
  static const std::pair<G4int, const G4double*> Z49[N49]={Z49N64, Z49N66};
  //==> Sn(Z=50)
  static const G4int N50=10;
  static const G4double pZ50N62[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N62(62,pZ50N62);
  static const G4double pZ50N64[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N64(64,pZ50N64);
  static const G4double pZ50N65[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N65(65,pZ50N65);
  static const G4double pZ50N66[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N66(66,pZ50N66);
  static const G4double pZ50N67[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N67(67,pZ50N67);
  static const G4double pZ50N68[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N68(68,pZ50N68);
  static const G4double pZ50N69[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N69(69,pZ50N69);
  static const G4double pZ50N70[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N70(70,pZ50N70);
  static const G4double pZ50N72[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N72(72,pZ50N72);
  static const G4double pZ50N74[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z50N74(74,pZ50N74);
  static const std::pair<G4int, const G4double*> Z50[N50]={Z50N62, Z50N64, Z50N65, Z50N66,
                                                           Z50N67, Z50N68, Z50N69, Z50N70,
                                                           Z50N72, Z50N74};
  //==> Sb(Z=51)
  static const G4int N51=2;
  static const G4double pZ51N70[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z51N70(70,pZ51N70);
  static const G4double pZ51N72[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z51N72(72,pZ51N72);
  static const std::pair<G4int, const G4double*> Z51[N51]={Z51N70, Z51N72};
  //==> Te(Z=52)
  static const G4int N52=8;
  static const G4double pZ52N68[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N68(68,pZ52N68);
  static const G4double pZ52N70[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N70(70,pZ52N70);
  static const G4double pZ52N71[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N71(71,pZ52N71);
  static const G4double pZ52N72[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N72(72,pZ52N72);
  static const G4double pZ52N73[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N73(73,pZ52N73);
  static const G4double pZ52N74[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N74(74,pZ52N74);
  static const G4double pZ52N76[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N76(76,pZ52N76);
  static const G4double pZ52N78[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z52N78(78,pZ52N78);
  static const std::pair<G4int, const G4double*> Z52[N52]={Z52N68, Z52N70, Z52N71, Z52N72,
                                                           Z52N73, Z52N74, Z52N76, Z52N78};
  //==> I (Z=53)
  static const G4int N53=1;
  static const G4double pZ53N74[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z53N74(74,pZ53N74);
  static const std::pair<G4int, const G4double*> Z53[N53]={Z53N74};
  //==> Xe(Z=54)
  static const G4int N54=9;
  static const G4double pZ54N70[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N70(70,pZ54N70);
  static const G4double pZ54N72[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N72(72,pZ54N72);
  static const G4double pZ54N74[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N74(74,pZ54N74);
  static const G4double pZ54N75[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N75(75,pZ54N75);
  static const G4double pZ54N76[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N76(76,pZ54N76);
  static const G4double pZ54N77[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N77(77,pZ54N77);
  static const G4double pZ54N78[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N78(78,pZ54N78);
  static const G4double pZ54N80[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N80(80,pZ54N80);
  static const G4double pZ54N82[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z54N82(82,pZ54N82);
  static const std::pair<G4int, const G4double*> Z54[N54]={Z54N70, Z54N72, Z54N74,
                                                           Z54N75, Z54N76, Z54N77,
                                                           Z54N78, Z54N80, Z54N82};
  //==> Cs(Z=55)
  static const G4int N55=1;
  static const G4double pZ55N78[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z55N78(78,pZ55N78);
  static const std::pair<G4int, const G4double*> Z55[N55]={Z55N78};
  //==> Ba(Z=56)
  static const G4int N56=7;
  static const G4double pZ56N74[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N74(70,pZ56N74);
  static const G4double pZ56N76[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N76(71,pZ56N76);
  static const G4double pZ56N78[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N78(72,pZ56N78);
  static const G4double pZ56N79[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N79(73,pZ56N79);
  static const G4double pZ56N80[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N80(74,pZ56N80);
  static const G4double pZ56N81[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N81(76,pZ56N81);
  static const G4double pZ56N82[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z56N82(78,pZ56N82);
  static const std::pair<G4int, const G4double*> Z56[N56]={Z56N74, Z56N76, Z56N78, Z56N79,
                                                           Z56N80, Z56N81, Z56N82};
  //==> La(Z=57)
  static const G4int N57=2;
  static const G4double pZ57N81[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z57N81(81,pZ57N81);
  static const G4double pZ57N82[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z57N82(82,pZ57N82);
  static const std::pair<G4int, const G4double*> Z57[N57]={Z57N81, Z57N82};
  //==> Ce(Z=58)
  static const G4int N58=4;
  static const G4double pZ58N78[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z58N78(78,pZ58N78);
  static const G4double pZ58N80[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z58N80(80,pZ58N80);
  static const G4double pZ58N82[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z58N82(82,pZ58N82);
  static const G4double pZ58N84[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z58N84(84,pZ58N84);
  static const std::pair<G4int, const G4double*> Z58[N58]={Z58N78, Z58N80, Z58N82, Z58N84};
  //==> Pr(Z=59)
  static const G4int N59=1;
  static const G4double pZ59N82[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z59N82(82,pZ59N82);
  static const std::pair<G4int, const G4double*> Z59[N59]={Z59N82};
  //==> Nd(Z=60)
  static const G4int N60=7;
  static const G4double pZ60N82[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N82(82,pZ60N82);
  static const G4double pZ60N83[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N83(83,pZ60N83);
  static const G4double pZ60N84[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N84(84,pZ60N84);
  static const G4double pZ60N85[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N85(85,pZ60N85);
  static const G4double pZ60N86[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N86(86,pZ60N86);
  static const G4double pZ60N88[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N88(88,pZ60N88);
  static const G4double pZ60N90[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z60N90(90,pZ60N90);
  static const std::pair<G4int, const G4double*> Z60[N60]={Z60N82, Z60N83, Z60N84, Z60N85,
                                                           Z60N86, Z60N88, Z60N90};
  //==> Mo(Z=61)
  static const G4int N61=1;
  static const G4double pZ61N0[5]={.018, 5., .1, .03, .002}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z61N0(0,pZ61N0);
  static const std::pair<G4int, const G4double*> Z61[N61]={Z61N0};
  //==> Sm(Z=62)
  static const G4int N62=7;
  static const G4double pZ62N82[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N82(82,pZ62N82);
  static const G4double pZ62N85[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N85(85,pZ62N85);
  static const G4double pZ62N86[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N86(86,pZ62N86);
  static const G4double pZ62N87[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N87(87,pZ62N87);
  static const G4double pZ62N88[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N88(88,pZ62N88);
  static const G4double pZ62N90[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N90(90,pZ62N90);
  static const G4double pZ62N92[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z62N92(92,pZ62N92);
  static const std::pair<G4int, const G4double*> Z62[N62]={Z62N82, Z62N85, Z62N86, Z62N87,
                                                           Z62N88, Z62N90, Z62N92};
  //==> Eu(Z=63)
  static const G4int N63=2;
  static const G4double pZ63N88[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z63N88(88,pZ63N88);
  static const G4double pZ63N90[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z63N90(90,pZ63N90);
  static const std::pair<G4int, const G4double*> Z63[N63]={Z63N88, Z63N90};
  //==> Gd(Z=64)
  static const G4int N64=7;
  static const G4double pZ64N88[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N88(88,pZ64N88);
  static const G4double pZ64N90[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N90(90,pZ64N90);
  static const G4double pZ64N91[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N91(91,pZ64N91);
  static const G4double pZ64N92[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N92(92,pZ64N92);
  static const G4double pZ64N93[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N93(93,pZ64N93);
  static const G4double pZ64N94[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N94(94,pZ64N94);
  static const G4double pZ64N96[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z64N96(96,pZ64N96);
  static const std::pair<G4int, const G4double*> Z64[N64]={Z64N88, Z64N90, Z64N91, Z64N92,
                                                           Z64N93, Z64N94, Z64N96};
  //==> Tb(Z=65)
  static const G4int N65=1;
  static const G4double pZ65N94[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z65N94(82,pZ65N94);
  static const std::pair<G4int, const G4double*> Z65[N65]={Z65N94};
  //==> Dy(Z=66)
  static const G4int N66=7;
  static const G4double pZ66N90[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N90(90,pZ66N90);
  static const G4double pZ66N92[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N92(92,pZ66N92);
  static const G4double pZ66N94[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N94(94,pZ66N94);
  static const G4double pZ66N95[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N95(95,pZ66N95);
  static const G4double pZ66N96[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N96(96,pZ66N96);
  static const G4double pZ66N97[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N97(97,pZ66N97);
  static const G4double pZ66N98[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z66N98(98,pZ66N98);
  static const std::pair<G4int, const G4double*> Z66[N66]={Z66N90, Z66N92, Z66N94, Z66N95,
                                                           Z66N96, Z66N97, Z66N98};
  //==> Ho(Z=67)
  static const G4int N67=1;
  static const G4double pZ67N98[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z67N98(98,pZ67N98);
  static const std::pair<G4int, const G4double*> Z67[N67]={Z67N98};
  //==> Er(Z=68)
  static const G4int N68=6;
  static const G4double pZ68N94[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N94(94,pZ68N94);
  static const G4double pZ68N96[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N96(96,pZ68N96);
  static const G4double pZ68N98[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N98(98,pZ68N98);
  static const G4double pZ68N99[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N99(99,pZ68N99);
  static const G4double pZ68N100[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N100(100,pZ68N100);
  static const G4double pZ68N102[5]={.01, 27., 0., 0., 1.}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z68N102(102,pZ68N102);
  static const std::pair<G4int, const G4double*> Z68[N68]={Z68N94, Z68N96, Z68N98,
                                                           Z68N99, Z68N100, Z68N102};
  //==> Tm(Z=69)
  static const G4int N69=1;
  static const G4double pZ69N100[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z69N100(100,pZ69N100);
  static const std::pair<G4int, const G4double*> Z69[N69]={Z69N100};
  //==> Yb(Z=70)
  static const G4int N70=7;
  static const G4double pZ70N98[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N98(98,pZ70N98);
  static const G4double pZ70N100[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N100(100,pZ70N100);
  static const G4double pZ70N101[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N101(101,pZ70N101);
  static const G4double pZ70N102[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N102(102,pZ70N102);
  static const G4double pZ70N103[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N103(103,pZ70N103);
  static const G4double pZ70N104[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N104(104,pZ70N104);
  static const G4double pZ70N106[5]={.01, 27., 0., 0., 1.}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z70N106(106,pZ70N106);
  static const std::pair<G4int, const G4double*> Z70[N70]={Z70N98, Z70N100, Z70N101,
                                                           Z70N102, Z70N103, Z70N104,
                                                           Z70N106};
  //==> Lu(Z=71)
  static const G4int N71=2;
  static const G4double pZ71N104[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z71N104(104,pZ71N104);
  static const G4double pZ71N105[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z71N105(105,pZ71N105);
  static const std::pair<G4int, const G4double*> Z71[N71]={Z71N104, Z71N105};
  //==> Hf(Z=72)
  static const G4int N72=6;
  static const G4double pZ72N102[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N102(102,pZ72N102);
  static const G4double pZ72N104[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N104(104,pZ72N104);
  static const G4double pZ72N105[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N105(105,pZ72N105);
  static const G4double pZ72N106[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N106(106,pZ72N106);
  static const G4double pZ72N107[5]={.018, 5., .1, .03, .002}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N107(107,pZ72N107);
  static const G4double pZ72N108[5]={.01, 27., 0., 0., 1.}; // *** NotImplemented ***
  static const std::pair<G4int, const G4double*> Z72N108(108,pZ72N108);
  static const std::pair<G4int, const G4double*> Z72[N72]={Z72N102, Z72N104, Z72N105,
                                                           Z72N106, Z72N107, Z72N108};
  //==> Ta(Z=73)
  static const G4int N73=1;
  static const G4double pZ73N108[5]={.0065, 2., .15, .001, .0012};
  static const std::pair<G4int, const G4double*> Z73N108(108,pZ73N108);
  static const std::pair<G4int, const G4double*> Z73[N73]={Z73N108};
  //==> W (Z=74)
  static const G4int N74=5;
  static const G4double pZ74N106[5]={.014, 5., .03, .04, .001}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z74N106(106,pZ74N106);
  static const G4double pZ74N108[5]={.013, 4.5, .35, .045, .001};
  static const std::pair<G4int, const G4double*> Z74N108(108,pZ74N108);
  static const G4double pZ74N109[5]={.012, 4., .02, .04, .001};
  static const std::pair<G4int, const G4double*> Z74N109(109,pZ74N109);
  static const G4double pZ74N110[5]={.014, 6., .02, .03, .0015};
  static const std::pair<G4int, const G4double*> Z74N110(110,pZ74N110);
  static const G4double pZ74N112[5]={.014, 6., .02, .03, .0015};
  static const std::pair<G4int, const G4double*> Z74N112(112,pZ74N112);
  static const std::pair<G4int, const G4double*> Z74[N74]={Z74N106, Z74N108, Z74N109,
                                                           Z74N110, Z74N112};
  //==> Re(Z=75)
  static const G4int N75=2;
  static const G4double pZ75N110[5]={.015, 4., .2, .0001, .0021};
  static const std::pair<G4int, const G4double*> Z75N110(110,pZ75N110);
  static const G4double pZ75N112[5]={.015, 4., .1, .0001, .002};
  static const std::pair<G4int, const G4double*> Z75N112(112,pZ75N112);
  static const std::pair<G4int, const G4double*> Z75[N75]={Z75N110, Z75N112};
  //==> Os(Z=76)
  static const G4int N76=7;
  static const G4double pZ76N108[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N108(108,pZ76N108);
  static const G4double pZ76N110[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N110(110,pZ76N110);
  static const G4double pZ76N111[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N111(111,pZ76N111);
  static const G4double pZ76N112[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N112(112,pZ76N112);
  static const G4double pZ76N113[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N113(113,pZ76N113);
  static const G4double pZ76N114[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N114(114,pZ76N114);
  static const G4double pZ76N116[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z76N116(116,pZ76N116);
  static const std::pair<G4int, const G4double*> Z76[N76]={Z76N108, Z76N110, Z76N111,
                                                           Z76N112, Z76N113, Z76N114,
                                                           Z76N116};
  //==> Ir(Z=77)
  static const G4int N77=2;
  static const G4double pZ77N114[5]={.012, 3., .1, .0001, .003};
  static const std::pair<G4int, const G4double*> Z77N114(114,pZ77N114);
  static const G4double pZ77N116[5]={.012, 3., .08, .0001, .002};
  static const std::pair<G4int, const G4double*> Z77N116(116,pZ77N116);
  static const std::pair<G4int, const G4double*> Z77[N77]={Z77N114, Z77N116};
  //==> Pt(Z=78)
  static const G4int N78=6;
  static const G4double pZ78N112[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N112(112,pZ78N112);
  static const G4double pZ78N114[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N114(114,pZ78N114);
  static const G4double pZ78N116[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N116(116,pZ78N116);
  static const G4double pZ78N117[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N117(117,pZ78N117);
  static const G4double pZ78N118[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N118(118,pZ78N118);
  static const G4double pZ78N120[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z78N120(120,pZ78N120);
  static const std::pair<G4int, const G4double*> Z78[N78]={Z78N112, Z78N114, Z78N116,
                                                           Z78N117, Z78N118, Z78N120};
  //==> Au(Z=79)
  static const G4int N79=1;
  static const G4double pZ79N118[5]={.012, 4., .2, .0001, .002};
  static const std::pair<G4int, const G4double*> Z79N118(118,pZ79N118);
  static const std::pair<G4int, const G4double*> Z79[N79]={Z79N118};
  //==> Hg(Z=80)
  static const G4int N80=7;
  static const G4double pZ80N116[5]={.022, 11., .006, .044, .001};
  static const std::pair<G4int, const G4double*> Z80N116(116,pZ80N116);
  static const G4double pZ80N118[5]={.024, 10., .04, .018, .0016};
  static const std::pair<G4int, const G4double*> Z80N118(118,pZ80N118);
  static const G4double pZ80N119[5]={.017, 8., .02, .02, .0016};
  static const std::pair<G4int, const G4double*> Z80N119(119,pZ80N119);
  static const G4double pZ80N120[5]={.021, 9., .03, .02, .0016};
  static const std::pair<G4int, const G4double*> Z80N120(120,pZ80N120);
  static const G4double pZ80N121[5]={.01, 7., .025, .025, .0016};
  static const std::pair<G4int, const G4double*> Z80N121(121,pZ80N121);
  static const G4double pZ80N122[5]={.023, 9., .008, .045, .0012};
  static const std::pair<G4int, const G4double*> Z80N122(122,pZ80N122);
  static const G4double pZ80N124[5]={.018, 8., .007, .04, .0013};
  static const std::pair<G4int, const G4double*> Z80N124(124,pZ80N124);
  static const std::pair<G4int, const G4double*> Z80[N80]={Z80N116, Z80N118, Z80N119,
                                                           Z80N120, Z80N121, Z80N122,
                                                           Z80N124};
  //==> Tl(Z=81)
  static const G4int N81=2;
  static const G4double pZ81N122[5]={.018, 5., .1, .03, .002}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z81N122(122,pZ81N122);
  static const G4double pZ81N124[5]={.01, 27., 0., 0., 1.}; // *** No DATA ***
  static const std::pair<G4int, const G4double*> Z81N124(124,pZ81N124);
  static const std::pair<G4int, const G4double*> Z81[N81]={Z81N122, Z81N124};
  //==> Pb(Z=82)
  static const G4int N82=4;
  static const G4double pZ82N122[5]={.032, 9., .001, .13, .004};
  static const std::pair<G4int, const G4double*> Z82N122(122,pZ82N122);
  static const G4double pZ82N124[5]={.02, 7., .005, .14, .002};
  static const std::pair<G4int, const G4double*> Z82N124(124,pZ82N124);
  static const G4double pZ82N125[5]={.021, 9., .0012, .05, .01};
  static const std::pair<G4int, const G4double*> Z82N125(125,pZ82N125);
  static const G4double pZ82N126[5]={.049, 14., .0007, .145, .001};
  static const std::pair<G4int, const G4double*> Z82N126(126,pZ82N126);
  static const std::pair<G4int, const G4double*> Z82[N82]={Z82N122, Z82N124, Z82N125,
                                                           Z82N126};
  //==> Bi(Z=83)
  static const G4int N83=1;
  static const G4double pZ83N126[5]={.033, 10., .001, .13, .006};
  static const std::pair<G4int, const G4double*> Z83N126(126,pZ83N126);
  static const std::pair<G4int, const G4double*> Z83[N83]={Z83N126};
  //==> Po(Z=84)
  static const G4int N84=1;
  static const G4double pZ84N0[5]={.01, 27., 0., 0., 1.}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z84N0(0,pZ84N0);
  static const std::pair<G4int, const G4double*> Z84[N84]={Z84N0};
  //==> At(Z=85)
  static const G4int N85=1;
  static const G4double pZ85N0[5]={.018, 5., .1, .03, .002}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z85N0(0,pZ85N0);
  static const std::pair<G4int, const G4double*> Z85[N85]={Z85N0};
  //==> Rn(Z=86)
  static const G4int N86=1;
  static const G4double pZ86N0[5]={.018, 5., .1, .03, .002}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z86N0(0,pZ86N0);
  static const std::pair<G4int, const G4double*> Z86[N86]={Z86N0};
  //==> Fr(Z=87)
  static const G4int N87=1;
  static const G4double pZ87N0[5]={.018, 5., .1, .03, .002}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z87N0(0,pZ87N0);
  static const std::pair<G4int, const G4double*> Z87[N87]={Z87N0};
  //==> Ra(Z=88)
  static const G4int N88=1;
  static const G4double pZ88N138[5]={.012, 4.4, .068, .051, .0008};
  static const std::pair<G4int, const G4double*> Z88N138(138,pZ88N138);
  static const std::pair<G4int, const G4double*> Z88[N88]={Z88N138};
  //==> Ac(Z=89)
  static const G4int N89=1;
  static const G4double pZ89N0[5]={.018, 5., .1, .03, .002}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z89N0(0,pZ89N0);
  static const std::pair<G4int, const G4double*> Z89[N89]={Z89N0};
  //==> Th(Z=90)
  static const G4int N90=1;
  static const G4double pZ90N142[5]={.01, 3.6, .07, .029, .0009};
  static const std::pair<G4int, const G4double*> Z90N142(142,pZ90N142);
  static const std::pair<G4int, const G4double*> Z90[N90]={Z90N142};
  //==> Pa(Z=91)
  static const G4int N91=1;
  static const G4double pZ91N0[5]={.01, 27., 0., 0., 1.}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z91N0(0,pZ91N0);
  static const std::pair<G4int, const G4double*> Z91[N91]={Z91N0};
  //==> U (Z=92)
  static const G4int N92=2;
  static const G4double pZ92N143[5]={.005, 2.5, .06, .008, .002};
  static const std::pair<G4int, const G4double*> Z92N143(143,pZ92N143);
  static const G4double pZ92N146[5]={.009, 3.1, .04, .025, .001};
  static const std::pair<G4int, const G4double*> Z92N146(146,pZ92N146);
  static const std::pair<G4int, const G4double*> Z92[N92]={Z92N143, Z92N146};
  //==> Np(Z=93)
  static const G4int N93=1;
  static const G4double pZ93N144[5]={.009, 2.35, .35, .003, .0007};
  static const std::pair<G4int, const G4double*> Z93N144(144,pZ93N144);
  static const std::pair<G4int, const G4double*> Z93[N93]={Z93N144};
  //==> Pu(Z=94)
  static const G4int N94=3;
  static const G4double pZ94N145[5]={.005, 2.6, .04, .023, .0002}; // *** Artificial ***
  static const std::pair<G4int, const G4double*> Z94N145(145,pZ94N145);
  static const G4double pZ94N148[5]={.009, 3.1, .06, .02, .001}; // *** Artificial ***
  static const std::pair<G4int, const G4double*> Z94N148(148,pZ94N148);
  static const G4double pZ94N150[5]={.01, 27., 0., 0., 1.};
  static const std::pair<G4int, const G4double*> Z94N150(150,pZ94N150);
  static const std::pair<G4int, const G4double*> Z94[N94]={Z94N145, Z94N148, Z94N150};
  //==> Am(Z=95)
  static const G4int N95=1;
  static const G4double pZ95N0[5]={.018, 5., .1, .03, .002}; // *** NoStableIsotopes ***
  static const std::pair<G4int, const G4double*> Z95N0(0,pZ95N0);
  static const std::pair<G4int, const G4double*> Z95[N95]={Z95N0};
  //==> Cm(Z=96)
  static const G4int N96=1;
  static const G4double pZ96N151[5]={.005, 2.5, .07, .027, .0009};
  static const std::pair<G4int, const G4double*> Z96N151(151,pZ96N151);
  static const std::pair<G4int, const G4double*> Z96[N96]={Z96N151};
 
  static const G4int NZ=97; // #of Elements covered by CHIPS
  static const std::pair<G4int, const G4double*>* Pars[NZ]={Z0,Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8,Z9,
    Z10,Z11,Z12,Z13,Z14,Z15,Z16,Z17,Z18,Z19,Z20,Z21,Z22,Z23,Z24,Z25,Z26,Z27,Z28,Z29,Z30,
    Z31,Z32,Z33,Z34,Z35,Z36,Z37,Z38,Z39,Z40,Z41,Z42,Z43,Z44,Z45,Z46,Z47,Z48,Z49,Z50,Z51,
    Z52,Z53,Z54,Z55,Z56,Z57,Z58,Z59,Z60,Z61,Z62,Z63,Z64,Z65,Z66,Z67,Z68,Z69,Z70,Z71,Z72,
    Z73,Z74,Z75,Z76,Z77,Z78,Z79,Z80,Z81,Z82,Z83,Z84,Z85,Z86,Z87,Z88,Z89,Z90,Z91,Z92,Z93,
    Z94,Z95,Z96};
  static const G4int NIso[NZ]={N0,N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,
    N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27,N28,N29,N30,N31,N32,N33,N34,N35,N36,N37,
    N38,N39,N40,N41,N42,N43,N44,N45,N46,N47,N48,N49,N50,N51,N52,N53,N54,N55,N56,N57,N58,
    N59,N60,N61,N62,N63,N64,N65,N66,N67,N68,N69,N70,N71,N72,N73,N74,N75,N76,N77,N78,N79,
    N80,N81,N82,N83,N84,N85,N86,N87,N88,N89,N90,N91,N92,N93,N94,N95,N96};
  static G4int    mZ=0;
  static G4int    mN=0;
  static G4double P=.1;  // s=SQRT(M_N^2+M_h^2+2*E_h*M_N)
  static G4double R=0.;  // Prototype of the result
  static G4double p0=0.;
  static G4double p1=0.;
  static G4double p2=0.;
  static G4double p3=0.;
  static G4double p4=0.;
  G4int A=Z+N;
  if(A<=1) return 0.;
  if(p<=0.001) return 1.;
  if(Z!=mZ || N!=mN)     // Recalculate the parameters for different isotope
  {
    mZ=Z;
    mN=N;
    P=p;
    if(Z<97 && N<248)                // General solution (*** Z/A limits ***)
    {
      p0=.01; // Default guess if there is no data
      p1=0.;
      p2=27.;
      p3=0.;
      p4=1.;
      G4int nn=NIso[Z];
      G4bool nfound=true;
      if(nn) for (G4int in=0; in<nn; in++)
      {
        std::pair<G4int, const G4double*> curIs=Pars[Z][in];
        if(curIs.first == N)
        {
          const G4double* curT=curIs.second;
          p0 = curT[0];
          if(p < p0)
          {
            R=1.;
            return R;
          }
          p1 = curT[1];
          p2 = curT[2];
          p3 = curT[3];
          p4 = curT[4];
          nfound = false;
          break;
        }
      }
      if(nfound) G4cout<<"-Warning-G4QNeutronCaptureRatio::CSLin: Z="<<Z<<", N="
                       <<N<<" isotope is not implemented in CHIPS"<<G4endl;
      R=std::pow(p0/p,p1);
      if(p2>0.)
      {
        G4double dp=p-p3;
        R+=p2*std::exp(-dp*dp/p4);
      }
      if(R>1.) R=1.;
    }
    else G4cerr<<"-Warning-G4QNeutronCaptureRatio::CalcR:*Bad A* Z="<<Z<<", N="<<N<<G4endl;
  }
  else if(std::fabs(p-P)/P<.0001) return R; // Normally used at high energies (direct) only
  return R;
} // End of CalcCap2IN_Ratio
