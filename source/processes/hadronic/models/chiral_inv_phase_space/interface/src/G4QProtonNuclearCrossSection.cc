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
// The lust update: M.V. Kossov, CERN/ITEP(Moscow) 17-June-02
// GEANT4 tag $Name: geant4-09-00 $
//
//
// G4 Physics class: G4QProtonNuclearCrossSection for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-03
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Feb-04
// --------------------------------------------------------------------------------
// ****************************************************************************************
// ***** This HEADER is a property of the CHIPS hadronic package in Geant4 (M. Kosov) *****
// *********** DO NOT MAKE ANY CHANGE without approval of Mikhail.Kossov@cern.ch **********
// ****************************************************************************************
//
//#define debug
//#define pdebug
//#define debug3
//#define debugn
//#define debugs

#include "G4QProtonNuclearCrossSection.hh"

// Initialization of the
G4double* G4QProtonNuclearCrossSection::lastLEN=0; // Pointer to the lastArray of LowEn CS
G4double* G4QProtonNuclearCrossSection::lastHEN=0; // Pointer to the lastArray of HighEn CS
G4double  G4QProtonNuclearCrossSection::lastSP=0.; // Last value of ShadowingPomeron(A-dep)
G4int     G4QProtonNuclearCrossSection::lastPDG=0; // The last PDG code of the projectile
G4int     G4QProtonNuclearCrossSection::lastN=0;   // The last N of calculated nucleus
G4int     G4QProtonNuclearCrossSection::lastZ=0;   // The last Z of calculated nucleus
G4double  G4QProtonNuclearCrossSection::lastP=0.;  // Last used in cross section Momentum
G4double  G4QProtonNuclearCrossSection::lastTH=0.; // Last threshold momentum
G4double  G4QProtonNuclearCrossSection::lastCS=0.; // Last value of the Cross Section
G4int     G4QProtonNuclearCrossSection::lastI=0;   // The last position in the DAMDB

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QProtonNuclearCrossSection::GetPointer()
{
  static G4QProtonNuclearCrossSection theCrossSection; //**Static body of Cross Section**
  return &theCrossSection;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4QProtonNuclearCrossSection::GetCrossSection(G4bool fCS, G4double pMom,
                                                       G4int tgZ, G4int tgN, G4int pPDG)
{
  static G4int j;                      // A#0f records found in DB for this projectile
  static std::vector <G4int>    colPDG;// Vector of the projectile PDG code
  static std::vector <G4int>    colN;  // Vector of N for calculated nuclei (isotops)
  static std::vector <G4int>    colZ;  // Vector of Z for calculated nuclei (isotops)
  static std::vector <G4double> colP;  // Vector of last momenta for the reaction
  static std::vector <G4double> colTH; // Vector of energy thresholds for the reaction
  static std::vector <G4double> colCS; // Vector of last cross sections for the reaction
  // ***---*** End of the mandatory Static Definitions of the Associative Memory ***---***
  G4double pEn=pMom;
#ifdef pdebug
  G4cout<<"G4QPrCS::GetCS:>>> f="<<fCS<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="<<tgN
        <<"("<<lastN<<"),PDG="<<pPDG<<"("<<lastPDG<<"), T="<<pEn<<"("<<lastTH<<")"<<",Sz="
        <<colN.size()<<G4endl;
		//CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  if(!pPDG)
  {
#ifdef pdebug
    G4cout<<"G4QPrCS::GetCS: *** Found pPDG="<<pPDG<<" ====> CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return 0.;                         // projectile PDG=0 is a mistake (?!) @@
  }
  G4bool in=false;                     // By default the isotope must be found in the AMDB
  if(tgN!=lastN || tgZ!=lastZ || pPDG!=lastPDG)// The nucleus was not the last used isotope
  {
    in = false;                        // By default the isotope haven't be found in AMDB  
    lastP   = 0.;                      // New momentum history (nothing to compare with)
    lastPDG = pPDG;                    // The last PDG of the projectile
    lastN   = tgN;                     // The last N of the calculated nucleus
    lastZ   = tgZ;                     // The last Z of the calculated nucleus
    lastI   = colN.size();             // Size of the Associative Memory DB in the heap
    j  = 0;                            // A#0f records found in DB for this projectile
    if(lastI) for(G4int i=0; i<lastI; i++) if(colPDG[i]==pPDG) // The partType is found
	   {                                  // The nucleus with projPDG is found in AMDB
      if(colN[i]==tgN && colZ[i]==tgZ)
						{
        lastI=i;
        lastTH =colTH[i];                // Last THreshold (A-dependent)
#ifdef pdebug
        G4cout<<"G4QPrCS::GetCS:*Found* P="<<pMom<<",Threshold="<<lastTH<<",j="<<j<<G4endl;
        //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(pEn<=lastTH)
        {
#ifdef pdebug
          G4cout<<"G4QPrCS::GetCS:Found T="<<pEn<<" < Threshold="<<lastTH<<",CS=0"<<G4endl;
          //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
          return 0.;                     // Energy is below the Threshold value
        }
        lastP  =colP [i];                // Last Momentum  (A-dependent)
        lastCS =colCS[i];                // Last CrossSect (A-dependent)
        if(std::fabs(lastP/pMom-1.)<tolerance)
        {
#ifdef pdebug
          G4cout<<"G4QPrCS::GetCS:P="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
#endif
          CalculateCrossSection(fCS,-1,j,lastPDG,lastZ,lastN,pMom); // Update param's only
          return lastCS*millibarn;     // Use theLastCS
        }
        in = true;                       // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
#ifdef pdebug
        G4cout<<"G4QPrCS::G:UpdatDB P="<<pMom<<",f="<<fCS<<",lI="<<lastI<<",j="<<j<<G4endl;
#endif
        lastCS=CalculateCrossSection(fCS,-1,j,lastPDG,lastZ,lastN,pMom); // read & update
#ifdef pdebug
        G4cout<<"G4QPrCS::GetCrosSec: *****> New (inDB) Calculated CS="<<lastCS<<G4endl;
        //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(lastCS<=0. && pEn>lastTH)    // Correct the threshold
        {
#ifdef pdebug
          G4cout<<"G4QPrCS::GetCS: New T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
        break;                           // Go out of the LOOP
      }
#ifdef pdebug
      G4cout<<"---G4QPrCrossSec::GetCrosSec:pPDG="<<pPDG<<",j="<<j<<",N="<<colN[i]
            <<",Z["<<i<<"]="<<colZ[i]<<",cPDG="<<colPDG[i]<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      j++;                             // Increment a#0f records found in DB for this pPDG
	   }
	   if(!in)                            // This nucleus has not been calculated previously
	   {
#ifdef pdebug
      G4cout<<"G4QPrCS::GetCrosSec:CalcNew P="<<pMom<<",f="<<fCS<<",lastI="<<lastI<<G4endl;
#endif
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastCS=CalculateCrossSection(fCS,0,j,lastPDG,lastZ,lastN,pMom); //calculate & create
      if(lastCS<=0.)
						{
        lastTH = ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
#ifdef pdebug
        G4cout<<"G4QPrCrossSection::GetCrossSect: NewThresh="<<lastTH<<",T="<<pEn<<G4endl;
#endif
        if(pEn>lastTH)
        {
#ifdef pdebug
          G4cout<<"G4QPrCS::GetCS: First T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
						}
#ifdef pdebug
      G4cout<<"G4QPrCS::GetCrosSec: New CS="<<lastCS<<",lZ="<<lastN<<",lN="<<lastZ<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      colN.push_back(tgN);
      colZ.push_back(tgZ);
      colPDG.push_back(pPDG);
      colP.push_back(pMom);
      colTH.push_back(lastTH);
      colCS.push_back(lastCS);
#ifdef pdebug
      G4cout<<"G4QPrCS::GetCS:1st,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      return lastCS*millibarn;
	   } // End of creation of the new set of parameters
    else
				{
#ifdef pdebug
      G4cout<<"G4QPrCS::GetCS: Update lastI="<<lastI<<",j="<<j<<G4endl;
#endif
      colP[lastI]=pMom;
      colPDG[lastI]=pPDG;
      colCS[lastI]=lastCS;
    }
  } // End of parameters udate
  else if(pEn<=lastTH)
  {
#ifdef pdebug
    G4cout<<"G4QPrCS::GetCS: Current T="<<pEn<<" < Threshold="<<lastTH<<", CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return 0.;                         // Momentum is below the Threshold Value -> CS=0
  }
  else if(std::fabs(lastP/pMom-1.)<tolerance)
  {
#ifdef pdebug
    G4cout<<"G4QPrCS::GetCS:OldCur P="<<pMom<<"="<<pMom<<", CS="<<lastCS*millibarn<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return lastCS*millibarn;     // Use theLastCS
  }
  else
  {
#ifdef pdebug
    G4cout<<"G4QPrCS::GetCS:UpdatCur P="<<pMom<<",f="<<fCS<<",I="<<lastI<<",j="<<j<<G4endl;
#endif
    lastCS=CalculateCrossSection(fCS,1,j,lastPDG,lastZ,lastN,pMom); // Only UpdateDB
    lastP=pMom;
  }
#ifdef pdebug
  G4cout<<"G4QPrCS::GetCroSec:End,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
  //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  return lastCS*millibarn;
}

// The main member function giving the gamma-A cross section (E in GeV, CS in mb)
G4double G4QProtonNuclearCrossSection::CalculateCrossSection(G4bool, G4int F, G4int I,
                                        G4int, G4int targZ, G4int targN, G4double Momentum)
{
  static const G4double THmin=0.;  // minimum Energy Threshold
  //static const G4double dP=1.;     // step for the LEN table
  static const G4int    nL=105;    // A#of LENesonance points in E (each MeV from 2 to 106)
  //static const G4double Pmin=THmin+(nL-1)*dP; // minE for the HighE part
  //static const G4double Pmax=50000.;       // maxE for the HighE part
  static const G4int    nH=224;            // A#of HResonance points in lnE
  //static const G4double milP=log(Pmin);    // Low logarithm energy for the HighE part
  //static const G4double malP=log(Pmax);    // High logarithm energy (each 2.75 percent)
  //static const G4double dlP=(malP-milP)/(nH-1); // Step in log energy in the HighE part
  //
  // Associative memory for acceleration
  static std::vector <G4double> spA;     // shadowing coefficients (A-dependent)
  static std::vector <G4double*> LEN;    // Vector of pointers to LowEnProtonCrossSection
  static std::vector <G4double*> HEN;    // Vector of pointers to HighEnProtonCrossSection
#ifdef debug
  G4cout<<"G4QProtonNuclearCrossSection::CalcCS: N="<<tN<<",Z="<<tZ<<",E="<<Energy<<G4endl;
#endif
  if (Momentum<THmin) return 0.;      // @@ This can be dangerouse for the heaviest nuc.!
  G4double sigma=0.;
  G4double A=targN+targZ;
  if(F<=0)                           // This isotope was not the last used isotop
  {
    if(F<0)                          // This isotope was found in DAMDB =========> RETRIEVE
				{
      lastLEN=LEN[I];                // Pointer to prepared LowEnergy cross sections
      lastHEN=HEN[I];                // Pointer to prepared High Energy cross sections
      lastSP =spA[I];                // Shadowing coefficient for UHE
    }
	   else                             // This isotope wasn't calculated previously => CREATE
	   {
      G4double lnA=std::log(A);           // The nucleus is not found in DB. It is new.
      if(A==1.) lastSP=1.;           // @@ The Reggeon shadowing (A=1)
      else      lastSP=lnA;          // @@ The Reggeon shadowing
#ifdef debug
      G4cout<<"G4QProtonNuclearCrossSection::CalcCS:lnA="<<lnA<<",lastSP="<<lastSP<<G4endl;
#endif
#ifdef debug3
      if(A==3) G4cout<<"G4QProtonNuclearCrossSection::CalcCS: lastSP="<<lastSP<<G4endl;
#endif
      lastLEN = new G4double[nL];        // Allocate memory for the new LEN cross sections
      lastHEN = new G4double[nH];        // Allocate memory for the new HEN cross sections
      G4int er=GetFunctions(A,lastLEN,lastHEN);// set newZeroPosition and fill theFunctions
	     if(er<1) G4cerr<<"***G4QProtNucCroSec::CalcCrossSection: A="<<A<<" failed"<<G4endl;
#ifdef debug
      G4cout<<"G4QProtonNuclearCrossSection::CalcCS: GetFunctions er="<<er<<G4endl;
#endif
      // *** The synchronization check ***
      G4int sync=LEN.size();
      if(sync!=I) G4cerr<<"***G4PhoronNuclCS::CalcCrossSect: Sync="<<sync<<"#"<<I<<G4endl;
      LEN.push_back(lastLEN);          // added LEN, found by AH 10/7/02
      HEN.push_back(lastHEN);          // added HEN, found by AH 10/7/02
      spA.push_back(lastSP);           // Pomeron Shadowing
	   } // End of creation of the new set of parameters
  } // End of parameters udate
  // ============================== NOW the Magic Formula =================================
  if (Momentum<lastTH) return 0.;      // It must be already checked in the interface class
		// else if (Momentum<Pmin)              // LEN region (approximated in E, not in lnE)
  //{
#ifdef debug
	 //  G4cout<<"G4QPrNCS::CalcCS:bLEN A="<<A<<", nL="<<nL<<",TH="<<THmin<<",dP="<<dP<<G4endl;
#endif
  //  if(A<=1.) sigma=0.;
  //  else      sigma=EquLinearFit(Momentum,nL,THmin,dP,lastLEN);
#ifdef debugn
	 //  if(sigma<0.)
  //    G4cout<<"G4QPrNuCS::CalcCS:A="<<A<<",E="<<Momentum<<",T="<<THmin<<",dP="<<dP<<G4endl;
#endif
  //}
  //else if (Momentum<Pmax)                     // High Energy region
  //{
  //  G4double lP=log(Momentum);
#ifdef debug
  //  G4cout<<"G4QProtNucCS::CalcCS: before HEN nH="<<nH<<",iE="<<milP<<",dlP="<<dlP<<G4endl;
#endif
  //  sigma=EquLinearFit(lP,nH,milP,dlP,lastHEN);
  //}
  else                                      // UHE region (calculation, not frequent)
  {
    G4double P=0.001*Momentum;              // Approximation formula is for P in GeV/c
    G4double lP=std::log(P);
    if(targZ==1&&!targN)                    // At present only for p, n, and d targets
    {
      G4double ds=lP-4.;
      sigma=5.3844/(.0018886+P*P)+(.3*ds*ds+39.+4.85/P)/(1+std::exp(-(.064+lP)/.27));
    }
    else if(!targZ&&targN==1)               // At present only for p, n, and d targets
    {
      G4double ds=lP-4.;
      sigma=18.045/(.00210946+P*P)+(.3*ds*ds+39.+4.85/P)/(1+std::exp((.1812-lP)/.3655));
    }
    else if(targZ==1&&targN==1)             // At present only for p, n, and d targets
    {
      G4double ds=lP-4.;
      sigma=7.4818/(.000656+P*P)+(.6*ds*ds+72.8+9.7/P)/(1+std::exp(-(.22+lP)/.299));
    }
    else
    {
      G4cerr<<"G4ProtonNucCroSect::CalcCS:only pp, pn, pd;Z="<<targZ<<",N="<<targN<<G4endl;
      sigma=0.;
    }
  }
#ifdef debug
  G4cout<<"G4ProtonNuclearCrossSection::CalcCS: sigma="<<sigma<<G4endl;
#endif
  if(sigma<0.) return 0.;
  return sigma;
}

// Linear fit for YN[N] tabulated (from X0 with fixed step DX) function to X point

// Calculate the functions for the log(A)
G4int G4QProtonNuclearCrossSection::GetFunctions(G4double a, G4double* y, G4double* z)
{
  static const G4int nLA=1;           // A#of Low Energy basic nuclei
  static const G4double LA[nLA]={1.};     
  static const G4int nL=105;          // A#of LE points in P (each MeV/C from 0 to 104)
  static const G4int nHA=1;           // A#of High Energy basic nuclei
  static const G4double HA[nHA]={1.};
  static const G4int nH=224;          // A#of HE points in lnE
  // If the cross section approximation formula is changed - replace from file.
  static const G4double SL0[nL]={
    7.094260e-1,1.532987e+0,2.449381e+0,2.785790e+0,2.525673e+0,2.128172e+0,1.780549e+0,
    1.506934e+0,1.294560e+0,1.128048e+0,9.953850e-1,8.879274e-1,7.995356e-1,7.258111e-1,
    6.635555e-1,6.104038e-1,5.645786e-1,5.247229e-1,4.897864e-1,4.589445e-1,4.315429e-1,
    4.070560e-1,3.850576e-1,3.651990e-1,3.471920e-1,3.307971e-1,3.158133e-1,3.020711e-1,
    2.894266e-1,2.777569e-1,2.669563e-1,2.569336e-1,2.476099e-1,2.389161e-1,2.307920e-1,
    2.231848e-1,2.160475e-1,2.093390e-1,2.030225e-1,1.970653e-1,1.914383e-1,1.861152e-1,
    1.810725e-1,1.762891e-1,1.717459e-1,1.674254e-1,1.633120e-1,1.593914e-1,1.556505e-1,
    1.520775e-1,1.486616e-1,1.453926e-1,1.422615e-1,1.392599e-1,1.363800e-1,1.336147e-1,
    1.309573e-1,1.284017e-1,1.259423e-1,1.235738e-1,1.212914e-1,1.190904e-1,1.169666e-1,
    1.149161e-1,1.129353e-1,1.110206e-1,1.091688e-1,1.073770e-1,1.056423e-1,1.039619e-1,
    1.023336e-1,1.007548e-1,9.922335e-2,9.773724e-2,9.629446e-2,9.489316e-2,9.353161e-2,
    9.220814e-2,9.092120e-2,8.966931e-2,8.845106e-2,8.726514e-2,8.611027e-2,8.498527e-2,
    8.388900e-2,8.282039e-2,8.177841e-2,8.076208e-2,7.977047e-2,7.880271e-2,7.785794e-2,
    7.693536e-2,7.603421e-2,7.515376e-2,7.429330e-2,7.345216e-2,7.262971e-2,7.182534e-2,
    7.103847e-2,7.026852e-2,6.951498e-2,6.877732e-2,6.805505e-2,6.734772e-2,6.665486e-2};
  static const G4double SH0[nH]={
    1.718841e-5,1.912141e-5,2.128656e-5,2.372770e-5,2.651339e-5,2.976162e-5,3.369201e-5,
    3.873597e-5,4.577051e-5,5.661516e-5,7.508997e-5,1.092699e-4,1.762839e-4,3.124886e-4,
    5.948094e-4,1.184449e-3,2.411855e-3,4.923726e-3,9.871386e-3,1.894320e-2,3.373152e-2,
    5.419455e-2,7.777948e-2,1.011811e-1,1.227807e-1,1.428966e-1,1.626818e-1,1.833195e-1,
    2.057743e-1,2.307930e-1,2.589428e-1,2.906090e-1,3.259289e-1,3.646554e-1,4.059556e-1,
    4.481828e-1,4.887166e-1,5.240358e-1,5.501959e-1,5.637401e-1,5.627614e-1,5.475832e-1,
    5.206446e-1,4.856647e-1,4.465759e-1,4.067172e-1,3.684796e-1,3.333189e-1,3.019524e-1,
    2.745971e-1,2.511726e-1,2.314485e-1,2.151395e-1,2.019637e-1,1.916740e-1,1.840748e-1,
    1.790291e-1,1.764601e-1,1.763488e-1,1.787259e-1,1.836564e-1,1.912090e-1,2.014025e-1,
    2.141163e-1,2.289594e-1,2.451064e-1,2.611598e-1,2.751583e-1,2.848795e-1,2.884723e-1,
    2.851743e-1,2.756664e-1,2.618121e-1,2.459864e-1,2.304469e-1,2.170242e-1,2.071089e-1,
    2.017331e-1,2.014838e-1,2.059886e-1,2.130499e-1,2.185478e-1,2.186039e-1,2.124513e-1,
    2.023557e-1,1.911989e-1,1.808918e-1,1.722630e-1,1.654744e-1,1.603770e-1,1.567046e-1,
    1.541608e-1,1.524546e-1,1.513189e-1,1.505256e-1,1.498980e-1,1.493175e-1,1.487199e-1,
    1.480828e-1,1.474096e-1,1.467148e-1,1.460147e-1,1.453221e-1,1.446452e-1,1.439881e-1,
    1.433514e-1,1.427339e-1,1.421336e-1,1.415477e-1,1.409739e-1,1.404099e-1,1.398539e-1,
    1.393046e-1,1.387609e-1,1.382221e-1,1.376879e-1,1.371581e-1,1.366326e-1,1.361116e-1,
    1.355952e-1,1.350837e-1,1.345775e-1,1.340767e-1,1.335816e-1,1.330926e-1,1.326099e-1,
    1.321338e-1,1.316644e-1,1.312019e-1,1.307465e-1,1.302983e-1,1.298574e-1,1.294239e-1,
    1.289978e-1,1.285792e-1,1.281681e-1,1.277645e-1,1.273684e-1,1.269797e-1,1.265984e-1,
    1.262246e-1,1.258580e-1,1.254987e-1,1.251465e-1,1.248015e-1,1.244635e-1,1.241324e-1,
    1.238082e-1,1.234908e-1,1.231801e-1,1.228760e-1,1.225784e-1,1.222872e-1,1.220024e-1,
    1.217239e-1,1.214515e-1,1.211852e-1,1.209249e-1,1.206706e-1,1.204221e-1,1.201793e-1,
    1.199423e-1,1.197109e-1,1.194850e-1,1.192646e-1,1.190497e-1,1.188400e-1,1.186357e-1,
    1.184365e-1,1.182425e-1,1.180536e-1,1.178697e-1,1.176908e-1,1.175169e-1,1.173477e-1,
    1.171834e-1,1.170239e-1,1.168690e-1,1.167189e-1,1.165733e-1,1.164323e-1,1.162959e-1,
    1.161639e-1,1.160364e-1,1.159132e-1,1.157944e-1,1.156800e-1,1.155698e-1,1.154639e-1,
    1.153622e-1,1.152646e-1,1.151712e-1,1.150819e-1,1.149967e-1,1.149155e-1,1.148384e-1,
    1.147652e-1,1.146960e-1,1.146307e-1,1.145693e-1,1.145118e-1,1.144581e-1,1.144082e-1,
    1.143621e-1,1.143198e-1,1.142812e-1,1.142464e-1,1.142152e-1,1.141877e-1,1.141639e-1,
    1.141437e-1,1.141271e-1,1.141140e-1,1.141046e-1,1.140986e-1,1.140962e-1,1.140973e-1,
    1.141019e-1,1.141099e-1,1.141214e-1,1.141363e-1,1.141546e-1,1.141763e-1,1.142013e-1};
  static const G4double* SL[nLA]={SL0};
  static const G4double* SH[nHA]={SH0};
  if(a<=.9)
  {
    G4cout<<"***G4QProtonNuclearCS::GetFunctions: A="<<a<<"(?). No CS returned!"<<G4endl;
    return -1;
  }
  G4int r=0;                            // Low channel for LEN (filling-flag for LEN)
  for(G4int i=0; i<nLA; i++) if(std::fabs(a-LA[i])<.0005)
  {
    for(G4int k=0; k<nL; k++) y[k]=SL[i][k];
    r=1;                                // Flag of filled LEN part 
  }
  G4int h=0;
  for(G4int j=0; j<nHA; j++) if(std::fabs(a-HA[j])<.0005)
  {
    for(G4int k=0; k<nH; k++) z[k]=SH[j][k];
    h=1;                                // Flag of filled LEN part 
  }
  if(!r)                                // LEN part is not filled
  {
    G4int k=0;                          // !! To be good for different compilers !!
    for(k=1; k<nLA; k++) if(a<LA[k]) break;
    if(k<1) k=1;                        // Extrapolation from the first bin (D/He)
    if(k>=nLA) k=nLA-1;                 // Extrapolation from the last bin (U)
    G4int     k1=k-1;
    G4double  xi=LA[k1];
    G4double   b=(a-xi)/(LA[k]-xi);
    for(G4int m=0; m<nL; m++)
    {
      if(a>1.5)
      {
        G4double yi=SL[k1][m];
        y[m]=yi+(SL[k][m]-yi)*b;
#ifdef debugs
        if(y[m]<0.)G4cout<<"G4QProtonNucleCS::GetF: y="<<y[m]<<",k="<<k<<",yi="<<yi<<",ya="
                         <<SL[k][m]<<",b="<<b<<",xi="<<xi<<",xa="<<LA[k]<<",a="<<a<<G4endl;
#endif
	  }
      else y[m]=0.;
    }
    r=1;
  }
  if(!h)                                // High Energy part is not filled
  {
    G4int k=0;
    for(k=1; k<nHA; k++) if(a<HA[k]) break;
    if(k<1) k=1;                        // Extrapolation from the first bin (D/He)
    if(k>=nHA) k=nHA-1;                 // Extrapolation from the last bin (Pu)
    G4int     k1=k-1;
    G4double  xi=HA[k1];
    G4double   b=(a-xi)/(HA[k]-xi);
    for(G4int m=0; m<nH; m++)
    {
      G4double zi=SH[k1][m];
      z[m]=zi+(SH[k][m]-zi)*b;
    }
    h=1;
  }
  return r*h;
}
