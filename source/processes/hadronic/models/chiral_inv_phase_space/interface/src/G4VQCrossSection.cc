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
// $Id: G4VQCrossSection.cc,v 1.10 2006/06/29 20:08:53 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//
// CHIPS virtual class: G4VQCrossSection for the collision cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-04
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 27-Nov-04
// 
//=========================================================================================
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************

//#define debug
#define edebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4VQCrossSection.hh"

// Initialization of the
G4int     G4VQCrossSection::lastPDG=0; // The last PDG code of the projectile
G4int     G4VQCrossSection::lastN=0;   // The last N of calculated nucleus
G4int     G4VQCrossSection::lastZ=0;   // The last Z of calculated nucleus
G4double  G4VQCrossSection::lastTH=0.; // Last threshold momentum
G4double  G4VQCrossSection::lastCS=0.;// Last value of the Cross Section
G4double  G4VQCrossSection::lastP=0.;  // Last used in cross section TheMomentum
//G4int     G4VQCrossSection::lastN1=0;  // Last used in cross section TheNumOfBin1
//G4int     G4VQCrossSection::lastF1=0;  // Last used in cross section TheFirstBin1
//G4int     G4VQCrossSection::lastL1=0;  // Last used in cross section TheLastBin1
//G4int     G4VQCrossSection::lastN2=0;  // Last used in cross section TheNumOfBin1
//G4int     G4VQCrossSection::lastF2=0;  // Last used in cross section TheFirstBin1
//G4int     G4VQCrossSection::lastL2=0;  // Last used in cross section TheLastBin1
//G4double  G4VQCrossSection::lastBP=0.; // Last value of the Boundary Momentum
//G4double  G4VQCrossSection::lastMP=0.; // Last value of the Maximum Momentum

G4int     G4VQCrossSection::lastI=0;          // The last position in the DAMDB
G4double  G4VQCrossSection::tolerance=.001;   // The last position in the DAMDB

// Set the new tolerance (abs(p_old/p_new-1)<tolerance)
void G4VQCrossSection::setTolerance(G4double tol)
//   ============================================
{
		tolerance=tol;
}

// Gives the threshold energy for different isotopes (can be improved in the derived class)
G4double G4VQCrossSection::ThresholdEnergy(G4int , G4int, G4int) {return 0.;} // Fake use

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4VQCrossSection::GetCrossSection(G4bool fCS, G4double pMom, G4int tgZ, G4int tgN,
                                                                                G4int pPDG)
{
  static const G4double mtu=1777.;     // Mass of a tau lepton in MeV
  static const G4double mtu2=mtu*mtu;  // Squared Mass of a tau-lepton in MeV^2
  static const G4double mmu=105.65839; // Mass of the muon in MeV
  static const G4double mmu2=mmu*mmu;  // Squared Mass of muon in MeV^2
  static const G4double mel=0.5109989; // Mass of the electron in MeV
  static const G4double mel2=mel*mel;  // Squared Mass of the electron in MeV
  static G4int j;                      // A#0f records found in DB for this projectile
  static std::vector <G4int>    colPDG;// Vector of the projectile PDG code
  static std::vector <G4int>    colN;  // Vector of N for calculated nuclei (isotops)
  static std::vector <G4int>    colZ;  // Vector of Z for calculated nuclei (isotops)
  static std::vector <G4double> colP;  // Vector of last momenta for the reaction
  static std::vector <G4double> colTH; // Vector of energy thresholds for the reaction
  static std::vector <G4double> colCS; // Vector of last cross sections for the reaction
  // ***---*** End of the mandatory Static Definitions of the Associative Memory ***---***
  G4double pEn=pMom;
  G4int apPDG=std::abs(pPDG);
  // @@ if the threshold exists for other particles, then p->T must be genergal (p=0->T=o)
  if     (apPDG==11) pEn=std::sqrt(pMom*pMom+mel2)-mel; // ==> electron/positron kinEnergy
  else if(apPDG==13) pEn=std::sqrt(pMom*pMom+mmu2)-mmu; // ==> mu-/mu+ kinEnergy
  else if(apPDG==15) pEn=std::sqrt(pMom*pMom+mtu2)-mtu; // ==> tau-/tau+ kinEnergy
#ifdef pdebug
  G4cout<<"G4VQCS::GetCS:>>>> f="<<fCS<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="<<tgN
        <<"("<<lastN<<"),PDG="<<pPDG<<"("<<lastPDG<<"), T="<<pEn<<"("<<lastTH<<")"<<",Sz="
        <<colN.size()<<G4endl;
		//CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  if(!pPDG)
  {
#ifdef pdebug
    G4cout<<"G4VQCS::GetCS: *** Found pPDG="<<pPDG<<" ====> CS=0"<<G4endl;
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
        G4cout<<"G4VQCS::GetCS: *Found* P="<<pMom<<",Threshold="<<lastTH<<",j="<<j<<G4endl;
        //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(pEn<=lastTH)
        {
#ifdef pdebug
          G4cout<<"G4VQCS::GetCS: Found T="<<pEn<<" < Threshold="<<lastTH<<",CS=0"<<G4endl;
          //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
          return 0.;                     // Energy is below the Threshold value
        }
        lastP  =colP [i];                // Last Momentum  (A-dependent)
        lastCS =colCS[i];                // Last CrossSect (A-dependent)
        if(std::fabs(lastP/pMom-1.)<tolerance)
        {
#ifdef pdebug
          G4cout<<"G4VQCS::GetCS:P="<<pMom<<"=Po="<<pMom<<",CS="<<lastCS*millibarn<<G4endl;
          //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
          return lastCS*millibarn;     // Use theLastCS
        }
        in = true;                       // This is the case when the isotop is found in DB
        // Momentum pMom is in IU ! @@ Units
#ifdef pdebug
        G4cout<<"G4VQCS::G:UpdateDB P="<<pMom<<",f="<<fCS<<",lI="<<lastI<<",j="<<j<<G4endl;
#endif
        lastCS=CalculateCrossSection(fCS,-1,j,lastPDG,lastZ,lastN,pMom); // read & update
#ifdef pdebug
        G4cout<<"G4VQCS::GetCrosSec: *****> New (inDB) Calculated CS="<<lastCS<<G4endl;
        //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
        if(lastCS<=0. && pEn>lastTH)    // Correct the threshold
        {
#ifdef pdebug
          G4cout<<"G4VQCS::GetCS: New T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
        break;                           // Go out of the LOOP
      }
#ifdef pdebug
      G4cout<<"---G4VQCrossSec::GetCrosSec:pPDG="<<pPDG<<",j="<<j<<",N="<<colN[i]
            <<",Z["<<i<<"]="<<colZ[i]<<",cPDG="<<colPDG[i]<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      j++;                             // Increment a#0f records found in DB for this pPDG
	   }
	   if(!in)                            // This nucleus has not been calculated previously
	   {
#ifdef pdebug
      G4cout<<"G4VQCS::GetCrosSec: CalcNew P="<<pMom<<",f="<<fCS<<",lastI="<<lastI<<G4endl;
#endif
      //!!The slave functions must provide cross-sections in millibarns (mb) !! (not in IU)
      lastCS=CalculateCrossSection(fCS,0,j,lastPDG,lastZ,lastN,pMom); //calculate & create
      if(lastCS<=0.)
						{
        lastTH = ThresholdEnergy(tgZ, tgN); // The Threshold Energy which is now the last
#ifdef pdebug
        G4cout<<"G4VQCrossSection::GetCrossSection:NewThresh="<<lastTH<<",T="<<pEn<<G4endl;
#endif
        if(pEn>lastTH)
        {
#ifdef pdebug
          G4cout<<"G4VQCS::GetCS: First T="<<pEn<<"(CS=0) > Threshold="<<lastTH<<G4endl;
#endif
          lastTH=pEn;
        }
						}
#ifdef pdebug
      G4cout<<"G4VQCS::GetCrosSec: New CS="<<lastCS<<",lZ="<<lastN<<",lN="<<lastZ<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      colN.push_back(tgN);
      colZ.push_back(tgZ);
      colPDG.push_back(pPDG);
      colP.push_back(pMom);
      colTH.push_back(lastTH);
      colCS.push_back(lastCS);
#ifdef pdebug
      G4cout<<"G4VQCS::GetCS:1st, P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
      //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
      return lastCS*millibarn;
	   } // End of creation of the new set of parameters
    else
				{
#ifdef pdebug
      G4cout<<"G4VQCS::GetCS: Update lastI="<<lastI<<",j="<<j<<G4endl;
#endif
      colP[lastI]=pMom;
      colPDG[lastI]=pPDG;
      colCS[lastI]=lastCS;
    }
  } // End of parameters udate
  else if(pEn<=lastTH)
  {
#ifdef pdebug
    G4cout<<"G4VQCS::GetCS: Current T="<<pEn<<" < Threshold="<<lastTH<<", CS=0"<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return 0.;                         // Momentum is below the Threshold Value -> CS=0
  }
  else if(std::fabs(lastP/pMom-1.)<tolerance)
  {
#ifdef pdebug
    G4cout<<"G4VQCS::GetCS: OldCur P="<<pMom<<"="<<pMom<<", CS="<<lastCS*millibarn<<G4endl;
    //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
    return lastCS*millibarn;     // Use theLastCS
  }
  else
  {
#ifdef pdebug
    G4cout<<"G4VQCS::GetCS:UpdateCur P="<<pMom<<",f="<<fCS<<",I="<<lastI<<",j="<<j<<G4endl;
#endif
    lastCS=CalculateCrossSection(fCS,1,j,lastPDG,lastZ,lastN,pMom); // Only UpdateDB
    lastP=pMom;
  }
#ifdef pdebug
  G4cout<<"G4VQCS::GetCrosSec:End,P="<<pMom<<"(MeV),CS="<<lastCS*millibarn<<"(mb)"<<G4endl;
  //CalculateCrossSection(fCS,-27,j,lastPDG,lastZ,lastN,pMom); // DUMMY TEST
#endif
  return lastCS*millibarn;
}

G4double G4VQCrossSection::GetDirectPart(G4double) {return 0.;} // Direct interaction

G4double G4VQCrossSection::GetNPartons(G4double) {return 3.;} // Direct interaction

G4double G4VQCrossSection::GetLastTOTCS() {return 0.;} // Get the last total CS

G4double G4VQCrossSection::GetLastQELCS() {return 0.;} // Get the last quasi-elast CS

G4double G4VQCrossSection::GetExchangeEnergy() {return 0.;}

G4double G4VQCrossSection::GetExchangeQ2(G4double) {return 0.;}

G4double G4VQCrossSection::GetExchangeT(G4int,G4int,G4int) {return 0.;}

G4double G4VQCrossSection::GetQEL_ExchangeQ2() {return 0.;}

G4double G4VQCrossSection::GetNQE_ExchangeQ2() {return 0.;}

G4int G4VQCrossSection::GetExchangePDGCode() {return 0;}

G4double G4VQCrossSection::GetVirtualFactor(G4double nu, G4double Q2) {return 0.*nu*Q2;}

// This function finds  the linear approximation Y-point for the XN(N), YN(N) table
G4double G4VQCrossSection::LinearFit(G4double X, G4int N, G4double* XN, G4double* YN)
{
  G4double Xj=XN[0];
  G4double Xh=XN[N-1];
  if(X<=Xj) return YN[0]; 
  else if(X>=Xh) return YN[N-1];
  G4double Xp=0.; G4int j=0; while (X>Xj && j<N) {j++; Xp=Xj; Xj=XN[j];}
  return YN[j]-(Xj-X)*(YN[j]-YN[j-1])/(Xj-Xp);
}

// This function finds the linear approximation Y-point for equidistant bins: XI=X0+I*DX
G4double G4VQCrossSection::EquLinearFit(G4double X, G4int N, G4double X0, G4double DX,
                                        G4double* Y)
{
#ifdef pdebug
		G4cout<<"G4VQCrossSection::EquLinearFit: ***Called*** X="<<DX<<", N="<<N<<", X0="<<X0
        <<", DX="<<DX<<G4endl;
		G4cout<<"G4VQCrossSection::EquLinearFit: Y[0]="<<Y[0]<<", Y[N-1]="<<Y[N-1]<<G4endl;
#endif
  if(DX<=0. || N<2)
  {
    G4cerr<<"***G4VQCrossSection::EquLinearFit: DX="<<DX<<", N="<<N<<G4endl;
    return Y[0];
  }
  G4int    N2=N-2;
  G4double d=(X-X0)/DX;
  G4int         j=static_cast<int>(d);
  if     (j<0)  j=0;
  else if(j>N2) j=N2;
  d-=j; // excess
  G4double yi=Y[j];
  G4double sigma=yi+(Y[j+1]-yi)*d;
  return sigma;
}
