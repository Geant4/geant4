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
// $Id: G4QElectroNuclearCrossSection.cc,v 1.1 2004-03-05 13:24:26 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: G4QElectroNuclearCrossSection for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-Oct-03
// 
//======================================== ================================================

///#define debug
#define edebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4QElectroNuclearCrossSection.hh"

// Initialization of the
G4double  G4QElectroNuclearCrossSection::lastE=0.;// Last used in cross section TheEnergy
G4int     G4QElectroNuclearCrossSection::lastF=0; // Last used in cross section TheFirstBin
G4double  G4QElectroNuclearCrossSection::lastG=0.;// Last used in cross section TheGamma
G4double  G4QElectroNuclearCrossSection::lastH=0.;// LastValue ofTheHighEnergy A-dependence
G4double* G4QElectroNuclearCrossSection::lastJ1=0;// Pointer to theLastArray of J1 function
G4double* G4QElectroNuclearCrossSection::lastJ2=0;// Pointer to theLastArray of J2 function
G4double* G4QElectroNuclearCrossSection::lastJ3=0;// Pointer to theLastArray of J3 function
G4int     G4QElectroNuclearCrossSection::lastL=0; // Last used in cross section TheLastBin
G4int     G4QElectroNuclearCrossSection::lastN=0; // The last N of calculated nucleus
G4int     G4QElectroNuclearCrossSection::lastZ=0; // The last Z of calculated nucleus
G4double  G4QElectroNuclearCrossSection::lastTH=0.; // Last energy threshold
G4double  G4QElectroNuclearCrossSection::lastSig=0.;// Last value of the Cross Section
G4double  G4QElectroNuclearCrossSection::ml=G4Electron::Electron()->GetPDGMass();// M_elect
G4double  G4QElectroNuclearCrossSection::ml2=
                       G4QElectroNuclearCrossSection::ml*G4QElectroNuclearCrossSection::ml;
G4double  G4QElectroNuclearCrossSection::lml=log(G4QElectroNuclearCrossSection::ml);//ln(M)

// The main member function giving the gamma-A cross section (E_kin in MeV, CS in mb)
G4double G4QElectroNuclearCrossSection::GetCrossSection(G4double Energy,
                                                     G4int targZ, G4int targN )
{
  static const G4int nE=336;    // !! If change this, change it in GetFunctions() (*.hh) !!
  static const G4int mL=nE-1;
  static const G4double EMi=2.0612;          // Minimum tabulated KineticEnergy of Electron
  static const G4double EMa=50000.;          // Maximum tabulated Energy of the Electron 
  static const G4double lEMi=log(EMi);       // Minimum tabulatedLogarithmKinEnergy of El.
  static const G4double lEMa=log(EMa);       // Maximum tabulatedLogarithmKinEnergy of El
  static const G4double dlnE=(lEMa-lEMi)/mL; // Logarithmic table step for ElectroKinEnergy
  static const G4double alop=1./137.036/3.14159265; //coeffitient for  E>50000 calculations
  // *** Begin of the Associative memory for acceleration of the cross section calculations
  static std::vector <G4int> colN;     // Vector of N for calculated nucleus (isotop)
  static std::vector <G4int> colZ;     // Vector of Z for calculated nucleus (isotop)
  static std::vector <G4int> colF;     // Vector of LastStartPosition in Ji-function tables
  static std::vector <G4double> colTH; // Vector of energy thresholds for eA->eX reactions
  static std::vector <G4double> colH;  // Vector of HighEnergyCoefficients (functional)
  static std::vector <G4double*> J1;   // Vector of pointers to the J1 tabulated functions
  static std::vector <G4double*> J2;   // Vector of pointers to the J2 tabulated functions
  static std::vector <G4double*> J3;   // Vector of pointers to the J3 tabulated functions
  // *** End of Static Definitions (Associative Memory) ***
  if (Energy<=EMi) return 0.;          // Energy is below the minimum energy in the table
  G4double A=targN+targZ;            // New A (can be different from targetAtomicNumber)
  if(targN!=lastN || targZ!=lastZ)   // This nucleus was not the last used isotop
  {
    lastE    = 0.;                   // New history in the electron Energy
    lastG    = 0.;                   // New history in the photon Energy
    lastN    = targN;                // The last N of calculated nucleus
    lastZ    = targZ;                // The last Z of calculated nucleus
    G4int n=colN.size();             // Size of the Associative Memory DB in the heap
    G4bool in=false;                 // "Found in AMDB" flag
    if(n) for(G4int i=0; i<n; i++) if(colN[i]==targN && colZ[i]==targZ) // Calcul. nuclei
	{ // The nucleus is found in AMDB
      in=true;                       // Rais the "Found in AMDB" flag
      lastTH =colTH[i];              // Last Energy threshold (A-dependent)
      lastF  =colF[i];               // Last ZeroPosition in the J-functions
      lastH  =colH[i];               // Last High Energy Coefficient (A-dependent)
      lastJ1 =J1[i];                 // Pointer to the prepared J1 function
      lastJ2 =J2[i];                 // Pointer to the prepared J2 function
      lastJ3 =J3[i];                 // Pointer to the prepared J3 function
	}
	if(!in)                          // This nucleus has not been calculated previously
	{
      lastJ1 = new G4double[nE];     // Allocate memory for the new J1 function
      lastJ2 = new G4double[nE];     // Allocate memory for the new J2 function
      lastJ3 = new G4double[nE];     // Allocate memory for the new J3 function
      lastF   = GetFunctions(A,lastJ1,lastJ2,lastJ3);//newZeroPos and J-functions filling
      lastH   = alop*A*(1.-.072*log(A));//similar to lastSP of G4PhotonuclearCrossSection
      lastTH  = G4QPhotoNuclearCrossSection().ThresholdEnergy(targZ,targN);//Last Threshold
#ifdef pdebug
      G4cout<<"G4QElectroNuclearCS::GetCrossSection: lastH="<<lastH<<",A="<<A<<G4endl;
#endif
      colN.push_back(targN);
      colZ.push_back(targZ);
      colF.push_back(lastF);
      J1.push_back(lastJ1);
      J2.push_back(lastJ2);
      J3.push_back(lastJ3);
      colH.push_back(lastH);
      colTH.push_back(lastTH);
	} // End of creation of the new set of parameters
  } // End of parameters udate
  else if(abs((lastE-Energy)/Energy)<.001) return lastSig*millibarn; // Came CS are calc.
  // ============================== NOW Calculate the Cross Section =====================
  lastE=Energy;                      // lastE - the electron energy
  if (Energy<=lastTH)                // Once more check that elKinE is higher than ThreshE
  {
    lastSig=0.;
    return 0.;
  }
  G4double lE=log(Energy);           // log(elE) (it is necessary for the fit)
  lastG=lE-lml;                     // Gamma of the electron (used to recover log(elE))
  G4double dlg1=lastG+lastG-1.;
  G4double lgoe=lastG/lastE;
  if(lE<lEMa) // Linear fit is made explicitly to fix the last bin for the randomization
  {
    G4double shift=(lE-lEMi)/dlnE;
    G4int    blast=static_cast<int>(shift);
    if(blast<0)   blast=0;
    if(blast>=mL) blast=mL-1;
    shift-=blast;
    lastL=blast+1;
    G4double YNi=dlg1*lastJ1[blast]
                -lgoe*(lastJ2[blast]+lastJ2[blast]-lastJ3[blast]/lastE);
    G4double YNj=dlg1*lastJ1[lastL]
                -lgoe*(lastJ2[lastL]+lastJ2[lastL]-lastJ3[lastL]/lastE);
    lastSig= YNi+shift*(YNj-YNi);
    if(lastSig>YNj)lastSig=YNj;
#ifdef pdebug
    G4cout<<"G4QElNCS::GCS:S="<<lastSig<<",E="<<lE<<",Yi="<<YNi<<",Yj="<<YNj<<",M="<<lEMa
          <<G4endl;
    G4cout<<"G4QElNCS::GCS:s="<<shift<<",Jb="<<lastJ1[blast]<<",J="<<lastJ1[lastL]<<",b="
          <<blast<<G4endl;
#endif
  }
  else
  {
    lastL=mL;
    G4double term1=lastJ1[mL]+lastH*HighEnergyJ1(lE);
    G4double term2=lastJ2[mL]+lastH*HighEnergyJ2(lE);
    G4double term3=lastJ3[mL]+lastH*HighEnergyJ3(lE);
    lastSig=dlg1*term1-lgoe*(term2+term2-term3/lastE);
#ifdef pdebug
    G4cout<<"G4QElNucCS::GetCrossSec:S="<<lastSig<<",lE="<<lE<<",J1="
          <<lastH*HighEnergyJ1(lE)<<",Pm="<<lastJ1[mL]<<",Fm="<<lastJ2[mL]<<",Fh="
          <<lastH*HighEnergyJ2(lE)<<",EM="<<lEMa<<G4endl;
#endif
  }
  if(lastSig<0.) lastSig = 0.;
  lastE=Energy;
  return lastSig*millibarn;
}

G4double G4QElectroNuclearCrossSection::GetEquivalentPhotonEnergy()
{
  // @@ All constants are copy of that from GetCrossSection funct. => Make them general.
  static const G4int nE=336; // !!  If change this, change it in GetFunctions() (*.hh) !!
  static const G4int mL=nE-1;
  static const G4double EMi=2.0612;          // Minimum Energy
  static const G4double EMa=50000.;          // Maximum Energy
  static const G4double lEMi=log(EMi);       // Minimum logarithmic Energy
  static const G4double lEMa=log(EMa);       // Maximum logarithmic Energy
  static const G4double dlnE=(lEMa-lEMi)/mL; // Logarithmic step in Energy
  G4double phLE=0.;                     // Prototype of the log(nu=E_gamma)
  G4double Y[nE];                       // Prepare the array for randomization
#ifdef debug
  G4cout<<"G4QElectroNuclCrossSec::GetEguPhE:B="<<lastF<<",l="<<lastL<<",1="<<lastJ1[lastL]
        <<",2="<<lastJ2[lastL]<<",3="<<lastJ3[lastL]<<",S="<<lastSig<<",E="<<lastE<<G4endl;
#endif
  G4double lastLE=lastG+lml;           // recover log(eE) from the gamma (lastG)
  G4double dlg1=lastG+lastG-1.;
  G4double lgoe=lastG/lastE;
  for(G4int i=lastF;i<=lastL;i++)
    Y[i]=dlg1*lastJ1[i]-lgoe*(lastJ2[i]+lastJ2[i]-lastJ3[i]/lastE);
  G4double ris=lastSig*G4UniformRand(); // Sig can be > Y[lastL=mL], then it is funct. reg.
#ifdef debug
  G4cout<<"G4QElectroNuclCrossSec::GetEquivalentPhotonEn: "<<ris<<",Y="<<Y[lastL]<<G4endl;
#endif
  if(ris<Y[lastL])                      // Search in the table
  {
	G4int j=lastF;
    G4double Yj=Y[j];                   // It mast be 0 (some times just very small)
    while (ris>Yj && j<lastL)           // Associative search
	{
      j++;
      Yj=Y[j];                          // High value
	}
    G4int j1=j-1;
    G4double Yi=Y[j1];                  // Low value
    phLE=lEMi+(j1+(ris-Yi)/(Yj-Yi))*dlnE;
#ifdef debug
	G4cout<<"G4QElNuclearCS::E="<<phLE<<",l="<<lEMi<<",j="<<j<<",ris="<<ris<<",Yi="<<Yi
          <<",Y="<<Yj<<G4endl;
#endif
  }
  else                                  // Search with the function
  {
    if(lastL<mL)
      G4cerr<<"**G4QEleNucCS::GetEfPhE:L="<<lastL<<",S="<<lastSig<<",Y="<<Y[lastL]<<G4endl;
    G4double f=(ris-Y[lastL])/lastH;    // ScaledResidualValue of the cross-sec. integral
#ifdef pdebug
	G4cout<<"G4QElNuCS::GetEfPhE:HighEnergy f="<<f<<",ris="<<ris<<",lastH="<<lastH<<G4endl;
#endif
    phLE=SolveTheEquation(f);           // Solve equation to find theLog(phE) (comp lastLE)
#ifdef pdebug
	G4cout<<"G4QEleNucCS::GetEfPhE:HighEnergy lphE="<<phLE<<G4endl;
#endif
  }
  if(phLE>lastLE)
  {
    G4cerr<<"***G4QElNuclearCS::GetEquPhotE:N="<<lastN<<",Z="<<lastZ<<",lpE"<<phLE<<">leE"
          <<lastLE<<",Sig="<<lastSig<<",rndSig="<<ris<<",Beg="<<lastF<<",End="<<lastL
          <<",Y="<<Y[lastL]<<G4endl;
    if(lastLE<7.2) phLE=log(exp(lastLE)-105.6584);
    else phLE=7.;
  }
  return exp(phLE);
}

G4double G4QElectroNuclearCrossSection::SolveTheEquation(G4double f)
{
  // This parameters must correspond to the G4PhotonuclearCrossSec::GetCrossSec parameters
  static const G4double shd=1.0734;                    // HE PomShadowing(D)
  static const G4double poc=0.0375;                    // HE Pomeron coefficient
  static const G4double pos=16.5;                      // HE Pomeron shift
  static const G4double reg=.11;                       // HE Reggeon slope
  static const G4double EMa=50000.;                    // Maximum Energy
  static const G4double z=log(EMa);                    // Initial argument
  static const G4double p=poc*(z-pos)+shd*exp(-reg*z); // CrossX on theHighTabEdge (small)
  static const G4int    imax=27;   // Not more than "imax" steps to find the solution
  static const G4double eps=0.001; // Accuracy which satisfies the search
  G4double lastLE=lastG+lml;                 // recover log(eE) from the gamma (lastG)
  G4double topLim=lastLE-.001;                // maximum log(phE) for equivalent photons
  G4double rE=EMa/exp(lastLE);                // r=EMa/Eel to make the firs guess
  G4double x=z+f/p/(lastG*(2.-rE*(2.-rE))-1.);// First guess (the first step from the edge)
#ifdef pdebug
  G4cout<<"G4QElNucCS::SolveTheEq: e="<<eps<<",f="<<f<<",z="<<z<<",p="<<p<<",lastG="<<lastG
        <<",x="<<x<<G4endl;
#endif
  if(x>topLim) x=topLim;
  for(G4int i=0; i<imax; i++)
  {
    G4double fx=Fun(x);
    G4double df=DFun(x);
    G4double d=(f-fx)/df;
    x=x+d;
#ifdef pdebug
    G4cout<<"G4QElNCS::SolveTheE:#"<<i<<",d="<<d<<",x="<<x<<",fx="<<fx<<",df="<<df<<G4endl;
#endif
    if(x>=lastLE)
	{
      G4cerr<<"**G4QElNCS::SolveTheEq:*Correction*"<<i<<",d="<<d<<",x="<<x<<">lE="<<lastLE
            <<",f="<<f<<",fx="<<fx<<",df="<<df<<",A(Z="<<lastZ<<",N="<<lastN<<")"<<G4endl;
      x=topLim;
      if(i)G4Exception("G4QElectroNuclearCrSec::SolveTheEq()","009",FatalException,"E>eE");
    }
    if(abs(d)<eps) break;
    if(i+1>=imax) G4cerr<<"G4QElNucCS::SolveTheE:"<<i+2<<">"<<imax<<"->Use bigMax. ln(eE)="
                        <<lastLE<<",Z="<<lastZ<<", N="<<lastN<<G4endl;
  }
  return x;
}

// Randomize Q2 for the scattered electron
G4double G4QElectroNuclearCrossSection::GetEquivalentPhotonQ2(G4double nu)
{
  G4double y=nu/lastE;                  // Part of energy carried by the equivalent pfoton
  if(y>=1.-1./(lastG+lastG)) return 0.; // The region where the method does not work
  G4double y2=y*y;                      // Squared photonic part of energy
  G4double ye=1.-y;                     // Part of energy carried by the secondary electron
  G4double Qi2=ml2*y2/ye;               // Minimum Q2
  G4double Qa2=4*lastE*lastE*ye;        // Maximum Q2
  G4double iar=Qi2/Qa2;                 // Q2min/Q2max ratio
  G4double Dy=ye+.5*y2;                 // D(y) function
  G4double Py=ye/Dy;                    // P(y) function
  G4double ePy=1.-exp(Py);              // 1-exp(P(y)) part
  G4double Uy=Py*(1.-iar);              // U(y) function
  G4double Fy=(ye+ye)*(1.+ye)*iar/y2;   // F(y) function
  G4double fr=iar/(1.-ePy*iar);         // Q-fraction
  if(Fy<=-fr)
  {
#ifdef edebug
    G4cerr<<"***G4QElectroNucCS::GetEqPhQ2:F="<<Fy<<"+fr="<<fr<<" <0"<<",ia="<<iar<<G4endl;
#endif
    return 0.;
  }    
  G4double LyQa2=log(Fy+fr);            // L(y,Q2max) function
  G4bool cond=true;
  G4int maxTry=3;
  G4int cntTry=0;
  G4double Q2=Qi2;
  while(cond&&cntTry<maxTry)            // The loop to avoid x>1.
  {
    G4double R=G4UniformRand();         // Random number (0,1)
    Q2=Qi2*(ePy+1./(exp(R*LyQa2-(1.-R)*Uy)-Fy));
    cntTry++;
    cond = Q2>1878.*nu;
  }
  if(Q2<Qi2)
  {
#ifdef edebug
    G4cerr<<"***G4QElectroNucCrossSec::GetEquPhQ2:Q2="<<Q2<<" < Q2min="<<Qi2<<G4endl;
#endif
    return Qi2;
  }  
  if(Q2>Qa2)
  {
#ifdef edebug
    G4cerr<<"***G4QElectroNucCrossSec::GetEquPhQ2:Q2="<<Q2<<" > Q2max="<<Qi2<<G4endl;
#endif
    return Qa2;
  }  
  return Q2;
}
