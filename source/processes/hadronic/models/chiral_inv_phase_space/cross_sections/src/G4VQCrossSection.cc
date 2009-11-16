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
// $Id: G4VQCrossSection.cc,v 1.1 2009-11-16 18:15:43 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// Short description: a basic class for all CHIPS reaction cross-sections.
// -----------------------------------------------------------------------

//#define debug
#define edebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4VQCrossSection.hh"

// Initialization of the
G4double  G4VQCrossSection::tolerance=.001;   // The relative tolarence for the same CrSec

// Gives the threshold energy for different isotopes (can be improved in the derived class)
G4double G4VQCrossSection::ThresholdEnergy(G4int , G4int, G4int) {return 0.;} // Fake use

G4double G4VQCrossSection::GetDirectPart(G4double) {return 0.;} // Direct interaction

G4double G4VQCrossSection::GetNPartons(G4double) {return 3.;} // Direct interaction

G4double G4VQCrossSection::GetLastTOTCS() {return 0.;} // Get the last total CS

G4double G4VQCrossSection::GetLastQELCS() {return 0.;} // Get the last quasi-elast CS

G4double G4VQCrossSection::GetExchangeEnergy() {return 0.;}

G4double G4VQCrossSection::GetExchangeQ2(G4double) {return 0.;}

G4double G4VQCrossSection::GetSlope(G4int,G4int,G4int) {return 0.;}

G4double G4VQCrossSection::GetExchangeT(G4int,G4int,G4int) {return 0.;}

G4double G4VQCrossSection::GetHMaxT() {return 0.;}

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
  G4cout<<"G4VQCrossSection::EquLinearFit: ***Called*** X="<<X<<", N="<<N<<", X0="<<X0
        <<", DX="<<DX<<G4endl;
  G4cout<<"G4VQCrossSection::EquLinearFit: Y[0]="<<Y[0]<<", Y[N-1]="<<Y[N-1]<<G4endl;
  //for(G4int i=1; i<N; i++)G4cout<<"-----G4VQCS::EquLinearFit: Y["<<i<<"]="<<Y[i]<<G4endl;
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
#ifdef pdebug
  G4cout<<"G4VQCrossSection::EquLinearFit: CS="<<sigma<<G4endl;
#endif
  return sigma;
}
