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
// $Id: G4VCrossSection.cc,v 1.1 2009-11-16 18:15:43 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// CHIPS virtual class: G4VCrossSection for the collision cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-04
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 27-Nov-04
// 
// This class has been extracted from the CHIPS model. 
// All the dependencies on CHIPS classes have been removed.
//
// Short description: a basic class for all CHIPS reaction cross-sections.
// -----------------------------------------------------------------------

//#define debug
#define edebug
//#define pdebug
//#define ppdebug
//#define tdebug
//#define sdebug

#include "G4VCrossSection.hh"

// Initialization of the
G4double  G4VCrossSection::tolerance=.001;   // The relative tolarence for the same CrSec

// Gives the threshold energy for different isotopes (can be improved in the derived class)
G4double G4VCrossSection::ThresholdEnergy(G4int , G4int, G4int) {return 0.;} // Fake use

G4double G4VCrossSection::GetDirectPart(G4double) {return 0.;} // Direct interaction

G4double G4VCrossSection::GetNPartons(G4double) {return 3.;} // Direct interaction

G4double G4VCrossSection::GetLastTOTCS() {return 0.;} // Get the last total CS

G4double G4VCrossSection::GetLastQELCS() {return 0.;} // Get the last quasi-elast CS

G4double G4VCrossSection::GetExchangeEnergy() {return 0.;}

G4double G4VCrossSection::GetExchangeQ2(G4double) {return 0.;}

G4double G4VCrossSection::GetSlope(G4int,G4int,G4int) {return 0.;}

G4double G4VCrossSection::GetExchangeT(G4int,G4int,G4int) {return 0.;}

G4double G4VCrossSection::GetHMaxT() {return 0.;}

G4double G4VCrossSection::GetQEL_ExchangeQ2() {return 0.;}

G4double G4VCrossSection::GetNQE_ExchangeQ2() {return 0.;}

G4int G4VCrossSection::GetExchangePDGCode() {return 0;}

G4double G4VCrossSection::GetVirtualFactor(G4double nu, G4double Q2) {return 0.*nu*Q2;}

// This function finds  the linear approximation Y-point for the XN(N), YN(N) table
G4double G4VCrossSection::LinearFit(G4double X, G4int N, G4double* XN, G4double* YN)
{
  G4double Xj=XN[0];
  G4double Xh=XN[N-1];
  if(X<=Xj) return YN[0]; 
  else if(X>=Xh) return YN[N-1];
  G4double Xp=0.; G4int j=0; while (X>Xj && j<N) {j++; Xp=Xj; Xj=XN[j];}
  return YN[j]-(Xj-X)*(YN[j]-YN[j-1])/(Xj-Xp);
}

// This function finds the linear approximation Y-point for equidistant bins: XI=X0+I*DX
G4double G4VCrossSection::EquLinearFit(G4double X, G4int N, G4double X0, G4double DX,
                                        G4double* Y)
{
  if(DX<=0. || N<2)
  {
    G4cerr<<"***G4VCrossSection::EquLinearFit: DX="<<DX<<", N="<<N<<G4endl;
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
