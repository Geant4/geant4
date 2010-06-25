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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Physics class: G4QKaonZeroNuclearCrossSection for gamma+A cross sections
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-03
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Feb-04
// --------------------------------------------------------------------------------
// ****************************************************************************************
// This Header is a part of the CHIPS physics package (author: M. Kosov)
// ****************************************************************************************
// Short description: CHIPS cross-sections for K0-nuclear interactions [K0=(K+ + K-)/2]
// -------------------------------------------------------------------------------------
//
//#define debug
//#define pdebug
//#define debug3
//#define debugn
//#define debugs

#include "G4QKaonZeroNuclearCrossSection.hh"

// Initialization of the
G4double* G4QKaonZeroNuclearCrossSection::lastLEN=0;// Pointer to the lastArray of LowEn CS
G4double* G4QKaonZeroNuclearCrossSection::lastHEN=0;// Pointer to the lastArray of HighEnCS
G4int     G4QKaonZeroNuclearCrossSection::lastN=0;  // The last N of calculated nucleus
G4int     G4QKaonZeroNuclearCrossSection::lastZ=0;  // The last Z of calculated nucleus
G4double  G4QKaonZeroNuclearCrossSection::lastP=0.; // Last used in cross section Momentum
G4double  G4QKaonZeroNuclearCrossSection::lastTH=0.;// Last threshold momentum
G4double  G4QKaonZeroNuclearCrossSection::lastCS=0.;// Last value of the Cross Section
G4int     G4QKaonZeroNuclearCrossSection::lastI=0;  // The last position in the DAMDB
G4VQCrossSection* G4QKaonZeroNuclearCrossSection::theKMinusCS =
                                             G4QKaonMinusNuclearCrossSection::GetPointer();
G4VQCrossSection*  G4QKaonZeroNuclearCrossSection::theKPlusCS  =
                                              G4QKaonPlusNuclearCrossSection::GetPointer();

// Returns Pointer to the G4VQCrossSection class
G4VQCrossSection* G4QKaonZeroNuclearCrossSection::GetPointer()
{
  static G4QKaonZeroNuclearCrossSection theCrossSection; //**Static body of Cross Section**
  return &theCrossSection;
}

// The main member function giving the collision cross section (P is in IU, CS is in mb)
// Make pMom in independent units ! (Now it is MeV)
G4double G4QKaonZeroNuclearCrossSection::GetCrossSection(G4bool fCS, G4double pMom,
                                                       G4int tgZ, G4int tgN, G4int PDG)
{
#ifdef debug
  G4cout<<"G4QKZCS::GetCS:>>> f="<<fCS<<", p="<<pMom<<", Z="<<tgZ<<"("<<lastZ<<") ,N="<<tgN
        <<"("<<lastN<<"), PDG=130/310, thresh="<<lastTH<<",Sz="<<colN.size()<<G4endl;
#endif
  if(PDG!=130 && PDG!=310 && PDG!=311 && PDG!=-311)
                 G4cout<<"-Warning-G4QKaonZeroCS::GetCS:***Not a K0***, PDG="<<PDG<<G4endl;
  G4double CS=(theKMinusCS->GetCrossSection(fCS,pMom,tgZ,tgN,-321)
              +theKPlusCS->GetCrossSection(fCS,pMom,tgZ,tgN,321))/2;
#ifdef debug
  G4cout<<"==>G4QKZCS::GetCroSec: P="<<pMom<<"(MeV),CS="<<CS<<"(mb)"<<G4endl;
#endif
  return CS;
}

// A fake function (never called) giving the K0-A cross section (Mom in GeV, CS in mb)
G4double G4QKaonZeroNuclearCrossSection::CalculateCrossSection(G4bool, G4int, G4int, G4int,
                                                               G4int, G4int, G4double)
{
  G4cout<<"-Warning-G4QKaonZeroCS::CalcCS:*A fake function is called, returns 0**"<<G4endl;
  return 0.;
} // It is kept because this is a pure virtual function of the G4VQCrossSection interface
