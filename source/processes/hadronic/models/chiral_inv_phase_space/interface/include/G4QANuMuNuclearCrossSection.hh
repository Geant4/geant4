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
// $Id: G4QANuMuNuclearCrossSection.hh,v 1.5 2006/06/29 20:07:58 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//
// GEANT4 physics class: G4QANuMuNuclearCrossSection -- header file
// M.V. Kossov, CERN-ITEP(Moscow), 20-DEC-2005
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 20-DEC-2005
//
// Short description: this G4 singletone class calculates (nu_mu,mu) Nuclear cross section
// (Energy limit: E<320GeV->badExtrapolation) for a particular isotope (proportional to A)
// ****************************************************************************************

#ifndef G4QANuMuNuclearCrossSection_h
#define G4QANuMuNuclearCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include <vector>
#include "Randomize.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4VQCrossSection.hh"

class G4QANuMuNuclearCrossSection : public G4VQCrossSection
{
protected:

  G4QANuMuNuclearCrossSection()  {};

public:

  ~G4QANuMuNuclearCrossSection()  {};

  static G4VQCrossSection* GetPointer(); // Gives a pointer to this singletone

  G4double ThresholdEnergy(G4int Z, G4int N, G4int PDG=-14);

  G4double CalculateCrossSection(G4bool CS, G4int F, G4int I, G4int PDG, G4int Z,
                                                               G4int N, G4double Momentum);

  G4int    GetExchangePDGCode();

  G4double GetDirectPart(G4double Q2);

  G4double GetNPartons(G4double Q2);

  G4double GetQEL_ExchangeQ2();

  G4double GetNQE_ExchangeQ2();

  // Get static members
  G4double GetLastTOTCS() {return lastSig;}
  G4double GetLastQELCS() {return lastQEL;}

private:
  G4int    GetFunctions(G4int z, G4int n, G4double* t, G4double* q, G4double* e);
  G4double HighEnergyTX(G4double lE);
  G4double HighEnergyQE(G4double lE);

// Body
private:
  static G4bool    onlyCS;   // flag to calculate only CS (not QE)
  static G4double  lastSig;  // Last calculated total cross section
  static G4double  lastQEL;  // Last calculated quasi-elastic cross section
  static G4int     lastL;    // Last bin used in the cross section TheLastBin
  static G4double  lastE;    // Last energy used in the cross section Energy
  static G4double* lastEN;   // Pointer to the last array of the energy axis
  static G4double* lastTX;   // Pointer to the last array of the total CS function
  static G4double* lastQE;   // Pointer to the last array of the QE CS function
};

#endif
