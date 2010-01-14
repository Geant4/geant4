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
// Short description: CHIPS cross-sectons for Ion-Ion interactions
// ---------------------------------------------------------------
//
//
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4QIonIonCrossSection -- header file
// Created: M.V. Kossov, CERN/ITEP(Moscow), 19-Aug-07
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 19-Aug-07
//-----------------------------------------------------------------------------------------
// In this class, just as an experiment the real PDG codes (not CHIPS) are used
//-----------------------------------------------------------------------------------------
//
// ****************************************************************************************
// ********* This HEADER is temporary moved from the chips/interface directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************

#ifndef G4QIonIonCrossSection_h
#define G4QIonIonCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include <vector>
#include "G4VQCrossSection.hh"
#include "G4QPDGCode.hh"
#include "G4QProtonElasticCrossSection.hh"
#include "G4QNeutronElasticCrossSection.hh"
#include "G4QProtonNuclearCrossSection.hh"
#include "G4QNeutronNuclearCrossSection.hh"


class G4QIonIonCrossSection : public G4VQCrossSection
{
protected:

  G4QIonIonCrossSection()  {}

public:

  ~G4QIonIonCrossSection() {}

  static G4VQCrossSection* GetPointer(); // Gives a pointer to this singletone

  // At present momentum (pMom) must be in GeV(@@ Units),fCS=true:Inelastic, =false:Elastic
  virtual G4double GetCrossSection(G4bool fCS, G4double pMom, G4int Z, G4int N, G4int PDG);

  // Momentum=p/A (MeV/c)
  G4double CalculateCrossSection(G4bool fCS, G4int F, G4int I, G4int PDG,
                                 G4int tZ, G4int tN, G4double Momentum);

  G4double ThresholdMomentum(G4int pZ, G4int pN, G4int tZ, G4int tN);// P-Threshold(MeV/c)

private:
  G4int    GetFunctions(G4int pZ, G4int pN, G4int tZ, G4int tN,
                        G4double* LI, G4double* HI,  // Inelastic
                        G4double* LE, G4double* HE); // Elastic

  G4double CalculateTotal(G4double pA, G4double tA, G4double Momentum);
  G4double CalculateElTot(G4double pA, G4double tA, G4double Momentum);
  // Momentum=p/A (MeV/c), first=InelasticCS, second=elasticCS (mb)
  std::pair<G4double,G4double> CalculateXS(G4int pZ,G4int pN,G4int tZ,G4int tN,G4double P);

// Body
private:
  static G4double* lastLENI;// Pointer to the last array of LowEnergy Inel cross sections
  static G4double* lastHENI;// Pointer to the last array of HighEnergy Inel cross sections
  static G4double* lastLENE;// Pointer to the last array of LowEnergy Elast cross sections
  static G4double* lastHENE;// Pointer to the last array of HighEnergy Elast cross sections
  static G4double  lastE;   // Last used in the cross section Energy
  static G4int     lastPDG; // The last projectile PDG
  static G4int     lastN;   // The last N of calculated nucleus
  static G4int     lastZ;   // The last Z of calculated nucleus
  static G4double  lastP;   // Last used in the cross section Momentum
  static G4double  lastTH;  // Last value of the Momentum Threshold
  static G4double  lastICS; // Last value of the Inelastic Cross Section
  static G4double  lastECS; // Last value of the Elastic Cross Section
  static G4int     lastI;   // The last position in the DAMDB
};

#endif
