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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4QPhotonNuclearCrossSection -- header file
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-May-02
//
// ****************************************************************************************
// ********* This HEADER is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************
// Short description: This is an original CHIPS process for photo-nuclear
// interactions, which does not include "fast and dirty" corrections for
// reactions near threshold, with respect to the GHAD application of CHIPS.
// ------------------------------------------------------------------------

#ifndef G4QPhotonNuclearCrossSection_h
#define G4QPhotonNuclearCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include <vector>
#include "G4VQCrossSection.hh"

class G4QPhotonNuclearCrossSection : public G4VQCrossSection
{
protected:

  G4QPhotonNuclearCrossSection()  {}

public:

  ~G4QPhotonNuclearCrossSection();

  static G4VQCrossSection* GetPointer(); // Gives a pointer to this singletone

  // At present momentum (pMom) must be in GeV (@@ Units)
  virtual G4double GetCrossSection(G4bool fCS, G4double pMom, G4int tgZ, G4int tgN,
                                                                             G4int pPDG=0);

  G4double CalculateCrossSection(G4bool CS, G4int F, G4int I, G4int PDG, G4int Z, G4int N,
                                                                        G4double Momentum);

  G4double ThresholdEnergy(G4int Z, G4int N, G4int PDG=22);

private:
  G4int    GetFunctions(G4double A, G4double* y, G4double* z);// y&z are pointers to arrays

// Body
private:
  static G4bool    onlyCS;  // flag to calculate only CrossSection
  static G4double  lastSig; // Last value of the Cross Section
  static G4double* lastGDR; // Pointer to the last array of GDR cross sections
  static G4double* lastHEN; // Pointer to the last array of HEn cross sections
  static G4double  lastE;   // Last used in the cross section Energy
  static G4double  lastSP;  // Last value of the ShadowingPomeron (A-dependent)
  static G4int     lastPDG;  // The last projectile PDG
  static G4int     lastN;    // The last N of calculated nucleus
  static G4int     lastZ;    // The last Z of calculated nucleus
  static G4double  lastP;    // Last used in the cross section Momentum
  static G4double  lastTH;   // Last value of the Momentum Threshold
  static G4double  lastCS;   // Last value of the Cross Section
  static G4int     lastI;    // The last position in the DAMDB
  static std::vector <G4double*>* GDR; // Vector of pointers to GDRPhotonuclearCrossSection
  static std::vector <G4double*>* HEN; // Vector of pointers to HighEnPhotonuclearCrossSect
};

#endif
