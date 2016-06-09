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
// GEANT4 tag $Name: geant4-08-01 $
//
//
// GEANT4 physics class: G4QProtonNuclearCrossSection -- header file
// Created: M.V. Kossov, CERN/ITEP(Moscow), 20-Dec-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 17-May-02
//
// ****************************************************************************************
// ********* This HEADER is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************

#ifndef G4QProtonNuclearCrossSection_h
#define G4QProtonNuclearCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include <vector>
#include "G4VQCrossSection.hh"

class G4QProtonNuclearCrossSection : public G4VQCrossSection
{
protected:

  G4QProtonNuclearCrossSection()  {}

public:

  ~G4QProtonNuclearCrossSection() {}

  static G4VQCrossSection* GetPointer(); // Gives a pointer to this singletone

  G4double CalculateCrossSection(G4bool CS, G4int F, G4int I, G4int PDG, G4int Z,
                                 G4int N, G4double Momentum);

private:
  G4int    GetFunctions(G4double A, G4double* y, G4double* z);// y&z are pointers to arrays

// Body
private:
  static G4double  lastSig; // Last value of the Cross Section
  static G4double* lastLEN; // Pointer to the last array of LowEnergy cross sections
  static G4double* lastHEN; // Pointer to the last array of HighEnergy cross sections
  static G4double  lastE;   // Last used in the cross section Energy
  static G4double  lastSP;  // Last value of the ShadowingPomeron (A-dependent)
};

#endif
