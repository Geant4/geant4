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
// GEANT4 tag $Name: geant4-07-00-patch-01 $
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

#ifndef G4QPhotonNuclearCrossSection_h
#define G4QPhotonNuclearCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include <vector>
#include "G4VQCrossSection.hh"

class G4QPhotonNuclearCrossSection : public G4VQCrossSection
{
protected:

  G4QPhotonNuclearCrossSection()  {}

public:

  ~G4QPhotonNuclearCrossSection() {}

  static G4VQCrossSection* GetPointer(); // Gives a pointer to this singletone

  G4double CalculateCrossSection(G4int F, G4int I, G4int Z, G4int N, G4double Momentum);

  G4double ThresholdEnergy(G4int Z, G4int N);

private:
  G4int    GetFunctions(G4double A, G4double* y, G4double* z);// y&z are pointers to arrays

// Body
private:
  static G4double  lastSig; // Last value of the Cross Section
  static G4double* lastGDR; // Pointer to the last array of GDR cross sections
  static G4double* lastHEN; // Pointer to the last array of HEn cross sections
  static G4double  lastE;   // Last used in the cross section Energy
  static G4double  lastSP;  // Last value of the ShadowingPomeron (A-dependent)
};

#endif
