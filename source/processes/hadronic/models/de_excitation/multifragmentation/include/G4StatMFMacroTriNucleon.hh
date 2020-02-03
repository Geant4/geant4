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
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMacroTriNucleon_h
#define G4StatMFMacroTriNucleon_h 1

#include "G4VStatMFMacroCluster.hh"
#include "G4NucleiProperties.hh"

class G4StatMFMacroTriNucleon : public G4VStatMFMacroCluster {

public:

  G4StatMFMacroTriNucleon();

  ~G4StatMFMacroTriNucleon();

  G4double CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
				const G4double nu, const G4double T);

  inline G4double CalcZARatio(const G4double ) {return theZARatio = 0.5;}

  G4double CalcEnergy(const G4double T);
	
  G4double CalcEntropy(const G4double T, const G4double FreeVol);
	
private:

  // Copy constructor
  G4StatMFMacroTriNucleon(const G4StatMFMacroTriNucleon & right);

  // operators
  G4StatMFMacroTriNucleon & operator=(const G4StatMFMacroTriNucleon & right);
  G4bool operator==(const G4StatMFMacroTriNucleon & right) const;
  G4bool operator!=(const G4StatMFMacroTriNucleon & right) const;

};

#endif
