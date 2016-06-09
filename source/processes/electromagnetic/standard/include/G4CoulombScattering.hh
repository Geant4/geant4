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
// $Id: G4CoulombScattering.hh,v 1.8 2007/07/31 17:24:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4CoulombScattering
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 12.03.2006
//
// Modifications:
//
// Class Description:
//
// This class manages the process of Coulomb elastic scattering
//

// -------------------------------------------------------------------
//

#ifndef G4CoulombScattering_h
#define G4CoulombScattering_h 1

#include "G4VEmProcess.hh"
#include "G4VEmModel.hh"

class G4CoulombScattering : public G4VEmProcess
{

public:

  G4CoulombScattering(const G4String& name = "eCoulombScat");

  virtual ~G4CoulombScattering();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  void SetThetaMin(G4double);

  void SetThetaMax(G4double);

  void SetQ2Max(G4double);

  // obsolete method to be removed
  void SetBuildTableFlag(G4bool);

  // Print out of the class parameters
  virtual void PrintInfo();

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);

private:

 // hide assignment operator
  G4CoulombScattering & operator=(const G4CoulombScattering &right);
  G4CoulombScattering(const G4CoulombScattering&);
  
  G4double thetaMin;
  G4double thetaMax;
  G4double q2Max;
  G4double thEnergy;
  G4double thEnergyElec;
  G4bool isInitialised;
  G4bool buildElmTableFlag;
  const G4ParticleDefinition* aParticle;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4CoulombScattering::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4CoulombScattering::SetThetaMin(G4double val)
{
  thetaMin = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4CoulombScattering::SetThetaMax(G4double val)
{
  thetaMax = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4CoulombScattering::SetQ2Max(G4double val)
{
  q2Max = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4CoulombScattering::SetBuildTableFlag(G4bool)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
