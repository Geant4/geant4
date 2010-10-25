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
// $Id: G4CoulombScattering.hh,v 1.15 2010-10-25 19:13:23 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  G4CoulombScattering(const G4String& name = "CoulombScat");

  virtual ~G4CoulombScattering();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  inline void SetThetaMin(G4double);

  inline void SetThetaMax(G4double);

  inline void SetQ2Max(G4double);

  // Set energy above which high energy model will be used
  inline void SetHEModelLimit(G4double);

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
  const G4ParticleDefinition* aParticle;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
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

inline void G4CoulombScattering::SetHEModelLimit(G4double val)
{
  thEnergy = val;
  thEnergyElec = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
