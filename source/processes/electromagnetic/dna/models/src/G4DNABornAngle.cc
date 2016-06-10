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
// $Id: G4DNABornAngle.cc 68380 2013-03-22 18:39:29Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4DNABornAngle
//
// Author:        Vladimir Ivantcheko
// 
// Creation date: 23 August 2013
//
// Modifications: 
//
// Class Description: 
//
// Delta-electron Angular Distribution Generation 
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#include "G4DNABornAngle.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

G4DNABornAngle::G4DNABornAngle(const G4String&)
  : G4VEmAngularDistribution("deltaBorn")
{
  fElectron = G4Electron::Electron();
}    

G4DNABornAngle::~G4DNABornAngle() 
{}

G4ThreeVector& 
G4DNABornAngle::SampleDirectionForShell(const G4DynamicParticle* dp,
					G4double secKinetic, G4int, 
					G4int, 
					const G4Material*)
{
  G4double k = dp->GetKineticEnergy();
  G4double cosTheta = 1.0;
  if (dp->GetDefinition() == fElectron)
    {
      if (secKinetic < 50.*eV) cosTheta = (2.*G4UniformRand())-1.;
      else if (secKinetic <= 200.*eV)
	{
	  if (G4UniformRand() <= 0.1) cosTheta = (2.*G4UniformRand())-1.;
	  else cosTheta = G4UniformRand()*(std::sqrt(2.)/2);
	}
      else
	{
	  G4double sin2O = 
	    (1.-secKinetic/k)/(1.+secKinetic/(2.*electron_mass_c2));
	  cosTheta = std::sqrt(1.-sin2O);
	}
    }
  else 
    {
      G4double mass = dp->GetDefinition()->GetPDGMass();
      G4double maxSecKinetic = 4.* (electron_mass_c2 / mass) * k;

      // cosTheta = std::sqrt(secKinetic / maxSecKinetic);

      // Restriction below 100 eV from Emfietzoglou (2000)

      if (secKinetic>100*eV) cosTheta = std::sqrt(secKinetic / maxSecKinetic);
      else cosTheta = (2.*G4UniformRand())-1.;
    }

  G4double sint = sqrt((1.0 - cosTheta)*(1.0 + cosTheta));
  G4double phi  = twopi*G4UniformRand(); 

  fLocalDirection.set(sint*cos(phi), sint*sin(phi), cosTheta);
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;
}

G4ThreeVector& 
G4DNABornAngle::SampleDirection(const G4DynamicParticle* dp,
				G4double secEkin, G4int Z, 
				const G4Material* mat)
{
  return SampleDirectionForShell(dp, secEkin, Z, 0, mat);
}

void G4DNABornAngle::PrintGeneratorInformation() const
{} 
