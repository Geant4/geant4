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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4BremsstrahlungCrossSectionHandler
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 17 September 2001
//
// Modified: 
//
// -------------------------------------------------------------------

// Class Description: 
//
// Provides build cross sections with cut for LowEnergyBremsstrahlung 
//
// Class Description: End 

// -------------------------------------------------------------------
//

#ifndef G4BremsstrahlungCrossSectionHandler_h
#define G4BremsstrahlungCrossSectionHandler_h 1

#include "G4VCrossSectionHandler.hh"
#include "G4VEnergySpectrum.hh"
#include "globals.hh"

class G4DataVector;
class G4VEMDataSet;
class G4VDataSetAlgorithm;

class G4BremsstrahlungCrossSectionHandler : public G4VCrossSectionHandler
{

public:

  G4BremsstrahlungCrossSectionHandler(const G4VEnergySpectrum* spectrum,
				      G4VDataSetAlgorithm* interpolation);

  ~G4BremsstrahlungCrossSectionHandler();
 
  G4double GetCrossSectionAboveThresholdForElement(G4double energy,
                                                   G4double cutEnergy,
                                                   G4int Z);

protected:

  std::vector<G4VEMDataSet*>* BuildCrossSectionsForMaterials(const G4DataVector& energyVector, 
							       const G4DataVector* energyCuts);

private:

  // Hide copy constructor and assignment operator 
  G4BremsstrahlungCrossSectionHandler& operator=(const G4BremsstrahlungCrossSectionHandler& right);
  G4BremsstrahlungCrossSectionHandler(const G4BremsstrahlungCrossSectionHandler&);

  const G4VEnergySpectrum* theBR;

  G4VDataSetAlgorithm* interp;
};
 
#endif

