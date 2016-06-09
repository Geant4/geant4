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
// File name:     G4RDBremsstrahlungCrossSectionHandler
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

#ifndef G4RDBremsstrahlungCrossSectionHandler_h
#define G4RDBremsstrahlungCrossSectionHandler_h 1

#include "G4RDVCrossSectionHandler.hh"
#include "G4RDVEnergySpectrum.hh"
#include "globals.hh"

class G4DataVector;
class G4RDVEMDataSet;
class G4RDVDataSetAlgorithm;

class G4RDBremsstrahlungCrossSectionHandler : public G4RDVCrossSectionHandler
{

public:

  G4RDBremsstrahlungCrossSectionHandler(const G4RDVEnergySpectrum* spectrum,
				      G4RDVDataSetAlgorithm* interpolation);

  ~G4RDBremsstrahlungCrossSectionHandler();
 
protected:

  std::vector<G4RDVEMDataSet*>* BuildCrossSectionsForMaterials(const G4DataVector& energyVector, 
							       const G4DataVector* energyCuts);

private:

  // Hide copy constructor and assignment operator 
  G4RDBremsstrahlungCrossSectionHandler& operator=(const G4RDBremsstrahlungCrossSectionHandler& right);
  G4RDBremsstrahlungCrossSectionHandler(const G4RDBremsstrahlungCrossSectionHandler&);

  const G4RDVEnergySpectrum* theBR;

  G4RDVDataSetAlgorithm* interp;
};
 
#endif

