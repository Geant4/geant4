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
 
protected:

  G4std::vector<G4VEMDataSet*>* BuildCrossSectionsForMaterials(const G4DataVector& energyVector, 
							       const G4DataVector* energyCuts);

private:

  // Hide copy constructor and assignment operator 
  G4BremsstrahlungCrossSectionHandler& operator=(const G4BremsstrahlungCrossSectionHandler& right);
  G4BremsstrahlungCrossSectionHandler(const G4BremsstrahlungCrossSectionHandler&);

  const G4VEnergySpectrum* theBR;

  G4VDataSetAlgorithm* interp;
};
 
#endif

