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
// File name:     G4eIonisationCrossSectionHandler
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
// Provides build cross sections with cut for LowEnergyIonisation 
//
// Class Description: End 

// -------------------------------------------------------------------
//

#ifndef G4eIonisationCrossSectionHandler_h
#define G4eIonisationCrossSectionHandler_h 1

#include "G4VCrossSectionHandler.hh"
#include "globals.hh"

class G4VEnergySpectrum;
class G4DataVector;
class G4VEMDataSet;
class G4VDataSetAlgorithm;

class G4eIonisationCrossSectionHandler : public G4VCrossSectionHandler
{

public:

  G4eIonisationCrossSectionHandler(const G4VEnergySpectrum* spec,
                                         G4VDataSetAlgorithm* alg,
                                         G4double emin, 
                                         G4double emax, 
                                         G4int nbin);

  ~G4eIonisationCrossSectionHandler();
 
protected:

  G4std::vector<G4VEMDataSet*>* BuildCrossSectionsForMaterials(
                                const G4DataVector& energyVector, 
				const G4DataVector* energyCuts);


private:

  const G4VEnergySpectrum* theParam;

  G4VDataSetAlgorithm* interp;
};
 
#endif

