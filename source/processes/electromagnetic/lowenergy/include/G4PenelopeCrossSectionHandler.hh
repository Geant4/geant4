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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PenelopeCrossSectionHandler
//
// Author:        Luciano Pandola
// 
// Creation date: 04 July 2003
//
// -------------------------------------------------------------------
// Class description: 
// Provides cross sections with cut for Penelope Ionisation 
// -------------------------------------------------------------------
//

#ifndef G4PENELOPECROSSSECTIONHANDLER_HH
#define G4PENELOPECROSSSECTIONHANDLER_HH 1

#include "G4VCrossSectionHandler.hh"
#include "globals.hh"

class G4PenelopeIonisation;
class G4DataVector;
class G4VEMDataSet;
class G4VDataSetAlgorithm;
class G4ParticleDefinition;

class G4PenelopeCrossSectionHandler : public G4VCrossSectionHandler
{
public:

  G4PenelopeCrossSectionHandler(G4PenelopeIonisation* theProcess,
				const G4ParticleDefinition& aPartycleType,
				G4VDataSetAlgorithm* alg,
				G4double emin, 
				G4double emax, 
				G4int nbin);

  ~G4PenelopeCrossSectionHandler();
 
protected:

  std::vector<G4VEMDataSet*>* BuildCrossSectionsForMaterials(
                                const G4DataVector& energyVector, 
				const G4DataVector* energyCuts);


private:

  // Hide copy constructor and assignment operator 
  G4PenelopeCrossSectionHandler& operator=(const G4PenelopeCrossSectionHandler& right);
  G4PenelopeCrossSectionHandler(const G4PenelopeCrossSectionHandler&);
  
  G4PenelopeIonisation* process;
  const G4ParticleDefinition& particle;

  G4VDataSetAlgorithm* interp;

};
#endif

