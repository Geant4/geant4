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
// File name:     G4eIonisationCrossSectionHandler
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 17 September 2001
//
// Modified: 
// 10 Oct 2001  M.G. Pia        Revision to improve code quality and consistency with design
// 28 Jan 2009  L.Pandola    Added public method to make a easier migration of
//                           G4LowEnergyIonisation to G4LivermoreIonisationModel
// -------------------------------------------------------------------

// Class description: 
// Provides cross sections with cut for LowEnergyIonisation 
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4EIONISATIONCROSSSECTIONHANDLER_HH
#define G4EIONISATIONCROSSSECTIONHANDLER_HH 1

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

  G4double GetCrossSectionAboveThresholdForElement(G4double energy,
						   G4double cutEnergy,
						   G4int Z);
 
protected:

  std::vector<G4VEMDataSet*>* BuildCrossSectionsForMaterials(
                                const G4DataVector& energyVector, 
				const G4DataVector* energyCuts);


private:

  // Hide copy constructor and assignment operator 
  G4eIonisationCrossSectionHandler& operator=(const G4eIonisationCrossSectionHandler& right);
  G4eIonisationCrossSectionHandler(const G4eIonisationCrossSectionHandler&);

 const G4VEnergySpectrum* theParam;

  G4VDataSetAlgorithm* interp;
  G4int verbose;
  
};
 
#endif

