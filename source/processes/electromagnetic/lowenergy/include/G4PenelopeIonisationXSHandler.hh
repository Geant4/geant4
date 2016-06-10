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
// $Id: G4PenelopeIonisationXSHandler.hh 74626 2013-10-17 07:00:59Z gcosmo $
//
// Author: Luciano Pandola
//
// History:
// -----------
// 09 Mar 2012   L. Pandola   1st implementation. 
//
// -------------------------------------------------------------------
//
//! Class description:
//!  This class is meant to calculate, store and provide the shell-per-shell
//!  and the total ionisation cross section calculated by the Penelope model.
//! The information is provided to other physics models, notably: Penelope 
//! ionisation model and Penelope "PIXE" model. Cross sections are calculated 
//! per material-cut couple and stored as G4PenelopeCrossSection objects.
//!
// -------------------------------------------------------------------

#ifndef G4PENELOPEIONISATIONXSHANDLER_HH
#define G4PENELOPEIONISATIONXSHANDLER_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4Material.hh"
#include <map>

class G4PhysicsFreeVector;
class G4PhysicsLogVector;
class G4ParticleDefinition;
class G4PenelopeOscillatorManager;
class G4PenelopeOscillator;
class G4PenelopeCrossSection;

class G4PenelopeIonisationXSHandler 
{

public:  
  //! Constructor. nBins is the number of intervals in the 
  //! energy grid. By default the energy grid goes from 100 eV
  //! to 100 GeV.
  G4PenelopeIonisationXSHandler(size_t nBins=200);
  
  //! Destructor. Clean all tables.
  virtual ~G4PenelopeIonisationXSHandler();

  //!Returns the density coeection for the material at the given energy
  G4double GetDensityCorrection(const G4Material*,const G4double energy) const;
  //! Returns the table of cross sections for the given particle, given 
  //! material and given cut as a G4PenelopeCrossSection* pointer.
  const G4PenelopeCrossSection* GetCrossSectionTableForCouple(const G4ParticleDefinition*,
							      const G4Material*,const G4double cut) const;
  //!Setter for the verbosity level
  void SetVerboseLevel(G4int vl){verboseLevel = vl;};

  //! This can be inkoved only by the master
  void BuildXSTable(const G4Material*,G4double cut,
		    const G4ParticleDefinition*,G4bool isMaster=true);

private:
  G4PenelopeIonisationXSHandler & operator=(const G4PenelopeIonisationXSHandler &right);
  G4PenelopeIonisationXSHandler(const G4PenelopeIonisationXSHandler&);

  void BuildDeltaTable(const G4Material*);

  G4DataVector* ComputeShellCrossSectionsElectron(G4PenelopeOscillator* ,
						  G4double energy,G4double cut,
					          G4double delta);
    
  G4DataVector* ComputeShellCrossSectionsPositron(G4PenelopeOscillator* ,
						  G4double energy,G4double cut,
		                                  G4double delta);

  //Oscillator manager
  G4PenelopeOscillatorManager* oscManager;

  //G4PenelopeCrossSection takes care of the logs
  std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*> *XSTableElectron;
  std::map< std::pair<const G4Material*,G4double>, G4PenelopeCrossSection*> *XSTablePositron;
  
  //delta vs. log(energy)
  std::map<const G4Material*,G4PhysicsFreeVector*> *theDeltaTable;

  //energy grid
  G4PhysicsLogVector* energyGrid;
  size_t nBins;

  G4int verboseLevel;

};

#endif

