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
//
// Authors: Luciano Pandola (luciano.pandola at lngs.infn.it)
//
// History:
// -----------
//  
//  03 Dec 2009  First implementation, Luciano Pandola
//  16 Feb 2010  Added methods to calculate and store also A and Z 
//               for molecules. Luciano Pandola
//  16 Mar 2010  Added methods to calculate and store mean exc energy
//               and plasma energy (used for Ionisation). L Pandola
//  18 Mar 2010  Added method to retrieve number of atoms per 
//               molecule. L. Pandola
//  15 Mar 2012  Added method to retrieve number of atom of given Z per 
//               molecule, L. Pandola
//
// -------------------------------------------------------------------
//
//! Class description:
//! Fills and manages G4PenelopeOscillator objects and takes care of 
//! building and managing the G4PenelopeOscillatorTables for the materials 
//! in the geometry. G4PenelopeOscillatorManager is a singleton.
//! This is compliant with Penelope2008, so different tables (with different 
//! grouping factors) are created for Ionisation and Compton processes.
//

// -------------------------------------------------------------------

#ifndef G4PenelopeOscillatorManager_h
#define G4PenelopeOscillatorManager_h 1

#include "globals.hh"
#include "G4PenelopeOscillator.hh"
#include <vector>
#include <map>

class G4Material;

typedef std::vector<G4PenelopeOscillator*> G4PenelopeOscillatorTable ;

// This class is a singleton
class G4PenelopeOscillatorManager {

public: 

  // The only way to get an instance of this class is to call the 
  // function GetOscillatorManager() 
  static G4PenelopeOscillatorManager* GetOscillatorManager();
 
  //Clear() is invoked by Initialise() of the processes, if required
  void Clear();
  void Dump(const G4Material*);
  
  //For ionisation
  G4PenelopeOscillatorTable* GetOscillatorTableIonisation(const G4Material*);  
  G4PenelopeOscillator* GetOscillatorIonisation(const G4Material*,G4int);
  
  //For Compton
  G4PenelopeOscillatorTable* GetOscillatorTableCompton(const G4Material*);  
  G4PenelopeOscillator* GetOscillatorCompton(const G4Material*,G4int);

  void SetVerbosityLevel(G4int vl){verbosityLevel = vl;};
  G4int GetVerbosityLevel(){return verbosityLevel;};
  
  //!These are cumulative for the molecule
  //! Returns the total Z for the molecule
  G4double GetTotalZ(const G4Material*);
  //!Returns the total A for the molecule
  G4double GetTotalA(const G4Material*);
  //! Returns the mean excitation energy
  G4double GetMeanExcitationEnergy(const G4Material*);
  //! Returns the squared plasma energy
  G4double GetPlasmaEnergySquared(const G4Material*);
  //! Returns the total number of atoms per molecule
  G4double GetAtomsPerMolecule(const G4Material*);

  //Components of each molecule
  G4double GetNumberOfZAtomsPerMolecule(const G4Material*,G4int Z);
  
protected:
  G4PenelopeOscillatorManager();
  ~G4PenelopeOscillatorManager();

private:
  // Hide copy constructor and assignment operator 
  G4PenelopeOscillatorManager& operator=(const 
					 G4PenelopeOscillatorManager& right);
  G4PenelopeOscillatorManager(const G4PenelopeOscillatorManager&);
 
  static G4ThreadLocal G4PenelopeOscillatorManager* instance;
  
  //In Penelope2008, the Ionisation and Compton oscillator tables are 
  //slightly different!
  std::map<const G4Material*,G4PenelopeOscillatorTable*> 
  *oscillatorStoreIonisation;

  std::map<const G4Material*,G4PenelopeOscillatorTable*> 
  *oscillatorStoreCompton;

  std::map<const G4Material*,G4double> *atomicNumber;
  std::map<const G4Material*,G4double> *atomicMass;

  std::map<const G4Material*,G4double> *excitationEnergy;
  std::map<const G4Material*,G4double> *plasmaSquared;

  std::map<const G4Material*,G4double> *atomsPerMolecule;

  std::map< std::pair<const G4Material*,G4int>, G4double> *atomTablePerMolecule;

  //create both tables simultaneously
  void CheckForTablesCreated();

  void ReadElementData();
  G4double elementData[5][2000];
  G4bool fReadElementData;

  void BuildOscillatorTable(const G4Material*);

  G4int verbosityLevel;

};

#endif
