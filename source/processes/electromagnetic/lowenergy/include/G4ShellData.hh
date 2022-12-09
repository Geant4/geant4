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
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
//  6 Aug 2001   MGP        Created
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Shell data set: shell identifiers and binding energies

// -------------------------------------------------------------------

#ifndef G4SHELLDATA_HH
#define G4SHELLDATA_HH 1

#include "globals.hh"
#include <vector>
#include <map>

class G4DataVector;

class G4ShellData 
{ 
public:
  explicit G4ShellData(G4int minZ = 1, G4int maxZ = 100, G4bool isOccupancy = false);
  ~G4ShellData();
 
  std::size_t NumberOfShells(G4int Z) const;
  G4int ShellId(G4int Z, G4int shellIndex) const;
  G4double ShellOccupancyProbability(G4int Z, G4int shellIndex) const;
  const std::vector<G4double>& ShellIdVector(G4int Z) const;
  G4double BindingEnergy(G4int Z, G4int shellIndex) const;
  void SetOccupancyData() { occupancyData = true; }
  void LoadData(const G4String& fileName);
  void PrintData() const;

  // Randomly select a shell based on shell occupancy
  G4int SelectRandomShell(G4int Z) const;

  // Hide copy constructor and assignment operator 
  G4ShellData& operator=(const G4ShellData& right) = delete;
  G4ShellData(const G4ShellData&) = delete;

private: 
  const std::vector<G4double>& ShellVector(G4int Z) const;

  std::map<G4int,std::vector<G4double>*,std::less<G4int> > idMap;
  std::map<G4int,G4DataVector*,std::less<G4int> > bindingMap;
  std::vector<G4int> nShells;
  std::map<G4int,std::vector<G4double>*,std::less<G4int> > occupancyPdfMap;

  G4int zMin;
  G4int zMax; 
  G4bool occupancyData;
};
 
#endif









