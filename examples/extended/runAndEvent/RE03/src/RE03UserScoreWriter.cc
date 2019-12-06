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
/// \file runAndEvent/RE03/src/RE03UserScoreWriter.cc
/// \brief Implementation of the RE03UserScoreWriter class
//
//

#include "RE03UserScoreWriter.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VScoringMesh.hh"

#include <map>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE03UserScoreWriter::RE03UserScoreWriter()
  : G4VScoreWriter() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE03UserScoreWriter::~RE03UserScoreWriter() 
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE03UserScoreWriter::DumpQuantityToFile(const G4String & psName, 
                                             const G4String & fileName, 
                                             const G4String & option) {
  using MeshScoreMap = G4VScoringMesh::MeshScoreMap;
  //
  if(verboseLevel > 0) {
    G4cout << "User-defined DumpQuantityToFile() method is invoked."
           << G4endl;
    G4cout << "  -- to obtain a projection of the quantity <"
           << psName
           << "> onto the x-y plane --" << G4endl;
  }

  // change the option string into lowercase to the case-insensitive.
  G4String opt = option;
  std::transform(opt.begin(), opt.end(), opt.begin(), (int (*)(int))(tolower));

  // confirm the option
  if(opt.size() == 0) opt = "csv";

  // open the file
  std::ofstream ofile(fileName);
  if(!ofile) {
    G4cerr << "ERROR : DumpToFile : File open error -> "
           << fileName << G4endl;
    return;
  }
  ofile << "# mesh name: " << fScoringMesh->GetWorldName() << G4endl;

  
  // retrieve the map
  MeshScoreMap scMap = fScoringMesh->GetScoreMap();

  MeshScoreMap::const_iterator msMapItr = scMap.find(psName);
  if(msMapItr == scMap.end()) {
    G4cerr << "ERROR : DumpToFile : Unknown quantity, \""
           << psName << "\"." << G4endl;
    return;
  }
  std::map<G4int, G4StatDouble*> * score = msMapItr->second->GetMap();
  ofile << "# primitive scorer name: " << msMapItr->first << G4endl;

  // write header info 
  ofile << "# xy projection" << G4endl;
  ofile << fNMeshSegments[0] << " " << fNMeshSegments[1] << " " << G4endl;

  // declare xy array
  std::vector<double> projy;
  for(int y = 0; y < fNMeshSegments[1]; y++) projy.push_back(0.);
  std::vector<std::vector<double> > projxy;
  for(int x = 0; x < fNMeshSegments[0]; x++) projxy.push_back(projy);
  // accumulate
  ofile << std::setprecision(16); // for double value with 8 bytes
  for(int x = 0; x < fNMeshSegments[0]; x++) {
    for(int y = 0; y < fNMeshSegments[1]; y++) {
      for(int z = 0; z < fNMeshSegments[2]; z++) {

        G4int idx = GetIndex(x, y, z);
        
        std::map<G4int, G4StatDouble*>::iterator value = score->find(idx);
        if(value != score->end()) projxy[x][y] += value->second->sum_wx();

      } // z
    } // y
  } // x

  // write quantity
  ofile << std::setprecision(16); // for double value with 8 bytes
  for(int x = 0; x < fNMeshSegments[0]; x++) {
    for(int y = 0; y < fNMeshSegments[1]; y++) {

      ofile << x << "," << y << ",";
      ofile << projxy[x][y] << G4endl;

    } // y
  } // x
  ofile << std::setprecision(6);

  // close the file
  ofile.close();
  
}

