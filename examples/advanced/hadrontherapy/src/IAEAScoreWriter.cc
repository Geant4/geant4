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
// This is the *basic* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// To obtain the full version visit the pages: http://sites.google.com/site/hadrontherapy/

#include "IAEAScoreWriter.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VScoringMesh.hh"

#include "HadrontherapyAnalysisManager.hh"

#include <map>
#include <fstream>

IAEAScoreWriter::IAEAScoreWriter()
  : G4VScoreWriter() {
  ;
}

IAEAScoreWriter::~IAEAScoreWriter() {
  ;
}

void IAEAScoreWriter::DumpQuantityToFiles(G4String & psName, G4String & option) {
  //
  if(verboseLevel > 0) {
    G4cout << "User-defined DumpQuantityToFile() method is invoked."
           << G4endl;
    G4cout << "Will write energy deposits along phantom"
	   << G4endl;
  }

  // change the option string into lowercase to the case-insensitive.
  G4String opt = option;
  std::transform(opt.begin(), opt.end(), opt.begin(), (int (*)(int))(tolower));

  // confirm the option
  if(opt.size() == 0) opt = "csv";

  // retrieve the map
  MeshScoreMap fSMap = fScoringMesh->GetScoreMap();
  

  MeshScoreMap::const_iterator msMapItr = fSMap.find(psName);
  if(msMapItr == fSMap.end()) {
    G4cerr << "ERROR : DumpToFile : Unknown quantity, \""
	   << psName << "\"." << G4endl;
    return;
  }
  std::map<G4int, G4double*> * score = msMapItr->second->GetMap();
  
  // declare xy array
  std::vector<double> projy;
  for(int y = 0; y < fNMeshSegments[1]; y++) projy.push_back(0.);
  std::vector<std::vector<double> > projxy;
  for(int x = 0; x < fNMeshSegments[0]; x++) projxy.push_back(projy);
  // accumulate
  // ofile << std::setprecision(16); // for double value with 8 bytes
  for(int x = 0; x < fNMeshSegments[0]; x++) {
    for(int y = 0; y < fNMeshSegments[1]; y++) {
      for(int z = 0; z < fNMeshSegments[2]; z++) {

	G4int idx = GetIndex(x, y, z);
	
	std::map<G4int, G4double*>::iterator value = score->find(idx);
	if(value != score->end()) projxy[x][y] += *(value->second);

      } // z
    } // y
  } // x

  // write quantity
  
#ifdef G4ANALYSIS_USE_ROOT // If we are using ROOT or AIDA analysis
  HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::GetInstance();
#endif
  for(int x = 0; x < fNMeshSegments[0]; x++) {
    for(int y = 0; y < fNMeshSegments[1]; y++) {
      /* There is one unused mashdimension here, but for now I've decided to ignore it */
	
      if(verboseLevel > 0) {
	std::cout << x << "\t" << projxy[x][y] << G4endl;
      }
#ifdef G4ANALYSIS_USE_ROOT
      analysis->BraggPeak(x, projxy[x][y]);
#endif
    } // y
  } // x
  // ofile << std::setprecision(6);

  // close the file
  //ofile.close();
  
}

