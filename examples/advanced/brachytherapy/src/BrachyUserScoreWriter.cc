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
/*
// Code developed by:
// S.Guatelli, susanna@uow.edu.au
//
Original code from geant4/examples/extended/runAndEvent/RE03, by M. Asai
*/
#include <map>
#include <fstream>
#include "BrachyUserScoreWriter.hh"

#include <CLHEP/Units/SystemOfUnits.h>

#ifdef ANALYSIS_USE  
#include "BrachyAnalysisManager.hh"
#endif

#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VScoringMesh.hh"
// The default output is
// voxelX, voxelY, voxelZ, edep
// The BrachyUserScoreWriter allows to change the format of the output file.
// in the specific case:
// xx (mm)  yy(mm) zz(mm) edep(keV)
// The same information is stored in a ntuple, in the 
// brachytherapy.root file

BrachyUserScoreWriter::BrachyUserScoreWriter():
G4VScoreWriter() 
{
}

BrachyUserScoreWriter::~BrachyUserScoreWriter() 
{;}

void BrachyUserScoreWriter::DumpQuantityToFile(const G4String & psName,
                                               const G4String & fileName, 
                                               const G4String & option) 
{
using MeshScoreMap = G4VScoringMesh::MeshScoreMap;

if(verboseLevel > 0) 
  {G4cout << "BrachyUserScorer-defined DumpQuantityToFile() method is invoked."
  << G4endl; 
  }

// change the option string into lowercase to the case-insensitive.
G4String opt = option;
std::transform(opt.begin(), opt.end(), opt.begin(), (int (*)(int))(tolower));

// confirm the option
if(opt.size() == 0) opt = "csv";

// open the file
std::ofstream ofile(fileName);
  
if(!ofile) 
{
   G4cerr << "ERROR : DumpToFile : File open error -> " << fileName << G4endl;
   return;
}
  ofile << "# mesh name: " << fScoringMesh->GetWorldName() << G4endl;

// retrieve the map
MeshScoreMap fSMap = fScoringMesh -> GetScoreMap();
  
MeshScoreMap::const_iterator msMapItr = fSMap.find(psName);
  
if(msMapItr == fSMap.end()) 
  {
   G4cerr << "ERROR : DumpToFile : Unknown quantity, \""<< psName 
   << "\"." << G4endl;
   return;
  }

std::map<G4int, G4StatDouble*> * score = msMapItr -> second-> GetMap();
  
ofile << "# primitive scorer name: " << msMapItr -> first << G4endl;
//
// Write quantity in the ASCII output file and in brachytherapy.root
//
ofile << std::setprecision(16); // for double value with 8 bytes
  
for(int x = 0; x < fNMeshSegments[0]; x++) {
   for(int y = 0; y < fNMeshSegments[1]; y++) {
     for(int z = 0; z < fNMeshSegments[2]; z++){
        G4int numberOfVoxel_x = fNMeshSegments[0];
        G4int numberOfVoxel_y = fNMeshSegments[1];
        G4int numberOfVoxel_z =fNMeshSegments[2];
        // If the voxel width is changed in the macro file, 
        // the voxel width variable must be updated
        G4double voxelWidth = 0.25 *CLHEP::mm;
        //
        G4double xx = ( - numberOfVoxel_x + 1+ 2*x )* voxelWidth/2;
        G4double yy = ( - numberOfVoxel_y + 1+ 2*y )* voxelWidth/2;
        G4double zz = ( - numberOfVoxel_z + 1+ 2*z )* voxelWidth/2;
        G4int idx = GetIndex(x, y, z);
        std::map<G4int, G4StatDouble*>::iterator value = score -> find(idx);

       if (value != score -> end()) 
        {
         // Print in the ASCII output file the information
 
         ofile << xx << "  " << yy << "  " << zz <<"  " 
               <<(value->second->sum_wx())/CLHEP::keV << G4endl;

#ifdef ANALYSIS_USE          
        // Save the same information in the output analysis file
       BrachyAnalysisManager* analysis = BrachyAnalysisManager::GetInstance();
   
       if(zz> -0.125 *CLHEP::mm && zz < 0.125*CLHEP::mm)
         analysis -> FillH2WithEnergyDeposition(xx,yy,
                       (value->second->sum_wx())/CLHEP::keV);
#endif
}}}} 

ofile << std::setprecision(6);

// Close the output ASCII file
ofile.close();
}
