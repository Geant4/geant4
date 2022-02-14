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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRScoreWriter.hh
//   Defines the printout format of primitive scorer
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#include "GRScoreWriter.hh"

#include <map>
#include <fstream>

#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VScoringMesh.hh"

GRScoreWriter::GRScoreWriter()
{;}

GRScoreWriter::~GRScoreWriter()
{;}

void GRScoreWriter::DumpQuantityToFile(const G4String& psName,
                                        const G4String& fileName,
                                        const G4String& option) {

  // change the option string into lowercase to the case-insensitive.
  G4String opt = option;
  std::transform(opt.begin(), opt.end(), opt.begin(), (int (*)(int))(tolower));

  // confirm the option
  if(opt.size() == 0) opt = "csv";
  if(opt.find("csv") == std::string::npos &&
     opt.find("sequence") == std::string::npos) {
    G4cerr << "ERROR : DumpToFile : Unknown option -> "
           << option << G4endl;
    return;
  }

  // open the file
  std::ofstream ofile(fileName);
  if(!ofile) {
    G4cerr << "ERROR : DumpToFile : File open error -> "
           << fileName << G4endl;
    return;
  }
  ofile << "# Mesh or volume name: " << fScoringMesh->GetWorldName() << G4endl;

  using MeshScoreMap = G4VScoringMesh::MeshScoreMap;
  // retrieve the map
  MeshScoreMap fSMap = fScoringMesh->GetScoreMap();


  MeshScoreMap::const_iterator msMapItr = fSMap.find(psName);
  if(msMapItr == fSMap.end()) {
    G4cerr << "ERROR : DumpToFile : Unknown quantity, \""
           << psName << "\"." << G4endl;
    return;
  }

  std::map<G4int, G4StatDouble*> * score = msMapItr->second->GetMap();
  ofile << "# Primitive scorer name: " << msMapItr->first << G4endl;
  if(fact!=1.0)
  { ofile << "# Multiplication factor : " << fact << G4endl; }

  G4double unitValue = fScoringMesh->GetPSUnitValue(psName);
  G4String unit = fScoringMesh->GetPSUnit(psName);
  G4String divisionAxisNames[3];
  fScoringMesh->GetDivisionAxisNames(divisionAxisNames);
  ofile << "# First three integer entries: index of a cell of a mesh, just one cell for a volume tally." << G4endl;
  ofile << "# Forth entry: sum of scores." << G4endl;
  ofile << "# Fifth entry: sum of squared scores." << G4endl;
  ofile << "# Sixth entry: number of events with non-zero effect." << G4endl;
  ofile << "# Seventh entry: relative statistical error in %." << G4endl << G4endl;
  // index of the cell
  ofile << "# i" << divisionAxisNames[0]
        << ", i" << divisionAxisNames[1]
        << ", i" << divisionAxisNames[2];
  // unit of scored value
  ofile << ", total(value) ";
  if(unit.size() > 0) ofile << "[" << unit << "]";
  ofile << ", total(value^2), number of entries, relative error (%)" << G4endl;

  // "sequence" option: write header info
  if(opt.find("sequence") != std::string::npos) {
    ofile << fNMeshSegments[0] << " " << fNMeshSegments[1] << " " << fNMeshSegments[2]
          << G4endl;
  }

  // write quantity
  long count = 0;
  ofile << std::setprecision(16); // for double value with 8 bytes
  for(int x = 0; x < fNMeshSegments[0]; x++) {
    for(int y = 0; y < fNMeshSegments[1]; y++) {
      for(int z = 0; z < fNMeshSegments[2]; z++) {
        G4int idx = GetIndex(x, y, z);

        if(opt.find("csv") != std::string::npos)
          ofile << x << "," << y << "," << z << ",";

        std::map<G4int, G4StatDouble*>::iterator value = score->find(idx);
        if(value == score->end()) {
          ofile << 0. << "," << 0. << "," << 0;
        } else {
          G4double x1 = value->second->sum_wx()/unitValue*fact;
          G4double x2 = value->second->sum_wx2()/unitValue/unitValue*fact*fact;
          G4int n = value->second->n();
          // rms  is sigma = sqrt((x2-x1^2/n)/(n-1))
          // var = mean +/- rms/sqrt(n): means that the relative error of the sum = rms*sqrt(n)/x1
          G4double rms = value->second->rms()/unitValue*fact;
          // Relative error in %
          G4double relError = rms*std::sqrt(n)*100./x1;
          ofile << x1 << ", " << x2 << ", " << n << ", " << relError;
          ofile << G4endl;
          G4double factor = relError*relError/100.;
          if (factor > 1.)
          {
             G4cout << "# Mesh or volume name: " << fScoringMesh->GetWorldName() 
               << "  -- # Primitive scorer name: " << msMapItr->first << G4endl
               << "  bin " << x << "," << y << "," << z << " : statistical error " << relError << "(%)" << G4endl
               << "   to reduce the statistical error below 10%, increase number of events approximately "
               << factor << " times." << G4endl;
          }
        }

        if(opt.find("csv") != std::string::npos) {
          ofile << G4endl;
        } else if(opt.find("sequence") != std::string::npos) {
          ofile << " ";
          if(count++%5 == 4) ofile << G4endl;
        }

      } // z
    } // y
  } // x
  ofile << std::setprecision(6);

  // close the file
  ofile.close();
}

void GRScoreWriter::DumpAllQuantitiesToFile(const G4String& fileName,
                                             const G4String& option) {

  // change the option string into lowercase to the case-insensitive.
  G4String opt = option;
  std::transform(opt.begin(), opt.end(), opt.begin(), (int (*)(int))(tolower));

  // confirm the option
  if(opt.size() == 0) opt = "csv";
  if(opt.find("csv") == std::string::npos &&
     opt.find("sequence") == std::string::npos) {
    G4cerr << "ERROR : DumpToFile : Unknown option -> "
           << option << G4endl;
    return;
  }

  // open the file
  std::ofstream ofile(fileName);
  if(!ofile) {
    G4cerr << "ERROR : DumpToFile : File open error -> "
           << fileName << G4endl;
    return;
  }
  ofile << "# Mesh or volume name: " << fScoringMesh->GetWorldName() << G4endl;
  if(fact!=1.0)
  { ofile << "# Multiplication factor : " << fact << G4endl; }
  ofile << "# First three integer entries: index of a cell of a mesh, just one cell for a volume tally." << G4endl;
  ofile << "# Forth entry: sum of scores." << G4endl;
  ofile << "# Fifth entry: sum of squared scores." << G4endl;
  ofile << "# Sixth entry: number of events with non-zero effect." << G4endl;
  ofile << "# Seventh entry: relative statistical error in %." << G4endl << G4endl;

  // retrieve the map
  using MeshScoreMap = G4VScoringMesh::MeshScoreMap;
  MeshScoreMap fSMap = fScoringMesh->GetScoreMap();
  MeshScoreMap::const_iterator msMapItr = fSMap.begin();
  std::map<G4int, G4StatDouble*> * score;
  for(; msMapItr != fSMap.end(); msMapItr++) {

    G4String psname = msMapItr->first;

    score = msMapItr->second->GetMap();
    ofile << "# Primitive scorer name: " << msMapItr->first << G4endl;

    G4double unitValue = fScoringMesh->GetPSUnitValue(psname);
    G4String unit = fScoringMesh->GetPSUnit(psname);
    G4String divisionAxisNames[3];
    fScoringMesh->GetDivisionAxisNames(divisionAxisNames);
    // index order
    ofile << "# i" << divisionAxisNames[0]
          << ", i" << divisionAxisNames[1]
          << ", i" << divisionAxisNames[2];
    // unit of scored value
    ofile << ", total(value) ";
    if(unit.size() > 0) ofile << "[" << unit << "]";
    ofile << ", total(value^2), number of entries, relative error (%)" << G4endl;


    // "sequence" option: write header info
    if(opt.find("sequence") != std::string::npos) {
      ofile << fNMeshSegments[0] << " " << fNMeshSegments[1] << " " << fNMeshSegments[2]
            << G4endl;
    }

    // write quantity
    long count = 0;
    ofile << std::setprecision(16); // for double value with 8 bytes
    for(int x = 0; x < fNMeshSegments[0]; x++) {
      for(int y = 0; y < fNMeshSegments[1]; y++) {
        for(int z = 0; z < fNMeshSegments[2]; z++) {
          G4int idx = GetIndex(x, y, z);

          if(opt.find("csv") != std::string::npos)
            ofile << x << "," << y << "," << z << ",";

          std::map<G4int, G4StatDouble*>::iterator value = score->find(idx);
          if(value == score->end()) {
            ofile << 0. << "," << 0. << "," << 0;
          } else {
            G4double x1 = value->second->sum_wx()/unitValue*fact;
            G4double x2 = value->second->sum_wx2()/unitValue/unitValue*fact*fact;
            G4int n = value->second->n();
            // rms  is sigma = sqrt((x2-x1^2/n)/(n-1))
            // var = mean +/- rms/sqrt(n): means that the relative error of the sum = rms*sqrt(n)/x1
            G4double rms = value->second->rms()/unitValue*fact;
            // Relative error in %
            G4double relError = rms*std::sqrt(n)*100./x1;
            ofile << x1 << ", " << x2 << ", " << n << ", " << relError;
            ofile << G4endl;
            G4double factor = relError*relError/100.;
            if (factor > 1.)
            {
               G4cout << "# Mesh or volume name: " << fScoringMesh->GetWorldName() 
                 << "  -- # Primitive scorer name: " << msMapItr->first << G4endl
                 << "  bin " << x << "," << y << "," << z << " : statistical error " << relError << "(%)" << G4endl
                 << "   to reduce the statistical error below 10%, increase number of events approximately "
                 << factor << " times." << G4endl;
            }
          }

          if(opt.find("csv") != std::string::npos) {
            ofile << G4endl;
          } else if(opt.find("sequence") != std::string::npos) {
            ofile << " ";
            if(count++%5 == 4) ofile << G4endl;
          }

        } // z
      } // y
    } // x
    ofile << std::setprecision(6);

  } // for(; msMapItr ....)

  // close the file
  ofile.close();

}


