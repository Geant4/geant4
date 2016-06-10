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
// $Id: G4LevelReader.cc 77025 2013-11-20 16:11:51Z gcosmo $
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4NucLevel
//
//      Author:        V.Ivanchenko (M.Kelsey reading method is used)
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//      
// -------------------------------------------------------------------

#include "G4LevelReader.hh"
#include "G4NucLevel.hh"
#include "G4SystemOfUnits.hh"

G4LevelReader::G4LevelReader() 
  : nLevels(0),nLevelMax(50),fVerbose(0),fMinProbability(1.e-10)
{
  fLevelEnergy = fNewEnergy = fDeltaEnergy = fNewTime 
    = fHalfLifeTime = fProbability = fICC = fx = 0.0;
  eGamma.resize(nLevelMax,0.0);
  wGamma.resize(nLevelMax,0.0);
  kICC.resize(nLevelMax,0.0);
  for(G4int i=0; i<30; ++i) { buffer[i] = 0; }
}

G4LevelReader::~G4LevelReader() 
{}

void G4LevelReader::FillLevels(G4int Z, G4int A,
			       std::vector<G4NucLevel*>* levels,
			       const G4String& filename) 
{ 
  std::ifstream inFile(filename);
  if (!inFile.is_open()) {
    if (fVerbose > 0) {
      G4cout << " G4LevelReader: nuclide (" 
	     << Z << "," << A 
	     << ") does not have a gamma levels file" << G4endl;
    }
    return;
  }

  // Read file with gamma data and fill levels
  fLevelEnergy = 0.0;
  nLevels = 0;
     
  // read next line
  while(Read(inFile)) {

    // create new level and start fill the next
    if(fNewEnergy != fLevelEnergy) {
      if(0 < nLevels) { MakeNewLevel(levels); }
      fLevelEnergy = fNewEnergy;
      fHalfLifeTime = fNewTime;
      nLevels = 0;
    }

    // fill data on a new daughter level
    eGamma[nLevels] = fDeltaEnergy*keV;
    wGamma[nLevels] = std::max(fProbability*0.01,fMinProbability);
    kICC[nLevels]   = fICC;
    ++nLevels;

    // check buffer size - should never happen
    if(nLevels > nLevelMax) {
      nLevelMax += 10;
      eGamma.resize(nLevelMax);
      wGamma.resize(nLevelMax);
      kICC.resize(nLevelMax);
    }
  }
  // end of reading
  if(0 < nLevels) {
    MakeNewLevel(levels);
    inFile.close();
  }
}

G4bool G4LevelReader::Read(std::ifstream& dataFile) 
{  
  // Each item will return iostream status
  return (ReadDataItem(dataFile, fNewEnergy) &&
	  ReadDataItem(dataFile, fDeltaEnergy) &&
	  ReadDataItem(dataFile, fProbability) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fNewTime) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fICC) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fx) &&
	  ReadDataItem(dataFile, fx) );
}

G4bool 
G4LevelReader::ReadDataItem(std::istream& dataFile, G4double& x) 
{
  // G4bool okay = (dataFile >> buffer) != 0;		// Get next token
  // if (okay) x = strtod(buffer, NULL);
  G4bool okay = true;
  dataFile >> buffer;
  if(dataFile.fail()) { okay = false; }
  else { x = strtod(buffer, NULL); }

  return okay;
}

void G4LevelReader::MakeNewLevel(std::vector<G4NucLevel*>* levels)
{
  // first normalize probabilities
  G4double norm = 0.0;
  for(size_t i=0; i<nLevels; ++i) { norm +=  wGamma[i]; }

  // should never happen
  if(norm <= 0.0) { return; }

  norm = 1.0/norm;
  for(size_t i=0; i<nLevels; ++i) { wGamma[i] *= norm; }

  // correct probabilities on ICC factor
  norm = 0.0;
  for(size_t i=0; i<nLevels; ++i) { 
    wGamma[i] /= (1.0 + kICC[i]);
    norm += wGamma[i]; 
  }
  norm = 1.0/norm;
  fHalfLifeTime *= norm*second;

  // cumulative sum
  if(1 == nLevels) {
    wGamma[0] = 1.0;
  } else if(2 == nLevels) {
    wGamma[0] *= norm;
    wGamma[1] = 1.0;
  } else {
    wGamma[0] *= norm;
    for(size_t i=1; i<nLevels-1; ++i) { 
      wGamma[i] = wGamma[i]*norm + wGamma[i-1];
    }
    wGamma[nLevels-1] = 1.0;
  }
  G4NucLevel* p = new G4NucLevel(fLevelEnergy, fHalfLifeTime,
				 eGamma, wGamma);
  levels->push_back(p);
  return;
}

