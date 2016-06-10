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
// $Id: G4LevelReader.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4NucLevel
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
// Helper class to read Geant4 nuclear level database ignoring 
// information on electron internal conversion probabilities
// 

#ifndef G4LEVELREADER_HH
#define G4LEVELREADER_HH 1

#include "globals.hh"
#include "G4NucLevel.hh"
#include <vector>
#include <fstream>

class G4LevelReader 
{

public:

  G4LevelReader();

  ~G4LevelReader();
  
  void FillLevels(G4int Z, G4int A,
		  std::vector<G4NucLevel*>* levels,
		  const G4String& filename); 

  inline void SetVerbose(G4int val);
  
private:

  G4bool Read(std::ifstream& aDataFile);

  G4bool ReadDataItem(std::istream& dataFile, G4double& x);

  void MakeNewLevel(std::vector<G4NucLevel*>* levels);
  
  G4LevelReader(const G4LevelReader & right);  
  const G4LevelReader& operator=(const G4LevelReader &right);
  G4bool operator==(const G4LevelReader &right) const;
  G4bool operator!=(const G4LevelReader &right) const;

  size_t   nLevels;
  size_t   nLevelMax;
  G4int    fVerbose;
  G4double fMinProbability;
  G4double fLevelEnergy;
  G4double fNewEnergy;
  G4double fDeltaEnergy;
  G4double fNewTime;
  G4double fHalfLifeTime;
  G4double fProbability;
  G4double fICC;
  G4double fx;

  std::vector<G4double> eGamma;
  std::vector<G4double> wGamma;
  std::vector<G4double> kICC;

  // Buffers for reading data file
  char buffer[30];	

};

inline void G4LevelReader::SetVerbose(G4int val)
{
  fVerbose = val;
}


#endif
