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
// $Id: G4LevelReader.hh 88382 2015-02-17 10:49:24Z vnivanch $
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4LevelReader
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
// Helper class to read Geant4 nuclear level data  
// 

#ifndef G4LEVELREADER_HH
#define G4LEVELREADER_HH 1

#include "globals.hh"
#include "G4LevelManager.hh"
#include <iosfwd>

class G4LevelReader 
{

public:

  G4LevelReader();

  ~G4LevelReader();

  // create run manager using G4LEVELGAMMADATA data 
  const G4LevelManager* CreateLevelManager(G4int Z, G4int A);

  // create run manager using whatever data
  const G4LevelManager* MakeLevelManager(G4int Z, G4int A,
					 const G4String& filename);
  
  inline void SetVerbose(G4int val);
  
private:

  G4bool ReadDataItem(std::istream& dataFile, G4double& x);

  G4bool ReadDataItem(std::istream& dataFile, G4String& x);
  
  const std::vector<G4float>* NormalizedICCProbability(G4int Z);

  G4LevelReader(const G4LevelReader & right);  
  const G4LevelReader& operator=(const G4LevelReader &right);
  G4bool operator==(const G4LevelReader &right) const;
  G4bool operator!=(const G4LevelReader &right) const;

  G4double fMinProbability;
  G4double fTimeFactor;

  G4double fEnergy;
  G4double fCurrEnergy;
  G4double fTrEnergy;
  G4double fProb;
  G4double fTime;
  G4double fSpin;
  G4double fAlpha;
  G4double fNorm1;
  G4double fNorm2;
  G4double fNorm3;
  G4double fICC[10];

  static G4String fTrans[10];
  G4String fDirectory;
  G4String fFile;
  G4String fPol;

  char buffer[20];
  char bufp[2];

  G4int    fVerbose;

  std::vector<G4float> vEnergy;
  std::vector<G4float> vTime;
  std::vector<G4float> vTimeg;
  std::vector<G4int>   vSpin;
  std::vector<const G4NucLevel*> vLevel;

  std::vector<G4float>  vTransEnergy;
  std::vector<G4float>  vGammaCumProbability;
  std::vector<G4float>  vGammaECumProbability;
  std::vector<G4float>  vGammaProbability;
  std::vector<G4int>    vTrans;
  std::vector<const std::vector<G4float>*> vShellProbability;
};

inline void G4LevelReader::SetVerbose(G4int val)
{
  fVerbose = val;
}

#endif
