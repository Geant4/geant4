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
// $Id: B4bRunData.hh 100946 2016-11-03 11:28:08Z gcosmo $
// 
/// \file B4bRunData.hh
/// \brief Definition of the B4bRunData class

#ifndef B4bRunData_h
#define B4bRunData_h 1

#include "G4Run.hh"
#include "globals.hh"

#include <array>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int kAbs = 0;
const G4int kGap = 1;
const G4int kDim = 2;

///  Run data class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers.
/// 
/// In order to reduce the number of data members a 2-dimensions array 
/// is introduced for each quantity:
/// - fEdep[], fTrackLength[].
///
/// The data are collected step by step in B4bSteppingAction, and
/// the accumulated values are filled in histograms and entuple
/// event by event in B4EventAction.

class B4bRunData : public G4Run
{
public:
  B4bRunData();
  virtual ~B4bRunData();

  void Add(G4int id, G4double de, G4double dl);
  void FillPerEvent();

  void Reset();

  // Get methods
  G4String  GetVolumeName(G4int id) const;
  G4double  GetEdep(G4int id) const;
  G4double  GetTrackLength(G4int id) const; 

private:
  std::array<G4String, kDim>  fVolumeNames;
  std::array<G4double, kDim>  fEdep;
  std::array<G4double, kDim>  fTrackLength; 
};

// inline functions

inline void B4bRunData::Add(G4int id, G4double de, G4double dl) {
  fEdep[id] += de; 
  fTrackLength[id] += dl;
}

inline G4String  B4bRunData::GetVolumeName(G4int id) const {
  return fVolumeNames[id];
}

inline G4double  B4bRunData::GetEdep(G4int id) const {
  return fEdep[id];
}   

inline G4double  B4bRunData::GetTrackLength(G4int id) const {
  return fTrackLength[id];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

