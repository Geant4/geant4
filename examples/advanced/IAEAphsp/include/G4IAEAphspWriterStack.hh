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
// Author: M.A. Cortes-Giraldo
//
// 2025-08-27: Its objects work in local runs, their mission is to store the
//   info of the particles to be written in the IAEAphsp output files at
//   the end of the run. In other words, they constitute a stack for
//   IAEAphsp particles until the run finishes, when the IAEAphsp file is
//   actually written.
//

#ifndef G4IAEAphspWriterStack_hh
#define G4IAEAphspWriterStack_hh 1


#include "globals.hh"
#include "G4ThreeVector.hh"

#include <set>
#include <vector>

class G4IAEAphspWriter;
class G4Step;


class G4IAEAphspWriterStack
{

public:

  G4IAEAphspWriterStack(const G4String filename);
  ~G4IAEAphspWriterStack();

  void AddZphsp(const G4double zphsp);
  void ClearZphspVec();
  void SetDataFromWriter(const G4IAEAphspWriter* );
  void PrepareRun();
  void PrepareNextEvent();
  void StoreParticleIfEligible(const G4Step*);
  void ClearRunVectors();

  void SetFileName(const G4String name)   { fFileName = name; }

  const G4String GetFileName() const                       {return fFileName;}
  const std::vector<G4double>* GetZphspVec() const         {return fZphspVec;}
  std::vector<std::vector<G4int>* >* GetPDGMtrx() const    {return fPDGMtrx;}
  std::vector<std::vector<G4ThreeVector>* >* GetPosMtrx() const
  {return fPosMtrx;}
  std::vector<std::vector<G4ThreeVector>* >* GetMomMtrx() const
  {return fMomMtrx;}
  std::vector<std::vector<G4double>* >* GetEneMtrx() const {return fEneMtrx;}
  std::vector<std::vector<G4double>* >* GetWtMtrx() const  {return fWtMtrx;}
  std::vector<std::vector<G4int>* >* GetNstatMtrx() const  {return fNstatMtrx;}

  
private:

  G4IAEAphspWriterStack() = default;
  void StoreIAEAParticle(const G4Step* aStep, const G4int zStopIdx,
			 const G4int pdgCode);


  // ------------
  // DATA MEMBERS
  // ------------
  
  // FILE PROPERTIES

  G4String fFileName;
  // Must include the path but not any of the IAEA extensions.
  // (This is set from G4IAEAphspWriter)

  std::vector<G4double>* fZphspVec = nullptr;
  // Vector storing the z-value of the phsp planes.

  // COUNTERS & TAGS

  std::vector<G4int>* fIncrNumberVec = nullptr;
  // Book-keeping of the number of previous events without having particles
  // crossing the phsp plane.
  // (i.e., the current incremental history number, or n_stat, of each phsp)

  std::vector< std::set<G4int>* >* fPassingTracksVec = nullptr;
  // Each set is meant to store the track ID of every particle
  // crossing one of the planes during an event.
  // This is done to avoid registering multiple crosses in the phsp file.

  // INFORMATION STORED DURING RUN

  // First component is the phsp plane, according to registration ordering.
  // Second components are the dynamic variables.
  std::vector< std::vector<G4int>* >*         fPDGMtrx = nullptr;
  std::vector< std::vector<G4ThreeVector>* >* fPosMtrx = nullptr;
  std::vector< std::vector<G4ThreeVector>* >* fMomMtrx = nullptr;
  std::vector< std::vector<G4double>* >*      fEneMtrx = nullptr;
  std::vector< std::vector<G4double>* >*      fWtMtrx = nullptr;
  std::vector< std::vector<G4int>* >*         fNstatMtrx = nullptr;

};

#endif
