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
// Author: A. Pfeiffer (Andreas.Pfeiffer@cern.ch) 
//         
//
// History:
// -----------
// 19 Mar 2013   LP   Migrated to G4AnalysisManager
//  7 Nov 2001   MGP  Implemented according to A. Pfeiffer's instructions
//
// -------------------------------------------------------------------
// Class description:
// Example of analysis in a simulation application (histograms, ntuples etc.)
// This class follows the singleton design pattern; 
// it is responsible for the analysis management and algorithms 
//
// -------------------------------------------------------------------

#ifndef G4PROCESSTESTANALYSIS_HH
#define G4PROCESSTESTANALYSIS_HH

#include "globals.hh"
#include "G4ios.hh"

// uncomment g4root.hh and comment g4xml.hh for a ROOT-based output file

#include "g4root.hh"
//#include "g4xml.hh"

class G4Track;

class XrayTelAnalysis
{
public:

  ~XrayTelAnalysis();

  void book(G4bool isMaster);
  
  void finish(G4bool isMaster);
  
  void analyseStepping(const G4Track& track, G4bool entering);

  static XrayTelAnalysis* getInstance();

  void Update(G4double energy,G4int threadID);

private:

  XrayTelAnalysis();

  static XrayTelAnalysis* instance;

  // Quantities for the ntuple
  G4double eKin;
  G4double x;
  G4double y;
  G4double z;
  G4double dirX;
  G4double dirY;
  G4double dirZ;

  G4String asciiFileName;
  G4String histFileName;
  
  std::ofstream *asciiFile;

  //global counters: log separately for each thread (or sequential)
  std::map<G4int,G4int> *nEnteringTracks;
  std::map<G4int,G4double> *totEnteringEnergy;

};

#endif 
