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
///////////////////////////////////////////////////////////////////////////////
// File: CCalAnalysis.hh
// Description: CCalAnalysis is a singleton class and interfaces all user
//              analysis code
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalAnalysis_h 
#define CCalAnalysis_h 1

#include "globals.hh"

//! Defines the format which is used for the output file
//! default is a ROOT file. Comment the g4root and un-comment
//! one of the others, to change the output format
#include "g4root.hh"
//! Alternative here is a XML file
//#include "g4xml.hh"

class CCalAnalysis {
public:
  virtual ~CCalAnalysis();
  
public:
  void BeginOfRun(G4int n);
  void EndOfRun(G4int n);
  void EndOfEvent(G4int flag);

  void Init();
  void Finish();

  int maxbin() {return numberOfTimeSlices;}

  void InsertEnergyHcal(float*);
  void InsertEnergyEcal(float*);
  void InsertEnergy(float v);
  void InsertLateralProfile(float*);
  void InsertTime(float*); 
  void InsertTimeProfile(int, double, double); 

  void setNtuple(float* hcalE, float* ecalE, float elab, float x, float y, 
		 float z, float edep, float edec, float edhc);

  static CCalAnalysis* getInstance();

private:
  CCalAnalysis();

private:
  static CCalAnalysis* instance;

  G4int fVerbosity;

  G4int numberOfTimeSlices;
  //enum {numberOfTimeSlices = 200}; 

  //
  //Id of the histograms (first of each block). Used internally for a 
  //quick mapping of the proper histoID
  //
  G4int energy;
  G4int hcalE; //  28 hadronic modules
  G4int ecalE; //  49 crystal towers
  G4int timeHist; // 200 nanoseconds time window
  G4int lateralProfile;//  70 centimeters lateral window
  // (indeed 64 should be enough)
  G4int timeProfile; // at step and in sensitive detector
};


#endif
