// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: XrayTelAnalysisManager.hh,v 1.1 2000-11-10 13:35:15 gbarrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// Example Analysis Manager implementing virtual function
//   RegisterAnalysisSystems.  Exploits C-pre-processor variables
//   G4ANALYSIS_USE_JAS, etc., which are set by the GNUmakefiles if
//   environment variables of the same name are set.

// So all you have to do is set environment variables and compile and
//   instantiate this in your main().

#ifndef XrayTelAnalysisManager_h
#define XrayTelAnalysisManager_h 1

#ifdef G4ANALYSIS_USE

#include "G4AnalysisManager.hh"

#include "g4std/vector"
#include "G4ThreeVector.hh"                                                                         

class G4Event;
class G4SteppingManager;

class IHistogram1D;
class IHistogram2D;

#ifdef G4ANALYSIS_NON_AIDA
class ICloud1D;
class ICloud2D;
#include <ITuple.h>
#endif

class XrayTelAnalysisManager: public G4AnalysisManager {
public:
  XrayTelAnalysisManager(const G4String&);
  ~XrayTelAnalysisManager();
public:
  void BeginOfRun(); 
  void EndOfRun(); 
  void BeginOfEvent(); 
  void EndOfEvent(const G4Event*,const G4String&);
  void Step(const G4SteppingManager*);
private:
  G4bool fDrawEvent;
  G4std::vector<G4double> fEnteringEnergy;
  G4std::vector<G4ThreeVector> fEnteringDirection;
  //
  IHistogram1D* fEnteringEnergyHistogram;
  IHistogram2D* fYZ_Histogram;
  //
  // Under study...
#ifdef G4ANALYSIS_NON_AIDA
  ICloud1D* fKinetic;
  ICloud2D* fYZ;
  ITuple* fTuple;
#endif
};

#endif

#endif
