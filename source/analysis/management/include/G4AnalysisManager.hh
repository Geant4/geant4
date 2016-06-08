// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4AnalysisManager.hh,v 1.8 2000/11/16 13:44:45 barrand Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// 
// Guy Barrand 20th Mai 2000
//
// Class description
//
//  Base class for the analysis manager proxy mechanism.
//

#ifndef G4ANALYSIS_HH
#define G4ANALYSIS_HH

#if defined(G4ANALYSIS_BUILD) || defined(G4ANALYSIS_USE)

#include "globals.hh"
#include "g4std/vector"
#include "G4VAnalysisManager.hh"

class G4AnalysisManager : public G4VAnalysisManager {
public: // I methods :
  G4bool RegisterAnalysisSystem(G4VAnalysisSystem*);
  IHistogramFactory* GetHistogramFactory(const G4String&);
  //
  virtual void Store(IHistogram* = 0,const G4String& = "");
  virtual void Plot(IHistogram*);
public:
  G4AnalysisManager();
  virtual ~G4AnalysisManager();
private:
  G4VAnalysisSystem* FindSystem(const G4String&);
  G4std::vector<G4VAnalysisSystem*> fSystems;
  G4VAnalysisSystem* fCurrentSystem;
  //
#ifdef G4ANALYSIS_NON_AIDA
public:
  ICloudFactory* GetCloudFactory(const G4String&);
  virtual void Plot(ICloud*);
  ITuple* CreateTuple(const G4String&,const G4String&,const G4String&);
#endif
};

#endif

#endif
