// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LizardSystem.hh,v 1.5 2000/11/16 13:44:41 barrand Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// 
// Guy Barrand 14 September 2000

#ifndef G4LIZARDSYSTEM_H
#define G4LIZARDSYSTEM_H

#if defined(G4ANALYSIS_BUILD_LIZARD) || defined(G4ANALYSIS_USE_LIZARD)

#include "G4VAnalysisSystem.hh"

class IHistogramFactory;
class IVectorFactory;
class IPlotter;

class G4LizardSystem : public G4VAnalysisSystem {
public: // I methods :
  virtual const G4String& GetName() const;
  virtual IHistogramFactory* GetHistogramFactory();
  //
  virtual void Store(IHistogram* = 0,const G4String& = "");
  virtual void Plot(IHistogram*);
public:
  G4LizardSystem(const G4String& name = "Lizard");
  virtual ~G4LizardSystem();
private:
  G4String fName;
  IHistogramFactory* fHistogramFactory;
  IVectorFactory* fVectorFactory;
  IPlotter* fPlotter;
};

#endif

#endif

