// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LabSystem.hh,v 1.2 2000/11/16 13:44:38 barrand Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// 
// Guy Barrand 14 September 2000

#ifndef G4LABSYSTEM_H
#define G4LABSYSTEM_H

#if defined(G4ANALYSIS_BUILD_LAB) || defined(G4ANALYSIS_USE_LAB)

#include "G4VAnalysisSystem.hh"

class IHistogramManager;
class IStorageSystemManager;
class ICloud;
class ICloudFactory;
class ICloudManager;
class ITuple;
class ITupleManager;
class ISession;

class G4LabSystem : public G4VAnalysisSystem {
public: // I methods :
  virtual const G4String& GetName() const;
  virtual IHistogramFactory* GetHistogramFactory();
  virtual ICloudFactory* GetCloudFactory();
  virtual ITuple* CreateTuple(const G4String&,const G4String&);
  //
  virtual void Store(IHistogram* = 0,const G4String& = "");
  virtual void Plot(IHistogram*);
  virtual void Plot(ICloud*);
public:
  G4LabSystem(const G4String& name = "Lab");
  virtual ~G4LabSystem();
private:
  G4String fName;
  IStorageSystemManager* fStorageManager;
  IHistogramManager* fHistogramManager;
  ICloudManager* fCloudManager;
  ITupleManager* fTupleManager;
  ISession* fSession;
};

#endif

#endif

