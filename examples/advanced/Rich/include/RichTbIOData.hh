// Rich advanced example for Geant4
// RichTbIOData.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbIOData_h
#define RichTbIOData_h 1
#include "globals.hh"
#include <iostream>
#include <fstream.h>
#include "RichTbRunConfig.hh"
#include "G4Event.hh"
class RichTbIOData {

public:

  RichTbIOData(RichTbRunConfig* );
  virtual ~RichTbIOData();
  void WriteOutEventHeaderData(const G4Event* );
  void WriteOutHitData(const G4Event* );
  
private:
  ofstream OutputDataFS;
  G4String aOutFileString;
};
#endif

