//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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

