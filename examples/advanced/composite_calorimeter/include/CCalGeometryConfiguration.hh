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
//              Use getInstance to retrieve it.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalGeometryConfiguration_h
#define CCalGeometryConfiguration_h 1

#include "g4std/map"
#include "CCalutils.hh"

class CCalGeometryConfiguration {

  struct GCInfo { 
    G4String FileName;
    int      ConstructFlag;
  };

  typedef G4std::map<G4String, GCInfo, G4std::less<G4String> > CCalGeometryConfTable;
  typedef G4std::map<G4String, GCInfo, G4std::less<G4String> >::iterator CCalGeometryConfIterator;
		    
  
public:
  ~CCalGeometryConfiguration() {}

  static CCalGeometryConfiguration* getInstance();

  int      getConstructFlag(const G4String& n) /*const*/;
  G4String getFileName(const G4String& n) /*const*/;

  CCalGeometryConfTable& getConfTable() { 
    return theConfiguration;
  }

private:
  CCalGeometryConfiguration();

private:
  static CCalGeometryConfiguration* instance;
  CCalGeometryConfTable theConfiguration;
};

#endif
