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
// File: CCalGeometryConfiguration.hh
// Description: This singleton holds the information in the file geometry.conf 
//              Use getInstance to retrieve it.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalGeometryConfiguration_h
#define CCalGeometryConfiguration_h 1

#include <map>
#include "CCalutils.hh"

class CCalGeometryConfiguration
{
  struct GCInfo { 
    G4String FileName;
    G4int      ConstructFlag;
  };

  typedef std::map<G4String, GCInfo, std::less<G4String> > CCalGeometryConfTable;
  
public:
  ~CCalGeometryConfiguration() {}

  static CCalGeometryConfiguration* getInstance();

  G4int getConstructFlag(const G4String& n) /*const*/;
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
