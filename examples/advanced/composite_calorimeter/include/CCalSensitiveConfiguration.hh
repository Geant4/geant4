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
//              g4geometry.conf 
//              Use getInstance to retrieve it.
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSensitiveConfiguration_h
#define CCalSensitiveConfiguration_h 1

#include "g4std/map"
#include "CCalutils.hh"

class CCalSensitiveConfiguration {

  struct GCInfo { 
    G4String FileName;
    int      SensitiveFlag;
  };

   typedef G4std::map<G4String, GCInfo, G4std::less<G4String> > CCalSensitiveConfTable;
   typedef G4std::map<G4String, GCInfo, G4std::less<G4String> >::iterator CCalSensitiveConfIterator;
		    
  
public:
  ~CCalSensitiveConfiguration() {}

  static CCalSensitiveConfiguration* getInstance();

  int      getSensitiveFlag(const G4String& n) /*const*/;
  G4String getFileName(const G4String& n) /*const*/;

private:
  CCalSensitiveConfiguration();

private:
  static CCalSensitiveConfiguration* instance;
  CCalSensitiveConfTable theConfiguration;
};

#endif
