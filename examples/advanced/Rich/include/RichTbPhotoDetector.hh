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
// RichTbPhotoDetector.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbPhotoDetector_h
#define RichTbPhotoDetector_h 1

#include "globals.hh"
#include "RichTbMaterial.hh"
#include "RichTbComponent.hh"
#include "RichTbHpd.hh"
#include "RichTbRunConfig.hh"
class RichTbPhotoDetector {
public:

  RichTbPhotoDetector();
  RichTbPhotoDetector(RichTbMaterial*, RichTbComponent*,  RichTbRunConfig*, G4bool);
  virtual ~RichTbPhotoDetector();
  G4int getNumberOfHpds()  {return NumOfHpds; }
  RichTbHpd* getRichHPD(G4int HpdNum){ return richHPD[HpdNum]; }
private:

  G4int NumOfHpds;
  RichTbHpd* richHPD[NumberOfHpds];
  G4bool ConstructTrackingGeometrySwitch;
};

#endif 

