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
//
// $Id: Tst33Geometry.hh,v 1.2 2002-10-29 16:37:09 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33Geometry
//
// Class description:
//
// ...

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33VGeometry_hh
#define Tst33VGeometry_hh Tst33VGeometry_hh

#include "globals.hh"

class G4VPhysicalVolume;

class Tst33VGeometry {
public:
  Tst33VGeometry();
  virtual ~Tst33VGeometry();

  const G4VPhysicalVolume &GetWorldVolume() const;
  const G4VPhysicalVolume &GetPhysicalVolumeByName(const G4String& name)
    const ;
  G4String ListPhysNamesAsG4String() const ;
  
}


#endif
