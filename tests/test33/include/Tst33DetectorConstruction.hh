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
// $Id: Tst33DetectorConstruction.hh,v 1.2 2002-11-20 09:38:25 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class Tst33DetectorConstruction
//
// Class description:
// This class only takes a pointer to the world volume and gives it
// to the Geant4 kernel when asked for it. The detector construction
// is done elsewhere.

#ifndef Tst33DetectorConstruction_hh
#define Tst33DetectorConstruction_hh Tst33DetectorConstruction_hh

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class Tst33DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  explicit Tst33DetectorConstruction();
  virtual ~Tst33DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();
  void SetWorldVolume(G4VPhysicalVolume *v);
  
private:
  Tst33DetectorConstruction(const Tst33DetectorConstruction &);
  Tst33DetectorConstruction &operator=(const Tst33DetectorConstruction &);
  G4VPhysicalVolume* fWorldVolume;
};

#endif
