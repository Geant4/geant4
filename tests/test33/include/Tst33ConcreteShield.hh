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
// $Id: Tst33ConcreteShield.hh,v 1.4 2002-11-20 13:09:15 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33ConcreteShield
//
// Class description:
//
// A block of concrete without cells to be used as mass geometry
// if the cells are given in the parallel geometry.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33ConcreteShield_hh
#define Tst33ConcreteShield_hh Tst33ConcreteShield_hh

#include "Tst33VGeometry.hh"
#include "Tst33PVolumeStore.hh"
#include "Tst33MaterialFactory.hh"


class Tst33ConcreteShield : public Tst33VGeometry {
public:
  Tst33ConcreteShield();
  virtual ~Tst33ConcreteShield();

  virtual G4VPhysicalVolume &GetWorldVolume() const;
  virtual G4GeometryCell GetGeometryCell(G4int i) const; 

  
private:
  Tst33ConcreteShield(const Tst33ConcreteShield &);
  Tst33ConcreteShield &operator=(const Tst33ConcreteShield &);
  void Construct();
  Tst33MaterialFactory fMaterialFactory;
  G4VPhysicalVolume *fWorldVolume;
  Tst33PVolumeStore fPVolumeStore;

  G4Material *fConcrete;
  G4Material *fGalactic;

};



#endif
