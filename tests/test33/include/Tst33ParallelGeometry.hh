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
// $Id: Tst33ParallelGeometry.hh,v 1.3 2002-11-20 09:38:25 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33ParallelGeometry
//
// Class description:
//
// ...

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33ParallelGeometry_hh
#define Tst33ParallelGeometry_hh Tst33ParallelGeometry_hh

#include "Tst33VGeometry.hh"
#include "Tst33PVolumeStore.hh"
#include "Tst33MaterialFactory.hh"


class Tst33ParallelGeometry : public Tst33VGeometry {
public:
  Tst33ParallelGeometry();
  virtual ~Tst33ParallelGeometry();

  virtual G4VPhysicalVolume &GetWorldVolume() const;

  virtual G4GeometryCell GetGeometryCell(G4int i) const; 

private:
  Tst33ParallelGeometry(const Tst33ParallelGeometry &);

  void Construct();

  Tst33ParallelGeometry &operator=(const Tst33ParallelGeometry &);

  Tst33MaterialFactory fMaterialFactory;
  G4VPhysicalVolume *fWorldVolume;
  Tst33PVolumeStore fPVolumeStore;

  G4Material *fGalactic;

};



#endif
