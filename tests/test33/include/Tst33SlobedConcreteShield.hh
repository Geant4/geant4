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
//
// $Id: Tst33SlobedConcreteShield.hh,v 1.8 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33SlobedConcreteShield
//
// Class description:
//
// The mass goeometry with cells.For usage as sampler geometry 
// without parallel geometry.
// 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33SlobedConcreteShield_hh
#define Tst33SlobedConcreteShield_hh Tst33SlobedConcreteShield_hh

#include "Tst33VGeometry.hh"
#include "Tst33PVolumeStore.hh"
#include "Tst33MaterialFactory.hh"

#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume; //ASO
class G4VIStore;
class G4VWeightWindowStore;


class Tst33SlobedConcreteShield : public Tst33VGeometry {
public:
  Tst33SlobedConcreteShield();
  virtual ~Tst33SlobedConcreteShield();

  virtual G4VPhysicalVolume &GetWorldVolumeAddress() const;
  virtual G4VPhysicalVolume *GetWorldVolume();

  virtual G4GeometryCell GetGeometryCell(G4int i, const G4String &) const; 

  void SetSensitive();  //ASO

  void Construct();

  G4String GetCellName(G4int i);
    
private:
  Tst33SlobedConcreteShield(const Tst33SlobedConcreteShield &);

  //xtest  void Construct();

  Tst33SlobedConcreteShield &operator=(const Tst33SlobedConcreteShield &);

  Tst33MaterialFactory fMaterialFactory;
  G4VPhysicalVolume *fWorldVolume;
  Tst33PVolumeStore fPVolumeStore;

  G4Material *fConcrete;
  G4Material *fGalactic;

  std::vector< G4VPhysicalVolume * > fPhysicalVolumeVector;
  std::vector< G4LogicalVolume * > fLogicalVolumeVector;  //ASO

};



#endif
