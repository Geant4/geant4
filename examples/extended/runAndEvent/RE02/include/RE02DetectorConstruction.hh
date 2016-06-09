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
// $Id: RE02DetectorConstruction.hh,v 1.1 2005/11/24 01:44:18 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#ifndef RE02DetectorConstruction_h
#define RE02DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MultiFunctionalDetector.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

//
class RE02DetectorConstruction : public G4VUserDetectorConstruction
{
public:
   // constructor and destructor.
   RE02DetectorConstruction();
   virtual ~RE02DetectorConstruction();

public:
  // virtual method from G4VUserDetectorCOnstruction.
  virtual G4VPhysicalVolume* Construct();

public:
  // Get/Set Access methods for data members
  // Size of Whater Phantom
  void SetPhantomSize(G4ThreeVector size) { fphantomSize=size; }
  const G4ThreeVector& GetPhantomSize() const { return fphantomSize; }
  // Number of segments of water phantom
  void SetNumberOfSegmentsInPhantom(G4int nx, G4int ny, G4int nz) 
      { fNx=nx; fNy=ny; fNz=nz; }
  void GetNumberOfSegmentsInPhantom(G4int& nx, G4int& ny, G4int& nz) 
      const{ nx = fNx; ny = fNy; nz = fNz; }

private:
  // Data members
  G4ThreeVector fphantomSize;   // Size of Water Phantom
  G4int         fNx,fNy,fNz;    // Number of segmentation of water phantom.
};
#endif
