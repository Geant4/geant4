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
// $Id: RE02DetectorConstruction.hh,v 1.3 2006-11-18 01:37:23 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
     const{ nx=fNx; ny = fNy; nz = fNz; }
  // Insert Lead plate in water or simple homogeneous water phantom
  void SetLeadSegment(G4bool flag=TRUE){ fInsertLead = flag; }
  G4bool IsLeadSegment(){ return fInsertLead; }

private:
  // Data members
  G4ThreeVector fphantomSize;   // Size of Water Phantom
  G4int         fNx,fNy,fNz;    // Number of segmentation of water phantom.
  G4bool        fInsertLead;    // Flag for inserting lead plate in water phantom




};
#endif
