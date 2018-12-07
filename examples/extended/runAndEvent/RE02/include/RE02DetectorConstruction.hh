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
/// \file runAndEvent/RE02/include/RE02DetectorConstruction.hh
/// \brief Definition of the RE02DetectorConstruction class
//
//
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
/// Uer detector construction class
///
/// (Description)
///
///    Detector construction for example RE02.
///   
///  [Geometry] 
///    The world volume is defined as 200 cm x 200 cm x 200 cm box with Air.
///  Water phantom is defined as  200 mm x 200 mm x 400 mm box with Water.
///  The water phantom is divided into 100 segments in x,y plane using
///  replication,
///  and then divided into 200 segments perpendicular to z axis using nested 
///  parameterised volume.  
///   These values are defined at constructor,
///   e.g. the size of water phantom (fPhantomSize), and number of segmentation
///  of water phantom (fNx, fNy, fNz).
///
///  By default, lead plates are inserted into the position of even order 
///  segments.
///  NIST database is used for materials.
///
///
///  [Scorer]
///   Assignment of G4MultiFunctionalDetector and G4PrimitiveScorer 
///  is demonstrated in this example.
///      -------------------------------------------------
///      The collection names of defined Primitives are
///       0       PhantomSD/totalEDep 
///       1       PhantomSD/protonEDep
///       2       PhantomSD/protonNStep
///       3       PhantomSD/chargedPassCellFlux
///       4       PhantomSD/chargedCellFlux 
///       5       PhantomSD/chargedSurfFlux 
///       6       PhantomSD/gammaSurfCurr000
///       7       PhantomSD/gammaSurfCurr001
///       9       PhantomSD/gammaSurdCurr002
///      10       PhantomSD/gammaSurdCurr003
///     -------------------------------------------------
///     Please see README for detail description.
///
///
/// - G4VPhysicalVolume* Construct()
///     retrieves material from NIST database,
///     constructs a water phantom "phantom" in the world volume "world" and
///     sets detector sensitivities with G4MultiFunctionalDetector
///
/// - void SetPhantomSize(G4ThreeVector size)
///     sets the water phantom size which is defined in G4Box
///
/// - const G4ThreeVector& GetPhantomSize() const
///     gets the water phantom size
///
/// - void SetNumberOfSegmentsInPhantom(G4int nx, G4int ny, G4int nz) 
///     sets the number of segments of the water phantom
///
/// - void GetNumberOfSegmentsInPhantom(G4int& nx, G4int& ny, G4int& nz) 
///     gets the number of segments of the water phantom
///
/// - void SetLeadSegment(G4bool flag=TRUE)
///     selects whether insert or not Lead plate in water or simple homogeneous
///     water phantom
///
/// - G4bool IsLeadSegment()
///     returns whether insert or not Lead plate
//
class RE02DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  // constructor and destructor.
  RE02DetectorConstruction();
  virtual ~RE02DetectorConstruction();

public:
  // virtual method from G4VUserDetectorConstruction.
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

public:
  // Get/Set Access methods for data members
  // Size of Whater Phantom
  void SetPhantomSize(G4ThreeVector size) { fPhantomSize=size; }
  const G4ThreeVector& GetPhantomSize() const { return fPhantomSize; }
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
  G4ThreeVector fPhantomSize;   // Size of Water Phantom
  G4int         fNx,fNy,fNz;    // Number of segmentation of water phantom.
  G4bool        fInsertLead;    // Flag for inserting lead plate in water phantom
  G4LogicalVolume* fLVPhantomSens;

};
#endif
