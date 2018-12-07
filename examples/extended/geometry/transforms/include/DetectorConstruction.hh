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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "globals.hh"

class DetectorMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;

/// Detector construction class to demonstrate various ways of placement
///
/// The geometry setup consists of two trapezoid volumes which are placed 
/// in a world so that their axial symmetry axis is in given theta and phi 
/// polar angles. The various ways of placement are implemented in the 
/// DetectorConstruction class in the following private functions:
///  - PlaceWithDirectMatrix()
///  - PlaceWithInverseMatrix()
///  - PlaceWithAxialRotations()
///  - PlaceWithEulerAngles()
///  - PlaceWithReflections()
///
/// which are then called from the Construct() function.
/// All method defines exactly same geometry except for the placement 
/// with reflection where trapezoids are placed with their symmetry axis 
/// in parallel with z-axis in order to make easier to check reflection 
/// visually. 

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    enum EMethod {
      kWithDirectMatrix,
      kWithInverseMatrix,
      kWithAxialRotations,
      kWithEulerAngles,
      kWithReflections
    };  

  public:
    DetectorConstruction();
   ~DetectorConstruction();

  public:
    // methods from base class 
    virtual G4VPhysicalVolume* Construct();
    
    // set methods
    void SetMethod(EMethod method);
                       
  private:
    // methods
    void PlaceWithDirectMatrix();
    void PlaceWithInverseMatrix();
    void PlaceWithAxialRotations();
    void PlaceWithEulerAngles();
    void PlaceWithReflections();

    // data members
    DetectorMessenger* fMessenger;
    EMethod  fMethod;
    G4LogicalVolume* fWorldVolume;
    G4LogicalVolume* fTrdVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

