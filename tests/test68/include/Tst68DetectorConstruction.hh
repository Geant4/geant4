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
#ifndef Tst68DetectorConstruction_H
#define Tst68DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"       

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4FieldManager;
class G4UniformMagField;
class G4Material;
class Tst68DetectorMessenger;


class Tst68DetectorConstruction : public G4VUserDetectorConstruction {

public:

  Tst68DetectorConstruction();
  ~Tst68DetectorConstruction();
  
  G4VPhysicalVolume* Construct();
  void ConstructSDandField();

  void SetMagField( const G4double fieldValue );
  void SetAbsorberMaterial( const G4String name );
  void SetActiveMaterial( const G4String name );
  // Use by the messenger.

  inline G4Material* GetAbsorberMaterial() const;
  inline G4Material* GetActiveMaterial() const;

  inline void SetIsCalHomogeneous( const G4bool choice );
  inline void SetIsUnitInLambda( const G4bool choice );
  inline void SetAbsorberTotalLength( const G4double value );
  inline void SetCalorimeterRadius( const G4double value );
  inline void SetActiveLayerNumber( const G4int value );
  inline void SetActiveLayerSize( const G4double value );
  inline void SetReadoutLayerNumber( const G4int value );
  // To define the calorimeter geometry.

  inline void SetIsRadiusUnitInLambda( const G4bool choice );
  inline void SetRadiusBinSize( const G4double value );
  inline void SetRadiusBinNumber( const G4int value );
  // To define the transverse shower analysis.
  
  void UpdateGeometry();

private:

  void DefineMaterials();
  // Define all the materials.

  G4VPhysicalVolume* ConstructCalorimeter();     
  // To be invoked each time the geometry needs to be updated.

  G4bool areParametersOK();
  // Return true if all the parameters are sensible, false otherwise.

  void PrintParameters();
  // Print the various parameters which define the calorimeter.

  G4Material* Vacuum;
  G4Material* Iron;
  G4Material* Copper;
  G4Material* Tungsten;
  G4Material* Lead;
  G4Material* Uranium;
  G4Material* PbWO4;
  G4Material* Polystyrene;
  G4Material* LiquidArgon;
  G4Material* Silicon;
  G4Material* Quartz;
  G4Material* Brass;
  G4Material* Aluminium;
  G4Material* Graphite;
  G4Material* theAbsorberMaterial;
  G4Material* theActiveMaterial;
  
  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;
  // World envelope. 
  
  G4LogicalVolume*  logicCalo;
  G4VPhysicalVolume* physiCalo;
  // "Calorimeter".
  
  G4LogicalVolume*  logicModule;
  G4VPhysicalVolume* physiModule;
  // Module of the "calorimeter".
  
  G4LogicalVolume*  logicAbsorber;
  G4VPhysicalVolume* physiAbsorber;
  // Absorber layer of the "calorimeter".
  
  G4LogicalVolume*  logicActive;
  G4VPhysicalVolume* physiActive;
  // Active layer of the "calorimeter".

  G4FieldManager* fieldMgr;
  // Pointer to the field manager.

  G4UniformMagField* uniformMagField; 
  // Pointer to the uniform magnetic field.
  
  Tst68DetectorMessenger* detectorMessenger;
  // Pointer to the Messenger.

  G4bool theIsCalHomogeneous; 
  // If false then Sampling calorimeter;
  // If true  then Homogeneous calorimeter.

  G4bool theIsUnitInLambda;
  // If false then normal unit of length to express the absorber total length.
  // If true  then lambda (interaction length) to express the absorber total length.

  G4double theAbsorberTotalLength;
  // This is the total length of the absorber material, expressed 
  // in unit of length (e.g. m, cm, mm) if theIsUnitInLambda is false, 
  // otherwise in number of lambdas (interaction lengths).
  // Notice that in the case of a sampling calorimeter (i.e. 
  // theIsCalHomogeneous is false), the active layers are not counted;
  // in the case of an homogenous calorimeter, this length account
  // for the overall dimension of the calorimeter.

  G4double theCalorimeterRadius;
  // This is the radius of the calorimeter which is a cylinder, expressed 
  // in unit of length (e.g. m, cm, mm) if theIsUnitInLambda is false, 
  // otherwise in number of lambdas (interaction lengths) of the absorber.

  G4int theActiveLayerNumber;
  G4double theActiveLayerSize;
  G4int theReadoutLayerNumber;
  // Number of active layers and length of each of them (in normal unit
  // of length, e.g. mm): in the case of sampling calorimeter
  // (i.e. theIsCalHomogeneous is false) the medium is theActiveMaterial;
  // in the case of an homogeneous calorimeter, the "active layers" are
  // only a fictitious way to sample the longitudinal energy deposits,
  // but they are actually made of the same absorber material, and their
  // thickness is taken into account in theAbsorberTotalLength.
  // Usually not all the active layers are read independently, but they
  // are combined: the number of actual readout layers must be specified
  // as well. If such value is not explicitly specified by the user, then
  // the number of readout layers is assumed to be the same as the number
  // of active layer. A warning is printed out if the specified number
  // of readout layers is not a divisor of the number of active layers.
  // and then the readout layer number is forced to be the same as the
  // active layer number.

  G4bool theIsRadiusUnitInLambda;
  // If false then normal unit of length to express the radius bin size.
  // If true  then lambda (interaction length of the absorber) to express 
  // the radius bin size.

  G4double theRadiusBinSize;
  G4int theRadiusBinNumber;
  // Bin size and number of bins in the radius, for the transverse shower
  // analysis. The bin size is expressed in unit of either lambda 
  // (interaction length of the absorber) or normal [mm] unit, according
  // to the value specified in theIsRadiusUnitInLambda.
  // The radius bin  i-th , where  0 <= i <= theRadiusBinNumber-1 ,
  // has the following content:
  // for all bins except the last one, i.e. 0 <= i < theRadiusBinNumber-1 
  // the content is the sum of all the energy deposited at a radius :  
  //     i * theRadiusBinSize  <=  radius  <  (i+1) * theRadiusBinSize  
  // with respect to the beam axis (i.e. y = 0.0, z = 0.0, because there
  // is no smearing around the nominal beam axis, which is the x-axis).
  // The last bin,  i = theRadiusBinNumber-1 , the content is the sum
  // of all the energy deposited at a radius :
  //            radius  >=  (theRadiusBinNumber - 1) * theRadiusBinSize
  // In other words, the last bin collects all the energy deposited 
  // (inside the calorimeter) at a radius greater than the ones so far
  // considered.

};


inline G4Material* Tst68DetectorConstruction::
GetAbsorberMaterial() const {
  return theAbsorberMaterial;
}


inline G4Material* Tst68DetectorConstruction::
GetActiveMaterial() const {
  return theActiveMaterial;
}


inline void Tst68DetectorConstruction::
SetIsCalHomogeneous( const G4bool choice ) {
  theIsCalHomogeneous = choice;
}

inline void Tst68DetectorConstruction::
SetIsUnitInLambda( const G4bool choice ) {
  theIsUnitInLambda = choice;
}

inline void Tst68DetectorConstruction::
SetAbsorberTotalLength( const G4double value ) {
  theAbsorberTotalLength = value;
}

inline void Tst68DetectorConstruction::
SetCalorimeterRadius( const G4double value ) {
  theCalorimeterRadius = value;
}

inline void Tst68DetectorConstruction::
SetActiveLayerNumber( const G4int value ) {
  theActiveLayerNumber = value;
  // By default, if the user does not set explicitly the number of
  // readout layers, this is assumed to be the number of active layers.
  theReadoutLayerNumber = theActiveLayerNumber; 
}

inline void Tst68DetectorConstruction::
SetActiveLayerSize( const G4double value ) {
  theActiveLayerSize = value;
}

inline void Tst68DetectorConstruction::
SetReadoutLayerNumber( const G4int value ) {
  theReadoutLayerNumber = value;
  if ( theActiveLayerNumber % theReadoutLayerNumber != 0 ) {
    G4cout << "***WARNING*** Tst68DetectorConstruction::SetReadoutLayerNumber"
           << G4endl << "\t theReadoutLayerNumber = " 
           << theReadoutLayerNumber
           << " is NOT compatible with theActiveLayerNumber = " 
           << theActiveLayerNumber 
	   << G4endl
           << "\t Forced theReadoutLayerNumber = theActiveLayerNumber" 
           << G4endl;
    theReadoutLayerNumber = theActiveLayerNumber;
  }
}


inline void Tst68DetectorConstruction::
SetIsRadiusUnitInLambda( const G4bool choice ) {
  theIsRadiusUnitInLambda = choice;
}

inline void Tst68DetectorConstruction::
SetRadiusBinSize( const G4double value ) {
  theRadiusBinSize = value;
}

inline void Tst68DetectorConstruction::
SetRadiusBinNumber( const G4int value ) {
  theRadiusBinNumber = value;
}


#endif
