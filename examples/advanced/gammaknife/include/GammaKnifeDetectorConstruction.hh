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

#ifndef GammaKnifeDetectorConstruction_H
#define GammaKnifeDetectorConstruction_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Tubs;
class G4Cons;
class GammaKnifeBeamLine;
class GammaKnifeDetectorMessenger;


class GammaKnifeDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  GammaKnifeDetectorConstruction();

  ~GammaKnifeDetectorConstruction();

  G4VPhysicalVolume* Construct();  

private: 

  void ConstructBeamLine();
  // This method allows to define the beam line geometry in the
  // experimental set-up

 void UpdateHelmet();
 // Called after the change of helmetSize to properly set the
 // the collimator radii. Can be called before initialization
 // as well as between the runs.

public: 

 // Set helmet size and (if necessary) resize the collimator volume
 // (valid values: 4, 8, 14, 18)
 void SetHelmetSize(G4int);

 
private:

  G4VPhysicalVolume* physicalTreatmentRoom;

  G4VPhysicalVolume* patientPhysicalVolume;
  G4LogicalVolume* patientLogicalVolume;
  G4Cons* solidColl_helmet;          


  G4Tubs*            solidTube_source;    
     G4LogicalVolume*   logicTube_source;     
     G4VPhysicalVolume* physiTube_source;    

     G4Tubs*            solidTube;     
     G4LogicalVolume*   logicTube;     
     G4VPhysicalVolume* physiTube;    

     G4Tubs*            solidTube_Al;    
     G4LogicalVolume*   logicTube_Al;    
     G4VPhysicalVolume* physiTube_Al;     

     G4Tubs*            solidTube_Fe;     
     G4LogicalVolume*   logicTube_Fe;    
     G4VPhysicalVolume* physiTube_Fe;     

     G4int helmetSize; // helmet size (4, 8, 14, 18)

     G4Tubs* solidTube_post;           
     G4LogicalVolume* logicTube_post;  
     G4VPhysicalVolume* physiTube_post;

     G4Tubs* solidTube_coll;          
     G4LogicalVolume* logicTube_coll;  
     G4VPhysicalVolume* physiTube_coll;

     G4Tubs* solidTube_coll_Fe;           
     G4LogicalVolume* logicTube_coll_Fe;  
     G4VPhysicalVolume* physiTube_coll_Fe;

    
     G4Cons* solidColl_fixed;          
     G4LogicalVolume* logicColl_fixed;  
     G4VPhysicalVolume* physiColl_fixed;

     G4Cons* solidColl_fixed_Fe;          
     G4LogicalVolume* logicColl_fixed_Fe;  
     G4VPhysicalVolume* physiColl_fixed_Fe;

    
     G4LogicalVolume* logicColl_helmet;  
     G4VPhysicalVolume* physiColl_helmet;

     G4Cons* solidColl_helmet_Fe;           
     G4LogicalVolume* logicColl_helmet_Fe;  
     G4VPhysicalVolume* physiColl_helmet_Fe;


  GammaKnifeDetectorMessenger* detectorMessenger;

};
#endif
