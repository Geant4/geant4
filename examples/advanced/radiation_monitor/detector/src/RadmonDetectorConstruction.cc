//
// File name:     RadmonDetectorConstruction.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorConstruction.cc,v 1.2 2005-09-14 12:28:31 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorConstruction.hh"

#include "RadmonVDetectorLayout.hh"
#include "RadmonVDetectorEntitiesConstructorsFactory.hh"

#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"



                                                RadmonDetectorConstruction :: RadmonDetectorConstruction(RadmonVDetectorLayout * layout, RadmonVDetectorEntitiesConstructorsFactory * factory)
:
 detectorLayout(layout),
 constructorsFactory(factory),
 motherPhysicalVolume(0),
 motherLogicalVolume(0)
{
 if (detectorLayout==0)
  G4Exception("RadmonDetectorConstruction::RadmonDetectorConstruction: layout==0.");
 
 detectorLayout->AttachObserver(this);
}



                                                RadmonDetectorConstruction :: ~RadmonDetectorConstruction()
{
 detectorLayout->DetachObserver(this);

 Destruct();
 
 delete constructorsFactory;

 // TO BE DONE
 G4cout << "RadmonDetectorConstruction::~RadmonDetectorConstruction(): PLEASE CHECK" << G4endl;
}





G4VPhysicalVolume *                             RadmonDetectorConstruction :: Construct(void)
{
 // TO BE DONE
 G4cout << "RadmonDetectorConstruction::Construct(): NOT IMPLEMENTED YET" << G4endl;
 
 G4Material * al = new G4Material("Aluminum", 13., 26.98*g/mole, 2.700*g/cm3);
 
 G4Box * worldBox = new G4Box("WBox", 5.*m, 5.*m, 5.*m);
 G4LogicalVolume * worldLogicalVolume  = new G4LogicalVolume(worldBox, al, "WLog", 0, 0, 0);
 G4VPhysicalVolume * worldPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(), "WPhys", worldLogicalVolume, 0, false, 0);
          
 return worldPhysicalVolume;
}





void                                            RadmonDetectorConstruction :: OnLayoutChange(void)
{
 Destruct();
 G4RunManager::GetRunManager()->DefineWorldVolume(Construct(), true);
}




void                                            RadmonDetectorConstruction :: Destruct(void)
{
 // TO BE DONE
 G4cout << "RadmonDetectorConstruction::Destruct(): NOT IMPLEMENTED YET" << G4endl;
}
