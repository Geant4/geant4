//
// File name:     RadmonDetectorConstruction.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorConstruction.cc,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorConstruction.hh"

#include "RadmonVDetectorLayout.hh"
#include "RadmonVDetectorEntitiesConstructorsFactory.hh"

#include "G4RunManager.hh"



                                                RadmonDetectorConstruction :: RadmonDetectorConstruction(RadmonVDetectorLayout * layout, RadmonVDetectorEntitiesConstructorsFactory * factory)
:
 detectorLayout(layout),
 constructorsFactory(factory),
 motherPhysicalVolume(0),
 motherLogicalVolume(0)
{
}



                                                RadmonDetectorConstruction :: ~RadmonDetectorConstruction()
{
 Destruct();
 
 delete constructorsFactory;

 // TO BE DONE
 G4cout << "RadmonDetectorConstruction::~RadmonDetectorConstruction(): PLEASE CHECK" << G4endl;
}





G4VPhysicalVolume *                             RadmonDetectorConstruction :: Construct(void)
{
 // TO BE DONE
 G4cout << "RadmonDetectorConstruction::Construct(): NOT IMPLEMENTED YED" << G4endl;
 
 return 0;
}





void                                            RadmonDetectorConstruction :: OnLayoutChange(void)
{
 Destruct();
 G4RunManager::GetRunManager()->DefineWorldVolume(Construct(), true);
}




void                                            RadmonDetectorConstruction :: Destruct(void)
{
 // TO BE DONE
 G4cout << "RadmonDetectorConstruction::Destruct(): NOT IMPLEMENTED YED" << G4endl;
}
