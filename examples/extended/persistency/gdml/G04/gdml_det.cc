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
/// \file persistency/gdml/G04/gdml_det.cc
/// \brief Main program of the persistency/gdml/G04 example
//
//
// $Id: gdml_det.cc 82285 2014-06-13 14:56:02Z gcosmo $
//
//
// --------------------------------------------------------------
//      GEANT 4 - gdml_det
//
// --------------------------------------------------------------

#include <vector>

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"

#include "G04PrimaryGeneratorAction.hh"
#include "G04DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "G04SensitiveDetector.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4GDMLParser.hh"

int main(int argc,char **argv)
{
   G4cout << G4endl;
   G4cout << "Usage: gdml_det <intput_gdml_file:mandatory>"
          << G4endl;
   G4cout << G4endl;

   if (argc<2)
   {
      G4cout << "Error! Mandatory input file is not specified!" << G4endl;
      G4cout << G4endl;
      return -1;
   }

   G4GDMLParser parser;
   parser.Read(argv[1]);

   
   G4RunManager* runManager = new G4RunManager;

   runManager->SetUserInitialization(new G04DetectorConstruction(
                                     parser.GetWorldVolume()));
   runManager->SetUserInitialization(new FTFP_BERT);
   runManager->SetUserAction(new G04PrimaryGeneratorAction);

   runManager->Initialize();

   G4UImanager* UImanager = G4UImanager::GetUIpointer();

   //------------------------------------------------ 
   // Sensitive detectors
   //------------------------------------------------ 
   
   G4SDManager* SDman = G4SDManager::GetSDMpointer();
   
   G4String trackerChamberSDname = "Tracker";
   G04SensitiveDetector* aTrackerSD = 
     new G04SensitiveDetector(trackerChamberSDname);
   SDman->AddNewDetector( aTrackerSD );
 
   ///////////////////////////////////////////////////////////////////////
   //
   // Example how to retrieve Auxiliary Information for sensitive detector
   //
   const G4GDMLAuxMapType* auxmap = parser.GetAuxMap();
   G4cout << "Found " << auxmap->size()
             << " volume(s) with auxiliary information."
             << G4endl << G4endl;
   for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin();
       iter!=auxmap->end(); iter++) 
   {
     G4cout << "Volume " << ((*iter).first)->GetName()
            << " has the following list of auxiliary information: "
            << G4endl << G4endl;
     for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();
          vit!=(*iter).second.end(); vit++)
     {
       G4cout << "--> Type: " << (*vit).type
                 << " Value: " << (*vit).value << G4endl;
     }
   }
   G4cout << G4endl;

   // The same as above, but now we are looking for
   // sensitive detectors setting them for the volumes

   for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin();
       iter!=auxmap->end(); iter++) 
   {
     G4cout << "Volume " << ((*iter).first)->GetName()
            << " has the following list of auxiliary information: "
            << G4endl << G4endl;
     for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();
          vit!=(*iter).second.end();vit++)
     {
       if ((*vit).type=="SensDet")
       {
         G4cout << "Attaching sensitive detector " << (*vit).value
                << " to volume " << ((*iter).first)->GetName()
                <<  G4endl << G4endl;

         G4VSensitiveDetector* mydet = 
           SDman->FindSensitiveDetector((*vit).value);
         if(mydet) 
         {
           G4LogicalVolume* myvol = (*iter).first;
           myvol->SetSensitiveDetector(mydet);
         }
         else
         {
           G4cout << (*vit).value << " detector not found" << G4endl;
         }
       }
     }
   }
   //
   // End of Auxiliary Information block
   //
   ////////////////////////////////////////////////////////////////////////

#ifdef G4VIS_USE
     G4VisManager* visManager = new G4VisExecutive;
     visManager->Initialize();
#endif

   if (argc==3)   // batch mode  
   {
     G4String command = "/control/execute ";
     G4String fileName = argv[2];
     UImanager->ApplyCommand(command+fileName);
   }
   else           // interactive mode
   {
#ifdef G4UI_USE
     G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
     UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
     ui->SessionStart();
     delete ui;
#endif
   }

#ifdef G4VIS_USE
   delete visManager;
#endif
   delete runManager;

   return 0;
}
