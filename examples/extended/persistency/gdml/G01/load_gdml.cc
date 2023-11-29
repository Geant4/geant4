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
/// \file persistency/gdml/G01/load_gdml.cc
/// \brief Main program of the persistency/gdml/G01 example
//
//
//
//
// --------------------------------------------------------------
//      GEANT 4 - load_gdml
//
// --------------------------------------------------------------

#include <vector>

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4TransportationManager.hh"

#include "G01PrimaryGeneratorAction.hh"
#include "G01DetectorConstruction.hh"
#include "G01ActionInitialization.hh"

#include "FTFP_BERT.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4GDMLParser.hh"

void print_aux(const G4GDMLAuxListType* auxInfoList, G4String prepend="|")
{
  for(std::vector<G4GDMLAuxStructType>::const_iterator
      iaux = auxInfoList->begin(); iaux != auxInfoList->end(); iaux++ )
    {
      G4String str=iaux->type;
      G4String val=iaux->value;
      G4String unit=iaux->unit;

      G4cout << prepend << str << " : " << val  << " " << unit << G4endl;

      if (iaux->auxList) print_aux(iaux->auxList, prepend + "|");
    }
  return;
}

// --------------------------------------------------------------

int main(int argc,char **argv)
{
   G4cout << G4endl;
   G4cout << "Usage: load_gdml <intput_gdml_file:mandatory>"
          << " <output_gdml_file:optional>" << G4endl;
   G4cout << G4endl;

   if (argc<2)
   {
      G4cout << "Error! Mandatory input file is not specified!" << G4endl;
      G4cout << G4endl;
      return -1;
   }

   G4GDMLParser parser;

// Uncomment the following if wish to avoid names stripping
// parser.SetStripFlag(false);

// Uncomment the following and set a string with proper absolute path and
// schema filename if wishing to use alternative schema for parsing validation
// parser.SetImportSchema("");

   parser.SetOverlapCheck(true);
   parser.Read(argv[1]);

   if (argc>4)
   {
      G4cout << "Error! Too many arguments!" << G4endl;
      G4cout << G4endl;
      return -1;
   }

   auto* runManager = G4RunManagerFactory::CreateRunManager();

   runManager->SetUserInitialization(new G01DetectorConstruction(
                                     parser.GetWorldVolume()));
   runManager->SetUserInitialization(new FTFP_BERT);
   runManager->SetUserInitialization(new G01ActionInitialization());

   runManager->Initialize();

   // Initialize visualization
   G4VisManager* visManager = new G4VisExecutive;
   visManager->Initialize();

   // Get the pointer to the User Interface manager
   G4UImanager* UImanager = G4UImanager::GetUIpointer();

   ///////////////////////////////////////////////////////////////////////
   //
   // Example how to retrieve Auxiliary Information
   //

   G4cout << std::endl;

   const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
   std::vector<G4LogicalVolume*>::const_iterator lvciter;
   for( lvciter = lvs->begin(); lvciter != lvs->end(); lvciter++ )
   {
     G4GDMLAuxListType auxInfo = parser.GetVolumeAuxiliaryInformation(*lvciter);

     if (auxInfo.size()>0)
       G4cout << "Auxiliary Information is found for Logical Volume :  "
              << (*lvciter)->GetName() << G4endl;

     print_aux(&auxInfo);
   }

   // now the 'global' auxiliary info
   G4cout << std::endl;
   G4cout << "Global auxiliary info:" << std::endl;
   G4cout << std::endl;

   print_aux(parser.GetAuxList());

   G4cout << std::endl;

   //
   // End of Auxiliary Information block
   //
   ////////////////////////////////////////////////////////////////////////


   runManager->BeamOn(0);

   // example of writing out

   if (argc>=3)
   {
/*
     G4GDMLAuxStructType mysubaux = {"mysubtype", "mysubvalue", "mysubunit", 0};
     G4GDMLAuxListType* myauxlist = new G4GDMLAuxListType();
     myauxlist->push_back(mysubaux);

     G4GDMLAuxStructType myaux = {"mytype", "myvalue", "myunit", myauxlist};
     parser.AddAuxiliary(myaux);


     // example of setting auxiliary info for world volume
     // (can be set for any volume)

     G4GDMLAuxStructType mylocalaux = {"sometype", "somevalue", "someunit", 0};

     parser.AddVolumeAuxiliary(mylocalaux,
       G4TransportationManager::GetTransportationManager()
       ->GetNavigatorForTracking()->GetWorldVolume()->GetLogicalVolume());
*/

     parser.SetRegionExport(true);
     //     parser.SetEnergyCutsExport(true);
     //     parser.SetOutputFileOverwrite(true);
     parser.Write(argv[2], G4TransportationManager::GetTransportationManager()
      ->GetNavigatorForTracking()->GetWorldVolume()->GetLogicalVolume());
   }


   if (argc==4)   // batch mode
   {
     G4String command = "/control/execute ";
     G4String fileName = argv[3];
     UImanager->ApplyCommand(command+fileName);
   }
   else           // interactive mode
   {
     G4UIExecutive* ui = new G4UIExecutive(argc, argv);
     UImanager->ApplyCommand("/control/execute vis.mac");
     ui->SessionStart();
     delete ui;
   }

   delete visManager;
   delete runManager;

   return 0;
}
