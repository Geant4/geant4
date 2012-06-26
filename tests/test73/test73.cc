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

// $Id:$
// GEANT4 tag $Name:$
//
//
#include "BeamTestDetectorConstruction.hh"
// User defined event action
#include "BeamTestEventAction.hh"
#include "BeamTestPhysicsList.hh"
#include "BeamTestPrimaryGeneratorAction.hh"
#include "BeamTestStackingAction.hh"
// User defined run action
#include "BeamTestRunAction.hh"
#include "BeamTestFileReader.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"
#include "G4PhysListFactory.hh"

#include <cstring>
#include <sstream>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4Version.hh"
#if G4VERSION_NUMBER>=930
#include <G4UIExecutive.hh> //For 9.3 or above
#else
#include "G4UIterminal.hh"
#ifdef G4UI_USE_TCSH
#include "G4UItcsh.hh"
#endif //G4UI_USE_TCSH
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif //G4UI_USE_XM
#endif //G4VERSION_NUMBER


int main(int argc,char** argv) 
{
#ifndef G4ANALYSIS_USEROOT
  G4Exception(__FILE__,"MSCTest001",FatalException,"ROOT analysis interface not compiled, test produces no"\
	      " results, set G4ANALYSIS_USEROOT env variable to 1 and recompile.");
#endif
	//choose the Random engine
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

	// Run manager
	G4RunManager * runManager = new G4RunManager;
    
	// Mandatory initialization classes
	BeamTestDetectorConstruction* detector = new BeamTestDetectorConstruction(/*parameter*/);
	runManager->SetUserInitialization(detector);
    const char* plname = 0;
    if ( argc > 2 ) { //Physics List via command line
        G4cout<<argv[2]<<G4endl;
        plname = argv[2];
    }
    else { //Search physics list via env.variable
        plname = getenv("PHYSLIST");
        if ( !plname ) plname = "local";//No command line and no env.
    }
    if ( strcmp(plname,"local")==0 ) {
        runManager->SetUserInitialization(new BeamTestPhysicsList);
    }
    else {
      G4PhysListFactory factory;
      runManager->SetUserInitialization( factory.GetReferencePhysList(plname) );
    }
	// User action classes
	BeamTestRunAction* run_action = new BeamTestRunAction(/*parameter*/);  
	runManager->SetUserAction(run_action);
	
	BeamTestEventAction* event_action = new BeamTestEventAction(/*parameter,*/ run_action);
	runManager->SetUserAction(event_action);
	
	runManager->SetUserAction(new BeamTestPrimaryGeneratorAction(/*parameter*/));
    
    //Stacking Action
    runManager->SetUserAction(new BeamTestStackingAction());
	//G4UserRunAction* run_action = new BeamTestRunAction;
	//runManager->SetUserAction(run_action);
    // User defined event action
    
#ifdef G4VIS_USE
	// Visualization manager
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
#endif

	// Initialize G4 kernel
	runManager->Initialize();


    G4UImanager* UI = G4UImanager::GetUIpointer();  
    if (argc==1)   // Define UI session for interactive mode.
	{
#if G4VERSION_NUMBER>=930
        G4UIExecutive* session = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
        UI->ApplyCommand("/control/execute init_vis.mac");
#else
        UI->ApplyCommand("/control/execute init.mac");
#endif
#else //G4VERSION_NUMBER<930
        G4UIsession* session = 0;
#ifdef G4UI_USE_XM
        session = new G4UIXm(argc,argv);
#else
#ifdef G4UI_USE_TCSH
        session = new G4UIterminal(new G4UItcsh);
#else
        session = new G4UIterminal();
#endif //G4UI_USE_TCSH
#endif //G4UI_USE_XM
#endif //G4VERSION_NUMBER>=930
        if (session)   // Define UI session for interactive mode.
        {
            session->SessionStart();
            delete session;
        }
    }
    else           // Batch mode
    { 
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        //G4String fileName = /*parameter->filename*/;
        UI->ApplyCommand(command+fileName);
    }
    
	// Job termination
#ifdef G4VIS_USE
	delete visManager;
#endif

	delete runManager;

	/*
	   double mydouble = parameter->pTransverse;
	   std::stringstream s; 
	   s << mydouble;
	   std::string result;
	   result = s.str();
	   TString Name = "IP_"+result+".root";

	   TFile* newFile = new TFile(Name, "recreate");
	   TTree* newTree = new TTree("newTree", "newTree");
	   Double_t bx1(0), by1(0), bz1(0), pTi(0);//, pTe(0);
	   newTree->Branch("pTi", &pTi, "pTi/D");
	//newTree->Branch("pTe", &pTe, "pTe/D");
	newTree->Branch("bx1", &bx1, "bx1/D");
	//newTree->Branch("by1", &by1, "by1/D");
	//newTree->Branch("bz1", &bz1, "bz1/D");
	// Now do the analysis
	// Need to return the variables saved to Data.dat


	//FileReader filetxt("Data"+result+".dat"); 
	G4int i(0);
	//std::vector<G4double> bx(0.), by(0.), bz(0.), ipT(0.);
	G4double sum(0);
	if (filetxt.isValid())
	{
	while (filetxt.nextLine())
	{
	if (filetxt.inputFailed()) break; // breaks the loop as soon as the next line fails to provide input

	sum += filetxt.getFieldAsDouble(1); 
	//bx.push_back(bx1);
	//by1 = filetxt.getFieldAsDouble(2); 
	//by.push_back(by1);
	//bz1 = filetxt.getFieldAsDouble(3); 
	//bz.push_back(bz1);
	//pTe = filetxt.getFieldAsDouble(4);
	//ipT.push_back(pTe); // number 'n' in getFieldAsDouble(n) refers to the column of data 

	// Actually know this so not needed as will now have inverse pT as initial input as particle travels at some angle with pT.
	//ipT.push_back(parameter->pTransverse); // number 'n' in getFieldAsDouble(n) refers to the column of data 
	//newTree->Fill();
	//std::cout << "ipT= " << ipT[i] << "   bx= " << bx[i] << "   by= " << by[i] << "   bz= " << bz[i] << std::endl;
	//ntuple->Fill(x,y,z);
	++i;
	}
	}

	//G4double bx_av, by_av, bz_av;
	//pTi = (parameter->pTransverse); 
	//Average(ipT,pTe);
	//bx1 = sum/i;
	//Average(bx,bx1);
	//Average(by,by1);
	//Average(bz,bz1);
	//StandardDev(ipT,pTe);
	//	StandardDev(bx,bx1);
	//	StandardDev(by,by1);
	//	StandardDev(bz,bz1);

	//newTree->Fill();

	G4cout << G4endl;
	//G4cout << "Average b= " << bx1 << " over " << i << " events." << G4endl;
	//G4cout << "Average b= " << bx1 << " over " << bx.size() << " events." << G4endl;
	//G4cout << "Average b_x= " << bx1 << " over " << bx.size() << " events." << G4endl;
	//G4cout << "Average b_y= " << by1 << " over " << by.size() << " events." << G4endl;
	//G4cout << "Average b_z= " << bz1 << " over " << bz.size() << " events." << G4endl;
	//	G4cout << "SD b_x= " << bx1 << " over " << bx.size() << " events." << G4endl;
	//	G4cout << "SD b_y= " << by1 << " over " << by.size() << " events." << G4endl;
	//	G4cout << "SD b_z= " << bz1 << " over " << bz.size() << " events." << G4endl;
	G4cout << G4endl;
	//newFile->Write();
	//newFile->Close();
	*/	
	//	delete parameter;

	return 0;
}

