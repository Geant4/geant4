
//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: test50.cc,v 1.19 2003-03-12 17:21:22 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//#include <iostream.h>
#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iostream"
#include "g4std/iomanip" 
#include "Tst50DetectorConstruction.hh"
#include "Tst50PhysicsList.hh"
#include "Tst50PrimaryGeneratorAction.hh"
#include "Tst50RunAction.hh"
#include "Tst50EventAction.hh"
#include "Tst50SteppingAction.hh"
#include "Tst50SteppingVerbose.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
#include "Tst50VisManager.hh"
#endif


int main(int argc,char** argv) {
       
  HepRandom::setTheEngine(new RanecuEngine);

  G4bool lowE = false;
  G4bool rangeOn = false;
  G4bool maxStep = false;
  G4bool radiationYield = false;
  G4bool end = true;
  G4bool stoppingPower = false;
  G4bool foil = false;
  G4bool hadronic = false;
  G4bool penelope = false;
  G4bool back = false;
  G4String filename = "test50.txt";
  G4cout << argc << ":argc" << G4endl;
  G4cout.setf(G4std::ios::scientific, G4std::ios::floatfield);
   

  if (argc == 1) { G4cout << "Input file is not specified!" << G4endl; }
  G4std::ifstream* fin=new G4std::ifstream();
  G4String fname=argv[1];
  fin->open(fname.c_str());
  if ( !fin->is_open() )
    { 
      G4cout << "InputFile<" << fname<<">doesn't exist!Exit" << G4endl;
      exit(1);
    }
  // Read input file//
  G4cout << "Available commands are: " << G4endl;
  G4cout << "#processes (LowE/Standard/Penelope)" << G4endl; 
  G4cout << "#range (on/off)" << G4endl;  
  G4cout << "#setMaxStep (on/off)" << G4endl;
  G4cout << "#StoppingPower (on/off)" << G4endl; 
  G4cout << "#RadiationYield (on/off)" << G4endl; 
  G4cout << "#Foil (on/off)" << G4endl;
  G4cout << "#Back (on/off)" << G4endl;
  G4cout << "#Hadronic (on/off)" << G4endl; 

  G4String line, line1, line2;
 
  do
    {
      (*fin) >> line;
      G4cout << "Next line" << line << G4endl; 
 
     if (line == "#processes")
	{ 
          line1="";
	  (*fin) >> line1;
	  G4cout << "Next line" << line1 << G4endl; 
	  if (line1 == "LowE") 
	    { 
	      lowE=true;
	      G4cout << "arrivo ai LowEnergy" << G4endl;
	    }
	  else if (line1 == "Penelope")
	    {
	      penelope = true; 
	      G4cout << penelope <<" :gamma Penelope processes switched on" << G4endl;
	    }
	}
 
     if (line == "#range") 
	{ 
	  line1 = ""; 
	  (*fin) >> line1;
	  if (line1 == "on") 
	    { 
	      rangeOn = true; 
	      G4cout << rangeOn << "range" << G4endl;
	    } 
	}       

      if (line == "#setMaxStep") { line1 =" "; (*fin) >> line1;
      if (line1 == "on") { maxStep = true; G4cout << maxStep << "maxStep" << G4endl;} }       

      if (line == "#StoppingPower") { line1 = ""; (*fin) >> line1;
      if (line1 == "on") { stoppingPower = true; G4cout << stoppingPower << "StoppingPower" << G4endl;} }       
   
      if (line == "#RadiationYield") { line1=""; (*fin) >> line1;
      if (line1 == "on") { radiationYield = true; G4cout << radiationYield << "RadiationYield" << G4endl;} } 
 
      if (line == "#Foil") { line1=""; (*fin) >> line1;
      if (line1 == "on") { foil = true; G4cout << foil << " :foil configuration set" << G4endl;} }   
 
      if (line == "#Back") { line1=""; (*fin) >> line1;
      if (line1 == "on") { back = true; G4cout << back << " :backscattering experiment set" << G4endl;} }   

      if (line == "#Hadronic") { line1=""; (*fin) >> line1;
      if (line1 == "on") { hadronic = true; G4cout << hadronic << " :Hadronic proton processes switched on" << G4endl;}}

      if (line == "end") { end = false;}

    } while(end);           
     
  if (rangeOn == true)        filename = "Range.txt";
  if (stoppingPower == true)  filename = "StoppingPower.txt";
  if (radiationYield == true) filename = "RadiationYeld.txt";
  if (foil == true)           filename = "Transmission.txt";
  
  G4int seed=time(NULL);
  HepRandom ::setTheSeed(seed);
  // my Verbose output class
  G4VSteppingVerbose::SetInstance(new Tst50SteppingVerbose);
  
  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  Tst50DetectorConstruction* Tst50detector = new Tst50DetectorConstruction(maxStep);
  runManager->SetUserInitialization(Tst50detector);

  Tst50PhysicsList* fisica = new Tst50PhysicsList(lowE,rangeOn,stoppingPower,
						  radiationYield,hadronic,penelope,back);
  runManager->SetUserInitialization(fisica);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new Tst50VisManager;
  visManager->Initialize();
#endif
   
  // UserAction classes
  Tst50PrimaryGeneratorAction* p_Primary = new Tst50PrimaryGeneratorAction(); 
  runManager->SetUserAction(p_Primary);
  Tst50RunAction* p_run=new Tst50RunAction(foil); 
  runManager->SetUserAction(p_run);  

  Tst50EventAction *pEventAction = new Tst50EventAction(p_Primary,radiationYield,filename,foil);
 
  runManager->SetUserAction(pEventAction);
     
  Tst50SteppingAction* steppingaction = new Tst50SteppingAction(pEventAction, p_Primary, p_run,
								Tst50detector, filename, stoppingPower,
								rangeOn, radiationYield, foil);
  runManager->SetUserAction(steppingaction);

  //Initialize G4 kernel
  runManager->Initialize();
    
  if (rangeOn)
    {
      G4std::ofstream ofs;
      ofs.open(filename);
      {
	ofs << "range(g/cm2)(y)" 
	    << '\t'
	    << "e- energy (MeV)(x)" 
	    <<'\t' 
	    << G4endl;
      }	       
      ofs.close();                     		
    }
 
  if (stoppingPower)
    {
      G4std::ofstream ofs;
      ofs.open(filename);
      {
	ofs << "StoppingPower(MeV*cm2/g)(y)" 
	    << '\t' 
	    << "e- energy (MeV)(x)" 
	    << '\t' 
	    << G4endl;
      }       
      ofs.close();                     		
    }
 
  if (radiationYield)
    {
      G4std::ofstream ofs;
      ofs.open(filename);
      {
	ofs << "RadiationYield"
	    << '\t'
	    << "e- energy (MeV)(x)"
	    << '\t' 
	    << G4endl;
      }       
      ofs.close();                     	
    }

  if (foil)
    {
      G4std::ofstream ofs;
      ofs.open(filename);
      {
	ofs << "Global information about primary particles (e+ or e- or proton or gamma)" << G4endl;
      } 
      ofs.close();                     
		
    }
  
  //get the pointer to the User Interface manager 
  G4UImanager * UI = G4UImanager::GetUIpointer();  
  UI->ApplyCommand("/run/verbose 0");
  UI->ApplyCommand("/event/verbose 0");
  UI->ApplyCommand("/tracking/verbose 0");
  
  if (argc == 2)
    // Define (G)UI terminal for interactive mode  
    { 
      // G4UIterminal is a (dumb) terminal.
      G4UIsession * session = 0;

#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    

      UI->ApplyCommand("/control/execute ");
      session->SessionStart();
      delete session;
    }
  else
    // Batch mode
    {     
      G4String command =("/control/execute ");
      G4String fileName = argv[2];
      G4cout << fileName << G4endl;
      UI->ApplyCommand(command+fileName);
    }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
