
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
// $Id: test50.cc,v 1.15 2003-02-10 11:21:04 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include <iostream.h>
#include <iomanip.h>
#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {
       
HepRandom::setTheEngine(new RanecuEngine);

 G4bool lowE=false;
 G4bool RangeOn=false;
 G4bool MaxStep=false;
 G4bool RY=false;
 G4bool end=true;
 G4bool SP=false;
 G4bool Foil= false;
 G4bool Adronic= false;
 G4String filename="test50.txt";
 G4cout<<argc<<":argc"<<G4endl;
 G4cout.setf(ios::scientific, ios::floatfield);
   

if (argc==1){G4cout <<"Input file is not specified!"<<G4endl;}
 ifstream* fin=new ifstream();
 string fname=argv[1];
 fin->open(fname.c_str());
 if (!fin->is_open())
{G4cout<<"InputFile<"<<fname<<">doesn't exist!Exit"<<G4endl;
 exit(1);}
 // Read input file//
 G4cout<<"Available commands are: "<<G4endl;
 G4cout<<"#processes (LowE/Standard)"<<G4endl; 
  G4cout<<"#range (on/off)"<<G4endl;  
 G4cout<<"#setMaxStep (on/off)"<<G4endl;
 G4cout<<"#StoppingPower(on/off)"<<G4endl; 
 G4cout<<"#RadiationYield(on/off)"<<G4endl; 
 G4cout<<"#Foil(on/off)"<<G4endl;
 G4cout<<"#Hadronic(on/off)"<<G4endl; 

 G4String line, line1, line2;
 
   do
     {
       (*fin)>>line;
       G4cout<<"Next line"<<line<<G4endl; 
       if (line=="#processes")
	 { 
          line1="";
	 (*fin)>>line1;
         G4cout<<"Next line"<<line1<<G4endl; 
	 if(line1=="LowE") {lowE=true;G4cout<<"arrivo ai LowEnergy"<<G4endl;}
	 else {lowE=false;G4cout<<"non arrivo ai LowEnergy"<<G4endl;}
         }
       if (line=="#range"){line1="";(*fin)>>line1;

       if(line1=="on"){RangeOn=true; G4cout<<RangeOn<<"range"<<G4endl;}}       

 if (line=="#setMaxStep"){line1="";(*fin)>>line1;

       if(line1=="on"){MaxStep=true; G4cout<<MaxStep<<"maxStep"<<G4endl;}}       
 if (line=="#StoppingPower"){line1="";(*fin)>>line1;

       if(line1=="on"){SP=true; G4cout<<SP<<"StoppingPower"<<G4endl;}}       
 if (line=="#RadiationYield"){line1="";(*fin)>>line1;

       if(line1=="on"){RY=true; G4cout<<RY<<"RadiationYield"<<G4endl;}} 
 
 if (line=="#Foil"){line1="";(*fin)>>line1;

       if(line1=="on"){Foil=true; G4cout<<Foil<<" :Foil configuration set"<<G4endl;}}    
 if (line=="#Hadronic"){line1="";(*fin)>>line1;

       if(line1=="on"){Adronic=true; G4cout<<Adronic<<" :Hadronic proton processes switched on"<<G4endl;}}

if (line=="end"){end=false;}



     }while(end);           
     
   if (RangeOn==true) filename="Range.txt";
   if (SP==true) filename="StoppingPower.txt"; 
   if (RY==true) filename="RadiationYeld.txt";
   if (Foil==true) filename="Transmission.txt";      
 G4int seed=time(NULL);
   HepRandom ::setTheSeed(seed);
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new Tst50SteppingVerbose);
  
  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  Tst50DetectorConstruction* Tst50detector = new Tst50DetectorConstruction( MaxStep);
  runManager->SetUserInitialization(Tst50detector);

  Tst50PhysicsList* fisica = new Tst50PhysicsList(lowE,RangeOn,SP,RY,Adronic);
  runManager->SetUserInitialization(fisica);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new Tst50VisManager;
  visManager->Initialize();
#endif
   
  // UserAction classes
  Tst50PrimaryGeneratorAction* p_Primary=new Tst50PrimaryGeneratorAction(); 
  runManager->SetUserAction(p_Primary);
  Tst50RunAction* p_run=new Tst50RunAction(Foil); 
  runManager->SetUserAction(p_run);  

  Tst50EventAction *pEventAction=new Tst50EventAction(p_Primary,RY,filename, Foil);
 
   runManager->SetUserAction(pEventAction );
     
 Tst50SteppingAction* steppingaction =new Tst50SteppingAction(pEventAction,p_Primary,p_run, Tst50detector, filename,SP,RangeOn,RY,Foil);
 runManager->SetUserAction(steppingaction);

  //Initialize G4 kernel
  runManager->Initialize();
    
  if(RangeOn)
    {
G4std::ofstream ofs;

	ofs.open(filename);
		{

		  ofs<<"range(g/cm2)(y)"<<'\t'<<"e- energy (MeV)(x)"<<'\t'<<G4endl;}
	       
       ofs.close();                     
		
    }
  if(SP)
    {
G4std::ofstream ofs;

	ofs.open(filename);
		{

		  ofs<<"StoppingPower(MeV*cm2/g)(y)"<<'\t'<<"e- energy (MeV)(x)"<<'\t'<<G4endl;}
	       
       ofs.close();                     
		
    }
  if(RY)
    {
G4std::ofstream ofs;

	ofs.open(filename);
		{

		  ofs<<"RadiationYield"<<'\t'<<"e- energy (MeV)(x)"<<'\t'<<G4endl;}


	       
       ofs.close();                     
		
    }

 if(Foil)
    {
G4std::ofstream ofs;

	ofs.open(filename);
		{

		  ofs<<"Global information about primary particles (e+ or e- or proton or gamma)"<<G4endl;}
	       
       ofs.close();                     
		
    }
   
  //get the pointer to the User Interface manager 
  G4UImanager * UI = G4UImanager::GetUIpointer();  
 UI->ApplyCommand("/run/verbose 0");
  UI->ApplyCommand("/event/verbose 0");
  UI->ApplyCommand("/tracking/verbose 0");
  
  if(argc==2)
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
    G4cout<<fileName<<G4endl;
    UI->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

