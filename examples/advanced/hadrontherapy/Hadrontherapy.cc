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

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "HadrontherapyEventAction.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyPhysicsList.hh"
#include "HadrontherapyPhantomSD.hh"
#include "HadrontherapyPrimaryGeneratorAction.hh"
#include "G4SDManager.hh"
#include "HadrontherapyRunAction.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "HadrontherapySteppingAction.hh"
#include "globals.hh"
#ifdef  G4ANALYSIS_USE
#include "HadrontherapyAnalysisManager.hh"
#endif
// ----------------------------------------------------------------
int main(int argc ,char ** argv)
{
  //Output matrix 
  G4int numberVoxelX = 80;
  G4int numberVoxelY = 80;
  G4int numberVoxelZ = 80;
 

  G4double* matrix = new G4double[numberVoxelX*numberVoxelY*numberVoxelZ];

  // Initialization of the matrix elemts to zero
  for(G4int i = 0; i < numberVoxelX; i++)
  {
      for(G4int j = 0; j < numberVoxelY; j++)
      {
	  for(G4int k = 0; k < numberVoxelZ; k++)

	matrix[(i*numberVoxelY+j)*numberVoxelZ+k] = 0.;
    }
  }

#ifdef G4ANALYSIS_USE
  HadrontherapyAnalysisManager* analysis = 
                          HadrontherapyAnalysisManager::getInstance();
  analysis -> book();
#endif
  

 //  G4double matrix[80][80][80]; // dimensions of the output matrix

  G4RunManager* pRunManager = new G4RunManager;

  // Initialize the geometry
  pRunManager -> SetUserInitialization(new HadrontherapyDetectorConstruction());
  
  //Initialize the physics 
  pRunManager -> SetUserInitialization(new HadrontherapyPhysicsList());
  
  // Initialize the primary particles  
  pRunManager -> SetUserAction(new HadrontherapyPrimaryGeneratorAction());

  // Optional UserActions: run, event, stepping
  pRunManager->SetUserAction(new HadrontherapyRunAction());
  HadrontherapyEventAction *pEventAction = new HadrontherapyEventAction( matrix, numberVoxelX,
                                                                         numberVoxelY, numberVoxelZ);
  pRunManager->SetUserAction(pEventAction );


  HadrontherapySteppingAction* steppingaction = new HadrontherapySteppingAction(); 
  pRunManager -> SetUserAction(steppingaction);    


#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  
  G4UIsession* session = 0;
  if (argc == 1)   // Define UI session for interactive mode.
    {
      session = new G4UIterminal();
    } 

  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  if (session)   // Define UI session for interactive mode.
    { 
      G4cout<<" UI session starts ..."<< G4endl;
      UI->ApplyCommand("/control/execute VisualisationMacro.mac");    
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }  


 if(matrix)
	{
	 std::ofstream ofs;

	// Output voxel data to text file
	// Format = i  <tab> j  <tab> k <tab> edep [MeV] <eol>
	ofs.open("EnergyDeposit.out");
		{

		  ofs<<" i "<<'\t'<<" j "<<'\t'<<"k"<<'\t' 
                     <<"released energy(Mev)"
		     <<G4endl;
		 
                  G4int k;
                  G4int j;
                  G4int i;             
                   
		    for(G4int l = 0; l < numberVoxelZ; l++) 
                      {
                        k = l;
                        
                       for(G4int m = 0; m < numberVoxelY; m++) 
                       { 
			 j = m * numberVoxelZ + k; 
                         
                        for(G4int n = 0; n <  numberVoxelX; n++)
			  {
			    i =  n* numberVoxelZ * numberVoxelY + j;
			     if(matrix[i]!=0)
			      {
			        ofs << n <<'\t'<< m <<'\t'<<
				k<<'\t'<<matrix[i]<<G4endl;
				// ofs<< i <<'\t'<<j<<'\t'<<
				//k<<'\t'<<matrix[i]<<G4endl;
#ifdef G4ANALYSIS_USE 
				//HadrontherapyAnalysisManager* analysis = 
				//HadrontherapyAnalysisManager::getInstance();
				analysis -> Energy_Dep(n, m, k, matrix[i]);
                                analysis -> BraggPeak(n, matrix[i]);
#endif
                             
			     }
			  }   
		       }
		      }
	       
		ofs.close();
		}
	}

 delete[] matrix;   

#ifdef G4ANALYSIS_USE
 //HadrontherapyAnalysisManager* analysis = 
 //                        HadrontherapyAnalysisManager::getInstance();
  analysis -> finish();
#endif
  
  // Job termination
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete pRunManager;

  return 0;
}
