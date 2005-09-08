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
// $Id: remsim.cc,v 1.11 2005-09-08 06:57:56 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "RemSimDetectorConstruction.hh"
#include "RemSimPhysicsList.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#include "RemSimEventAction.hh"
#include "RemSimRunAction.hh"
#include "RemSimSteppingAction.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#ifdef G4ANALYSIS_USE
#include "RemSimAnalysisManager.hh"
#endif 

bool messaggi(false);

int main(int argc,char** argv)
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  //  HepRandom :: setTheSeed(0);
  // Set mandatory initialization classes
  
  // Geometry
  RemSimDetectorConstruction* detector = new RemSimDetectorConstruction();
  runManager -> SetUserInitialization(detector);
  
  // Physics
  runManager->SetUserInitialization(new RemSimPhysicsList);

  // Set mandatory user action class

  // Primary particles
  RemSimPrimaryGeneratorAction* primary = new RemSimPrimaryGeneratorAction();
  runManager -> SetUserAction(primary);
 
  // Set optional user action class
  RemSimEventAction* event = new RemSimEventAction();
  runManager -> SetUserAction(event);

  RemSimRunAction* run = new RemSimRunAction();
  runManager -> SetUserAction(run);

  runManager -> SetUserAction(new RemSimSteppingAction(primary));

#ifdef G4VIS_USE
   // Visualisation
   G4VisManager* visManager = new G4VisExecutive;
   visManager -> Initialize();
#endif

#ifdef G4ANALYSIS_USE
 RemSimAnalysisManager* analysis = RemSimAnalysisManager::getInstance();
 analysis -> SetFormat("hbook");
#endif
 
  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
   
  if(argc == 1)
    // Define (G)UI terminal for interactive mode  
    { 
      // G4UIterminal is a (dumb) terminal.
      G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    

      UI -> ApplyCommand("/control/execute vis.mac");    
      session -> SessionStart();
      delete session;
    }
  else
    // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI -> ApplyCommand(command+fileName);
    }

#ifdef G4ANALYSIS_USE
 analysis -> finish();
#endif

#ifdef G4VIS_USE
 delete visManager;
#endif

  // job termination
  delete runManager;
  messaggi=1;
  return 0;
}

// // RICCARDO
// #include <stdlib.h>

// #define MARGINI 10
// #define VALORE 0x3F
// #define VALORE2 0xFA

// void Dump(void *ptr, size_t size)
// {
//   char *ptrC=(char *)ptr;

//   for (size_t i=0; i<size; i++)
//     {
//       if (i%32==0)
// 	printf("%04x   ", i);

//       printf("%02lx ", (((long)ptrC[i])&0xFF));

//       if (i%32==31)
// 	printf("\n");
//     }

//   printf("\n");
// }

// void *operator new(size_t size) throw(std::bad_alloc)
// {
//  size_t sizeFull;
//  sizeFull=size+MARGINI*2+sizeof(size_t);

//  void * ptr;

//  ptr=malloc(sizeFull);

//  size_t * ptrL=(size_t *)ptr;
//  ptrL[0]=size;
//  char * ptrC=(char *)(ptrL+1);

//  int i(MARGINI);

//  while (i>0)
//    {
//      i--;
//      ptrC[i]=VALORE;
//      ptrC[i+size+MARGINI]=VALORE;
//    }

//  i=size;
//  while (i>0)
//    {
//      i--;
//      ptrC[i+MARGINI]=VALORE2;
//    }

//  // -4- -100- ----------------------  -100-
//  //           |

//  if (messaggi)
//    {
//   printf("MALLOC: %p %p %p %ld\n", ptrL, ptrC, ptrC+MARGINI, size);

//   Dump(ptr, sizeFull);
//    }
//  return (void*)(ptrC+MARGINI);
// }

// void operator delete(void *data) throw()
// {
//  char *ptrC=(char *)data;
//  ptrC-=MARGINI;

//  size_t size;
//  size_t * ptrL=(size_t *)ptrC;
//  ptrL--;

//  size=ptrL[0];

//  int i(MARGINI);

//  if (messaggi)
//    {
//   printf("FREE: %p %p %p %ld\n", ptrL, ptrC, data, size);

//   Dump((void *)ptrL, size+MARGINI*2+sizeof(size_t));
//    }
//  while (i>0)
//    {
//      i--;
//      if (ptrC[i]!=VALORE)
//        abort();
//      if (ptrC[i+size+MARGINI]!=VALORE)
//        abort();
//    }


//  free((void *) ptrL);
// }
 
