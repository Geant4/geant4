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
//
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TToolMenuGUI.h"


ClassImp(G4TToolMenuGUI)

using namespace std;
using namespace ROOT;
using namespace TMath;




//______________________________________________________________________________
G4TToolMenuGUI::G4TToolMenuGUI(const TGWindow *p, UInt_t w, UInt_t h)
  :TGMainFrame(p, w, h)
{

 G4TSimHelper::LoadLibraries();
 Initialize();

}

//______________________________________________________________________________
void G4TToolMenuGUI::Initialize()
{
    this->SetLayoutBroken(kTRUE);

    bRun1 = new TGTextButton(this,"Make Publication",1);
    bRun1->SetTextJustify(36);
    bRun1->SetMargins(0,0,0,0);
    bRun1->SetWrapLength(-1);
    bRun1->Resize(264,40);
    this->AddFrame(bRun1, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    bRun1->MoveResize(8,8,264,40);


    bRun2 = new TGTextButton(this,"Run Simulation",2);
    bRun2->SetTextJustify(36);
    bRun2->SetMargins(0,0,0,0);
    bRun2->SetWrapLength(-1);
    bRun2->Resize(264,40);
    this->AddFrame(bRun2, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    bRun2->MoveResize(8,56,264,40);


    bRun3 = new TGTextButton(this,"Run Analysis",3);
    bRun3->SetTextJustify(36);
    bRun3->SetMargins(0,0,0,0);
    bRun3->SetWrapLength(-1);
    bRun3->Resize(264,40);
    this->AddFrame(bRun3, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    bRun3->MoveResize(8,104,264,40);



    this->SetMWMHints(kMWMDecorAll,
                         kMWMFuncAll,
                         kMWMInputModeless);
    this->MapSubwindows();

    this->MapWindow();
    this->Resize(278,154);
}



//______________________________________________________________________________
G4TToolMenuGUI::~G4TToolMenuGUI()
{
  TGFrameElement *ptr;

  // delete all frames and layout hints
  if (fList) {
    TIter next(fList);
    while ((ptr = (TGFrameElement *) next())) {
      if (ptr->fLayout)
        delete ptr->fLayout;
      if (ptr->fFrame)
        delete ptr->fFrame;
    }
  }
}


//______________________________________________________________________________
void G4TToolMenuGUI::Run1_Click()
{
 G4TSimHelper::LoadLibraries();
 G4TPublicationGUI* Window = new G4TPublicationGUI(gClient->GetRoot(), 0, 0);
 Window->MapWindow();
 //this->CloseWindow();
}

//______________________________________________________________________________
void G4TToolMenuGUI::Run2_Click()
{
 G4TSimHelper::LoadLibraries();
 G4TSimulationGUI* Window = new G4TSimulationGUI(gClient->GetRoot(), 0, 0);
 Window->MapWindow();
 //this->CloseWindow();
}

//______________________________________________________________________________
void G4TToolMenuGUI::Run3_Click()
{
 G4TSimHelper::LoadLibraries();
 G4TAnalysisGUI* Window = new G4TAnalysisGUI(gClient->GetRoot(), 0, 0);
 Window->MapWindow();
 //this->CloseWindow();
}



//______________________________________________________________________________
void G4TToolMenuGUI::CloseWindow()
{
  //Terminate application
  gApplication->Terminate(0);
  //this->DestroyWindow();
}


//______________________________________________________________________________
Bool_t G4TToolMenuGUI::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
  switch (GET_MSG(msg))
  {
    case kC_COMMAND:
      switch (GET_SUBMSG(msg))
      {
        case kCM_BUTTON:
          // Only one button, therefore it's run button
         if     (parm1 == 1) Run1_Click();
         else if(parm1 == 2) Run2_Click();
         else if(parm1 == 3) Run3_Click();
         break;
       default: break;
     }
     default: break;
  }
  return kTRUE;
}

