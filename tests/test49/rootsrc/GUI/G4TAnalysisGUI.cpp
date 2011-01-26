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

#include "G4TAnalysisGUI.h"


ClassImp(G4TAnalysisGUI)

using namespace std;
using namespace ROOT;
using namespace TMath;





G4TAnalysisGUI::G4TAnalysisGUI(const TGWindow *p, UInt_t w, UInt_t h)
  :TGMainFrame(p, w, h)
{

 G4TSimHelper::LoadLibraries();
 Initialize();

}


void G4TAnalysisGUI::Initialize()
{
    // composite frame
    TGCompositeFrame *fMainFrame631 = new TGCompositeFrame(this,378,181,kVerticalFrame);
    fMainFrame631->SetLayoutBroken(kTRUE);

    // composite frame
    TGCompositeFrame *fCompositeFrame632 = new TGCompositeFrame(fMainFrame631,378,372,kVerticalFrame);
    fCompositeFrame632->SetLayoutBroken(kTRUE);

    // composite frame
    TGCompositeFrame *fCompositeFrame633 = new TGCompositeFrame(fCompositeFrame632,378,372,kVerticalFrame);
    fCompositeFrame633->SetLayoutBroken(kTRUE);

    // composite frame
    TGCompositeFrame *fCompositeFrame634 = new TGCompositeFrame(fCompositeFrame633,379,372,kVerticalFrame);
    fCompositeFrame634->SetLayoutBroken(kTRUE);

    // composite frame
    TGCompositeFrame *fCompositeFrame635 = new TGCompositeFrame(fCompositeFrame634,380,184,kVerticalFrame);
    fCompositeFrame635->SetLayoutBroken(kTRUE);

    // composite frame
    TGCompositeFrame *fCompositeFrame636 = new TGCompositeFrame(fCompositeFrame635,378,176,kVerticalFrame);
    fCompositeFrame636->SetLayoutBroken(kTRUE);

    TGFont *ufont;         // will reflect user font changes
    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    TGGC   *uGC;           // will reflect user GC changes
    // graphics context changes
    GCValues_t vall637;
    vall637.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall637.fForeground);
    gClient->GetColorByName("#c0c0c0",vall637.fBackground);
    vall637.fFillStyle = kFillSolid;
    vall637.fFont = ufont->GetFontHandle();
    vall637.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall637, kTRUE);
    TGLabel *fLabel637 = new TGLabel(fCompositeFrame636,"Publication file:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel637->SetTextJustify(9);
    fLabel637->SetMargins(0,0,0,0);
    fLabel637->SetWrapLength(-1);
    fCompositeFrame636->AddFrame(fLabel637, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel637->MoveResize(8,8,176,18);

    ULong_t ucolor;        // will reflect user color changes
    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cPublication = new TGComboBox(fCompositeFrame636,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);

  //Fill automatically
  vector<TString> publications = gCatalog->GetPublications();
  for(UInt_t i=0;i< publications.size(); ++i)
  {
   cPublication->AddEntry(publications[i].Data() , i);
  }


  cPublication->Resize(176,22);

  fCompositeFrame636->AddFrame(cPublication, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

  cPublication->MoveResize(192,8,176,22);
  // pre-select
  if(publications.size() > 0)
    cPublication->Select(0);

  bRun = new TGTextButton(this/*fCompositeFrame799*/,"Run",100);
  bRun->SetTextJustify(36);
  bRun->SetMargins(0,0,0,0);
  bRun->SetWrapLength(-1);
  bRun->Resize(99,24);
  bRun->SetState(kButtonUp);
  fCompositeFrame636->AddFrame(bRun, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  bRun->MoveResize(264,144,99,24);

  bClose = new TGTextButton(this/*fCompositeFrame799*/,"Close ROOT",200);
  bClose->SetTextJustify(36);
  bClose->SetMargins(0,0,0,0);
  bClose->SetWrapLength(-1);
  bClose->Resize(99,24);
  bClose->SetState(kButtonUp);
  fCompositeFrame636->AddFrame(bClose, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  bClose->MoveResize(160,144,100,24);

  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

  // graphics context changes
  GCValues_t vall657;
  vall657.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
  gClient->GetColorByName("#000000",vall657.fForeground);
  gClient->GetColorByName("#c0c0c0",vall657.fBackground);
  vall657.fFillStyle = kFillSolid;
  vall657.fFont = ufont->GetFontHandle();
  vall657.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall657, kTRUE);
  TGLabel *fLabel657 = new TGLabel(fCompositeFrame636,"Particle number (0 for all):",uGC->GetGC(),ufont->GetFontStruct());
  fLabel657->SetTextJustify(9);
  fLabel657->SetMargins(0,0,0,0);
  fLabel657->SetWrapLength(-1);
  fCompositeFrame636->AddFrame(fLabel657, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  fLabel657->MoveResize(8,32,176,18);
  nParticleIdx = new TGNumberEntry(fCompositeFrame636, (Double_t) 0,11,-1,(TGNumberFormat::EStyle) 5);
  fCompositeFrame636->AddFrame(nParticleIdx, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  nParticleIdx->MoveResize(192,32,96,22);

  fCompositeFrame635->AddFrame(fCompositeFrame636, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  fCompositeFrame636->MoveResize(0,0,378,176);

  fCompositeFrame634->AddFrame(fCompositeFrame635, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  fCompositeFrame635->MoveResize(0,0,380,184);

  fCompositeFrame633->AddFrame(fCompositeFrame634, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  fCompositeFrame634->MoveResize(0,0,379,372);

  fCompositeFrame632->AddFrame(fCompositeFrame633, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  fCompositeFrame633->MoveResize(0,0,378,372);

  fMainFrame631->AddFrame(fCompositeFrame632, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  fCompositeFrame632->MoveResize(0,0,378,372);

  this->AddFrame(fMainFrame631, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  fMainFrame631->MoveResize(0,0,378,181);

  this->SetMWMHints(kMWMDecorAll,
                       kMWMFuncAll,
                       kMWMInputModeless);
  this->MapSubwindows();

  this->MapWindow();
  this->Resize(379,182);
}



/* Destructor MainFrame */
G4TAnalysisGUI::~G4TAnalysisGUI()
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


void G4TAnalysisGUI::CloseWindow()
{
  //Terminate application
  //gApplication->Terminate(0);
  this->DestroyWindow();
}


//______________________________________________________________________________
void G4TAnalysisGUI::Run_Click()
{
  TString publication(cPublication->GetSelectedEntry()->GetTitle());
  Int_t idx = nParticleIdx->GetIntNumber();
  gAnalysisTool->Run(publication, idx);
}

//______________________________________________________________________________
void G4TAnalysisGUI::Close_Click()
{
  gApplication->Terminate(0);
  //this->DestroyWindow();
}



Bool_t G4TAnalysisGUI::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
  switch (GET_MSG(msg))
  {
    case kC_COMMAND:
      switch (GET_SUBMSG(msg))
      {
        case kCM_BUTTON:
          if     (parm1 == 100) Run_Click();
          else if(parm1 == 200) Close_Click();
          break;
        default: break;
      }
    default: break;
  }
  return kTRUE;
}

