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

#include "G4TSimulationGUI.h"


ClassImp(G4TSimulationGUI)

using namespace std;
using namespace ROOT;
using namespace TMath;

G4TSimulationGUI::G4TSimulationGUI(const TGWindow *p, UInt_t w, UInt_t h)
  :TGMainFrame(p, w, h)
{
  G4TSimHelper::LoadLibraries();
  Initialize();
}


void G4TSimulationGUI::Initialize()
{
  // composite frame
  TGCompositeFrame *fMainFrame795 = new TGCompositeFrame(this, 378, 372, kVerticalFrame);
  fMainFrame795->SetLayoutBroken(kTRUE);

  // composite frame
  TGCompositeFrame *fCompositeFrame796 = new TGCompositeFrame(fMainFrame795, 378, 372,
                                                              kVerticalFrame);
  fCompositeFrame796->SetLayoutBroken(kTRUE);

  // composite frame
  TGCompositeFrame *fCompositeFrame797 = new TGCompositeFrame(fCompositeFrame796, 379, 372,
                                                             kVerticalFrame);
  fCompositeFrame797->SetLayoutBroken(kTRUE);

  // composite frame
  TGCompositeFrame *fCompositeFrame798 = new TGCompositeFrame(fCompositeFrame797, 380, 184,
                                                             kVerticalFrame);
  fCompositeFrame798->SetLayoutBroken(kTRUE);

  // composite frame
  TGCompositeFrame *fCompositeFrame799 = new TGCompositeFrame(fCompositeFrame798, 378, 176,
                                                              kVerticalFrame);
  fCompositeFrame799->SetLayoutBroken(kTRUE);

  TGFont *ufont;         // will reflect user font changes
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

  TGGC   *uGC;           // will reflect user GC changes
  // graphics context changes
  GCValues_t vall800;
  vall800.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000",vall800.fForeground);
  gClient->GetColorByName("#c0c0c0",vall800.fBackground);
  vall800.fFillStyle = kFillSolid;
  vall800.fFont = ufont->GetFontHandle();
  vall800.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall800, kTRUE);
  TGLabel *fLabel800 = new TGLabel(fCompositeFrame799, "Publication:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel800->SetTextJustify(9);
  fLabel800->SetMargins(0,0,0,0);
  fLabel800->SetWrapLength(-1);
  fCompositeFrame799->AddFrame( fLabel800, new TGLayoutHints( kLHintsLeft | kLHintsTop,
                                                            2, 2, 2, 2 ) );
  fLabel800->MoveResize(8, 8, 80, 18);

  ULong_t ucolor;        // will reflect user color changes
  gClient->GetColorByName("#ffffff",ucolor);

  // combo box
  cPublication = new TGComboBox(fCompositeFrame799, -1, kHorizontalFrame | kSunkenFrame |
                                                        kDoubleBorder    | kOwnBackground);
 
  //Fill automatically (List of names of the existing publication files)
  vector<TString> publications = gCatalog->GetPublications();
  for(UInt_t i=0;i< publications.size(); ++i)
                                        cPublication->AddEntry(publications[i].Data() , i);
  cPublication->Resize(176,22);

  fCompositeFrame799->AddFrame( cPublication, new TGLayoutHints( kLHintsLeft | kLHintsTop,
                                                                 2, 2, 2, 2) );
  cPublication->MoveResize(96, 8, 272, 22);

  // pre-select
  if(publications.size() > 0) cPublication->Select(0);

  bRun = new TGTextButton(this/*fCompositeFrame799*/,"Run",100);
  bRun->SetTextJustify(36);
  bRun->SetMargins(0,0,0,0);
  bRun->SetWrapLength(-1);
  bRun->Resize(99,24);
  bRun->SetState(kButtonUp);
  fCompositeFrame799->AddFrame(bRun, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2,2,2,2));
  bRun->MoveResize(264, 144, 99, 24);

  bClose = new TGTextButton(this/*fCompositeFrame799*/,"Close ROOT",200);
  bClose->SetTextJustify(36);
  bClose->SetMargins(0,0,0,0);
  bClose->SetWrapLength(-1);
  bClose->Resize(99,24);
  bClose->SetState(kButtonUp);
  fCompositeFrame799->AddFrame(bClose, new TGLayoutHints(kLHintsLeft |
                                                         kLHintsTop, 2, 2, 2, 2) );
  bClose->MoveResize(160, 144, 100, 24);

  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

  // graphics context changes "Reaction Cross Section:"
  GCValues_t vall820;
  vall820.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000",vall820.fForeground);
  gClient->GetColorByName("#c0c0c0",vall820.fBackground);
  vall820.fFillStyle = kFillSolid;
  vall820.fFont = ufont->GetFontHandle();
  vall820.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall820, kTRUE);
  TGLabel *fLabel820 = new TGLabel(fCompositeFrame799, "Reaction Cross Section:",
                                   uGC->GetGC(), ufont->GetFontStruct());
  fLabel820->SetTextJustify(9);
  fLabel820->SetMargins(0,0,0,0);
  fLabel820->SetWrapLength(-1);
  fCompositeFrame799->AddFrame( fLabel820, new TGLayoutHints( kLHintsLeft | kLHintsTop,
                                                              2, 2, 2, 2 ) );
  fLabel820->MoveResize(8, 32, 144, 18);

  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

  // graphics context changes "Reaction Cross Section:"
  GCValues_t vall840;
  vall840.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000",vall840.fForeground);
  gClient->GetColorByName("#c0c0c0",vall840.fBackground);
  vall840.fFillStyle = kFillSolid;
  vall840.fFont = ufont->GetFontHandle();
  vall840.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall840, kTRUE);
  TGLabel *fLabel840 = new TGLabel(fCompositeFrame799, "Units:",
                                   uGC->GetGC(), ufont->GetFontStruct());
  fLabel840->SetTextJustify(9);
  fLabel840->SetMargins(0,0,0,0);
  fLabel840->SetWrapLength(-1);
  fCompositeFrame799->AddFrame( fLabel840, new TGLayoutHints( kLHintsLeft | kLHintsTop,
                                                             2, 2, 2, 2 ) );
  fLabel840->MoveResize(256, 32, 32, 18);

  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

  // graphics context changes for "MC Model"
  GCValues_t vall821;
  vall821.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000",vall821.fForeground);
  gClient->GetColorByName("#c0c0c0",vall821.fBackground);
  vall821.fFillStyle = kFillSolid;
  vall821.fFont = ufont->GetFontHandle();
  vall821.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall821, kTRUE);
  TGLabel *fLabel821 = new TGLabel( fCompositeFrame799, "MC Model:", uGC->GetGC(),
                                   ufont->GetFontStruct() );
  fLabel821->SetTextJustify(9);
  fLabel821->SetMargins(0,0,0,0);
  fLabel821->SetWrapLength(-1);
  fCompositeFrame799->AddFrame( fLabel821, new TGLayoutHints( kLHintsLeft | kLHintsTop,
                                                             2, 2, 2, 2 ) );
  fLabel821->MoveResize(8, 56, 48, 18);

  gClient->GetColorByName("#ffffff",ucolor);

  // combo box for the "Model Name"
  cModel = new TGComboBox(fCompositeFrame799, -1, kHorizontalFrame | kSunkenFrame |
                                                 kDoubleBorder | kOwnBackground);
  // fill with the models
  cModel->AddEntry("CHIPS",  0);
  cModel->AddEntry("BERTINI",1);
  cModel->AddEntry("PRECO",  2);
  cModel->AddEntry("LEP",    3);
  cModel->AddEntry("HEP",    4);
  cModel->AddEntry("BINARY", 5);
  cModel->AddEntry("QGSP",   6);
  cModel->AddEntry("QGSC",   7);
  cModel->AddEntry("FTFP",   8);
  cModel->AddEntry("FTFC",   9);
  cModel->Resize(176,22);
  fCompositeFrame799->AddFrame(cModel, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  cModel->MoveResize(64, 56, 120, 22);
  cModel->Select(0);

  // graphics context changes for "Postfix"
  GCValues_t vall822;
  vall822.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000",vall822.fForeground);
  gClient->GetColorByName("#c0c0c0",vall822.fBackground);
  vall822.fFillStyle = kFillSolid;
  vall822.fFont = ufont->GetFontHandle();
  vall822.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall822, kTRUE);
  TGLabel *fLabel822 = new TGLabel( fCompositeFrame799, "Postfix:", uGC->GetGC(),
                                   ufont->GetFontStruct() );
  fLabel822->SetTextJustify(9);
  fLabel822->SetMargins(0,0,0,0);
  fLabel822->SetWrapLength(-1);
  fCompositeFrame799->AddFrame( fLabel822, new TGLayoutHints( kLHintsLeft | kLHintsTop,
                                                             2, 2, 2, 2 ) );
  fLabel822->MoveResize(192, 56, 48, 18);

  gClient->GetColorByName("#ffffff",ucolor);

  // text box for the Postfix
  tPostfix = new TGTextEntry( fCompositeFrame799, new TGTextBuffer(15), -1, uGC->GetGC(),
                              ufont->GetFontStruct(), kSunkenFrame | kDoubleBorder |
                                                      kOwnBackground);
  tPostfix->SetMaxLength(32);
  tPostfix->SetAlignment(kTextLeft);
  tPostfix->SetText("g493p5");
  tPostfix->Resize(120,tPostfix->GetDefaultHeight());
  fCompositeFrame799->AddFrame(tPostfix, new TGLayoutHints(kLHintsLeft |
                                                           kLHintsTop, 2, 2, 2, 2) );
  tPostfix->MoveResize(248, 56, 120, 22);

  // Insert Reaction Cross-Section value
  nCrossSection = new TGNumberEntry( fCompositeFrame799, 0., 11, -1,
                                     (TGNumberFormat::EStyle) 5 );
  fCompositeFrame799->AddFrame( nCrossSection, new TGLayoutHints( kLHintsLeft | kLHintsTop,
                                                                  2, 2, 2, 2 ) );
  nCrossSection->MoveResize(160, 32, 80, 22);
  nCrossSection->SetIntNumber(450); // default CS

  // combo box for the Cross-section Units
  nUnits = new TGComboBox( fCompositeFrame799, -1, kHorizontalFrame | kSunkenFrame |
                                                   kDoubleBorder | kOwnBackground);
  nUnits->AddEntry("mb ",    mBarn);
  nUnits->AddEntry("Barn ",  Barn);
  nUnits->AddEntry("mkb ",   mkBarn);
  nUnits->Resize(176,22);
  fCompositeFrame799->AddFrame(nUnits, new TGLayoutHints( kLHintsLeft |
                                                          kLHintsTop, 2, 2, 2, 2) );
  nUnits->MoveResize(304, 32, 64, 22);
  nUnits->Select(mBarn);

  // Insert Number of Runs value
  nRuns = new TGNumberEntry( fCompositeFrame799, (Double_t) 0, 11, -1,
                             (TGNumberFormat::EStyle) 5 );
  fCompositeFrame799->AddFrame(nRuns, new TGLayoutHints(kLHintsLeft |
                                                        kLHintsTop, 2, 2, 2, 2) );
  nRuns->MoveResize(192,80,96,22);
  nRuns->SetIntNumber(25); // default #ofRuns

  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

  // graphics context changes
  GCValues_t vall848;
  vall848.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000",vall848.fForeground);
  gClient->GetColorByName("#c0c0c0",vall848.fBackground);
  vall848.fFillStyle = kFillSolid;
  vall848.fFont = ufont->GetFontHandle();
  vall848.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall848, kTRUE);
  TGLabel *fLabel848 = new TGLabel(fCompositeFrame799, "#ofEents/2000:", uGC->GetGC(),
                                  ufont->GetFontStruct());
  fLabel848->SetTextJustify(33);
  fLabel848->SetMargins(0,0,0,0);
  fLabel848->SetWrapLength(-1);
  fCompositeFrame799->AddFrame(fLabel848, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  fLabel848->MoveResize(8,80,176,18);
  cExisting = new TGCheckButton(fCompositeFrame799," ");
  cExisting->SetTextJustify(36);
  cExisting->SetMargins(0,0,0,0);
  cExisting->SetWrapLength(-1);
  fCompositeFrame799->AddFrame(cExisting, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  cExisting->MoveResize(192,104,109,19);

  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

  // graphics context changes
  GCValues_t vall850;
  vall850.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000",vall850.fForeground);
  gClient->GetColorByName("#c0c0c0",vall850.fBackground);
  vall850.fFillStyle = kFillSolid;
  vall850.fFont = ufont->GetFontHandle();
  vall850.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall850, kTRUE);
  TGLabel *fLabel850 = new TGLabel(fCompositeFrame799, "Use existing data:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel850->SetTextJustify(33);
  fLabel850->SetMargins(0,0,0,0);
  fLabel850->SetWrapLength(-1);
  fCompositeFrame799->AddFrame(fLabel850, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  fLabel850->MoveResize(8,104,176,18);

  fCompositeFrame798->AddFrame(fCompositeFrame799, new TGLayoutHints(kLHintsExpandX |
                                                                     kLHintsExpandY) );
  fCompositeFrame799->MoveResize(0,0,378,176);

  fCompositeFrame797->AddFrame(fCompositeFrame798, new TGLayoutHints(kLHintsExpandX |
                                                                    kLHintsExpandY) );
  fCompositeFrame798->MoveResize(0,0,380,184);

  fCompositeFrame796->AddFrame(fCompositeFrame797, new TGLayoutHints(kLHintsExpandX |
                                                                    kLHintsExpandY) );
  fCompositeFrame797->MoveResize(0,0,379,372);

  fMainFrame795->AddFrame(fCompositeFrame796, new TGLayoutHints(kLHintsExpandX |
                                                                kLHintsExpandY) );
  fCompositeFrame796->MoveResize(0,0,378,372);

  AddFrame(fMainFrame795, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  fMainFrame795->MoveResize(0,0,378,372);

  SetMWMHints(kMWMDecorAll,
              kMWMFuncAll,
              kMWMInputModeless);
              this->MapSubwindows();

  MapWindow();
  Resize(378,181);
}

/* Destructor MainFrame */
G4TSimulationGUI::~G4TSimulationGUI()
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

void G4TSimulationGUI::CloseWindow()
{
  //Terminate application
  //gApplication->Terminate(0);
  this->DestroyWindow();
}


//______________________________________________________________________________
void G4TSimulationGUI::Run_Click()
{
  TString publication(cPublication->GetSelectedEntry()->GetTitle());
  TString model(cModel->GetSelectedEntry()->GetTitle());
  TString postfix(tPostfix->GetText());
  Double_t cs = static_cast<Double_t>(nCrossSection->GetIntNumber());
  switch(static_cast<SigmaEnum>(nUnits->GetSelected()))
  {
  case mBarn:
    break;
  case Barn:
    cs*=1000.;
    break;
  case mkBarn:
    cs/=1000.;
    break;
  }
  Int_t runs = nRuns->GetIntNumber();
  Bool_t existing = (cExisting->GetState() == 0 ? false: true);
  model.ToLower(); // switch the model name to low case !
  postfix.ToLower(); // switch the postfix name to low case !
  // Only one button, therefore it's run button
  cout<<"G4TSimulationGUI::Run_Click:Publication="<<publication<<", Model="<<model<<", CS="
      <<cs<<", Postfix="<<postfix<<", Runs="<<runs<<", Exist="<<existing<<endl;
  // run the simulation
  gSimulationTool->Run(publication, cs, model, postfix, runs, existing);
}

//______________________________________________________________________________
void G4TSimulationGUI::Close_Click()
{
  gApplication->Terminate(0);
  //this->DestroyWindow();
}

Bool_t G4TSimulationGUI::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
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

        default:
          break;
      }
    default:
      break;
  }
  return kTRUE;
}

