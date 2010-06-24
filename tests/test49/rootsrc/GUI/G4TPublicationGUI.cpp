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

#include "G4TPublicationGUI.h"

ClassImp(G4TPublicationGUI)

using namespace std;
using namespace ROOT;
using namespace TMath;


//______________________________________________________________________________
G4TPublicationGUI::G4TPublicationGUI(const TGWindow *p, UInt_t w, UInt_t h)
  : TGMainFrame(p, w, h)
{
   G4TSimHelper::LoadLibraries();
   Initialize();
}

//______________________________________________________________________________
void G4TPublicationGUI::Initialize()
{
  // composite frame
  TGCompositeFrame *fMainFrame635 = new TGCompositeFrame(this, 549, 487, kVerticalFrame);
  fMainFrame635->SetLayoutBroken(kTRUE);
  // composite frame
  TGCompositeFrame *fCompositeFrame636 = new TGCompositeFrame(fMainFrame635, 551, 487,
                                                              kVerticalFrame);
  fCompositeFrame636->SetLayoutBroken(kTRUE);
  // composite frame
  TGCompositeFrame *fCompositeFrame637 = new TGCompositeFrame(fCompositeFrame636, 551, 489,
                                                              kVerticalFrame);
  fCompositeFrame637->SetLayoutBroken(kTRUE);
  TGFont* ufont;         // will reflect user font changes
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  TGGC*   uGC;           // will reflect user GC changes
  // graphics context changes for "Projectile PDG"
  GCValues_t vall638;
  vall638.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall638.fForeground);
  gClient->GetColorByName("#c0c0c0", vall638.fBackground);
  vall638.fFillStyle = kFillSolid;
  vall638.fFont = ufont->GetFontHandle();
  vall638.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall638, kTRUE);
  TGLabel* fLabel638 = new TGLabel(fCompositeFrame637, "Projectile PDG:",
                                   uGC->GetGC(), ufont->GetFontStruct());
  fLabel638->SetTextJustify(33);
  fLabel638->SetMargins(0, 0, 0, 0);
  fLabel638->SetWrapLength(-1);
  fCompositeFrame637->AddFrame(fLabel638, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  fLabel638->MoveResize(8, 8, 120, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes "A:":
  GCValues_t vall654;
  vall654.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall654.fForeground);
  gClient->GetColorByName("#c0c0c0", vall654.fBackground);
  vall654.fFillStyle = kFillSolid;
  vall654.fFont = ufont->GetFontHandle();
  vall654.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall654, kTRUE);
  TGLabel* fLabel654 = new TGLabel(fCompositeFrame637, "Comment:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel654->SetTextJustify(33);
  fLabel654->SetMargins(0, 0, 0, 0);
  fLabel654->SetWrapLength(-1);
  fCompositeFrame637->AddFrame(fLabel654, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  fLabel654->MoveResize(302, 8, 48, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes "Target(Z,A)":
  GCValues_t vall639;
  vall639.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall639.fForeground);
  gClient->GetColorByName("#c0c0c0", vall639.fBackground);
  vall639.fFillStyle = kFillSolid;
  vall639.fFont = ufont->GetFontHandle();
  vall639.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall639, kTRUE);
  TGLabel* fLabel639 = new TGLabel(fCompositeFrame637, "Target (Z, A):", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel639->SetTextJustify(33);
  fLabel639->SetMargins(0, 0, 0, 0);
  fLabel639->SetWrapLength(-1);
  fCompositeFrame637->AddFrame(fLabel639, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  fLabel639->MoveResize(8, 32, 120, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes "A:":
  GCValues_t vall651;
  vall651.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall651.fForeground);
  gClient->GetColorByName("#c0c0c0", vall651.fBackground);
  vall651.fFillStyle = kFillSolid;
  vall651.fFont = ufont->GetFontHandle();
  vall651.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall651, kTRUE);
  TGLabel* fLabel651 = new TGLabel(fCompositeFrame637, "A:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel651->SetTextJustify(33);
  fLabel651->SetMargins(0, 0, 0, 0);
  fLabel651->SetWrapLength(-1);
  fCompositeFrame637->AddFrame(fLabel651, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  fLabel651->MoveResize(302, 32, 48, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes "Argument Type"
  GCValues_t vall640;
  vall640.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall640.fForeground);
  gClient->GetColorByName("#c0c0c0", vall640.fBackground);
  vall640.fFillStyle = kFillSolid;
  vall640.fFont = ufont->GetFontHandle();
  vall640.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall640, kTRUE);
  TGLabel* fLabel640 = new TGLabel(fCompositeFrame637, "Argument Type:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel640->SetTextJustify(33);
  fLabel640->SetMargins(0, 0, 0, 0);
  fLabel640->SetWrapLength(-1);
  fCompositeFrame637->AddFrame(fLabel640, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  fLabel640->MoveResize(8, 56, 120, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes "Argument Value"
  GCValues_t vall641;
  vall641.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall641.fForeground);
  gClient->GetColorByName("#c0c0c0", vall641.fBackground);
  vall641.fFillStyle = kFillSolid;
  vall641.fFont = ufont->GetFontHandle();
  vall641.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall641, kTRUE);
  TGLabel* fLabel641 = new TGLabel(fCompositeFrame637, "Argument Value:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel641->SetTextJustify(33);
  fLabel641->SetMargins(0, 0, 0, 0);
  fLabel641->SetWrapLength(-1);
  fCompositeFrame637->AddFrame(fLabel641, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  fLabel641->MoveResize(8, 80, 120, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes "Argument Units"
  GCValues_t vall642;
  vall642.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall642.fForeground);
  gClient->GetColorByName("#c0c0c0", vall642.fBackground);
  vall642.fFillStyle = kFillSolid;
  vall642.fFont = ufont->GetFontHandle();
  vall642.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall642, kTRUE);
  TGLabel* fLabel642 = new TGLabel(fCompositeFrame637, "Argument Units:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel642->SetTextJustify(33);
  fLabel642->SetMargins(0, 0, 0, 0);
  fLabel642->SetWrapLength(-1);
  fCompositeFrame637->AddFrame(fLabel642, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2));
  fLabel642->MoveResize(240, 80, 120, 18);
  ULong_t ucolor;        // will reflect user color changes
  gClient->GetColorByName("#ffffff", ucolor);
  // graphics context changes "Argument Value"
  GCValues_t vall653;
  vall653.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall653.fForeground);
  gClient->GetColorByName("#c0c0c0", vall653.fBackground);
  vall653.fFillStyle = kFillSolid;
  vall653.fFont = ufont->GetFontHandle();
  vall653.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall653, kTRUE);
  TGLabel* fLabel653 = new TGLabel(fCompositeFrame637, "Sigma Value:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel653->SetTextJustify(33);
  fLabel653->SetMargins(0, 0, 0, 0);
  fLabel653->SetWrapLength(-1);
  fCompositeFrame637->AddFrame(fLabel653, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  fLabel653->MoveResize(8, 104, 120, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes "Argument Units"
  GCValues_t vall652;
  vall652.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall652.fForeground);
  gClient->GetColorByName("#c0c0c0", vall652.fBackground);
  vall652.fFillStyle = kFillSolid;
  vall652.fFont = ufont->GetFontHandle();
  vall652.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall652, kTRUE);
  TGLabel* fLabel652 = new TGLabel(fCompositeFrame637, "Sigma Units:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel652->SetTextJustify(33);
  fLabel652->SetMargins(0, 0, 0, 0);
  fLabel652->SetWrapLength(-1);
  fCompositeFrame637->AddFrame(fLabel652, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2));
  fLabel652->MoveResize(240, 104, 120, 18);
  gClient->GetColorByName("#ffffff", ucolor);

  // combo box "Projectile Type"
  cProjectilePDG = new TGComboBox(fCompositeFrame637, -1, kHorizontalFrame | kSunkenFrame |
                                                          kDoubleBorder | kOwnBackground);
  cProjectilePDG->AddEntry("p",   0);
  cProjectilePDG->AddEntry("n",   1);
  cProjectilePDG->AddEntry("pi-", 2);
  cProjectilePDG->AddEntry("pi+", 3);
  cProjectilePDG->AddEntry("K+",  4);
  cProjectilePDG->AddEntry("K-",  5); // Add other projectiles only in the end
  cProjectilePDG->Resize(160,22);
  fCompositeFrame637->AddFrame(cProjectilePDG, new TGLayoutHints(kLHintsLeft |
                                                                 kLHintsTop, 2, 2, 2, 2) );
  cProjectilePDG->Select(0);
  cProjectilePDG->MoveResize(136, 8, 160, 22);
  gClient->GetColorByName("#ffffff", ucolor);

  // text box for the Comment
  tComment = new TGTextEntry(fCompositeFrame637, new TGTextBuffer(15), -1, uGC->GetGC(),
                            ufont->GetFontStruct(), kSunkenFrame | kDoubleBorder |
                                                    kOwnBackground);
  tComment->SetMaxLength(32);
  tComment->SetAlignment(kTextLeft);
  tComment->SetText("prc80");
  tComment->Resize(120,tComment->GetDefaultHeight());
  fCompositeFrame637->AddFrame(tComment, new TGLayoutHints(kLHintsLeft |
                                                           kLHintsTop, 2, 2, 2, 2) );
  tComment->MoveResize(358, 8, 48, 22);

  // combo box "Target Name & A"
  cTargetPDG = new TGComboBox(fCompositeFrame637, -1, kHorizontalFrame | kSunkenFrame |
                              kDoubleBorder | kOwnBackground);
  for(Int_t i=0; i<111; ++i)
    cTargetPDG->AddEntry(TString::Format("%d - %s", i+1,
                                         gParticlesDAL->GetElementName(i+1).Data()).Data(),
                                         i );
  cTargetPDG->Resize(160, 22);
  fCompositeFrame637->AddFrame(cTargetPDG, new TGLayoutHints(kLHintsLeft |
                                                             kLHintsTop, 2, 2, 2, 2) );
  cTargetPDG->Select(0);
  cTargetPDG->MoveResize(136, 32, 160, 22);
  nTargetA = new TGNumberEntry(fCompositeFrame637, 0., 19, -1, (TGNumberFormat::EStyle) 5);
  fCompositeFrame637->AddFrame(nTargetA, new TGLayoutHints(kLHintsLeft |
                                                           kLHintsTop, 2, 2, 2, 2) );
  nTargetA->MoveResize(358, 32, 48, 22);
  gClient->GetColorByName("#ffffff", ucolor);
  // combo box "Projectile Energy"
  cArgType = new TGComboBox(fCompositeFrame637,-1,kHorizontalFrame | kSunkenFrame |
                            kDoubleBorder | kOwnBackground);
  cArgType->AddEntry("T_kin", T_kin);
  cArgType->AddEntry("p_mom", p_mom);
  cArgType->AddEntry("E_tot", E_tot);
  cArgType->Resize(160, 22);
  fCompositeFrame637->AddFrame(cArgType, new TGLayoutHints(kLHintsLeft |
                                                           kLHintsTop, 2, 2, 2, 2) );
  cArgType->Select(T_kin);
  cArgType->MoveResize(136, 56, 160, 22);
  // digi box with ArgValue
  nArgValue = new TGNumberEntry(fCompositeFrame637, 1., 12, -1,
                                (TGNumberFormat::EStyle) 5);
  fCompositeFrame637->AddFrame(nArgValue, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2));
  nArgValue->MoveResize(136, 80, 96, 22);
  gClient->GetColorByName("#ffffff", ucolor);
  // combo box "Projectile Units"
  cArgUnits = new TGComboBox(fCompositeFrame637, -1, kHorizontalFrame | kSunkenFrame |
                             kDoubleBorder | kOwnBackground);
  cArgUnits->AddEntry("MeV ", MeV);
  cArgUnits->AddEntry("GeV ", GeV);
  cArgUnits->AddEntry("TeV ", TeV);
  cArgUnits->AddEntry("keV ", keV);
  cArgUnits->Resize(160, 22);
  fCompositeFrame637->AddFrame(cArgUnits, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  cArgUnits->Select(MeV);
  cArgUnits->MoveResize(368, 80, 160, 22);
  // digi box with SigmaValue
  nSigValue = new TGNumberEntry(fCompositeFrame637, 1., 12, -1,
                                (TGNumberFormat::EStyle) 5);
  fCompositeFrame637->AddFrame(nSigValue, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2));
  nSigValue->MoveResize(136, 104, 96, 22);
  gClient->GetColorByName("#ffffff", ucolor);
  // combo box "Sigma Units"
  cSigUnits = new TGComboBox(fCompositeFrame637, -1, kHorizontalFrame | kSunkenFrame |
                             kDoubleBorder | kOwnBackground);
  cSigUnits->AddEntry("mb ",    mBarn);
  cSigUnits->AddEntry("Barn ",  Barn);
  cSigUnits->AddEntry("mkb ",   mkBarn);
  cSigUnits->Resize(160, 22);
  fCompositeFrame637->AddFrame(cSigUnits, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  cSigUnits->Select(mBarn);
  cSigUnits->MoveResize(368, 104, 160, 22);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall719;
  vall719.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall719.fForeground);
  gClient->GetColorByName("#c0c0c0", vall719.fBackground);
  vall719.fFillStyle = kFillSolid;
  vall719.fFont = ufont->GetFontHandle();
  vall719.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall719, kTRUE);
  TGLabel* fLabel719 = new TGLabel(fCompositeFrame637, "ASCII File:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel719->SetTextJustify(33);
  fLabel719->SetMargins(0, 0, 0, 0);
  fLabel719->SetWrapLength(-1);
  fCompositeFrame637->AddFrame(fLabel719, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2));
  fLabel719->MoveResize(8, 136, 144, 18);

  // "Secondary Particle " group frame
  gSecondary = new TGGroupFrame(fCompositeFrame637, "Secondary Particle Parameters");
  gSecondary->SetLayoutBroken(kTRUE);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall721;
  vall721.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000",vall721.fForeground);
  gClient->GetColorByName("#c0c0c0",vall721.fBackground);
  vall721.fFillStyle = kFillSolid;
  vall721.fFont = ufont->GetFontHandle();
  vall721.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall721, kTRUE);
  TGLabel* fLabel721 = new TGLabel(gSecondary,"*>Repeat step2 for each PublicationSubItem",
                                   uGC->GetGC(), ufont->GetFontStruct());
  fLabel721->SetTextJustify(33);
  fLabel721->SetMargins(0, 0, 0, 0);
  fLabel721->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel721, new TGLayoutHints(kLHintsLeft | kLHintsTop,2, 2, 2, 2) );
  fLabel721->MoveResize(8,254,520,26);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall722;
  vall722.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall722.fForeground);
  gClient->GetColorByName("#c0c0c0", vall722.fBackground);
  vall722.fFillStyle = kFillSolid;
  vall722.fFont = ufont->GetFontHandle();
  vall722.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall722, kTRUE);
  TGLabel* fLabel722 = new TGLabel(gSecondary, "Secondary Particle PDG:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel722->SetTextJustify(33);
  fLabel722->SetMargins(0, 0, 0, 0);
  fLabel722->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel722, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel722->MoveResize(8, 24, 248, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

  // graphics context changes
  GCValues_t vall723;
  vall723.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall723.fForeground);
  gClient->GetColorByName("#c0c0c0", vall723.fBackground);
  vall723.fFillStyle = kFillSolid;
  vall723.fFont = ufont->GetFontHandle();
  vall723.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall723, kTRUE);
  TGLabel* fLabel723 = new TGLabel(gSecondary, "Cut Type:", uGC->GetGC(),
                                   ufont->GetFontStruct() );
  fLabel723->SetTextJustify(33);
  fLabel723->SetMargins(0, 0, 0, 0);
  fLabel723->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel723, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel723->MoveResize(8,48,248,18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall724;
  vall724.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall724.fForeground);
  gClient->GetColorByName("#c0c0c0", vall724.fBackground);
  vall724.fFillStyle = kFillSolid;
  vall724.fFont = ufont->GetFontHandle();
  vall724.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall724, kTRUE);
  TGLabel* fLabel724 = new TGLabel(gSecondary, "Cut Value:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel724->SetTextJustify(33);
  fLabel724->SetMargins(0, 0, 0, 0);
  fLabel724->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel724, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel724->MoveResize(8, 72, 248, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall725;
  vall725.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall725.fForeground);
  gClient->GetColorByName("#c0c0c0", vall725.fBackground);
  vall725.fFillStyle = kFillSolid;
  vall725.fFont = ufont->GetFontHandle();
  vall725.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall725, kTRUE);
  TGLabel* fLabel725 = new TGLabel(gSecondary, "Cut Delta:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel725->SetTextJustify(33);
  fLabel725->SetMargins(0, 0, 0, 0);
  fLabel725->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel725, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel725->MoveResize(8, 96, 248, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall726;
  vall726.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall726.fForeground);
  gClient->GetColorByName("#c0c0c0", vall726.fBackground);
  vall726.fFillStyle = kFillSolid;
  vall726.fFont = ufont->GetFontHandle();
  vall726.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall726, kTRUE);
  TGLabel* fLabel726 = new TGLabel(gSecondary, "Cut Units:", uGC->GetGC(),
                                  ufont->GetFontStruct());
  fLabel726->SetTextJustify(33);
  fLabel726->SetMargins(0, 0, 0, 0);
  fLabel726->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel726, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel726->MoveResize(8,120,248,18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall727;
  vall727.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall727.fForeground);
  gClient->GetColorByName("#c0c0c0", vall727.fBackground);
  vall727.fFillStyle = kFillSolid;
  vall727.fFont = ufont->GetFontHandle();
  vall727.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall727, kTRUE);
  TGLabel* fLabel727 = new TGLabel(gSecondary, "Function Type:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel727->SetTextJustify(33);
  fLabel727->SetMargins(0, 0, 0, 0);
  fLabel727->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel727, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel727->MoveResize(8, 144, 248, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall728;
  vall728.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall728.fForeground);
  gClient->GetColorByName("#c0c0c0", vall728.fBackground);
  vall728.fFillStyle = kFillSolid;
  vall728.fFont = ufont->GetFontHandle();
  vall728.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall728, kTRUE);
  TGLabel* fLabel728 = new TGLabel(gSecondary, "Function Units:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel728->SetTextJustify(33);
  fLabel728->SetMargins(0, 0, 0, 0);
  fLabel728->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel728, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel728->MoveResize(8, 168, 248, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall729;
  vall729.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall729.fForeground);
  gClient->GetColorByName("#c0c0c0", vall729.fBackground);
  vall729.fFillStyle = kFillSolid;
  vall729.fFont = ufont->GetFontHandle();
  vall729.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall729, kTRUE);
  TGLabel* fLabel729 = new TGLabel(gSecondary, "FunctionErrorType:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel729->SetTextJustify(33);
  fLabel729->SetMargins(0, 0, 0, 0);
  fLabel729->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel729, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel729->MoveResize(8, 192, 248, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall730;
  vall730.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall730.fForeground);
  gClient->GetColorByName("#c0c0c0", vall730.fBackground);
  vall730.fFillStyle = kFillSolid;
  vall730.fFont = ufont->GetFontHandle();
  vall730.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall730, kTRUE);
  TGLabel* fLabel730 = new TGLabel(gSecondary, "Argument Type:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel730->SetTextJustify(33);
  fLabel730->SetMargins(0, 0, 0, 0);
  fLabel730->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel730, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel730->MoveResize(8, 216, 248, 18);
  ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");
  // graphics context changes
  GCValues_t vall731;
  vall731.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                  kGCGraphicsExposures;
  gClient->GetColorByName("#000000", vall731.fForeground);
  gClient->GetColorByName("#c0c0c0", vall731.fBackground);
  vall731.fFillStyle = kFillSolid;
  vall731.fFont = ufont->GetFontHandle();
  vall731.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&vall731, kTRUE);
  TGLabel* fLabel731 = new TGLabel(gSecondary,"Argument Units:", uGC->GetGC(),
                                   ufont->GetFontStruct());
  fLabel731->SetTextJustify(33);
  fLabel731->SetMargins(0, 0, 0, 0);
  fLabel731->SetWrapLength(-1);
  gSecondary->AddFrame(fLabel731, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  fLabel731->MoveResize(8, 240, 248, 18);

  gClient->GetColorByName("#ffffff",ucolor);

  // Text input for the file name
  tFilename = new TGTextEntry(fCompositeFrame637, new TGTextBuffer(15), -1, uGC->GetGC(),
                              ufont->GetFontStruct(), kSunkenFrame | kDoubleBorder |
                                                      kOwnBackground);
  tFilename->SetMaxLength(4096);
  tFilename->SetAlignment(kTextLeft);
  tFilename->SetText(TString::Format("%s/",  gSystem->WorkingDirectory()).Data());
  tFilename->Resize(384,tFilename->GetDefaultHeight());
  fCompositeFrame637->AddFrame(tFilename, new TGLayoutHints(kLHintsLeft |
                                                            kLHintsTop, 2, 2, 2, 2) );
  tFilename->MoveResize(160, 136, 384, 22);

  // combo box "Secondary Particle PDG"
  cSecondaryPDG = new TGComboBox(gSecondary, -1, kHorizontalFrame | kSunkenFrame |
                                                 kDoubleBorder | kOwnBackground);
  cSecondaryPDG->AddEntry("p",  0);
  cSecondaryPDG->AddEntry("n",  1);
  cSecondaryPDG->AddEntry("d",  2);
  cSecondaryPDG->AddEntry("t",  3);
  cSecondaryPDG->AddEntry("He3",4);
  cSecondaryPDG->AddEntry("He4",5);
  cSecondaryPDG->AddEntry("pi+",6);
  cSecondaryPDG->AddEntry("pi-",7);
  cSecondaryPDG->AddEntry("K+", 8);
  cSecondaryPDG->AddEntry("K-", 9);
  cSecondaryPDG->Resize(264, 22);
  gSecondary->AddFrame(cSecondaryPDG, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  cSecondaryPDG->Select(0);
  cSecondaryPDG->MoveResize(264, 24, 264, 22);
  gClient->GetColorByName("#ffffff", ucolor);
  // combo box "Cut Type"
  cCutType = new TGComboBox(gSecondary, -1, kHorizontalFrame | kSunkenFrame |
                                            kDoubleBorder | kOwnBackground );
  cCutType->AddEntry("Theta",    Theta);
  cCutType->AddEntry("CosTheta", CosTheta);
  cCutType->AddEntry("eta",      LogTgHalfTheta);
  cCutType->AddEntry("p_mom",    p_mom);
  cCutType->AddEntry("E_tot",    E_tot);
  cCutType->AddEntry("T_kin",    T_kin);
  cCutType->AddEntry("p_L",      p_L);
  cCutType->AddEntry("p_T",      p_T);
  cCutType->AddEntry("x",        x_F);
  cCutType->AddEntry("y",        y_R);
  cCutType->Resize(264, 22);
  gSecondary->AddFrame(cCutType, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  cCutType->Select(Theta);
  cCutType->MoveResize(264, 48, 264, 22);

  // combo box "Cut Value"
  nCutValue = new TGNumberEntry(gSecondary, 11., 19, -1, (TGNumberFormat::EStyle) 5);
  gSecondary->AddFrame(nCutValue, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  nCutValue->MoveResize(264, 72, 152, 22);

  // combo box "Cut Delta"
  nCutDelta = new TGNumberEntry(gSecondary, 5., 19, -1, (TGNumberFormat::EStyle) 5);
  gSecondary->AddFrame(nCutDelta, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  nCutDelta->MoveResize(264, 96, 152, 22);
  gClient->GetColorByName("#ffffff", ucolor);

  // combo box "Cut Units"
  cCutUnits = new TGComboBox(gSecondary, -1, kHorizontalFrame | kSunkenFrame |
                                             kDoubleBorder | kOwnBackground );
  cCutUnits->AddEntry("MeV",         MeV);
  cCutUnits->AddEntry("GeV",         GeV);
  cCutUnits->AddEntry("TeV",         TeV);
  cCutUnits->AddEntry("keV",         keV);
  cCutUnits->AddEntry("Degrees",     Degrees);
  cCutUnits->AddEntry("NoDimension", NoDim);
  cCutUnits->AddEntry("Radians",     Radians);

  cCutUnits->Resize(264, 22);
  gSecondary->AddFrame(cCutUnits, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
  cCutUnits->Select(Degrees);
  cCutUnits->MoveResize(264, 120, 264, 22);
  gClient->GetColorByName("#ffffff", ucolor);
  // combo box "Function Type"
  cFunctionType = new TGComboBox(gSecondary, -1, kHorizontalFrame | kSunkenFrame |
                                                 kDoubleBorder | kOwnBackground );
  cFunctionType->AddEntry("dN/dEdO",    dN_dEdO);
  cFunctionType->AddEntry("dS/dEdO",    dS_dEdO);
  cFunctionType->AddEntry("dN/dpdO",    dN_dpdO);
  cFunctionType->AddEntry("dS/dpdO",    dS_dpdO);
  cFunctionType->AddEntry("dN/pdEdO",   dN_pdEdO);
  cFunctionType->AddEntry("dS/pdEdO",   dS_pdEdO);
  cFunctionType->AddEntry("EdN/d3p",    EdN_d3p);
  cFunctionType->AddEntry("EdS/d3p",    EdS_d3p);
  cFunctionType->AddEntry("EdS/Ad3p",   EdS_Ad3p);
  cFunctionType->AddEntry("dN/dydpT",   dN_dydpT);
  cFunctionType->AddEntry("dS/dydpT",   dS_dydpT);
  cFunctionType->AddEntry("dN/dxdpT",   dN_dxdpT);
  cFunctionType->AddEntry("dS/dxdpT",   dS_dxdpT);
  cFunctionType->AddEntry("dN/dxdET",   dN_dxdET);
  cFunctionType->AddEntry("dS/dxdET",   dS_dxdET);
  cFunctionType->AddEntry("dN/dydET",   dN_dydET);
  cFunctionType->AddEntry("dS/dydET",   dS_dydET);
  cFunctionType->AddEntry("dN/detadET", dN_detadET);
  cFunctionType->AddEntry("dS/detadET", dS_detadET);
  cFunctionType->AddEntry("dS/dt",      dS_dt);
  cFunctionType->AddEntry("dS/du",      dS_du);
  cFunctionType->Resize(264, 22);
  gSecondary->AddFrame(cFunctionType, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  cFunctionType->Select(dS_dEdO);
  cFunctionType->MoveResize(264, 144, 264, 22);
  gClient->GetColorByName("#ffffff", ucolor);
  // combo box "Function Units"
  cFunctionUnits = new TGComboBox(gSecondary, -1, kHorizontalFrame | kSunkenFrame |
                                                   kDoubleBorder | kOwnBackground );
  cFunctionUnits->AddEntry("1/MeV/sr",   _MeV_sr);
  cFunctionUnits->AddEntry("mb/MeV/sr",  mb_MeV_sr);
  cFunctionUnits->AddEntry("1/MeV2/sr",  _MeV2_sr);
  cFunctionUnits->AddEntry("mb/MeV2/sr", mb_MeV2_sr);
  cFunctionUnits->AddEntry("1/MeV",      _MeV);
  cFunctionUnits->AddEntry("mb/MeV",     mb_MeV);
  cFunctionUnits->AddEntry("1/MeV2",     _MeV2);
  cFunctionUnits->AddEntry("mb/MeV2",    mb_MeV2);
  cFunctionUnits->AddEntry("1/GeV/sr",   _GeV_sr);
  cFunctionUnits->AddEntry("mb/GeV/sr",  mb_GeV_sr);
  cFunctionUnits->AddEntry("1/GeV2/sr",  _GeV2_sr);
  cFunctionUnits->AddEntry("mb/GeV2/sr", mb_GeV2_sr);
  cFunctionUnits->AddEntry("1/GeV",      _GeV);
  cFunctionUnits->AddEntry("mb/GeV",     mb_GeV);
  cFunctionUnits->AddEntry("1/GeV2",     _GeV2);
  cFunctionUnits->AddEntry("mb/GeV2",    mb_GeV2);
  cFunctionUnits->AddEntry("1/1",        _1);
  cFunctionUnits->AddEntry("mb/1",       mb_1);
  cFunctionUnits->Resize(264, 22);
  gSecondary->AddFrame(cFunctionUnits, new TGLayoutHints(kLHintsLeft |
                                                          kLHintsTop, 2, 2, 2, 2) );
  cFunctionUnits->Select(mb_MeV_sr);
  cFunctionUnits->MoveResize(264, 168, 264, 22);
  gClient->GetColorByName("#ffffff", ucolor);
  // combo box "FunctionErrorType"
  cFunctErrorType = new TGComboBox(gSecondary, -1, kHorizontalFrame | kSunkenFrame |
                                                  kDoubleBorder | kOwnBackground );
  cFunctErrorType->AddEntry("NoError",   NoError);
  cFunctErrorType->AddEntry("Absolute",  Absolute);
  cFunctErrorType->AddEntry("Relative",  Relative);
  cFunctErrorType->AddEntry("InPerCent", InPerCent);
  cFunctErrorType->Resize(264, 22);
  gSecondary->AddFrame(cFunctErrorType, new TGLayoutHints(kLHintsLeft |
                                                         kLHintsTop, 2, 2, 2, 2));
  cFunctErrorType->Select(Absolute);
  cFunctErrorType->MoveResize(264, 192, 264, 22);
  gClient->GetColorByName("#ffffff", ucolor);
  // combo box "Argument Type"
  cSecArgType = new TGComboBox(gSecondary, -1, kHorizontalFrame | kSunkenFrame |
                                               kDoubleBorder | kOwnBackground);
  cSecArgType->AddEntry("Theta",    Theta);
  cSecArgType->AddEntry("CosTheta", CosTheta);
  cSecArgType->AddEntry("eta",      LogTgHalfTheta);
  cSecArgType->AddEntry("p_mom",    p_mom);
  cSecArgType->AddEntry("E_tot",    E_tot);
  cSecArgType->AddEntry("T_kin",    T_kin);
  cSecArgType->AddEntry("p_L",      p_L);
  cSecArgType->AddEntry("p_T",      p_T);
  cSecArgType->AddEntry("x",        x_F);
  cSecArgType->AddEntry("y",        y_R);
  cSecArgType->Resize(264, 22);
  gSecondary->AddFrame(cSecArgType, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2,2,2,2) );
  cSecArgType->Select(T_kin);
  cSecArgType->MoveResize(264, 216, 264, 22);
  gClient->GetColorByName("#ffffff", ucolor);
  // combo box "Argument Units"
  cSecArgUnits = new TGComboBox(gSecondary, -1, kHorizontalFrame | kSunkenFrame |
                                                kDoubleBorder | kOwnBackground);
  cSecArgUnits->AddEntry("NoDimension", NoDim);
  cSecArgUnits->AddEntry("Degrees",     Degrees);
  cSecArgUnits->AddEntry("Radians",     Radians);
  cSecArgUnits->AddEntry("keV",         keV);
  cSecArgUnits->AddEntry("MeV",         MeV);
  cSecArgUnits->AddEntry("GeV",         GeV);
  cSecArgUnits->AddEntry("TeV",         TeV);
  cSecArgUnits->Resize(264, 22);
  gSecondary->AddFrame(cSecArgUnits, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2) );
  cSecArgUnits->Select(MeV);
  cSecArgUnits->MoveResize(264, 240, 264, 22);
  //---------------------------------------------------------------- Bottom commands ------
  gSecondary->SetLayoutManager(new TGVerticalLayout(gSecondary));
  gSecondary->Resize(536,288);
  fCompositeFrame637->AddFrame(gSecondary, new TGLayoutHints(kLHintsLeft |
                                                             kLHintsTop, 2, 2, 2, 2) );
  gSecondary->MoveResize(8, 160, 536, 288);
  bRun2 = new TGTextButton(this,"2. Add Secondary",2);
  bRun2->SetTextJustify(36);
  bRun2->SetMargins(0, 0, 0, 0);
  bRun2->SetWrapLength(-1);
  bRun2->Resize(140, 24);
  fCompositeFrame637->AddFrame(bRun2, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  bRun2->MoveResize(264, 456, 140, 24);
  
  bRun0 = new TGTextButton(this,"Close ROOT",4);
  bRun0->SetTextJustify(36);
  bRun0->SetMargins(0, 0, 0, 0);
  bRun0->SetWrapLength(-1);
  //bRun0->Resize(136,24);
  fCompositeFrame637->AddFrame(bRun0, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  bRun0->MoveResize(20, 456, 100, 24);

  bRun1 = new TGTextButton(this,"1. Create Publication",1);
  bRun1->SetTextJustify(36);
  bRun1->SetMargins(0, 0, 0, 0);
  bRun1->SetWrapLength(-1);
  bRun1->Resize(136, 24);
  fCompositeFrame637->AddFrame(bRun1, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  bRun1->MoveResize(124, 456, 136, 24);
  ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
  // graphics context changes
  GCValues_t valEntry867;
  valEntry867.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont |
                      kGCGraphicsExposures;
  gClient->GetColorByName("#000000", valEntry867.fForeground);
  gClient->GetColorByName("#c0c0c0", valEntry867.fBackground);
  valEntry867.fFillStyle = kFillSolid;
  valEntry867.fFont = ufont->GetFontHandle();
  valEntry867.fGraphicsExposures = kFALSE;
  uGC = gClient->GetGC(&valEntry867, kTRUE);

  bRun3 = new TGTextButton(this,"3. Create File",3);
  bRun3->SetTextJustify(36);
  bRun3->SetMargins(0, 0, 0, 0);
  bRun3->SetWrapLength(-1);
  bRun3->Resize(136, 24);
  fCompositeFrame637->AddFrame(bRun3, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
  bRun3->MoveResize(408, 456, 136, 24);
  fCompositeFrame636->AddFrame(fCompositeFrame637, new TGLayoutHints(kLHintsExpandX |
                                                                     kLHintsExpandY));
  fCompositeFrame637->MoveResize(0, 0, 551, 489);
  fMainFrame635->AddFrame(fCompositeFrame636, new TGLayoutHints(kLHintsExpandX |
                                                                kLHintsExpandY));
  fCompositeFrame636->MoveResize(0, 0, 551, 487);
  this->AddFrame(fMainFrame635, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  fMainFrame635->MoveResize(0, 0, 549, 487);

  // set the initial state
  ChangeState(0);

  SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputModeless);
  MapSubwindows();

  MapWindow();
  Resize(549, 487);
}

//______________________________________________________________________________
void G4TPublicationGUI::Run1_Click() // Create Publication
{
  // get values
  TString comment(tComment->GetText());
  Int_t projectilePDG =
            gParticlesDAL->GetPDG(TString(cProjectilePDG->GetSelectedEntry()->GetTitle()));
  Int_t targetZ       = cTargetPDG->GetSelected() + 1;
  Int_t targetA       = nTargetA->GetIntNumber();
  ArgEnum argType     = (ArgEnum) cArgType->GetSelected();
  Double_t argValue   = nArgValue->GetNumber();
  UnitsEnum argUnits  = (UnitsEnum)cArgUnits->GetSelected();
  Double_t sigValue   = nSigValue->GetNumber();
  SigmaEnum sigUnits  = (SigmaEnum)cSigUnits->GetSelected();
  Int_t targetPDG     = gParticlesDAL->GetPDG(targetZ, targetA);
  if(gParticlesDAL->GetParticle(targetPDG) == 0)
  {
    cout<<"Error>G4TPublicationGUI::Run1_Click: target with Z="<<targetZ<<", A="<<targetA
        <<" isn't found"<<endl;
    return;
  }
  cout<<"G4TPublicationGUI::Run1_Click:pPDG="<<projectilePDG<<",tPDG="<<targetPDG<<",comt="
      <<comment<<",ArgType="<<argType<<",ArgVal="<< argValue<<",ArgUnits="<<argUnits<<endl;
  //                                                  *publ*
  fPublication = new G4TData(projectilePDG, targetPDG, true, "data", comment, argType,
                             argValue, argUnits, sigValue, sigUnits, 6); // 6 = color type 
  // go to next state
  ChangeState(1);
}

//______________________________________________________________________________
void G4TPublicationGUI::Run2_Click() // Add Secondary
{
  TString filename(tFilename->GetText());
  ifstream asciiFile;
  asciiFile.open(filename.Data());
  if (!asciiFile)
  {
    cout<<"Error>G4TPublicationGUI::Run2_Click: No file="<<filename.Data()<<endl;
    return;
  }
  Int_t secondaryPDG      =
             gParticlesDAL->GetPDG(TString(cSecondaryPDG->GetSelectedEntry()->GetTitle()));
  ArgEnum   cutType       = (ArgEnum)   cCutType->GetSelected();
  UnitsEnum cutUnits      = (UnitsEnum) cCutUnits->GetSelected();
  Double_t  cutValue      = nCutValue->GetNumber();
  Double_t  cutDelta      = nCutDelta->GetNumber();
  FunctEnum functionType  = (FunctEnum) cFunctionType->GetSelected();
  FunUnEnum functionUnits = (FunUnEnum) cFunctionUnits->GetSelected();
  ErrorType functErrType  = (ErrorType) cFunctErrorType->GetSelected();
  ArgEnum   secArgType    = (ArgEnum)   cSecArgType->GetSelected();
  UnitsEnum secArgUnits   = (UnitsEnum) cSecArgUnits->GetSelected();

  cout<<"G4TPublicationGUI::Run2_Click: PDG="<<secondaryPDG<<",CutTp="<<cutType<<",CutUn="
      <<cutUnits<<",CutVl="<<cutValue<<",CutDv="<<cutDelta<<",FunTp="<<functionType
      <<",FunUn="<<functionUnits<<",FunEr="<<functErrType<<",ArgTp="<<secArgType<<",ArgUn="
      <<secArgUnits<<endl;

  cout<<"G4TPublicationGUI::Run2_Click: Loading ASCII file " << filename << "..." << endl;
  fPublication->AddItem(secondaryPDG,
                        cutType,
                        cutUnits,
                        cutValue,
                        cutDelta,
                        functionType,
                        functionUnits,
                        functErrType,
                        secArgType,
                        secArgUnits)->LoadFromASCII(filename);
}

//______________________________________________________________________________
void G4TPublicationGUI::Run3_Click() // Create File
{
  cout<<"G4TPublicationGUI::Run3_Click: Saving the publication file... "<<endl;
  gTestDB->SaveData(fPublication);
  cout<<"G4TPublicationGUI::Run3_Click: Data are saved. Done!"<<endl;

  // go to new publication state
  ChangeState(0);
}

//______________________________________________________________________________
void G4TPublicationGUI::Run4_Click() // Close ROOT
{
  gApplication->Terminate(0); // Stops ROOT application
  //this->DestroyWindow();    // Only destroyes the Publication window
}

//______________________________________________________________________________
void G4TPublicationGUI::ChangeState(Int_t state)
{
  if(!state)
  {
    bRun1->SetEnabled(true);
    bRun2->SetEnabled(false);
    bRun3->SetEnabled(false);
    cProjectilePDG->SetEnabled(true);
    tComment->SetEnabled(true);
    cTargetPDG->SetEnabled(true);
    nTargetA->SetState(true);
    cArgType->SetEnabled(true);
    nArgValue->SetState(true);
    cArgUnits->SetEnabled(true);
    nSigValue->SetState(true);
    cSigUnits->SetEnabled(true);
    tFilename->SetEnabled(false);
    cSecondaryPDG->SetEnabled(false);
    cCutType->SetEnabled(false);
    nCutValue->SetState(false);
    nCutDelta->SetState(false);
    cCutUnits->SetEnabled(false);
    cFunctionType->SetEnabled(false);
    cFunctionUnits->SetEnabled(false);
    cFunctErrorType->SetEnabled(false);
    cSecArgType->SetEnabled(false);
    cSecArgUnits->SetEnabled(false);
  }
  else if(state == 1)
  {
    bRun1->SetEnabled(false);
    bRun2->SetEnabled(true);
    bRun3->SetEnabled(true);
    cProjectilePDG->SetEnabled(false);
    tComment->SetEnabled(false);
    cTargetPDG->SetEnabled(false);
    nTargetA->SetState(false);
    cArgType->SetEnabled(false);
    nArgValue->SetState(false);
    cArgUnits->SetEnabled(false);
    nSigValue->SetState(false);
    cSigUnits->SetEnabled(false);
    tFilename->SetEnabled(true);
    cSecondaryPDG->SetEnabled(true);
    cCutType->SetEnabled(true);
    nCutValue->SetState(true);
    nCutDelta->SetState(true);
    cCutUnits->SetEnabled(true);
    cFunctionType->SetEnabled(true);
    cFunctionUnits->SetEnabled(true);
    cFunctErrorType->SetEnabled(true);
    cSecArgType->SetEnabled(true);
    cSecArgUnits->SetEnabled(true);
  }
}

//______________________________________________________________________________
G4TPublicationGUI::~G4TPublicationGUI()
{
  TGFrameElement *ptr;
  // delete all frames and layout hints
  if (fList)
  {
    TIter next(fList);
    while ((ptr = (TGFrameElement *) next()))
    {
      if (ptr->fLayout) delete ptr->fLayout;
      if (ptr->fFrame) delete ptr->fFrame;
    }
  }
}

//______________________________________________________________________________
void G4TPublicationGUI::CloseWindow() // for the x of the Window
{
  //Terminate application
  //gApplication->Terminate(0);
  this->DestroyWindow();
}


//______________________________________________________________________________
Bool_t G4TPublicationGUI::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
  switch (GET_MSG(msg))
  {
    case kC_COMMAND:
      switch (GET_SUBMSG(msg))
      {
        case kCM_BUTTON:
          if     (parm1 == 1) Run1_Click();
          else if(parm1 == 2) Run2_Click();
          else if(parm1 == 3) Run3_Click();
          else if(parm1 == 4) Run4_Click();
          break;
        default:
          break;
      }
    default:
      break;
  }
  return kTRUE;
}
