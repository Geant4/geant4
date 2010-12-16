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
  :TGMainFrame(p, w, h)
{

 G4TSimHelper::LoadLibraries();
 Initialize();

}

//______________________________________________________________________________
void G4TPublicationGUI::Initialize()
{
    // composite frame
    TGCompositeFrame *fMainFrame635 = new TGCompositeFrame(this,549,487,kVerticalFrame);
    fMainFrame635->SetLayoutBroken(kTRUE);

    // composite frame
    TGCompositeFrame *fCompositeFrame636 = new TGCompositeFrame(fMainFrame635,551,487,kVerticalFrame);
    fCompositeFrame636->SetLayoutBroken(kTRUE);

    // composite frame
    TGCompositeFrame *fCompositeFrame637 = new TGCompositeFrame(fCompositeFrame636,551,489,kVerticalFrame);
    fCompositeFrame637->SetLayoutBroken(kTRUE);

    TGFont *ufont;         // will reflect user font changes
    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    TGGC   *uGC;           // will reflect user GC changes
    // graphics context changes
    GCValues_t vall638;
    vall638.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall638.fForeground);
    gClient->GetColorByName("#c0c0c0",vall638.fBackground);
    vall638.fFillStyle = kFillSolid;
    vall638.fFont = ufont->GetFontHandle();
    vall638.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall638, kTRUE);
    TGLabel *fLabel638 = new TGLabel(fCompositeFrame637,"Projectile PDG:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel638->SetTextJustify(33);
    fLabel638->SetMargins(0,0,0,0);
    fLabel638->SetWrapLength(-1);
    fCompositeFrame637->AddFrame(fLabel638, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel638->MoveResize(8,8,144,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall639;
    vall639.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall639.fForeground);
    gClient->GetColorByName("#c0c0c0",vall639.fBackground);
    vall639.fFillStyle = kFillSolid;
    vall639.fFont = ufont->GetFontHandle();
    vall639.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall639, kTRUE);
    TGLabel *fLabel639 = new TGLabel(fCompositeFrame637,"Target (Z, A):",uGC->GetGC(),ufont->GetFontStruct());
    fLabel639->SetTextJustify(33);
    fLabel639->SetMargins(0,0,0,0);
    fLabel639->SetWrapLength(-1);
    fCompositeFrame637->AddFrame(fLabel639, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel639->MoveResize(8,32,144,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall640;
    vall640.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall640.fForeground);
    gClient->GetColorByName("#c0c0c0",vall640.fBackground);
    vall640.fFillStyle = kFillSolid;
    vall640.fFont = ufont->GetFontHandle();
    vall640.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall640, kTRUE);
    TGLabel *fLabel640 = new TGLabel(fCompositeFrame637,"Argument Type:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel640->SetTextJustify(33);
    fLabel640->SetMargins(0,0,0,0);
    fLabel640->SetWrapLength(-1);
    fCompositeFrame637->AddFrame(fLabel640, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel640->MoveResize(8,56,144,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall641;
    vall641.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall641.fForeground);
    gClient->GetColorByName("#c0c0c0",vall641.fBackground);
    vall641.fFillStyle = kFillSolid;
    vall641.fFont = ufont->GetFontHandle();
    vall641.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall641, kTRUE);
    TGLabel *fLabel641 = new TGLabel(fCompositeFrame637,"Argument Value:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel641->SetTextJustify(33);
    fLabel641->SetMargins(0,0,0,0);
    fLabel641->SetWrapLength(-1);
    fCompositeFrame637->AddFrame(fLabel641, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel641->MoveResize(8,80,144,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall642;
    vall642.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall642.fForeground);
    gClient->GetColorByName("#c0c0c0",vall642.fBackground);
    vall642.fFillStyle = kFillSolid;
    vall642.fFont = ufont->GetFontHandle();
    vall642.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall642, kTRUE);
    TGLabel *fLabel642 = new TGLabel(fCompositeFrame637,"Argument Units:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel642->SetTextJustify(33);
    fLabel642->SetMargins(0,0,0,0);
    fLabel642->SetWrapLength(-1);
    fCompositeFrame637->AddFrame(fLabel642, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel642->MoveResize(8,104,144,18);

    ULong_t ucolor;        // will reflect user color changes
    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cProjectilePDG = new TGComboBox(fCompositeFrame637,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    cProjectilePDG->AddEntry("p",0);
    cProjectilePDG->AddEntry("n",1);
    cProjectilePDG->Resize(192,22);
    fCompositeFrame637->AddFrame(cProjectilePDG, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cProjectilePDG->Select(0);
    cProjectilePDG->MoveResize(160,8,192,22);

    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cTargetPDG = new TGComboBox(fCompositeFrame637,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);

    for(Int_t i=0; i<111;i++)
    {
     cTargetPDG->AddEntry(TString::Format("%d - %s", i+1, gParticlesDAL->GetElementName(i+1).Data()).Data(),i);
    }

    cTargetPDG->Resize(192,22);
    fCompositeFrame637->AddFrame(cTargetPDG, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cTargetPDG->Select(0);
    cTargetPDG->MoveResize(160,32,192,22);


    nTargetA = new TGNumberEntry(fCompositeFrame637, (Double_t) 0,19,-1,(TGNumberFormat::EStyle) 5);
    fCompositeFrame637->AddFrame(nTargetA, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    nTargetA->MoveResize(357,32,50,22);



    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cArgType = new TGComboBox(fCompositeFrame637,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    cArgType->AddEntry("E_Kin",E_Kin);
    cArgType->Resize(192,22);
    fCompositeFrame637->AddFrame(cArgType, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cArgType->Select(E_Kin);
    cArgType->MoveResize(160,56,192,22);


    nArgValue = new TGNumberEntry(fCompositeFrame637, (Double_t) 29,12,-1,(TGNumberFormat::EStyle) 5);
    fCompositeFrame637->AddFrame(nArgValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    nArgValue->MoveResize(160,80,104,22);

    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cArgUnits = new TGComboBox(fCompositeFrame637,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    cArgUnits->AddEntry("MeV ",MeV);
    cArgUnits->Resize(192,22);
    fCompositeFrame637->AddFrame(cArgUnits, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cArgUnits->Select(MeV);
    cArgUnits->MoveResize(160,104,192,22);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall719;
    vall719.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall719.fForeground);
    gClient->GetColorByName("#c0c0c0",vall719.fBackground);
    vall719.fFillStyle = kFillSolid;
    vall719.fFont = ufont->GetFontHandle();
    vall719.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall719, kTRUE);
    TGLabel *fLabel719 = new TGLabel(fCompositeFrame637,"ASCII File:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel719->SetTextJustify(33);
    fLabel719->SetMargins(0,0,0,0);
    fLabel719->SetWrapLength(-1);
    fCompositeFrame637->AddFrame(fLabel719, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel719->MoveResize(8,136,144,18);

    // "Secondary Particle (using the angles)" group frame
    gSecondary = new TGGroupFrame(fCompositeFrame637,"Secondary Particle (using the angles)");
    gSecondary->SetLayoutBroken(kTRUE);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall721;
    vall721.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall721.fForeground);
    gClient->GetColorByName("#c0c0c0",vall721.fBackground);
    vall721.fFillStyle = kFillSolid;
    vall721.fFont = ufont->GetFontHandle();
    vall721.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall721, kTRUE);
    TGLabel *fLabel721 = new TGLabel(gSecondary,"NOTE: repeat step2 for each sub-item of the publication",uGC->GetGC(),ufont->GetFontStruct());
    fLabel721->SetTextJustify(33);
    fLabel721->SetMargins(0,0,0,0);
    fLabel721->SetWrapLength(-1);
    gSecondary->AddFrame(fLabel721, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel721->MoveResize(8,248,520,26);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall722;
    vall722.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall722.fForeground);
    gClient->GetColorByName("#c0c0c0",vall722.fBackground);
    vall722.fFillStyle = kFillSolid;
    vall722.fFont = ufont->GetFontHandle();
    vall722.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall722, kTRUE);
    TGLabel *fLabel722 = new TGLabel(gSecondary,"Secondary Particle PDG:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel722->SetTextJustify(33);
    fLabel722->SetMargins(0,0,0,0);
    fLabel722->SetWrapLength(-1);
    gSecondary->AddFrame(fLabel722, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel722->MoveResize(8,24,248,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall723;
    vall723.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall723.fForeground);
    gClient->GetColorByName("#c0c0c0",vall723.fBackground);
    vall723.fFillStyle = kFillSolid;
    vall723.fFont = ufont->GetFontHandle();
    vall723.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall723, kTRUE);
    TGLabel *fLabel723 = new TGLabel(gSecondary,"Cut Type:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel723->SetTextJustify(33);
    fLabel723->SetMargins(0,0,0,0);
    fLabel723->SetWrapLength(-1);
    gSecondary->AddFrame(fLabel723, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel723->MoveResize(8,48,248,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall724;
    vall724.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall724.fForeground);
    gClient->GetColorByName("#c0c0c0",vall724.fBackground);
    vall724.fFillStyle = kFillSolid;
    vall724.fFont = ufont->GetFontHandle();
    vall724.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall724, kTRUE);
    TGLabel *fLabel724 = new TGLabel(gSecondary,"Cut Value:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel724->SetTextJustify(33);
    fLabel724->SetMargins(0,0,0,0);
    fLabel724->SetWrapLength(-1);
    gSecondary->AddFrame(fLabel724, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel724->MoveResize(8,72,248,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall725;
    vall725.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall725.fForeground);
    gClient->GetColorByName("#c0c0c0",vall725.fBackground);
    vall725.fFillStyle = kFillSolid;
    vall725.fFont = ufont->GetFontHandle();
    vall725.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall725, kTRUE);
    TGLabel *fLabel725 = new TGLabel(gSecondary,"Cut Delta:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel725->SetTextJustify(33);
    fLabel725->SetMargins(0,0,0,0);
    fLabel725->SetWrapLength(-1);
    gSecondary->AddFrame(fLabel725, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel725->MoveResize(8,96,248,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall726;
    vall726.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall726.fForeground);
    gClient->GetColorByName("#c0c0c0",vall726.fBackground);
    vall726.fFillStyle = kFillSolid;
    vall726.fFont = ufont->GetFontHandle();
    vall726.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall726, kTRUE);
    TGLabel *fLabel726 = new TGLabel(gSecondary,"Cut Units:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel726->SetTextJustify(33);
    fLabel726->SetMargins(0,0,0,0);
    fLabel726->SetWrapLength(-1);
    gSecondary->AddFrame(fLabel726, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel726->MoveResize(8,120,248,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall727;
    vall727.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall727.fForeground);
    gClient->GetColorByName("#c0c0c0",vall727.fBackground);
    vall727.fFillStyle = kFillSolid;
    vall727.fFont = ufont->GetFontHandle();
    vall727.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall727, kTRUE);
    TGLabel *fLabel727 = new TGLabel(gSecondary,"Function Type:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel727->SetTextJustify(33);
    fLabel727->SetMargins(0,0,0,0);
    fLabel727->SetWrapLength(-1);
    gSecondary->AddFrame(fLabel727, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel727->MoveResize(8,144,248,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall728;
    vall728.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall728.fForeground);
    gClient->GetColorByName("#c0c0c0",vall728.fBackground);
    vall728.fFillStyle = kFillSolid;
    vall728.fFont = ufont->GetFontHandle();
    vall728.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall728, kTRUE);
    TGLabel *fLabel728 = new TGLabel(gSecondary,"Function Units:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel728->SetTextJustify(33);
    fLabel728->SetMargins(0,0,0,0);
    fLabel728->SetWrapLength(-1);
    gSecondary->AddFrame(fLabel728, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel728->MoveResize(8,168,248,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall729;
    vall729.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall729.fForeground);
    gClient->GetColorByName("#c0c0c0",vall729.fBackground);
    vall729.fFillStyle = kFillSolid;
    vall729.fFont = ufont->GetFontHandle();
    vall729.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall729, kTRUE);
    TGLabel *fLabel729 = new TGLabel(gSecondary,"Argument Type:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel729->SetTextJustify(33);
    fLabel729->SetMargins(0,0,0,0);
    fLabel729->SetWrapLength(-1);
    gSecondary->AddFrame(fLabel729, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel729->MoveResize(8,192,248,18);

    ufont = gClient->GetFont("-*-newspaper-(null)-*-*-0-*-*-*-*-*-*-*");

    // graphics context changes
    GCValues_t vall730;
    vall730.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",vall730.fForeground);
    gClient->GetColorByName("#c0c0c0",vall730.fBackground);
    vall730.fFillStyle = kFillSolid;
    vall730.fFont = ufont->GetFontHandle();
    vall730.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&vall730, kTRUE);
    TGLabel *fLabel730 = new TGLabel(gSecondary,"Argument Units:",uGC->GetGC(),ufont->GetFontStruct());
    fLabel730->SetTextJustify(33);
    fLabel730->SetMargins(0,0,0,0);
    fLabel730->SetWrapLength(-1);
    gSecondary->AddFrame(fLabel730, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    fLabel730->MoveResize(8,216,248,18);

    gClient->GetColorByName("#ffffff",ucolor);


    tFilename = new TGTextEntry(fCompositeFrame637, new TGTextBuffer(15),-1,uGC->GetGC(),ufont->GetFontStruct(),kSunkenFrame | kDoubleBorder | kOwnBackground);
    tFilename->SetMaxLength(4096);
    tFilename->SetAlignment(kTextLeft);
    tFilename->SetText(TString::Format("%s/",  gSystem->WorkingDirectory()).Data());
    tFilename->Resize(384,tFilename->GetDefaultHeight());
    fCompositeFrame637->AddFrame(tFilename, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    tFilename->MoveResize(160,136,384,22);


    // combo box
    cSecondaryPDG = new TGComboBox(gSecondary,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    cSecondaryPDG->AddEntry("p",0);
    cSecondaryPDG->AddEntry("n",1);
    cSecondaryPDG->AddEntry("H2",2);
    cSecondaryPDG->AddEntry("H3",3);
    cSecondaryPDG->AddEntry("He3",4);
    cSecondaryPDG->AddEntry("He4",5);
    cSecondaryPDG->Resize(264,22);
    gSecondary->AddFrame(cSecondaryPDG, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cSecondaryPDG->Select(0);
    cSecondaryPDG->MoveResize(264,24,264,22);

    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cCutType = new TGComboBox(gSecondary,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    cCutType->AddEntry("Theta",Theta);
    cCutType->Resize(264,22);
    gSecondary->AddFrame(cCutType, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cCutType->Select(Theta);
    cCutType->MoveResize(264,48,264,22);

    nCutValue = new TGNumberEntry(gSecondary, (Double_t) 11,19,-1,(TGNumberFormat::EStyle) 5);
    gSecondary->AddFrame(nCutValue, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    nCutValue->MoveResize(264,72,152,22);

    nCutDelta = new TGNumberEntry(gSecondary, (Double_t) 5,19,-1,(TGNumberFormat::EStyle) 5);
    gSecondary->AddFrame(nCutDelta, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    nCutDelta->MoveResize(264,96,152,22);

    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cCutUnits = new TGComboBox(gSecondary,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    cCutUnits->AddEntry("Degrees",Degrees);
    cCutUnits->Resize(264,22);
    gSecondary->AddFrame(cCutUnits, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cCutUnits->Select(Degrees);
    cCutUnits->MoveResize(264,120,264,22);

    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cFunctionType = new TGComboBox(gSecondary,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    cFunctionType->AddEntry("dS_over_dE",dS_over_dE);
    cFunctionType->Resize(264,22);
    gSecondary->AddFrame(cFunctionType, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cFunctionType->Select(dS_over_dE);
    cFunctionType->MoveResize(264,144,264,22);

    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cFunctionUnits = new TGComboBox(gSecondary,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    cFunctionUnits->AddEntry("MeV",MeV);
    cFunctionUnits->Resize(264,22);
    gSecondary->AddFrame(cFunctionUnits, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cFunctionUnits->Select(MeV);
    cFunctionUnits->MoveResize(264,168,264,22);

    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cSecArgType = new TGComboBox(gSecondary,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    cSecArgType->AddEntry("E_Kin",E_Kin);
    cSecArgType->Resize(264,22);
    gSecondary->AddFrame(cSecArgType, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cSecArgType->Select(E_Kin);
    cSecArgType->MoveResize(264,192,264,22);

    gClient->GetColorByName("#ffffff",ucolor);

    // combo box
    cSecArgUnits = new TGComboBox(gSecondary,-1,kHorizontalFrame | kSunkenFrame | kDoubleBorder | kOwnBackground);
    cSecArgUnits->AddEntry("MeV",MeV);
    cSecArgUnits->Resize(264,22);
    gSecondary->AddFrame(cSecArgUnits, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    cSecArgUnits->Select(MeV);
    cSecArgUnits->MoveResize(264,216,264,22);

    gSecondary->SetLayoutManager(new TGVerticalLayout(gSecondary));
    gSecondary->Resize(536,288);
    fCompositeFrame637->AddFrame(gSecondary, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    gSecondary->MoveResize(8,160,536,288);
    bRun2 = new TGTextButton(this,"2. Add Secondary",2);
    bRun2->SetTextJustify(36);
    bRun2->SetMargins(0,0,0,0);
    bRun2->SetWrapLength(-1);
    bRun2->Resize(140,24);
    fCompositeFrame637->AddFrame(bRun2, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    bRun2->MoveResize(264,456,140,24);
   
    bRun0 = new TGTextButton(this,"Close ROOT",4);
    bRun0->SetTextJustify(36);
    bRun0->SetMargins(0,0,0,0);
    bRun0->SetWrapLength(-1);
    //bRun0->Resize(136,24);
    fCompositeFrame637->AddFrame(bRun0, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    bRun0->MoveResize(20,456,100,24);

    bRun1 = new TGTextButton(this,"1. Create Publication",1);
    bRun1->SetTextJustify(36);
    bRun1->SetMargins(0,0,0,0);
    bRun1->SetWrapLength(-1);
    bRun1->Resize(136,24);
    fCompositeFrame637->AddFrame(bRun1, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    bRun1->MoveResize(124,456,136,24);

    ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");

    // graphics context changes
    GCValues_t valEntry867;
    valEntry867.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
    gClient->GetColorByName("#000000",valEntry867.fForeground);
    gClient->GetColorByName("#c0c0c0",valEntry867.fBackground);
    valEntry867.fFillStyle = kFillSolid;
    valEntry867.fFont = ufont->GetFontHandle();
    valEntry867.fGraphicsExposures = kFALSE;
    uGC = gClient->GetGC(&valEntry867, kTRUE);

    bRun3 = new TGTextButton(this,"3. Create File",3);
    bRun3->SetTextJustify(36);
    bRun3->SetMargins(0,0,0,0);
    bRun3->SetWrapLength(-1);
    bRun3->Resize(136,24);
    fCompositeFrame637->AddFrame(bRun3, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
    bRun3->MoveResize(408,456,136,24);

    fCompositeFrame636->AddFrame(fCompositeFrame637, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fCompositeFrame637->MoveResize(0,0,551,489);

    fMainFrame635->AddFrame(fCompositeFrame636, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fCompositeFrame636->MoveResize(0,0,551,487);

    this->AddFrame(fMainFrame635, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    fMainFrame635->MoveResize(0,0,549,487);


    // set the initial state
    ChangeState(0);


    this->SetMWMHints(kMWMDecorAll,
                         kMWMFuncAll,
                         kMWMInputModeless);
    this->MapSubwindows();

    this->MapWindow();
    this->Resize(549,487);
}

//______________________________________________________________________________
void G4TPublicationGUI::Run1_Click()
{
 // get values
 Int_t projectilePDG  = gParticlesDAL->GetPDG(TString(cProjectilePDG->GetSelectedEntry()->GetTitle()));
 //Int_t targetZ   = gParticlesDAL->GetPDG(TString(cTargetPDG->GetSelectedEntry()->GetTitle()));
 Int_t targetZ   = cTargetPDG->GetSelected() + 1;
 Int_t targetA   = (Int_t)nTargetA->GetNumber();
 ArgEnum argType   = (ArgEnum) cArgType->GetSelected();
 Double_t argValue  = nArgValue->GetNumber();
 UnitsEnum argUnits  = (UnitsEnum)cArgUnits->GetSelected();
 Int_t targetPDG   = gParticlesDAL->GetPDG(targetZ, targetA);


 if(gParticlesDAL->GetParticle(targetPDG) == 0)
 {
  cout << "Target Z="<< targetZ << ", A=" << targetA << " was not found, aborting..." <<  endl;
  return;
 }

 cout << "Projectile PDG = " << projectilePDG << endl;
 cout << "TargetPDG      = " << targetPDG << endl;
 cout << "Argument Type  = " << argType << endl;
 cout << "Argument Value = " << argValue << endl;
 cout << "Argument Units = " << argUnits << endl;

 fPublication = new G4TData(projectilePDG,targetPDG,true,"data",argType, argValue, argUnits, 6);

 // go to next state
 ChangeState(1);
}

//______________________________________________________________________________
void G4TPublicationGUI::Run2_Click()
{
 TString filename(tFilename->GetText());

 ifstream asciiFile;
 asciiFile.open(filename.Data());
 if (!asciiFile) {
  cout << "Error opening file "<< filename.Data() << ", please check the filename." <<  endl;
  return;
 }


 Int_t secondaryPDG  = gParticlesDAL->GetPDG(TString(cSecondaryPDG->GetSelectedEntry()->GetTitle()));
 CutEnum cutType   = (CutEnum) cCutType->GetSelected();
 UnitsEnum cutUnits  = (UnitsEnum) cCutUnits->GetSelected();
 Double_t cutValue  = nCutValue->GetNumber();
 Double_t cutDelta  = nCutDelta->GetNumber();
 FuncEnum functionType = (FuncEnum) cFunctionType->GetSelected();
 UnitsEnum functionUnits = (UnitsEnum) cFunctionUnits->GetSelected();
 ArgEnum secArgType  = (ArgEnum) cSecArgType->GetSelected();
 UnitsEnum secArgUnits = (UnitsEnum) cSecArgUnits->GetSelected();

 cout << "Secondary PDG  = " << secondaryPDG << endl;
 cout << "Cut Type       = " << cutType << endl;
 cout << "Cut Units      = " << cutUnits << endl;
 cout << "Cut Value      = " << cutValue << endl;
 cout << "Cut Delta      = " << cutDelta << endl;
 cout << "Function Type  = " << functionType << endl;
 cout << "Function Units = " << functionUnits << endl;
 cout << "Argument Type  = " << secArgType << endl;
 cout << "Argument Units = " << secArgUnits << endl;


 cout << "Loading ASCII file " << filename << "..." << endl;
 fPublication->AddItem(secondaryPDG,
   cutType,
   cutUnits,
   cutValue,
   cutDelta,
   functionType,
   functionUnits,
   secArgType,
   secArgUnits)->LoadFromASCII(filename);
}

//______________________________________________________________________________
void G4TPublicationGUI::Run3_Click()
{
 cout << "Saving the publication file... ";
 gTestDB->SaveData(fPublication);
 cout << "done!" << endl;
}

//______________________________________________________________________________
void G4TPublicationGUI::Run4_Click()
{
  gApplication->Terminate(0);
  //this->DestroyWindow();
}


//______________________________________________________________________________
void G4TPublicationGUI::ChangeState(Int_t state)
{
 if(state == 0)
 {

  bRun1->SetEnabled(true);
  bRun2->SetEnabled(false);
  bRun3->SetEnabled(false);

  cProjectilePDG->SetEnabled(true);
  cTargetPDG->SetEnabled(true);
  cArgType->SetEnabled(true);

  nTargetA->SetState(true);
  nArgValue->SetState(true);

  cArgUnits->SetEnabled(true);

  tFilename->SetEnabled(false);
  cSecondaryPDG->SetEnabled(false);
  cCutType->SetEnabled(false);
  cCutUnits->SetEnabled(false);
  nCutValue->SetState(false);
  nCutDelta->SetState(false);
  cFunctionType->SetEnabled(false);
  cFunctionUnits->SetEnabled(false);
  cSecArgType->SetEnabled(false);
  cSecArgUnits->SetEnabled(false);


 }
 else if(state == 1)
 {
  bRun1->SetEnabled(false);
  bRun2->SetEnabled(true);
  bRun3->SetEnabled(true);

  cProjectilePDG->SetEnabled(false);
  cTargetPDG->SetEnabled(false);
  cArgType->SetEnabled(false);
  nArgValue->SetState(false);
  nTargetA->SetState(false);
  cArgUnits->SetEnabled(false);

  tFilename->SetEnabled(true);
  cSecondaryPDG->SetEnabled(true);
  cCutType->SetEnabled(true);
  cCutUnits->SetEnabled(true);
  nCutValue->SetState(true);
  nCutDelta->SetState(true);
  cFunctionType->SetEnabled(true);
  cFunctionUnits->SetEnabled(true);
  cSecArgType->SetEnabled(true);
  cSecArgUnits->SetEnabled(true);
 }
}


//______________________________________________________________________________
G4TPublicationGUI::~G4TPublicationGUI()
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

