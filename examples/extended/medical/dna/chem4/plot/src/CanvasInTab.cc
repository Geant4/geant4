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
#include "CanvasInTab.hh"

#include "TGFileDialog.h"
#include "TGTab.h"

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGButton.h>
#include <TGClient.h>
#include <TRandom.h>
#include <TRootEmbeddedCanvas.h>
#include <iostream>
#include <map>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const char* SaveFileDialog()
{
  // Prompt for file to be saved. Depending on navigation in
  // dialog the current working directory can be changed.
  // The returned file name is always with respect to the
  // current directory.

  const char* gSaveAsTypes[] = {"Macro files",
                                "*.C",
                                "ROOT files",
                                "*.root",
                                "PostScript",
                                "*.ps",
                                "Encapsulated PostScript",
                                "*.eps",
                                "PDF files",
                                "*.pdf",
                                "Gif files",
                                "*.gif",
                                "PNG files",
                                "*.png",
                                "All files",
                                "*",
                                0,
                                0};

  static TGFileInfo fi;
  fi.fFileTypes = gSaveAsTypes;

  new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDSave, &fi);

  return fi.fFilename;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CanvasInTab::CanvasInTab(const TGWindow* p, UInt_t w, UInt_t h) : TGMainFrame(p, w, h)
{
  fpTab = new TGTab(this, 200, 200);

  fHintPlots = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 2);

  AddFrame(fpTab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1));

  fpTab->Resize();

  //-----
  TGHorizontalFrame* hframe = new TGHorizontalFrame(this, 200, 40);

  TGTextButton* save = new TGTextButton(hframe, "&Save as ...");
  save->Connect("Clicked()", "CanvasInTab", this, "SaveCanvas()");
  hframe->AddFrame(save, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  TGTextButton* exit = new TGTextButton(hframe, "&Exit ", "gApplication->Terminate()");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

  //-----
  // Sets window name and shows the main frame

  SetWindowName("PlotG");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CanvasInTab::~CanvasInTab()
{
  Cleanup();
  //  if(fpTab)
  //  {
  ////    fpTab->Cleanup();
  //    delete fpTab;
  //  }

  //  if(fHintPlots)
  //    delete fHintPlots;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

size_t CanvasInTab::AddCanvas(const char* name)
{
  size_t output = fEcanvas.size();
  auto compositeFrame = fpTab->AddTab(name);
  TRootEmbeddedCanvas* embeddedCanvas = new TRootEmbeddedCanvas(name, compositeFrame, 500, 300);
  embeddedCanvas->SetAutoFit();
  fEcanvas.push_back(embeddedCanvas);
  compositeFrame->AddFrame(embeddedCanvas, fHintPlots);
  embeddedCanvas->SetContainer(compositeFrame);
  fpTab->Resize();
  fpTab->MapSubwindows();
  //  fpTab->MapWindow();
  Resize();
  return output;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TCanvas* CanvasInTab::GetCanvas(int i)
{
  return fEcanvas[i]->GetCanvas();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CanvasInTab::SaveCanvas()
{
  if (fpTab->GetNumberOfTabs() == 0) return;

  const char* name = SaveFileDialog();

  if (name == 0 || strlen(name) == 0) return;

  int current = fpTab->GetCurrent();
  TCanvas* canvas = fEcanvas[current]->GetCanvas();
  canvas->SaveAs(name);
}
