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
// CommonHeaders.h file
//
// Class description:
//
// Contains the headers for the ROOT classes, included into each header
// file by default.
//
// History:
// Created by Roman Atachiants, 18/08/2009
//
// --------------------------------------------------------------------

// Standard/Core Libraries
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

// Root Includes
#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TNamed.h>
#include <TString.h>

// Files and DataStructures
#include <TTree.h>
#include <TDirectory.h>
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TFolder.h>
#include <TKey.h>



// Functions
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

#include <TPostScript.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TCanvas.h>

#include <TMath.h>
#include <TLatex.h>
#include <TClassTable.h>
#include <TBrowser.h>


#include <TGraph.h>
#include <TGraphErrors.h>
#include <TView.h>
#include <TFrame.h>
#include <TCut.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TPolyMarker.h>
#include <TBenchmark.h>

#include <TCollection.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include <TRef.h>
#include <TBits.h>
#include <TDatabasePDG.h>
#include <TGeoManager.h>


//GUI
#include <TApplication.h>
#include "TGDockableFrame.h"
#include "TGMenu.h"
#include "TGMdiDecorFrame.h"
#include "TG3DLine.h"
#include "TGMdiFrame.h"
#include "TGMdiMainFrame.h"
#include "TGuiBldHintsButton.h"
#include "TRootBrowserLite.h"
#include "TGMdiMenu.h"
#include "TGListBox.h"
#include "TGNumberEntry.h"
#include "TGScrollBar.h"
#include "TGComboBox.h"
#include "TGuiBldHintsEditor.h"
#include "TGFrame.h"
#include "TGFileDialog.h"
#include "TGShutter.h"
#include "TGButtonGroup.h"
#include "TGCanvas.h"
#include "TGFSContainer.h"
#include "TGFontDialog.h"
#include "TGuiBldEditor.h"
#include "TGColorSelect.h"
#include "TGButton.h"
#include "TGFSComboBox.h"
#include "TGLabel.h"
#include "TGMsgBox.h"
#include "TRootGuiBuilder.h"
#include "TGTab.h"
#include "TGListView.h"
#include "TGSplitter.h"
#include "TGStatusBar.h"
#include "TGListTree.h"
#include "TGToolTip.h"
#include "TGToolBar.h"
#include "TGuiBldDragManager.h"



