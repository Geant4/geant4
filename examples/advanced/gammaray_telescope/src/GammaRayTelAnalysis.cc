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
// $Id: GammaRayTelAnalysis.cc,v 1.8 2001-12-07 13:20:20 pfeiffer Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayAnalysisManager  ------
//           by R.Giannitrapani, F.Longo & G.Santin (03 dic 2000)
//
// 07.12.2001 A.Pfeiffer
// - integrated Guy's addition of the ntuple
//
// 06.12.2001 A.Pfeiffer
// - updating to new design (singleton)
//
// 22.11.2001 G.Barrand
// - Adaptation to AIDA
//
// ************************************************************


#include "G4RunManager.hh" 

#include "GammaRayTelAnalysis.hh"
#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelAnalysisMessenger.hh"

#ifdef  G4ANALYSIS_USE
#include <AIDA/IAnalysisFactory.h>
#include <AIDA/ITreeFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IPlotterFactory.h>
#include <AIDA/IPlotter.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ITuple.h>
#endif

GammaRayTelAnalysis* GammaRayTelAnalysis::instance = 0;
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysis::GammaRayTelAnalysis(int argc,char** argv)
:GammaRayTelDetector(0)
,analysisFactory(0), tree(0), plotter(0), tuple(0)
,energy(0), hits(0), posXZ(0), posYZ(0)
,histo1DDraw("enable"),histo1DSave("enable"),histo2DDraw("enable")
,histo2DSave("enable"),histo2DMode("strip")
{
  G4RunManager* runManager = G4RunManager::GetRunManager();
  GammaRayTelDetector =
(GammaRayTelDetectorConstruction*)(runManager->GetUserDetectorConstruction());

#ifdef  G4ANALYSIS_USE
  // Define the messenger and the analysis system
  analysisMessenger = new GammaRayTelAnalysisMessenger(this);

  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory) {

    ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if(treeFactory) {
      // Tree in memory :
      //tree = treeFactory->create();
      // Create a "tree" associated to a ROOT "store" (in RECREATE mode).
      tree = treeFactory->create("GammaRayTel.root",false,false,"ROOT");
      if(tree) {
	IHistogramFactory* histoFactory	= 
	  analysisFactory->createHistogramFactory(*tree);  
	if(histoFactory) {
	  // Create histos :
	  int Nplane = GammaRayTelDetector->GetNbOfTKRLayers();
	  int Nstrip = GammaRayTelDetector->GetNbOfTKRStrips();
	  int Ntile = GammaRayTelDetector->GetNbOfTKRTiles();
	  float sizexy = GammaRayTelDetector->GetTKRSizeXY();
	  float sizez = GammaRayTelDetector->GetTKRSizeZ();
	  int N = Nstrip*Ntile;      

	  // 1D histogram that store the energy deposition of the
	  // particle in the last (number 0) TKR X-plane
	  energy = histoFactory->create1D("1","Energy deposition in the last X plane (keV)", 100, 50, 200);
	  
	  // 1D histogram that store the hits distribution along the TKR X-planes
	  hits = histoFactory->create1D("2","Hits distribution in the TKR X planes",Nplane, 0, Nplane-1);
	  
	  // 2D histogram that store the position (mm) of the hits (XZ projection)
	  if (histo2DMode == "strip")
	    posXZ = histoFactory->create2D("3","Tracker Hits XZ (strip,plane)", 
					   N, 0, N-1, 
					   2*Nplane, 0, Nplane-1);
	  else
	    posXZ = histoFactory->create2D("3","Tracker Hits XZ (x,z) in mm", 
					   sizexy/5, -sizexy/2, sizexy/2, 
					   sizez/5, -sizez/2, sizez/2);
	  
	  // 2D histogram that store the position (mm) of the hits (YZ projection)
	  if(histo2DMode=="strip")
	    posYZ = histoFactory->create2D("4","Tracker Hits YZ (strip,plane)", 
					   N, 0, N-1, 
					   2*Nplane, 0, Nplane-1);
	  else
	    posYZ = histoFactory->create2D("4","Tracker Hits YZ (y,z) in mm", 
					   sizexy/5, -sizexy/2, sizexy/2, 
					   sizez/5, -sizez/2, sizez/2);
	  
	  delete histoFactory;
	}
    

	// Get a tuple factory :
	ITupleFactory* tupleFactory = analysisFactory->createTupleFactory(*tree);
	if(tupleFactory) {
	  // Create a tuple :
	  tuple = tupleFactory->create("tuple","tuple",
				       "energy plane x y z");
	  delete tupleFactory;
	}
      }
      delete treeFactory; // Will not delete the ITree.
    }

    IPlotterFactory* plotterFactory = 
      analysisFactory->createPlotterFactory(argc,argv);
    if(plotterFactory) {
      plotter = plotterFactory->create();
      if(plotter) {
	plotter->show();
	plotter->setParameter("pageTitle","Gamma Ray Tel");
      }
      delete plotterFactory;
    }

  }


#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysis::~GammaRayTelAnalysis() {
  Finish();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnalysis::Init()
{
}                       

void GammaRayTelAnalysis::Finish()
{
#ifdef  G4ANALYSIS_USE
  delete plotter;
  delete analysisFactory; // Will delete tree and histos.
  delete analysisMessenger;
  analysisMessenger = 0;
#endif
}             

GammaRayTelAnalysis* GammaRayTelAnalysis::getInstance(int argc,char** argv)
{
  if (instance == 0) instance = new GammaRayTelAnalysis(argc,argv);
  return instance;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 2d histogram of the XZ positions
void GammaRayTelAnalysis::InsertPositionXZ(double x, double z)
{
#ifdef  G4ANALYSIS_USE
  if(posXZ) posXZ->fill(x, z);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 2d histogram of the YZ positions
void GammaRayTelAnalysis::InsertPositionYZ(double y, double z)
{
#ifdef  G4ANALYSIS_USE
  if(posYZ) posYZ->fill(y, z);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 1d histogram of the energy released in the last Si plane
void GammaRayTelAnalysis::InsertEnergy(double en)
{
#ifdef  G4ANALYSIS_USE
  if(energy) energy->fill(en);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 1d histogram of the hits distribution along the TKR planes
void GammaRayTelAnalysis::InsertHits(int nplane)
{
#ifdef  G4ANALYSIS_USE
  if(hits) hits->fill(nplane);
#endif
}

void GammaRayTelAnalysis::setNtuple(float E, float p, float x, float y, float z)
{
#ifdef  G4ANALYSIS_USE
  if(tuple) {
    tuple->fill(0,E);
    tuple->fill(1,p);
    tuple->fill(2,x);
    tuple->fill(3,y);
    tuple->fill(4,z);
    tuple->addRow();
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/
void GammaRayTelAnalysis::BeginOfRun(G4int n) 
{ 
#ifdef  G4ANALYSIS_USE

  if(energy) energy->reset();
  if(hits) hits->reset();
  if(posXZ) posXZ->reset();
  if(posYZ) posYZ->reset();

  //  Set the plotter ; set the number of regions and attach histograms
  // to plot for each region.
  //  It is done here, since then EndOfRun set regions
  // for paper output.
  if(plotter) {
    if((histo2DDraw == "enable") && (histo1DDraw == "enable")) {
      plotter->createRegions(2,2);
      plotter->plot(*posXZ);
      plotter->next();
      plotter->plot(*posYZ);
      plotter->next();
      plotter->plot(*energy);
      plotter->next();
      plotter->plot(*hits);
    } else if((histo1DDraw == "enable") && (histo2DDraw != "enable")) {
      plotter->createRegions(1,2);
      plotter->plot(*energy);
      plotter->next();
      plotter->plot(*hits);
    } else if((histo1DDraw != "enable") && (histo2DDraw == "enable")) {
      plotter->createRegions(1,2);
      plotter->plot(*posXZ);
      plotter->next();
      plotter->plot(*posYZ);
    } else { // Nothing to plot.
      plotter->createRegions(1,1);
    }
    plotter->refresh();
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* 
   This member is called at the end of each run 
*/
void GammaRayTelAnalysis::EndOfRun(G4int n) 
{
#ifdef  G4ANALYSIS_USE
  if(tree) tree->commit();

  if(plotter) {
    // We set one single region for the plotter
    plotter->createRegions(1,1);
  
    // We now print the histograms, each one in a separate file
    if(histo2DSave == "enable") {
      char name[15];
      sprintf(name,"posxz_%d.ps", n);
      plotter->plot(*posXZ);
      plotter->write(name,"PS");

      sprintf(name,"posyz_%d.ps", n);
      plotter->clearRegion();
      plotter->plot(*posYZ);
      plotter->write(name,"PS");
    }

    if(histo1DSave == "enable") {
      char name[15];
      sprintf(name,"energy_%d.ps", n);
      plotter->plot(*energy);
      plotter->write(name,"PS");
      
      sprintf(name,"hits_%d.ps", n);
      plotter->clearRegion();
      plotter->plot(*hits);
      plotter->write(name,"PS");
    }
  }
#endif
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* This member is called at the end of every event */
void GammaRayTelAnalysis::EndOfEvent(G4int flag) 
{
  // The plotter is updated only if there is some
  // hits in the event
  if(!flag) return;
#ifdef  G4ANALYSIS_USE
  if(plotter) plotter->refresh();
#endif
}

