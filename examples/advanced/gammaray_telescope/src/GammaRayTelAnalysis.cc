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
#ifdef G4ANALYSIS_USE
//
// $Id$
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayAnalysisManager  ------
//           by R.Giannitrapani, F.Longo & G.Santin (03 dic 2000)
//
// 29.05.2003 F.Longo 
// - anaphe 5.0.5 compliant
//
// 18.06.2002 R.Giannitrapani, F.Longo & G.Santin
// - new release for Anaphe 4.0.3
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
#include <fstream>

#include "G4RunManager.hh" 

#include "GammaRayTelAnalysis.hh"
#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelAnalysisMessenger.hh"

GammaRayTelAnalysis* GammaRayTelAnalysis::instance = 0;

//-------------------------------------------------------------------------------- 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysis::GammaRayTelAnalysis()
  :GammaRayTelDetector(0),analysisFactory(0), tree(0)//, plotter(0), 
  ,tuple(0)
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

  analysisFactory = AIDA_createAnalysisFactory(); // create the Analysis Factory
  if(analysisFactory) {
     AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
  // create Tree Factory

    if(treeFactory) {
      // Tree in memory :
      // Create a "tree" associated to an xml file 

      //      tree = treeFactory->create("gammaraytel.hbook", "hbook", false, false);
      // (hbook implementation)
      tree = treeFactory->create("gammaraytel.aida","xml",false,true,"");
      if(tree) {
	// Get a tuple factory :
	ITupleFactory* tupleFactory = analysisFactory->createTupleFactory(*tree);
	if(tupleFactory) {
	  // Create a tuple :
	    tuple = tupleFactory->create("1","1", "float energy, plane, x, y, z");
	    assert(tuple);
	    
	    delete tupleFactory;
	}
	
	IHistogramFactory* histoFactory = analysisFactory->createHistogramFactory(*tree);

	if(histoFactory) {
	  // Create histos :

	  int Nplane = GammaRayTelDetector->GetNbOfTKRLayers();
	  int Nstrip = GammaRayTelDetector->GetNbOfTKRStrips();
	  int Ntile = GammaRayTelDetector->GetNbOfTKRTiles();
	  double sizexy = GammaRayTelDetector->GetTKRSizeXY();
	  double sizez = GammaRayTelDetector->GetTKRSizeZ();
	  int N = Nstrip*Ntile;      

	  // 1D histogram that store the energy deposition of the
	  // particle in the last (number 0) TKR X-plane

	  energy = histoFactory->createHistogram1D("10","Edep in the last X plane (keV)", 100, 50, 200);
	  
	  // 1D histogram that store the hits distribution along the TKR X-planes

	  hits = histoFactory->createHistogram1D("20","Hits dist in TKR X planes",Nplane, 0, Nplane-1);
	  // 2D histogram that store the position (mm) of the hits (XZ projection)

	  if (histo2DMode == "strip"){
	    posXZ = histoFactory->createHistogram2D("30","Tracker Hits XZ (strip,plane)", 
					   N, 0, N-1, 
					   2*Nplane, 0, Nplane-1); 
            }
	  else
	    {
	    posXZ = histoFactory->createHistogram2D("30","Tracker Hits XZ (x,z) in mm", 
					   int(sizexy/5), -sizexy/2, sizexy/2, 
						    int(sizez/5), -sizez/2, sizez/2);
	  }
	  // 2D histogram that store the position (mm) of the hits (YZ projection)

	  if(histo2DMode=="strip")
	    posYZ = histoFactory->createHistogram2D("40","Tracker Hits YZ (strip,plane)", 
						    N, 0, N-1, 
						    2*Nplane, 0, Nplane-1);
	  else
	    posYZ = histoFactory->createHistogram2D("40","Tracker Hits YZ (y,z) in mm", 
						    int(sizexy/5), -sizexy/2, sizexy/2, 
					   int(sizez/5), -sizez/2, sizez/2);
	  
	  delete histoFactory;
	}
    
      }
     }
      delete treeFactory; // Will not delete the ITree. 
     
    }
    
 //  IPlotterFactory* plotterFactory = analysisFactory->createPlotterFactory(0,0);    
//   if(plotterFactory) {
//     plotter  = plotterFactory->create();
//     if(plotter) {
//       plotter->show();
//       plotter->setParameter("pageTitle","Gamma Ray Tel");
        
//     }
//     delete plotterFactory;
//   }

#endif
}


GammaRayTelAnalysis::~GammaRayTelAnalysis() {
  Finish();
}


void GammaRayTelAnalysis::Init()
{
}                       

void GammaRayTelAnalysis::Finish()
{
#ifdef  G4ANALYSIS_USE
  delete tree;
  //delete plotter;
  //  delete analysisFactory; // Will delete tree and histos.
  delete analysisMessenger;
  
  analysisMessenger = 0;
#endif
}             

GammaRayTelAnalysis* GammaRayTelAnalysis::getInstance()
{
  if (instance == 0) instance = new GammaRayTelAnalysis();
  return instance;
}

// This function fill the 2d histogram of the XZ positions
void GammaRayTelAnalysis::InsertPositionXZ(double x, double z)
{
#ifdef  G4ANALYSIS_USE
  if(posXZ) posXZ->fill(x, z);
#endif
}

// This function fill the 2d histogram of the YZ positions
void GammaRayTelAnalysis::InsertPositionYZ(double y, double z)
{
#ifdef  G4ANALYSIS_USE
  if(posYZ) posYZ->fill(y, z);
#endif
}

// This function fill the 1d histogram of the energy released in the last Si plane
void GammaRayTelAnalysis::InsertEnergy(double en)
{
#ifdef  G4ANALYSIS_USE
  if(energy) energy->fill(en);
#endif
}

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
    tuple->fill(tuple->findColumn("energy"),E);
    tuple->fill(tuple->findColumn("plane"),p);
    tuple->fill(tuple->findColumn("x"),x);
    tuple->fill(tuple->findColumn("y"),y);
    tuple->fill(tuple->findColumn("z"),z);
    tuple->addRow();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/

//void GammaRayTelAnalysis::BeginOfRun(G4int n) // to be reintroduced
void GammaRayTelAnalysis::BeginOfRun() 
{ 
#ifdef  G4ANALYSIS_USE

  if(energy) energy->reset();
  if(hits) hits->reset();
  if(posXZ) posXZ->reset();
  if(posYZ) posYZ->reset();

#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* 
   This member is called at the end of each run 
*/
void GammaRayTelAnalysis::EndOfRun(G4int) 
{
#ifdef  G4ANALYSIS_USE
  if(tree) {
    tree->commit();
    tree->close();
  }

//   if(plotter) {

//     // We set one single region for the plotter
//     // We now print the histograms, each one in a separate file

//     if(histo2DSave == "enable") {
//       char name[15];
//       plotter->createRegions(1,1);
//       sprintf(name,"posxz_%d.ps", n);
//       plotter->currentRegion().plot(*posXZ);
//       plotter->refresh();
//       //      plotter->write(name,"ps"); // temporary unavailable
      
//       plotter->createRegions(1,1);
//       sprintf(name,"posyz_%d.ps", n);
//       plotter->currentRegion().plot(*posYZ);
//       plotter->next().plot(*posYZ);
//       plotter->refresh();
//       // plotter->write(name,"ps"); // temporary unavailable
//     }

//     if(histo1DSave == "enable") {
//       plotter->createRegions(1,1);
//       char name[15];
//       sprintf(name,"energy_%d.ps", n);
//       plotter->currentRegion().plot(*energy);
//       plotter->refresh();
//       //      plotter->write(name,"ps"); // temporary unavailable
//       plotter->createRegions(1,1);   
//       sprintf(name,"hits_%d.ps", n);
//       plotter->currentRegion().plot(*hits);
//       plotter->refresh();
//       //      plotter->write(name,"ps"); // temporary unavailable
//       plotter->createRegions(1,2);
//       plotter->currentRegion().plot(*energy);
//       plotter->next().plot(*hits);
//       plotter->refresh();
  //  }
#endif
}

/* This member is called at the end of every event */
void GammaRayTelAnalysis::EndOfEvent(G4int flag) 
{
  // The plotter is updated only if there is some
  // hits in the event
  if(!flag) return;
#ifdef  G4ANALYSIS_USE
  //  Set the plotter ; set the number of regions and attach histograms
  // to plot for each region.
  //  It is done here, since then EndOfRun set regions
  // for paper output.
  // if(plotter) {
//     if((histo2DDraw == "enable") && (histo1DDraw == "enable")) {
//       plotter->createRegions(1,2);
//       //plotter->currentRegion().plot(*posXZ); //temporary unavailable
//       plotter->currentRegion().plot(*hits);
//       //      plotter->next().plot(*posYZ); //temporary unavailable
//       plotter->next().plot(*energy);
//       //plotter->next().plot(*energy);
//       //      plotter->currentRegion().plot(*hits);
//       //plotter->next().plot(*hits);
//     } else if((histo1DDraw == "enable") && (histo2DDraw != "enable")) {
//       plotter->createRegions(1,2);
//       plotter->currentRegion().plot(*energy);
//       plotter->next().plot(*hits);
//     } else if((histo1DDraw != "enable") && (histo2DDraw == "enable")) {
//       /*      plotter->createRegions(1,2);
//       plotter->currentRegion().plot(*posXZ);
//       plotter->next().plot(*posYZ);*/
//       G4cout << "Temporary Unavailable " << G4endl;
//     } else { // Nothing to plot.
//       plotter->createRegions(1,1);
//     }
//     plotter->refresh();
//   }

#endif
}
#endif







