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
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelAnalysisManager.cc                       *     
// * -------                                                            *
// *                                                                    *
// * Version:           0.2                                             *
// * Date:              30/11/00                                        *
// * Author:            A. Pfeiffer, G. Barrand, MG Pia, R Nartallo     *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
//
// CHANGE HISTORY
// --------------
//
// 07.12.2001 A.Pfeiffer
// - merged with Guy Barrand's AIDA 2.2 port
//
// 22.11.2001 G.Barrand
// - Adaptation to AIDA
//
// 30.11.2000 M.G. Pia, R. Nartallo
// - Simplification of code: removal of non-Lizard specific code
// - Inheritance directly from the base class G4VAnalysisManager instead
//   of the derived class G4AnalysisManager
//
// 15.11.2000 A. Pfeiffer
// - Adaptation to Lizard 
//
// 16.10.2000 G. Barrand
// - First implementation of XrayAnalysisManager
// - Provision of code for various AIDA and non-AIDA systems
//
// **********************************************************************


#include "g4std/fstream"

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

#ifdef G4ANALYSIS_USE
#include <AIDA/IAnalysisFactory.h>
#include <AIDA/ITreeFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IPlotterFactory.h>
#include <AIDA/IPlotter.h>
#endif

#include "XrayTelAnalysis.hh"

XrayTelAnalysis* XrayTelAnalysis::instance = 0;

XrayTelAnalysis::XrayTelAnalysis(int argc,char** argv)
:analysisFactory(0)
,tree(0)
,enteringEnergyHistogram(0)
,yzHistogram(0)
{
#ifdef G4ANALYSIS_USE
  // Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory) {

    ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if(treeFactory) {
      // tree = treeFactory->create(); // Tree in memory.
      tree = treeFactory->create("XrayTel.root",false,false,"ROOT");
      if(tree) {
	IHistogramFactory* hFactory = 
	  analysisFactory->createHistogramFactory(*tree);  
	if(hFactory) {
	  // Histogram creation (the below has a name and a label) :
	  enteringEnergyHistogram = 
	    hFactory->create1D("ekin.vop","Entering energy",100,0,0.5);
	  // Book the histogram for the 2D position. 
	  // Instead of using a scatter plot, just book enough bins ...
	  yzHistogram = 
	    hFactory->create2D("position.vop","YZ",200, -50, 50, 200, -50, 50);
	  delete hFactory;
	}
      }
      delete treeFactory; // Will not delete the ITree.
    }


  /* 
     The following lines set the plotter that is
     needed in this example for a multiple histograms
     visualization. Please see the README for more information 
  */
    IPlotterFactory* plotterFactory = 
      analysisFactory->createPlotterFactory(argc,argv);
    if(plotterFactory) {
      plotter = plotterFactory->create();
      if(plotter) {
        // Map the plotter on screen :
	plotter->show();
        // Set the page title :
	plotter->setParameter("pageTitle","XrayTel");
	// Have two plotting regions (one column, two rows).
	plotter->createRegions(1,2); 
        // Attach histograms to plotter regions :
	if(enteringEnergyHistogram) {
	  plotter->plot(*enteringEnergyHistogram);
	  plotter->next();
	}
	if(yzHistogram) {
	  plotter->plot(*yzHistogram);
	}
      }
      delete plotterFactory;
    }
  }

#endif
}

XrayTelAnalysis::~XrayTelAnalysis()
{
#ifdef G4ANALYSIS_USE
  delete plotter;
  delete enteringEnergyHistogram;
  delete yzHistogram;
  delete analysisFactory;
#endif
}

XrayTelAnalysis* XrayTelAnalysis::getInstance(G4int argc, char** argv)
{
  if (instance == 0) instance = new XrayTelAnalysis(argc, argv);
  return instance;
}


void XrayTelAnalysis::book(){
#ifdef G4ANALYSIS_USE
  if(enteringEnergyHistogram) enteringEnergyHistogram->reset();
  if(yzHistogram) yzHistogram->reset();
#endif
}

void XrayTelAnalysis::finish(){
#ifdef G4ANALYSIS_USE
  // the following things cannot be done in Run::EndOfRun()
  // only now plot the energy of the particles

  if(tree) tree->commit(); // Write histos in file.

  // Give control to the GUI so that someone 
  // can play with the plotter. The plotter
  // should be equipped with an "escape" (button) 
  // to return of the "interact" method.
  G4cout << "Click the escape button to continue..." << G4endl;
  if(plotter) plotter->interact();
#endif
}

void XrayTelAnalysis::analyseStepping(const G4Track& track, G4bool entering) {
#ifdef G4ANALYSIS_USE

  G4ThreeVector pos = track.GetPosition();

  if(enteringEnergyHistogram) {
    enteringEnergyHistogram->fill(track.GetKineticEnergy());
  }
  if(yzHistogram) {
    yzHistogram->fill(pos.y(), pos.z());
  }
  if(plotter) plotter->refresh();
#endif
} 










