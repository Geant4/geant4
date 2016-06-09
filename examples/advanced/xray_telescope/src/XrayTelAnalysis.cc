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
// $Id: XrayTelAnalysis.cc,v 1.10 2003/08/13 13:06:02 santin Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Author:  A. Pfeiffer (Andreas.Pfeiffer@cern.ch) 
//         (copied from his UserAnalyser class)
//
// History:
// -----------
//  8 Nov 2002   GS         Migration to AIDA 3
//  7 Nov 2001   MGP        Implementation
//
// -------------------------------------------------------------------

#include "XrayTelAnalysis.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include "G4SteppingManager.hh"
#include "G4ThreeVector.hh"


XrayTelAnalysis* XrayTelAnalysis::instance = 0;

XrayTelAnalysis::XrayTelAnalysis()
#ifdef G4ANALYSIS_USE
  : analysisFactory(0)
  , tree(0)
  , histoFactory(0)
  , tupleFactory(0)
#ifdef G4ANALYSIS_USE_PLOTTER
  , plotterFactory(0)
  , plotter(0)
#endif
#endif
{
#ifdef G4ANALYSIS_USE
  histFileName = "xraytel";
  histFileType = "hbook";
#endif

  asciiFileName="xraytel.out";
  std::ofstream asciiFile(asciiFileName, std::ios::app);
  if(asciiFile.is_open()) {
    asciiFile << "Energy (keV)  x (mm)    y (mm)    z (mm)" << G4endl << G4endl;
  }
}

XrayTelAnalysis::~XrayTelAnalysis()
{ 
#ifdef G4ANALYSIS_USE

#ifdef G4ANALYSIS_USE_PLOTTER
  if (plotterFactory) delete plotterFactory;
  plotterFactory = 0;
#endif

  if (tupleFactory) delete tupleFactory;
  tupleFactory = 0;

  if (histoFactory) delete histoFactory;
  histoFactory = 0;

  if (tree) delete tree;
  tree = 0;

  if (analysisFactory) delete analysisFactory;
  analysisFactory = 0;
#endif
}

XrayTelAnalysis* XrayTelAnalysis::getInstance()
{
  if (instance == 0) instance = new XrayTelAnalysis;
  return instance;
}


void XrayTelAnalysis::book()
{
#ifdef G4ANALYSIS_USE
  //build up  the  factories
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory) {
    //parameters for the TreeFactory
    G4bool fileExists = false;
    G4bool readOnly   = false;
    AIDA::ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    if(treeFactory) {
      G4String histFileNameComplete; 
      histFileNameComplete = histFileName+".hbook";
      tree = treeFactory->create(histFileNameComplete, "hbook", readOnly, fileExists);
      G4cout << " Histogramfile: " << histFileNameComplete << G4endl;
      
      if (tree) {
	G4cout << "Tree store : " << tree->storeName() << G4endl;
	G4cout << "Booked Hbook File " << G4endl;

	//HistoFactory and TupleFactory depend on theTree
	histoFactory = analysisFactory->createHistogramFactory ( *tree );
	tupleFactory = analysisFactory->createTupleFactory     ( *tree );

	// Book histograms
	histoFactory->createHistogram1D("1","Energy, all /keV",  100,0.,100.);
	histoFactory->createHistogram2D("2","y-z, all /mm", 100,-500.,500.,100,-500.,500.);
	histoFactory->createHistogram1D("3","Energy, entering detector /keV", 500,0.,500.);
	histoFactory->createHistogram2D("4","y-z, entering detector /mm", 200,-50.,50.,200,-50.,50.);

	// Book ntuples
	AIDA::ITuple* ntuple10 = tupleFactory->create( "10", "Track ntuple", 
						       "double energy,x,y,z,dirx,diry,dirz" );
	assert(ntuple10);

      }
      delete treeFactory;
    }
  }
#endif

}

void XrayTelAnalysis::finish()
{
#ifdef G4ANALYSIS_USE
  if (tree) {
    // Committing the transaction with the tree
    std::cout << "Committing..." << std::endl;
    // write all histograms to file
    tree->commit();

    std::cout << "Closing the tree..." << std::endl;

    // close (will again commit)
    tree->close();
  }

  // extra delete as objects are created in book() method rather than during
  // initialisation of class

#ifdef G4ANALYSIS_USE_PLOTTER
  if (plotterFactory)  delete plotterFactory;
#endif

  if (tupleFactory)    delete tupleFactory;
  if (histoFactory)    delete histoFactory;
  if (tree)            delete tree;
  if (analysisFactory) delete analysisFactory;
#endif
}

void XrayTelAnalysis::analyseStepping(const G4Track& track, G4bool entering)
{
  eKin = track.GetKineticEnergy()/keV;
  G4ThreeVector pos = track.GetPosition()/mm;
  y = pos.y();
  z = pos.z();
  G4ThreeVector dir = track.GetMomentumDirection();
  dirX = dir.x();
  dirY = dir.y();
  dirZ = dir.z();

#ifdef G4ANALYSIS_USE
  // Fill histograms, all tracks
  AIDA::IHistogram1D* h1 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("1") );
  h1->fill(eKin);  // fill(x,y,weight)
  AIDA::IHistogram2D* h2 = dynamic_cast<AIDA::IHistogram2D *> ( tree->find("2") );
  h2->fill(y,z);

  // Fill histograms and ntuple, tracks entering the detector
  if (entering) {
    // Fill and plot histograms
    AIDA::IHistogram1D* h3 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("3") );
    h3->fill(eKin);
    AIDA::IHistogram2D* h4 = dynamic_cast<AIDA::IHistogram2D *> ( tree->find("4") );
    h4->fill(y,z);
#ifdef G4ANALYSIS_USE_PLOTTER
    plotAll();
#endif
  }

  // Fill ntuple
  if (entering) {
    
    AIDA::ITuple * ntuple = dynamic_cast<AIDA::ITuple *> ( tree->find("10") );
    if (ntuple) {
      // Fill the secondaries ntuple
      ntuple->fill( ntuple->findColumn( "energy" ), (G4double) eKin );
      ntuple->fill( ntuple->findColumn( "x"      ), (G4double) x    );
      ntuple->fill( ntuple->findColumn( "y"      ), (G4double) y    );
      ntuple->fill( ntuple->findColumn( "z"      ), (G4double) z    );
      ntuple->fill( ntuple->findColumn( "dirx"   ), (G4double) dirX );
      ntuple->fill( ntuple->findColumn( "diry"   ), (G4double) dirY );
      ntuple->fill( ntuple->findColumn( "dirz"   ), (G4double) dirZ );

      ntuple->addRow(); // check for returning true ...
    } else {
      G4cout << "Ntuple not found" << G4endl;
    }
  }

#endif

  // Write to file
  if (entering) {
    std::ofstream asciiFile(asciiFileName, std::ios::app);
    if(asciiFile.is_open()) {
      asciiFile << std::setiosflags(std::ios::fixed)
		<< std::setprecision(3)
		<< std::setiosflags(std::ios::right)
		<< std::setw(10);
      asciiFile << eKin;
      asciiFile << std::setiosflags(std::ios::fixed)
		<< std::setprecision(3)
		<< std::setiosflags(std::ios::right)
		<< std::setw(10);
      asciiFile << x;
      asciiFile << std::setiosflags(std::ios::fixed)
		<< std::setprecision(3)
		<< std::setiosflags(std::ios::right)
		<< std::setw(10);
      asciiFile << y;
      asciiFile << std::setiosflags(std::ios::fixed)
		<< std::setprecision(3)
		<< std::setiosflags(std::ios::right)
		<< std::setw(10);
      asciiFile << z
		<< G4endl;
      asciiFile.close();
    }
  }

}

#ifdef G4ANALYSIS_USE_PLOTTER
void XrayTelAnalysis::plotAll()
{
  if (!plotter) {
    AIDA::IPlotterFactory* plotterFactory = 
      analysisFactory->createPlotterFactory();
    if(plotterFactory) {
      G4cout << "Creating the Plotter" << G4endl;
      plotter = plotterFactory->create();
      if(plotter) {
	// Map the plotter on screen :
	G4cout << "Showing the Plotter on screen" << G4endl;
	plotter->show();
      } else {
	G4cout << "XrayTelAnalysis::plotAll: WARNING: Plotter not created" << G4endl;
      }
      delete plotterFactory;
    } else {
      G4cout << "XrayTelAnalysis::plotAll: WARNING: Plotter Factory not created" << G4endl;
    }
  }

  if (plotter) {
    plotter->createRegions(2,1,0);
    AIDA::IHistogram1D* hp = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("3") );
    AIDA::IHistogram1D& h  = *hp;  
    (plotter->currentRegion()).plot(h);
    plotter->refresh();
    plotter->setCurrentRegionNumber(1);
    AIDA::IHistogram1D* hp2 = dynamic_cast<AIDA::IHistogram1D *> ( tree->find("1") );
    AIDA::IHistogram1D& h2  = *hp2;  
    (plotter->currentRegion()).plot(h2);
    plotter->refresh();
  }
}
#endif

