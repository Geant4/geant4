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
// Author:  A. Pfeiffer (Andreas.Pfeiffer@cern.ch) 
//         (copied from his UserAnalyser class)
//
// History:
// -----------
//  7 Nov 2001   MGP        Implementation
// 03 dic 2000   AnalysisManager by R.Giannitrapani, F.Longo & G.Santin
// 18 Nov 2001   G.Santin GammaRayTel analysis management modified
//               according to the new design
//
// -------------------------------------------------------------------
#ifdef  G4ANALYSIS_USE

#include "GammaRayTelAnalysis.hh"
#include "GammaRayTelAnalysisMessenger.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include "G4ThreeVector.hh"
#include "GammaRayTelDetectorConstruction.hh"

GammaRayTelAnalysis* GammaRayTelAnalysis::instance = 0;

GammaRayTelAnalysis::GammaRayTelAnalysis() :
  histoManager(0), vectorFactory(0),plotter(0),
  histo1DDraw("enable"),histo1DSave("enable"),histo2DDraw("enable"),
  histo2DSave("enable"),histo2DMode("strip"),
  analysisMessenger(0)
{
  G4RunManager* runManager = G4RunManager::GetRunManager();
  GammaRayTelDetector = (GammaRayTelDetectorConstruction*)(runManager->GetUserDetectorConstruction());

  analysisMessenger = new GammaRayTelAnalysisMessenger(this);

  histoManager = createIHistoManager();
  ntFactory = Lizard::createNTupleFactory();
  vectorFactory = createIVectorFactory();
  plotter = createIPlotter();
}

GammaRayTelAnalysis::~GammaRayTelAnalysis()
{ 
  delete ntFactory;
  ntFactory = 0;

  delete histoManager;
  histoManager = 0;

  delete vectorFactory;
  vectorFactory = 0;

  delete plotter;
  plotter = 0;

  delete analysisMessenger;
  analysisMessenger = 0;
}

GammaRayTelAnalysis* GammaRayTelAnalysis::getInstance()
{
  if (instance == 0) instance = new GammaRayTelAnalysis;
  return instance;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 2d histogram of the XZ positions
void GammaRayTelAnalysis::InsertPositionXZ(double x, double z)
{
  IHistogram2D* posXZ = histoManager->retrieveHisto2D("3");
  posXZ->fill(x, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 2d histogram of the YZ positions
void GammaRayTelAnalysis::InsertPositionYZ(double y, double z)
{
  IHistogram2D* posYZ = histoManager->retrieveHisto2D("4");
  posYZ->fill(y, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 1d histogram of the energy released in the last Si plane
void GammaRayTelAnalysis::InsertEnergy(double en)
{
  IHistogram1D* energy = histoManager->retrieveHisto1D("1");
  energy->fill(en);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 1d histogram of the hits distribution along the TKR planes
void GammaRayTelAnalysis::InsertHits(int nplane)
{
  IHistogram1D* hits = histoManager->retrieveHisto1D("2");
  hits->fill(nplane);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/
void GammaRayTelAnalysis::BeginOfRun() 
{ 
  float sizexy, sizez;
  int Nstrip, Nplane, Ntile, N;

  // Relevant data from the detector to set the histograms dimensions
  Nplane = GammaRayTelDetector->GetNbOfTKRLayers();
  Nstrip = GammaRayTelDetector->GetNbOfTKRStrips();
  Ntile = GammaRayTelDetector->GetNbOfTKRTiles();
  sizexy = GammaRayTelDetector->GetTKRSizeXY();
  sizez = GammaRayTelDetector->GetTKRSizeZ();
  N = Nstrip*Ntile;      

  if (histoManager) {
    histoManager->selectStore("GammaRayTel.his");
    // 1D histogram that store the energy deposition of the
    // particle in the last (number 0) TKR X-plane

    // Check if deleteHisto is enough, maybe we need to delete the object
    histoManager->deleteHisto("1");
    histoManager->create1D("1","Energy deposition in the last X plane (keV)", 100, 50., 200.);
      
    // 1D histogram that store the hits distribution along the TKR X-planes
    histoManager->deleteHisto("2");
    histoManager->create1D("2","Hits distribution in the TKR X planes",
				  Nplane, 0, Nplane-1);
    
    // 2D histogram that store the position (mm) of the hits (XZ projection)
    histoManager->deleteHisto("3");
    if (histo2DMode == "strip")
      histoManager->create2D("3","Tracker Hits XZ (strip,plane)", 
				     N, 0, N-1, 
				     2*Nplane, 0, Nplane-1);
    else
      histoManager->create2D("3","Tracker Hits XZ (x,z) in mm", 
				     sizexy/5, -sizexy/2, sizexy/2, 
				     sizez/5, -sizez/2, sizez/2);
    
    // 2D histogram that store the position (mm) of the hits (YZ projection)
    histoManager->deleteHisto("4");
    if(histo2DMode=="strip")
      histoManager->create2D("4","Tracker Hits YZ (strip,plane)", 
				     N, 0, N-1, 
				     2*Nplane, 0, Nplane-1);
    else
      histoManager->create2D("4","Tracker Hits YZ (y,z) in mm", 
				     sizexy/5, -sizexy/2, sizexy/2, 
				     sizez/5, -sizez/2, sizez/2);
  }
  
  IHistogram2D* posXZ = histoManager->retrieveHisto2D("3");
  if (posXZ) 
    posXZ->reset();
  IHistogram2D* posYZ = histoManager->retrieveHisto2D("4");
  if(posYZ) 
    posYZ->reset();
  IHistogram1D* energy = histoManager->retrieveHisto1D("1");
  if(energy)
    energy->reset();
  IHistogram1D* hits = histoManager->retrieveHisto1D("2");
  if(hits)
    hits->reset();
  
  // We divide the plotter in the right nuber of zone depending
  // on which histograms the user want to draw
  if((histo2DDraw == "enable") && (histo1DDraw == "enable"))
    plotter->zone(2,2,0,0);
  else if((histo1DDraw == "enable") || (histo2DDraw == "enable"))
    plotter->zone(1,2,0,0);
  else
    plotter->zone(1,1,0,0);

  // Book ntuples
  ntuple = ntFactory->createC("GammaRayTel.his::1");
  
  //  Add and bind the attributes to the ntuple
  if ( !( ntuple->addAndBind( "energy", ntEnergy) &&
  	  ntuple->addAndBind( "plane" , ntPlane) &&
  	  ntuple->addAndBind( "x"     , ntX   ) &&
  	  ntuple->addAndBind( "y"     , ntY   ) &&
  	  ntuple->addAndBind( "z"     , ntZ   ) ) ) {
    delete ntuple;
    G4Exception(" GammaRayTelAnalysis::BeginOfRun - Could not addAndBind ntuple");
  }  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* 
   This member is called at the end of each run 
*/
void GammaRayTelAnalysis::EndOfRun(G4int n) 
{
  // This variable contains the names of the PS files
  char name[15];
  // We define some vectors
  IVector* vxz  = 0;
  IVector* vyz  = 0;
  IVector* ve   = 0;
  IVector* vhit = 0;


  // Temporary we set one single zone for the plotter
  plotter->zone(1,1,0,0);
  
  // We now print the histograms, each one in a separate file
  if(histo2DSave == "enable") {
    IHistogram2D* posXZ = histoManager->retrieveHisto2D("3");
    IHistogram2D* posYZ = histoManager->retrieveHisto2D("4");
    plot2D(posXZ);
    sprintf(name,"posxz_%d.ps", n);
    plotter->psPrint(name);
    
    plot2D(posYZ);
    sprintf(name,"posyz_%d.ps", n);
    plotter->psPrint(name);
  }

  if(histo1DSave == "enable") {
    IHistogram1D* energy = histoManager->retrieveHisto1D("1");
    IHistogram1D* hits   = histoManager->retrieveHisto1D("2");
    
    plot1D(energy);
    sprintf(name,"energy_%d.ps", n);
    plotter->psPrint(name);
    
    plot1D(hits);
    sprintf(name,"hits_%d.ps", n);
    plotter->psPrint(name);  
  }
  
  delete vxz;
  delete vyz;
  delete ve;
  delete vhit;

  // Store histograms

  // Because of a Lizard feature, ntuples must be deleted at this stage, 
  // not in the destructor (otherwise the ntuples are not stored)

  // In version 3.6.4 of Anaphe, no licence, storing other histograms after 
  // having stored a 2D one generates an error message "ERROR in HLNEXT"
  // Histograms are stored correctly, in spite of the error message
  // The same error message is issued if attempting to store 1D or 2D
  // histograms after deleting a ntuple

  histoManager->store("1");
  histoManager->store("3");
  histoManager->store("2");
  histoManager->store("4");

  delete ntuple;
  ntuple = 0;
  G4cout << "Deleted ntuple" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* This member is called at the end of every event */
void GammaRayTelAnalysis::EndOfEvent(G4int flag) 
{
  // The histograms are updated only if there is some
  // hits in the event
  if(flag) Plot();
}


void GammaRayTelAnalysis::Plot()
{
  // In a normal case we use the following line to use the standard plot
  // analysisSystem->Plot(histo);
  
  // We define some vectors
  IVector* vxz = 0;
  IVector* vyz = 0;
  IVector* ve  = 0;
  IVector* vhit= 0;

  // We fill them with the histograms and
  // draw them 
  if(histo2DDraw == "enable")
    {
      IHistogram2D* posXZ = histoManager->retrieveHisto2D("3");
      IHistogram2D* posYZ = histoManager->retrieveHisto2D("4");
      plot2D(posXZ);
      plot2D(posYZ);
      plotter->refresh();
    }

  if(histo1DDraw == "enable")
    {
      IHistogram1D* energy = histoManager->retrieveHisto1D("1");
      IHistogram1D* hits   = histoManager->retrieveHisto1D("2");
      plot1D(energy);
      plot1D(hits);
      plotter->refresh();
    }

  delete vxz;
  delete vyz;
  delete ve;
  delete vhit;
}

void GammaRayTelAnalysis::plot1D(IHistogram1D* histo)
{
  IVector* v = vectorFactory->from1D(histo);
  plotter->plot(v,0);
  plotter->refresh();
  delete v;	
}

void GammaRayTelAnalysis::plot2D(IHistogram2D* histo)
{
  IVector* v = vectorFactory->from2D(histo);
  plotter->plot(v,0);
  plotter->refresh();
  delete v;	
}

#endif
