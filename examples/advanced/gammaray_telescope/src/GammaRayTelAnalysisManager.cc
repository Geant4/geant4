// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelAnalysisManager.cc,v 1.1 2000-12-06 16:53:13 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayAnalysisManager  ------
//           by R.Giannitrapani, F.Longo & G.Santin (03 dic 2000)
//
// ************************************************************
#ifdef  G4ANALYSIS_USE

#include <stdlib.h>
#include "g4std/fstream"
#include "GammaRayTelAnalysisManager.hh"
#include "G4VAnalysisSystem.hh"
#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelAnalysisMessenger.hh"
#include <IHistogramFactory.h>
#include <IHistogram1D.h>
#include <IHistogram2D.h>
#include <IPlotter.h>
#include <IVector.h>
#include <IVectorFactory.h>

#include "G4LizardSystem.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysisManager::GammaRayTelAnalysisManager(GammaRayTelDetectorConstruction* GammaRayTelDC): 
  GammaRayTelDetector(GammaRayTelDC), 
  posXZ(0), posYZ(0), energy(0), hits(0), histoFactory(0), pl(0),
  histo1DDraw("enable"),histo1DSave("enable"),histo2DDraw("enable"),
  histo2DSave("enable"),histo2DMode("strip")
{
  // Define the messenger and the analysis system
  analysisMessenger = new GammaRayTelAnalysisMessenger(this);
  analysisSystem = new G4LizardSystem;

  histoFactory = analysisSystem->GetHistogramFactory();  

  /* 
     The following lines set the plotter and the vectorfactory that
     are needed in this example for a multiple histograms
     visualization. Please see the README for more information 
  */
  fVectorFactory = createIVectorFactory();          
  pl = createIPlotter();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysisManager::~GammaRayTelAnalysisManager() {
  delete posXZ;
  delete posYZ;
  delete energy;
  delete hits;
  delete analysisSystem;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool GammaRayTelAnalysisManager::RegisterAnalysisSystem(G4VAnalysisSystem*)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

IHistogramFactory* GammaRayTelAnalysisManager::GetHistogramFactory(const G4String& aSystem)
{
  return histoFactory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnalysisManager::Store(IHistogram* histo, const G4String& ID)
{
  analysisSystem->Store(histo, ID);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* 
   Since the Lizard interface and analysis classes in G4 are still
   experimental, they lack for now the possibility to show directly
   more than one histograms at the same time. In order to show 2 or 4 histo
   in a single view in our example, we decided to override this implementation 
   and use directly the IPlotter interface for our needs. This will change
   in future releases.

   Please note that for visualization purpouses the histograms result
   stretched along the x axis; when they are saved in a PostScript
   file the proportions are the right ones.
*/
void GammaRayTelAnalysisManager::Plot(IHistogram* histo = 0)
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
      vxz = fVectorFactory->from2D(dynamic_cast<IHistogram2D*>(posXZ));  
      vyz = fVectorFactory->from2D(dynamic_cast<IHistogram2D*>(posYZ));  
      pl->plot(vxz);
      pl->plot(vyz);
      pl->refresh();
    }

  if(histo1DDraw == "enable")
    {
      ve = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(energy));  
      vhit = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(hits));  
      pl->plot(ve);
      pl->plot(vhit);
      pl->refresh();
    }

  delete vxz;
  delete vyz;
  delete ve;
  delete vhit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 2d histogram of the XZ positions
void GammaRayTelAnalysisManager::InsertPositionXZ(double x, double z)
{
  posXZ->fill(x, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 2d histogram of the YZ positions
void GammaRayTelAnalysisManager::InsertPositionYZ(double y, double z)
{
  posYZ->fill(y, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 1d histogram of the energy released in the last Si plane
void GammaRayTelAnalysisManager::InsertEnergy(double en)
{
  energy->fill(en);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function fill the 1d histogram of the hits distribution along the TKR planes
void GammaRayTelAnalysisManager::InsertHits(int nplane)
{
  hits->fill(nplane);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/
void GammaRayTelAnalysisManager::BeginOfRun() 
{ 
  float sizexy, sizez;
  int nplane;
  int Nstrip, Nplane, Ntile, N;

  // Relevant data from the detector to set the histograms dimensions
  Nplane = GammaRayTelDetector->GetNbOfTKRLayers();
  Nstrip = GammaRayTelDetector->GetNbOfTKRStrips();
  Ntile = GammaRayTelDetector->GetNbOfTKRTiles();
  sizexy = GammaRayTelDetector->GetTKRSizeXY();
  sizez = GammaRayTelDetector->GetTKRSizeZ();
  N = Nstrip*Ntile;      

  if (histoFactory) 
    {
      // 1D histogram that store the energy deposition of the
      // particle in the last (number 0) TKR X-plane
      histoFactory->destroy(energy);
      energy = histoFactory->create1D("Energy deposition in the last X plane (keV)", 100, 50, 200);
      
      // 1D histogram that store the hits distribution along the TKR X-planes
      histoFactory->destroy(hits);
      hits = histoFactory->create1D("Hits distribution in the TKR X planes",
				     Nplane, 0, Nplane-1);

      // 2D histogram that store the position (mm) of the hits (XZ projection)
      histoFactory->destroy(posXZ);
      if (histo2DMode == "strip")
	posXZ = histoFactory->create2D("Tracker Hits XZ (strip,plane)", 
				       N, 0, N-1, 
				       2*Nplane, 0, Nplane-1);
      else
	posXZ = histoFactory->create2D("Tracker Hits XZ (x,z) in mm", 
				       sizexy/5, -sizexy/2, sizexy/2, 
				       sizez/5, -sizez/2, sizez/2);

      // 2D histogram that store the position (mm) of the hits (YZ projection)
      histoFactory->destroy(posYZ);
      if(histo2DMode=="strip")
	posYZ = histoFactory->create2D("Tracker Hits YZ (strip,plane)", 
				       N, 0, N-1, 
				       2*Nplane, 0, Nplane-1);
      else
	posYZ = histoFactory->create2D("Tracker Hits YZ (y,z) in mm", 
				       sizexy/5, -sizexy/2, sizexy/2, 
				       sizez/5, -sizez/2, sizez/2);
    }

  if(posXZ) 
    posXZ->reset();
  if(posYZ) 
    posYZ->reset();
  if(energy)
    energy->reset();
  if(hits)
    hits->reset();

  // We divide the plotter in the right nuber of zone depending
  // on which histograms the user want to draw
  if((histo2DDraw == "enable") && (histo1DDraw == "enable"))
    pl->zone(2,2);
  else if((histo1DDraw == "enable") || (histo2DDraw == "enable"))
    pl->zone(1,2);
  else
    pl->zone(1,1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* 
   This member is called at the end of each run 
*/
void GammaRayTelAnalysisManager::EndOfRun(G4int n) 
{
  // This variable contains the names of the PS files
  char name[15];
  // We define some vectors
  IVector* vxz  = 0;
  IVector* vyz  = 0;
  IVector* ve   = 0;
  IVector* vhit = 0;


  // Temporary we set one single zone for the plotter
  pl->zone(1,1);
  
  // We now print the histograms, each one in a separate file
  if(histo2DSave == "enable")
    {
      vxz = fVectorFactory->from2D(dynamic_cast<IHistogram2D*>(posXZ));  
      vyz = fVectorFactory->from2D(dynamic_cast<IHistogram2D*>(posYZ));  

      sprintf(name,"posxz_%d.ps", n);
      pl->plot(vxz);
      pl->psPrint(name);

      sprintf(name,"posyz_%d.ps", n);
      pl->plot(vyz);
      pl->psPrint(name);
    }

  if(histo1DSave == "enable")
    {
      ve = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(energy));  
      vhit = fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(hits));  

      sprintf(name,"energy_%d.ps", n);
      pl->plot(ve);
      pl->psPrint(name);
  
      sprintf(name,"hits_%d.ps", n);
      pl->plot(vhit);
      pl->psPrint(name);

    }
  
  delete vxz;
  delete vyz;
  delete ve;
  delete vhit;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* This member is called at the end of every event */
void GammaRayTelAnalysisManager::EndOfEvent(G4int flag) 
{
  // The histograms are updated only if there is some
  // hits in the event
  if(flag) Plot();
}

#endif











