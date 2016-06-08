// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
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
#include "G4SteppingManager.hh"
#include "G4VAnalysisSystem.hh"

#ifdef G4ANALYSIS_USE
#include <IHistogramFactory.h>
#endif

#ifdef G4ANALYSIS_USE_LIZARD
#include <IHistogram1D.h>
#include <IHistogram2D.h>
#include "G4LizardSystem.hh"
#endif

#include "XrayTelAnalysisManager.hh"

#ifndef G4ANALYSIS_USE
XrayTelAnalysisManager::XrayTelAnalysisManager(G4String const& aSystem)
{
}
XrayTelAnalysisManager::~XrayTelAnalysisManager()
{
}
#endif

#ifdef G4ANALYSIS_USE
XrayTelAnalysisManager::XrayTelAnalysisManager(G4String const& aSystem)
  :
  enteringEnergyHistogram(0),
  yzHistogram(0)
{
  // The Lizard analysis system is provided here as default
  // other analysis systems will be implemented at a later stage
  analysisSystem = new G4LizardSystem; 
  IHistogramFactory* hFactory = analysisSystem->GetHistogramFactory();

  if(hFactory) {
    // Histogram creation :
    enteringEnergyHistogram = hFactory->create1D("Entering energy",100,0,0.5);
    // Book the histogram for the 2D position. 
    // Instead of using a scatter plot, just book enough bins ...
    yzHistogram = hFactory->create2D("YZ",200, -50, 50, 200, -50, 50);
  }
}

XrayTelAnalysisManager::~XrayTelAnalysisManager()
{
  // delete hFactory;
  delete enteringEnergyHistogram;
  delete yzHistogram;
  delete analysisSystem;
}

G4bool XrayTelAnalysisManager::RegisterAnalysisSystem(
G4VAnalysisSystem* )
{
  return true;
}

IHistogramFactory* XrayTelAnalysisManager::GetHistogramFactory(
 const G4String& aSystem)
{
  return hFactory;
}

void XrayTelAnalysisManager::Store(IHistogram* aHistogram,const G4String& aSID)
{
  analysisSystem->Store(aHistogram,aSID);
}
void XrayTelAnalysisManager::Plot(IHistogram* aHistogram)
{
  analysisSystem->Plot(aHistogram);
}

void XrayTelAnalysisManager::BeginOfRun(){
  if(enteringEnergyHistogram) 
    enteringEnergyHistogram->reset();
  if(yzHistogram) 
    yzHistogram->reset();
}
void XrayTelAnalysisManager::EndOfRun(){
  // the following things cannot be done in Run::EndOfRun()
  // only now plot the energy of the particles
  if(enteringEnergyHistogram) {
    Plot(enteringEnergyHistogram);
  }
  // and store the histograms
  if(enteringEnergyHistogram) {
    Store(enteringEnergyHistogram,"ekin.vop");
  }
  if(yzHistogram) {
    Store(yzHistogram,"position.vop");
  }
}

void XrayTelAnalysisManager::Step(const G4SteppingManager* aSteppingManager) {
  if(!aSteppingManager) return;

  G4Track* track = aSteppingManager->GetTrack();
  G4ThreeVector pos = track->GetPosition();

  if(enteringEnergyHistogram) {
    enteringEnergyHistogram->fill(track->GetKineticEnergy());
  }
  if(yzHistogram) {
    yzHistogram->fill(pos.y(), pos.z());
    Plot(yzHistogram);
  }
} 

#endif








