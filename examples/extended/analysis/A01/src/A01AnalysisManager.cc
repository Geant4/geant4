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
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 tutorial 1                              *
// *                                                                    *
// * MODULE:            A01AnalysisManager.cc                           *
// * -------                                                            *
// *                                                                    *
// * Version:           0.1                                             *
// * Date:              January 28 2002                                 *
// * Author:            T.Johnson                                       *
// * Organisation:      SLAC                                            *
// *                                                                    *
// **********************************************************************
//
// CHANGE HISTORY
// --------------
//
// Nov 4 2002 -- Upgrade to AIDA 3.0
// **********************************************************************
#ifdef G4ANALYSIS_USE

#include <fstream>

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

#include <AIDA/AIDA.h>

#include "A01AnalysisManager.hh"

A01AnalysisManager* A01AnalysisManager::instance = 0;

A01AnalysisManager::A01AnalysisManager()
:analysisFactory(0), hFactory(0), tFactory(0), plotter(0)
{
  // Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory)
  {
    ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    tree = treeFactory->create("A01.aida","xml",false,true,"compress=yes");
    hFactory = analysisFactory->createHistogramFactory(*tree);
    tFactory = analysisFactory->createTupleFactory(*tree);
    IPlotterFactory* pf = analysisFactory->createPlotterFactory(0,0);
    if (pf) {
      plotter = pf->create("Plotter");
      delete pf;
    }
    delete treeFactory; // Will not delete the ITree.
  }
}

A01AnalysisManager::~A01AnalysisManager()
{
  if (analysisFactory)
  {
    if (!tree->commit()) G4cout << "Commit failed: no AIDA file produced!" << G4endl;
    delete tree;
    delete tFactory;
    delete hFactory;
    delete plotter;
    G4cout << "Warning: In case of working with JAS-AIDA, Geant4 will NOT exit unless you close the JAS-AIDA window." << G4endl;
    delete analysisFactory;
  }
}
IHistogramFactory* A01AnalysisManager::getHistogramFactory()
{
  return hFactory;
}
ITupleFactory* A01AnalysisManager::getTupleFactory()
{
  return tFactory;
}
IPlotter* A01AnalysisManager::getPlotter()
{
  return plotter;
}

A01AnalysisManager* A01AnalysisManager::getInstance()
{
  if (instance == 0) instance = new A01AnalysisManager();
  return instance;
}

void A01AnalysisManager::dispose()
{
  if (instance != 0)
  {
    delete instance;
    instance = 0;
  }
}

#endif // G4ANALYSIS_USE

