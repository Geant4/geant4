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
// * MODULE:            TargetAnalysisManager.cc                        *     
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
//
// **********************************************************************


#include <fstream>

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

#include <AIDA/AIDA.h>


#include "TargetAnalysisManager.hh"

TargetAnalysisManager* TargetAnalysisManager::instance = 0;

TargetAnalysisManager::TargetAnalysisManager()
:analysisFactory(0), hFactory(0), tFactory(0)
{
  // Hooking an AIDA compliant analysis system.
  analysisFactory = AIDA_createAnalysisFactory();
  if(analysisFactory) 
  {
    ITreeFactory* treeFactory = analysisFactory->createTreeFactory();
    tree = treeFactory->create("LEP.aida",false,false,"type=xml;compress=yes"); // Tree in memory.
	hFactory = analysisFactory->createHistogramFactory(*tree);  
	tFactory = analysisFactory->createTupleFactory(*tree);
	cFactory = analysisFactory->createCloudFactory(*tree);
    delete treeFactory; // Will not delete the ITree.
  }
}

TargetAnalysisManager::~TargetAnalysisManager()
{
  if (analysisFactory)
  {
    if (!tree->commit()) G4cout << "Commit failed: no AIDA file produced!" << G4endl;
    delete tree;
	delete cFactory;
	delete tFactory;
	delete hFactory;
    delete analysisFactory;
  }
}

IHistogramFactory* TargetAnalysisManager::getHistogramFactory()
{
  return hFactory;
}
ITupleFactory* TargetAnalysisManager::getTupleFactory()
{
  return tFactory;
}
ICloudFactory* TargetAnalysisManager::getCloudFactory()
{
  return cFactory;
}
IPlotter* TargetAnalysisManager::createPlotter()
{
  if (analysisFactory)
  {
    IPlotterFactory* pf = analysisFactory->createPlotterFactory(0,0);
	if (pf) return pf->create("Plotter");
  }
  return 0;
}

TargetAnalysisManager* TargetAnalysisManager::getInstance()
{
  if (instance == 0) instance = new TargetAnalysisManager();
  return instance;
}

void TargetAnalysisManager::dispose()
{
  if (instance != 0) 
  {
    delete instance;
	instance = 0; 
  }
}











