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
// $Id$
//
/// \file analysis/A01/src/A01AnalysisManager.cc
/// \brief Implementation of the A01AnalysisManager class

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

A01AnalysisManager* A01AnalysisManager::fInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

A01AnalysisManager::A01AnalysisManager()
:fHisfTupleFactory(0), fTupleFactory(0), fPlotter(0), fTree(0)
{
    // Hooking an AIDA compliant analysis system.
    fAnalysisFactory = AIDA_createAnalysisFactory();
    if(fAnalysisFactory)
    {
        ITreeFactory* treeFactory = fAnalysisFactory->createTreeFactory();
        fTree = treeFactory->create("A01.aida","xml",false,true,"compress=yes");
        fHisfTupleFactory = fAnalysisFactory->createHistogramFactory(*fTree);
        fTupleFactory = fAnalysisFactory->createTupleFactory(*fTree);
        IPlotterFactory* pf = fAnalysisFactory->createPlotterFactory(0,0);
        if (pf) {
            fPlotter = pf->create("Plotter");
            delete pf;
        }
        delete treeFactory; // Will not delete the ITree.
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

A01AnalysisManager::~A01AnalysisManager()
{
    if (fAnalysisFactory)
    {
        if (!fTree->commit()) G4cout << "Commit failed: no AIDA file produced!" << G4endl;
        delete fTree;
        delete fTupleFactory;
        delete fHisfTupleFactory;
        delete fPlotter;
        G4cout << "Warning: In case of working with JAS-AIDA:" << G4endl;
        G4cout << "Geant4 will NOT exit unless you close the JAS-AIDA window." << G4endl;
        delete fAnalysisFactory;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IHistogramFactory* A01AnalysisManager::GetHistogramFactory()
{
    return fHisfTupleFactory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ITupleFactory* A01AnalysisManager::GetTupleFactory()
{
    return fTupleFactory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IPlotter* A01AnalysisManager::GetPlotter()
{
    return fPlotter;
}

A01AnalysisManager* A01AnalysisManager::GetInstance()
{
    if (fInstance == 0) fInstance = new A01AnalysisManager();
    return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void A01AnalysisManager::Dispose()
{
    if (fInstance != 0)
    {
        delete fInstance;
        fInstance = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif // G4ANALYSIS_USE
