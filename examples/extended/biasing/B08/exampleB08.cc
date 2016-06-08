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
// $Id: exampleB08.cc,v 1.1 2002/06/04 11:14:51 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB08
//
// --------------------------------------------------------------
// Comments
// This example demonstrates how to use the importance sampling 
// and weight window samling together with scoring in a parallel 
// geometry.
// The example for weight window sampling is somewhat artificial 
// since no other weight changing
// process than the importance sampling is used.
// Therefore high particle energys are choosen (for which importance 
// sampling would not be neccessary). That way
// secondaries with enough energy to create neutrons are 
// produced. Since the secondaries are not importance sampled
// they may create neutrons with a weight not coresponding to 
// the importance of a cell. In this case the weight window
// algorithm can be demonstrated: If the upper and lower limit
// is set to 1 it corrects the weight to exactly the inverse
// of the importance. 
// The init.mac file is used to set the upper and lower limit and the 
// importances (See init.mac).
//  
// --------------------------------------------------------------

#include "g4std/iostream"

#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "B08DetectorConstruction.hh"
#include "B08PhysicsList.hh"
#include "B08ImportanceDetectorConstruction.hh"
#include "B08PrimaryGeneratorAction.hh"

#include "G4VIStore.hh"
#include "G4Sigma.hh"
#include "G4WeightWindowAlgorithm.hh"
#include "G4PImportanceWWindowScoreSampler.hh"

#include "G4ParallelScoreSampler.hh"

#include "B08Scorer.hh"
#include "B08ScorePrinter.hh"

#include "B08MainMessenger.hh"

int main(int argc, char **argv)
{  

  
  G4long myseed = 345354;

  HepRandom::setTheSeed(myseed);

  G4RunManager *runManager = new G4RunManager;
  
  // create the   "tracking"    detector      -----------------
  runManager->SetUserInitialization(new B08DetectorConstruction);
  //  ---------------------------------------------------
  runManager->SetUserInitialization(new B08PhysicsList);
  runManager->SetUserAction(new B08PrimaryGeneratorAction);
  runManager->Initialize();

  ////////////////////////////////////////////////////////////////
  // create a messanger reading an initialication script to specify:
  // - the umber of events to be processed
  // - name of the outptfile
  // - upper limit and
  // - lower limit for the weight window sampling
  G4int numberOfEvent = 10;

  B08MainMessenger mess;

  // create the detector for the importances and scoring
  // this object also creates a messanger for the importance values
  B08ImportanceDetectorConstruction importancedetector;

  // read the init.mac file 
  G4UImanager *UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/control/execute init.mac"); 

  numberOfEvent = mess.GetNumberOfEvents();

  G4String resultsfilename = mess.GetResultsFileName();
  G4std::ofstream of(resultsfilename);

  //////////////////////////////////////////////////////////////////////



  /////////////////////////////////////////////////////////////////
  // setting up importance asmpling and scoring for neutrons and gammas

  // the IStore is filled during detector construction
  G4VIStore &aIstore = *importancedetector.GetIStore();
  
  // create scorers for neutrons and gammas
  B08Scorer nScorer;
  B08Scorer gScorer;

  // create the weight window algorithm 
  G4WeightWindowAlgorithm wwalg;
  wwalg.SetUpperLimit(mess.GetUpperLimit());
  wwalg.SetLowerLimit(mess.GetLowerLimit());
  
  // create the sampler for importance nad weight window sampling 
  // of neutron
  G4PImportanceWWindowScoreSampler nwwsampler(aIstore, 
					      nScorer, 
					      "neutron",
					      wwalg);
  nwwsampler.Initialize();

  // create the sampler for scoring gammas
  G4ParallelScoreSampler gsamp(*importancedetector.GetWorldVolume(),
			       "gamma",
			       gScorer);
  
  gsamp.Initialize();

  /////////////////////////////////////////////////////////////////
  // running

  runManager->BeamOn(numberOfEvent);

  /////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////
  // create ouput

  of << "Results from running of: " << numberOfEvent << " events" << G4endl; 
  of << "results for neutron:" << G4endl;
  B08ScorePrinter nsp(&aIstore);
  nsp.PrintHeader(&of);
  nsp.PrintTable(nScorer.GetMapPtkTallys(),  &of);
  of << G4endl;

  of << "results for gamma:" << G4endl;
  B08ScorePrinter sp;
  sp.PrintHeader(&of);
  sp.PrintTable(gScorer.GetMapPtkTallys(), &of);
  

  return 0;
}




