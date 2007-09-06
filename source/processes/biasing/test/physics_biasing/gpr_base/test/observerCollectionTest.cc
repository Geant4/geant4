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
// $Id: observerCollectionTest.cc,v 1.3 2007-09-06 22:07:04 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, May 2007. ObserverCollectionT demonstration.
//
#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"
#include "G4GPRObserverCollectionT.hh"
// Regular G4VProcess
struct VProcess : public G4VDiscreteProcess
{
  VProcess():G4VDiscreteProcess("test"){}
  G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*) {return 0;} // Required to implement

  void StartTracking(G4Track*)
  {
    G4cout<<"Execute VProcess::StartTracking"<<G4endl;
  }
};

// Simple class
struct OtherProcess
{
  void Method(G4Track*)
  {
    G4cout<<"Execute OtherProcess::Method"<<G4endl;
  }

  void operator()(G4Track*)
  {
    G4cout<<"Execute OtherProcess::operator"<<G4endl;
  }
};

// Function
void MyFunc(G4Track*)
{
  G4cout<<"Execute MyFunc"<<G4endl;
}

int main(int argc, char** argv) {

  // Create demo processes
  G4VProcess* vProcess = new VProcess();
  OtherProcess* otherProcess = new OtherProcess();

  // Create observer collection
  G4GPRObserverCollectionT<G4GPRTypeList_1(G4Track*)> observerCollection;

  // Register demo observers
  observerCollection.RegisterObserver("VProcess", vProcess, &G4VProcess::StartTracking);
  observerCollection.RegisterObserver("OtherProcess::Method", otherProcess, &OtherProcess::Method);
  observerCollection.RegisterObserver("OtherProcess::Operator", otherProcess, &OtherProcess::operator());
  observerCollection.RegisterObserver("MyFunc", &MyFunc);

  G4Track* dummyTrk(0);

  // Notify observers for no particular reason - just a demo
  observerCollection(dummyTrk);

  // Cleanup
  delete vProcess;
  delete otherProcess;

  return 0;
}
