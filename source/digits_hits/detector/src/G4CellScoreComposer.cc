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
// $Id: G4CellScoreComposer.cc 94771 2015-12-09 09:44:05Z gcosmo $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4CellSCoreComposer.cc
//
// ----------------------------------------------------------------------

#include "G4CellScoreComposer.hh"
#include "G4Step.hh"

G4CellScoreComposer::G4CellScoreComposer():
  fSCScoreValues()
{}

G4CellScoreComposer::~G4CellScoreComposer()
{}

void G4CellScoreComposer::EstimatorCalculation(const G4Step &aStep){

  G4StepPoint *p = aStep.GetPreStepPoint();
  if (!p) {
    G4Exception("G4CellScoreComposer::EstimatorCalculation","Det0191",FatalException," no pointer to pre PreStepPoint!");
  }
  G4double sl = aStep.GetStepLength();
  G4double slw = sl * p->GetWeight();
  G4double slwe = slw * p->GetKineticEnergy();
  
  G4double v = p->GetVelocity();
  if (!(v>0.)) {
    v = 10e-9;
  }

  fSCScoreValues.fSumSL += sl;
  fSCScoreValues.fSumSLW += slw;
  fSCScoreValues.fSumSLW_v += slw / v;
  fSCScoreValues.fSumSLWE +=  slwe;
  fSCScoreValues.fSumSLWE_v += slwe / v;
 
}
void G4CellScoreComposer::TrackEnters(){
  fSCScoreValues.fSumTracksEntering++;
}
void G4CellScoreComposer::NewTrackPopedUp(){
  fSCScoreValues.fSumPopulation++;
}

void G4CellScoreComposer::SetCollisionWeight(G4double weight){
  fSCScoreValues.fSumCollisions++;
  fSCScoreValues.fSumCollisionsWeight+=weight;
}


const G4CellScoreValues &G4CellScoreComposer::
GetStandardCellScoreValues() const {
  if (fSCScoreValues.fSumSLW > 0.) {
    //divide by SumSLW or SumSLW_v ?
    fSCScoreValues.fNumberWeightedEnergy = 
      fSCScoreValues.fSumSLWE_v / fSCScoreValues.fSumSLW_v; 

    fSCScoreValues.fFluxWeightedEnergy = 
      fSCScoreValues.fSumSLWE / fSCScoreValues.fSumSLW;

    fSCScoreValues.fAverageTrackWeight = 
      fSCScoreValues.fSumSLW / fSCScoreValues.fSumSL;
  }
  return fSCScoreValues;
}

void G4CellScoreComposer::SetImportnace(G4double importance){
  fSCScoreValues.fImportance = importance;
}

std::ostream& operator<<(std::ostream &out, 
                           const G4CellScoreComposer &ps) {
  G4CellScoreValues scores =  ps.GetStandardCellScoreValues();
  out << "Tracks entering: " << scores.fSumTracksEntering << G4endl;
  out << "Population:      " << scores.fSumPopulation << G4endl;
  out << "Collisions:      " << scores.fSumCollisions << G4endl;
  out << "Collisions*Wgt:  " << scores.fSumCollisionsWeight << G4endl;
  out << "NumWGTedEnergy:  " << scores.fNumberWeightedEnergy << G4endl;
  out << "FluxWGTedEnergy: " << scores.fFluxWeightedEnergy << G4endl;
  out << "Aver.TrackWGT*I: " << scores.fAverageTrackWeight*
    scores.fImportance << G4endl;
  return out;
}

