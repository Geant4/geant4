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
// $Id: G4CellScoreComposer.cc,v 1.2 2004/07/01 09:19:12 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02-patch-01 $
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

  G4StepPoint *p = 0;
  p = aStep.GetPreStepPoint();
  if (!p) {
    G4Exception("G4CellScoreComposer::EstimatorCalculation: no pointer to pre PreStepPoint!");
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

