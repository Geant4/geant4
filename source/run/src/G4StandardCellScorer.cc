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
// $Id: G4StandardCellScorer.cc,v 1.1 2002-07-11 16:19:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4StandardCellScorer.cc
//
// ----------------------------------------------------------------------

#include "G4StandardCellScorer.hh"

G4StandardCellScorer::G4StandardCellScorer()
{}
G4StandardCellScorer::~G4StandardCellScorer()
{}

void G4StandardCellScorer::SetSLRawValues(G4SLRawValues slr){
  G4double sl = slr.fStepLength;
  G4double slw = sl * slr.fPre_Weight;
  G4double slwe = slw * slr.fPre_Energy;

  G4double v = slr.fPre_Velocity;
  if (v==0.) v = 10e-9;

  fSCScoreValues.fSumSL += sl;
  fSCScoreValues.fSumSLW += slw;
  fSCScoreValues.fSumSLWE +=  slwe;
  fSCScoreValues.fSumSLWE_v += slwe / v;
  
}
void G4StandardCellScorer::TrackEnters(){
  fSCScoreValues.fSumTracksEntering++;
}
void G4StandardCellScorer::NewTrackPopedUp(){
  fSCScoreValues.fSumPopulation++;
}

void G4StandardCellScorer::SetCollisionWeight(G4double weight){
  fSCScoreValues.fSumCollisions++;
  fSCScoreValues.fSumCollisionsWeight+=weight;
}


G4StandardCellScoreValues G4StandardCellScorer::GetStandardCellScoreValues() const {
  if (fSCScoreValues.fSumSLW > 0.) {
    fSCScoreValues.fNumberWeightedEnergy = 
      fSCScoreValues.fSumSLWE_v / fSCScoreValues.fSumSLW;
    fSCScoreValues.fFluxWeightedEnergy = 
      fSCScoreValues.fSumSLWE / fSCScoreValues.fSumSLW;
    fSCScoreValues.fAverageTrackWeight = 
      fSCScoreValues.fSumSLW / fSCScoreValues.fSumSL;
  }
  return fSCScoreValues;
}

void G4StandardCellScorer::SetImportnace(G4double importance){
  fSCScoreValues.fImportance = importance;
};

G4std::ostream& operator<<(G4std::ostream &out, 
                           const G4StandardCellScorer &ps) {
  G4StandardCellScoreValues scores =  ps.GetStandardCellScoreValues();
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
