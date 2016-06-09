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
// $Id: TiaraCellScorer.cc,v 1.3 2003/06/18 16:40:28 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// TiaraCellScorer.cc
//
// ----------------------------------------------------------------------

#include <cmath>
#include "TiaraCellScorer.hh"
#include "G4GeometryCell.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4CellScoreValues.hh"
#include "Randomize.hh"

#ifdef G4ANALYSIS_USE
TiaraCellScorer::TiaraCellScorer(AIDA::IHistogramFactory *hf,
				 const G4String &histBaseName,
				 const std::vector<double> 
				 &binEdgesScinti,
				 const std::vector<double> 
				 &binEdgesBonner,
				 const TiaraTally &tally) :
  fBaseName(histBaseName),
  fEnergyHisto(hf->createHistogram1D(histBaseName, 
				     histBaseName + ": flux vs. energy", 
				     binEdgesScinti)),
  fEnergyFluxHisto(hf->createHistogram1D(histBaseName + "Eflux", 
					 histBaseName + ": energy flux vs. energy", 
					 binEdgesScinti)),
  
  fEnergyHistoBonner(hf->createHistogram1D(histBaseName + "Bonner", 
					   histBaseName + ": flux vs. energy", 
					   binEdgesBonner)),
  fEnergyFluxHistoBonner(hf->createHistogram1D(histBaseName + "BonnerEflux", 
					       histBaseName + ": energy flux vs. energy", 
					       binEdgesBonner)),
  
  fTally(tally)
{}
#else
TiaraCellScorer::TiaraCellScorer(const G4String &histBaseName,
				 const TiaraTally &tally) :
  fBaseName(histBaseName),
  fTally(tally)
{}
#endif


TiaraCellScorer::~TiaraCellScorer()
{}

void TiaraCellScorer::ScoreAnExitingStep(const G4Step &aStep,
				       const G4GeometryCell &pre_gCell){
  fG4CellScorer.ScoreAnExitingStep(aStep, pre_gCell);
  FillHisto(aStep);
}

void TiaraCellScorer::ScoreAnEnteringStep(const G4Step &aStep,
					const G4GeometryCell &post_gCell){
  fG4CellScorer.ScoreAnEnteringStep(aStep, post_gCell);
  return;
}

void TiaraCellScorer::ScoreAnInVolumeStep(const G4Step &aStep,
					const G4GeometryCell &post_gCell){
  fG4CellScorer.ScoreAnInVolumeStep(aStep, post_gCell);
  FillHisto(aStep);
}

void TiaraCellScorer::FillHisto(const G4Step &aStep){

  G4StepPoint *p = 0;
  p = aStep.GetPreStepPoint();

  G4double sl = aStep.GetStepLength() / cm;  
  G4double w = p->GetWeight();
  G4double slw = sl * w;

  G4double e(p->GetKineticEnergy());

#ifdef G4ANALYSIS_USE
  fEnergyHisto->fill(e, slw);
  fEnergyFluxHisto->fill(e, slw*e);
  fEnergyHistoBonner->fill(e/eV, slw);
  fEnergyFluxHistoBonner->fill(e/eV, slw*e/eV);
#else
  G4cout << fBaseName << ": Energy: " << e << ", sl*w: " << slw << "\n";
#endif

  fTally.fill(e, slw);

}

void TiaraCellScorer::EndOfEventAction() {
#ifndef G4ANALYSIS_USE
  G4cout << fBaseName << ": end of event" << G4endl;
#endif
  fTally.EndOfEventAction();
}

const TiaraTally &TiaraCellScorer::GetTally() const {
  return fTally;
}

