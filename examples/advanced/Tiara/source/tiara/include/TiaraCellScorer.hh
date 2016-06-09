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
// $Id: TiaraCellScorer.hh,v 1.3 2003/06/18 16:40:23 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// ----------------------------------------------------------------------
// Class TiaraCellScorer
//
// Class description:
//
// This class is an example how to build a customized cell scorer
// derived from G4VCellScorer that also uses the G4CellScorer.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef TiaraCellScorer_hh
#define TiaraCellScorer_hh TiaraCellScorer_hh


#include "G4CellScoreComposer.hh"
#include "G4CellScorer.hh"
#include "G4VCellScorer.hh"

#include <vector>

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

#include "TiaraTally.hh"


class TiaraCellScorer : public G4VCellScorer{
public:  
#ifdef G4ANALYSIS_USE
  TiaraCellScorer(AIDA::IHistogramFactory *hf,
		  const G4String &histBaseName,
		  const std::vector<double>  & binEdgesScinti,
		  const std::vector<double>  & binEdgesBonner,
		  const TiaraTally &tally);
#else
  TiaraCellScorer(const G4String &histBaseName,
		  const TiaraTally &tally);
#endif

  virtual ~TiaraCellScorer();
  virtual void ScoreAnExitingStep(const G4Step &aStep, 
				  const G4GeometryCell &gCell);
  virtual void ScoreAnEnteringStep(const G4Step &aStep, 
				   const G4GeometryCell &gCell);
  virtual void ScoreAnInVolumeStep(const G4Step &aStep, 
				   const G4GeometryCell &gCell);

  const G4CellScorer &GetG4CellScorer() const {
    return fG4CellScorer;
  }
  const TiaraTally &GetTally() const ;

  void EndOfEventAction();

private:
  void FillHisto(const G4Step &aStep);

  G4CellScorer fG4CellScorer;
  G4String fBaseName;
#ifdef G4ANALYSIS_USE
  AIDA::IHistogram1D* fEnergyHisto;
  AIDA::IHistogram1D* fEnergyFluxHisto;
  AIDA::IHistogram1D* fEnergyHistoBonner;
  AIDA::IHistogram1D* fEnergyFluxHistoBonner;
#endif
  TiaraTally fTally;
};



#endif
