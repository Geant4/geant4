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
// $Id: G4CellScoreComposer.hh,v 1.1 2003/10/03 10:06:49 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// ----------------------------------------------------------------------
// Class G4CellScoreComposer
//
// Class description:
// This class will be created for every cell standard
// scoring should be applied. It does the actual scoring.
// GetStandardCellScoreValues() delivers the struct
// G4CellScoreValues  does calculations based on the 
// sums of scores and delivers the results in
// G4CellScoreValues.
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4CellScoreComposer_hh
#define G4CellScoreComposer_hh G4CellScoreComposer_hh

#include "G4CellScoreValues.hh"

class G4Step;

class G4CellScoreComposer {
public: // with description

  G4CellScoreComposer();

  ~G4CellScoreComposer();

  void EstimatorCalculation(const G4Step &step);
    // get values for estimators based on
    // track length

  void TrackEnters();
    // to be called if a track enters the cell 

  void NewTrackPopedUp();
    // to be caled if the cell popultion is increased 

  void SetCollisionWeight(G4double weight);
    // to be called for every collision
    // in the cell with the weight of the colliding particle 

  void SetImportnace(G4double importance);
    // informs G4CellScoreComposer about  the importance of the cell 

  const G4CellScoreValues &GetStandardCellScoreValues() const;
    // return scores in G4CellScoreValues

private:
  mutable G4CellScoreValues fSCScoreValues;
};

std::ostream& operator<<(std::ostream &out, 
                           const G4CellScoreComposer &ps);



#endif

