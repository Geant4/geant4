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
// $Id: G4StandardCellScorer.hh,v 1.1 2002-07-11 16:19:43 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4StandardCellScorer
//
// Class description:
// This class will be created for every cell standard
// scoring should be applied. It does the actual scoring.
// GetStandardCellScoreValues() delivers the struct
// G4StandardCellScoreValues  does calculations based on the 
// sums of scores and delivers the results in
// G4StandardCellScoreValues.
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4StandardCellScorer_hh
#define G4StandardCellScorer_hh G4StandardCellScorer_hh

#include "G4SLRawValues.hh"
#include "G4StandardCellScoreValues.hh"

class G4StandardCellScorer {
public:
  G4StandardCellScorer();
  ~G4StandardCellScorer();
  void SetSLRawValues(G4SLRawValues slr);
    // count in the raw scoe values for estimators based on
    // track length
  void TrackEnters();
    // to be called if a track enters the cell 
  void NewTrackPopedUp();
    // to be caled if the cell popultion is increased 
  void SetCollisionWeight(G4double weight);
    // to be called for every collision
    // in the cell with the weight of the colliding particle 
  void SetImportnace(G4double importance);
    // informs G4StandardCellScorer about  the importance of the cell 
  G4StandardCellScoreValues GetStandardCellScoreValues() const;
    // return scores in G4StandardCellScoreValues

private:
  mutable G4StandardCellScoreValues fSCScoreValues;
};

G4std::ostream& operator<<(G4std::ostream &out, 
                           const G4StandardCellScorer &ps);



#endif

