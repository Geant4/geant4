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
// $Id: G4PIScorer.hh,v 1.5 2002-04-10 13:14:16 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PIScorer
//
// Class description:
//
// A scorer delegating the scoring to G4PScorer. It performs
// a check for consistency between the importance value
// assigned  to a post G4PTouchableKey and the weight of a track.
// The relation should be importance value * weight = 1.
// CAUTION: This holds only in well defined situations:
// A track starts in a "cell" with importance 1 and has a weight of 1.
// The particle can not be produced by other not importance sampled
// particles.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4PIScorer_hh
#define G4PIScorer_hh G4PIScorer_hh

#include "g4std/iostream"

#include "G4VPScorer.hh"
#include "G4PMapPtkTallys.hh"
#include "G4PScorer.hh"

class G4Step;
class G4PStep;
class G4VIStore;

class G4PIScorer : public G4VPScorer
{

public:  // with description

  G4PIScorer(const G4VIStore &IStore);
    // simple construction

  ~G4PIScorer();
    // simple destruction

  void Score(const G4Step &aStep, const G4PStep &aPStep);
     // perform scoring for the G4Step and G4PStep
 
  const G4PMapPtkTallys &GetMapPtkTallys() const 
  {return fPScorer.GetMapPtkTallys();}
    // get a reference to the G4PMapPtkTallys containing 
    // the tallies
  
  G4bool CorrectWeight() const {return fCorrectWeight;}
    // returns false if the weight was not correct
    // a user might want to get this information 
    // to mark the event a incorrect weight occurred

  void ResetCorrectWeight() {fCorrectWeight = true;}
    // reset fCorrectWeight 
    // the user has to reset the flag fCorrectWeight
    // if he wants to use CorrectWeight()

private:

  const G4VIStore &fIStore; 
  G4PScorer fPScorer;
  G4int fNerr;
  G4bool fCorrectWeight;
};

G4std::ostream& operator<<(G4std::ostream &out, const G4PIScorer &ps);

#endif
