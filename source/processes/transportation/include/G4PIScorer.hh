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
// $Id: G4PIScorer.hh,v 1.4 2002-04-09 17:40:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PIScorer
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4PIScorer_hh
#define G4PIScorer_hh

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
  ~G4PIScorer();

  void Score(const G4Step &aStep, const G4PStep &aPStep);
  const G4PMapPtkTallys &GetMapPtkTallys() const {return fPScorer.GetMapPtkTallys();}
  G4bool CorrectWeight() const {return fCorrectWeight;}
  void ResetCorrectWeight() {fCorrectWeight = true;}

private:

  const G4VIStore &fIStore; 
  G4PScorer fPScorer;
  G4int fNerr;
  G4bool fCorrectWeight;
};

G4std::ostream& operator<<(G4std::ostream &out, const G4PIScorer &ps);

#endif
