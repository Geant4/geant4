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
// $Id: G4PScorer.hh,v 1.5 2002-04-10 13:14:16 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PScorer
//
// Class description:
//
// A scorer for standard tallies. Under development.
// See also G4VPScorer.
// It scores using informations form G4Step and G4PStep.
// The tallies are contained and related with "cells"
// in the container G4PMapPtkTallys (see G4PMapPtkTallys.hh)

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4Pscorer_hh
#define G4Pscorer_hh G4Pscorer_hh

#include "g4std/iostream"

#include "G4VPScorer.hh"
#include "G4PMapPtkTallys.hh"

class G4Step;
class G4PStep;

class G4PScorer : public G4VPScorer
{

public:  // with description

  G4PScorer();
    // simple construction

  ~G4PScorer();
    // simple destruction

  void Score(const G4Step &aStep, const G4PStep &aPStep);
     // perform scoring for the G4Step and G4PStep
 
  const G4PMapPtkTallys &GetMapPtkTallys() const {return fPtkTallys;}
    // get a reference to the G4PMapPtkTallys containing 
    // the tallies

private:

  G4PMapPtkTallys fPtkTallys;
};

G4std::ostream& operator<<(G4std::ostream &out, const G4PScorer &ps);

#endif
