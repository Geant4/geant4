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
// $Id: Tst33VEventAction.hh,v 1.3 2002-11-20 13:09:17 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Classs Tst33VEventAction
//
// Class description:
//
// Base class for the event actions used in this test.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
// 


#ifndef Tst33VEventAction_h
#define Tst33VEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"


class G4CellScorer;

class Tst33VEventAction : public G4UserEventAction
{
public:
  Tst33VEventAction();
  virtual ~Tst33VEventAction();
  virtual void SpecialCellScorer(const G4CellScorer *scorer) = 0;
  virtual void Clear() = 0;
};

#endif

    
