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
// $Id: G4VPrintableScorer.hh,v 1.1 2002-07-29 15:56:19 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VPrintableScorer
//
// Class description:
//
// This interface extends the G4VPScorer. A scorer with this
// interface can print it's values to an ostream.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VPrintableScorer_hh
#define G4VPrintableScorer_hh G4VPrintableScorer_hh 

#include "globals.hh"
#include "g4std/iostream"
#include "G4VPScorer.hh"

class G4VIStore;

class G4VPrintableScorer : public G4VPScorer {
public:
  virtual ~G4VPrintableScorer() {}
  virtual void Score(const G4Step &step, const G4PStep &pstep) = 0;
  virtual void PrintTable(G4std::ostream *out,
			  const G4VIStore *istore = 0) = 0;
  virtual G4VPScorer *GetPointerToScorer() = 0;
  // this pointer may be used to strip off one level of inheritance
};

#endif
