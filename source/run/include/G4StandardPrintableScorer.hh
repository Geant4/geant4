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
// $Id: G4StandardPrintableScorer.hh,v 1.1 2002-07-29 16:03:15 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4StandardPrintableScorer
//
// Class description:
//
// The standard implementation of a G4VPrintableScorer.
// It uses the G4StandadScorer and the G4StandradScoreTable.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4StandardPrintableScorer_hh
#define G4StandardPrintableScorer_hh G4StandardPrintableScorer_hh

#include "G4VPrintableScorer.hh"
class G4StandardScorer;

class G4StandardPrintableScorer: public G4VPrintableScorer {
public:
  G4StandardPrintableScorer();
  ~G4StandardPrintableScorer();
  void Score(const G4Step &step, const G4PStep &pstep);
  void PrintTable(G4std::ostream *out,
		  const G4VIStore *istore = 0);
  G4VPScorer *GetPointerToScorer();
private:
  G4StandardScorer *fStandardScorer;
};

#endif
