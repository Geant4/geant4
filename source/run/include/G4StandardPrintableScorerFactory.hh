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
// $Id: G4StandardPrintableScorerFactory.hh,v 1.1 2002-07-29 16:03:15 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4StandardPrintableScorer
//
// Class description:
//
// The standard implementation of a G4VPrintableScorerFactory.
// This class may be used to create G4StandradPrintableScorers.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4StandardPrintableScorerFactory_hh
#define G4StandardPrintableScorerFactory_hh G4StandardPrintableScorerFactory_hh

#include "G4VPrintableScorerFactory.hh"

class G4StandardPrintableScorerFactory : public G4VPrintableScorerFactory{
public:
  G4StandardPrintableScorerFactory(){}
  ~G4StandardPrintableScorerFactory(){}
  G4VPrintableScorer *
  CreatePrintableScorer(const G4String &particlename) const ;
};

#endif
