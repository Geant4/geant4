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
// $Id: G4VProcessPlacer.hh,v 1.2 2002-04-09 17:40:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VProcessPlacer
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VProcessPlacer_hh
#define G4VProcessPlacer_hh G4VProcessPlacer_hh

#include "globals.hh"

class G4VProcess;
 
class G4VProcessPlacer
{

public:  // with description

  G4VProcessPlacer(const G4String &particlename){}
  virtual ~G4VProcessPlacer(){}

  virtual void AddProcessAsLastDoIt(G4VProcess *process) = 0;
  virtual void AddProcessAsSecondDoIt(G4VProcess *process) = 0;
};

#endif
