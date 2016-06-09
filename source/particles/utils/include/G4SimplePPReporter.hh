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
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SimplePPReporter.hh,v 1.1 2003/09/21 19:38:50 kurasige Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// ---------------------------------------------------------------
#ifndef G4SimplePPReporter_h
#define G4SimplePPReporter_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4VParticlePropertyReporter.hh"

class G4SimplePPReporter: public G4VParticlePropertyReporter
{
 public:
  //constructors
  G4SimplePPReporter();
  
  //destructor
  virtual ~G4SimplePPReporter();
  
public:
  virtual void Print(const G4String& option="");
};


#endif
