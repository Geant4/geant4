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
// $Id: Test17SteppingVerbose.hh,v 1.2 2001-07-11 10:10:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//---------------------------------------------------------------
//
// Test17SteppingVerbose.hh
//
// Description:
//   This class manages the vervose outputs in G4SteppingManager. 
//   
//
//---------------------------------------------------------------

class Test17SteppingVerbose;

#ifndef Test17SteppingVerbose_h
#define Test17SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class Test17SteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  Test17SteppingVerbose();
 ~Test17SteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
