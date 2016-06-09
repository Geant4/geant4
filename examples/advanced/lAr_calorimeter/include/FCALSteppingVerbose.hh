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
// $Id: FCALSteppingVerbose.hh,v 1.5 2004/11/29 18:03:05 ribon Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//  
//---------------------------------------------------------------
//
// FCALSteppingVerbose.hh
//
// Description:
//   This class manages the verbose outputs in G4SteppingManager. 
//   It inherits from G4SteppingVerbose   
//
//---------------------------------------------------------------

class FCALSteppingVerbose;

#ifndef FCALSteppingVerbose_h
#define FCALSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class FCALSteppingVerbose : public G4SteppingVerbose {
public:   
// Constructor/Destructor
  FCALSteppingVerbose();
 ~FCALSteppingVerbose();
//
  void StepInfo();
  void TrackingStarted();
//


};

#endif
