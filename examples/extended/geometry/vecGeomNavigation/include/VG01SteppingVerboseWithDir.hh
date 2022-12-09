//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file VG01SteppingVerboseWithDir.hh
/// \brief Definition of the VG01SteppingVerboseWithDir class
//
// Alternative Stepping Verbose method that prints out:
//   - momentum direction
//   - safety  ( distance to nearest boundary at pre-step point )
//
// ********************************************************************
// Code developed by:
//    J. Apostolakis,   April 2021 
// ********************************************************************

#ifndef VG01SteppingVerboseWithDir_hh
#define VG01SteppingVerboseWithDir_hh 1

#include "G4SteppingVerbose.hh"

class VG01SteppingVerboseWithDir : public G4SteppingVerbose {

public:   
  
  //Constructor/Destructor
  VG01SteppingVerboseWithDir();
  ~VG01SteppingVerboseWithDir();
  
  void StepInfo() override;
  void TrackingStarted() override;

  // implementation methods - simplify/reduce code
  void Banner();
};

#endif
