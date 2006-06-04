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
// $Id: EzDetectorConstruction.hh,v 1.2 2006-06-04 21:36:33 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   EzDetectorConstruction.hh
//
//                                         2005 Q
// ====================================================================
#ifndef EZ_DETECTOR_CONSTRUCTION_H
#define EZ_DETECTOR_CONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class EzDetectorConstruction : public G4VUserDetectorConstruction {

public:
  EzDetectorConstruction();
  ~EzDetectorConstruction();

  virtual G4VPhysicalVolume* Construct();

};

#endif
