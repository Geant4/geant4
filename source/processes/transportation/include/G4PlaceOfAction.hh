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
// $Id: G4PlaceOfAction.hh,v 1.3 2003/11/26 14:51:48 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// ----------------------------------------------------------------------
// Class G4PlaceOfAction
//
// Class description:
// This is an enum to specify when weight window sampling shoul be
// applied: onBoundary, onCollision, onBoundaryAndCollision.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4PlaceOfAction_hh
#define G4PlaceOfAction_hh G4PlaceOfAction_hh

enum G4PlaceOfAction
{
  onBoundary = 1, 
  onCollision = 2, 
  onBoundaryAndCollision = 3 
};

#endif
