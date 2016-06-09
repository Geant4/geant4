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
// $Id: G4VTrajectoryPoint.cc,v 1.1 2004/07/05 17:08:16 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02-patch-02 $
//
//
//---------------------------------------------------------------
//
// G4VTrajectoryPoint.cc
//
// ---------------------------------------------------------------

#include "G4VTrajectoryPoint.hh"

G4VTrajectoryPoint::G4VTrajectoryPoint() {;}
G4VTrajectoryPoint::~G4VTrajectoryPoint() {;}

G4bool G4VTrajectoryPoint::operator==(const G4VTrajectoryPoint& right) const
{
  return (this==&right);
}
