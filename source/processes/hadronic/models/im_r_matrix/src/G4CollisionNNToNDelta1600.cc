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
// $Id: G4CollisionNNToNDelta1600.cc,v 1.2.2.1 2004/03/24 13:18:36 hpw Exp $ //

#include "globals.hh"
#include "G4CollisionNNToNDelta1600.hh"
#include "G4ConcreteNNToNDeltaStar.hh"

G4CollisionNNToNDelta1600::G4CollisionNNToNDelta1600()
{ 
  MakeNNToNDelta<Dm_1600PC, D0_1600PC, Dp_1600PC, Dpp_1600PC, G4ConcreteNNToNDeltaStar>::Make(this);
}
