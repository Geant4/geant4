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
// $Id: G4GeometryCellComp.cc,v 1.3 2002-12-13 11:54:19 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeometryCellCmp.cc
//
// ----------------------------------------------------------------------

#include "G4GeometryCellComp.hh"
#include "G4GeometryCell.hh"
//#include "G4VPhysicalVolume.hh"


G4GeometryCellComp::G4GeometryCellComp()
{}

G4bool G4GeometryCellComp::operator() (const G4GeometryCell &k1,
                              const G4GeometryCell &k2) const
{
  G4bool smaler=false;
  if (&(k1.GetPhysicalVolume()) != &(k2.GetPhysicalVolume())) {
    smaler = &(k1.GetPhysicalVolume()) < &(k2.GetPhysicalVolume());
  } else {
    smaler =  k1.GetReplicaNumber() < k2.GetReplicaNumber();
  }
  return smaler;
}
