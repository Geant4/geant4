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
// $Id: G4GeometryCell.cc,v 1.1 2002-08-29 15:31:38 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeometryCell.cc
//
// ----------------------------------------------------------------------

#include "G4GeometryCell.hh"

G4GeometryCell::G4GeometryCell(const G4VPhysicalVolume &aVolume,
                                 G4int RepNum)
 : fVPhysiclaVolume(&aVolume),
   fRepNum(RepNum)
{}

G4GeometryCell::~G4GeometryCell()
{}

G4bool G4GeometryCellComp::operator() (const G4GeometryCell &k1,
                              const G4GeometryCell &k2) const
{
  if (k1.fVPhysiclaVolume != k2.fVPhysiclaVolume) {
    return  k1.fVPhysiclaVolume < k2.fVPhysiclaVolume;
  } else {
    return k1.fRepNum < k2.fRepNum;
  }
}

G4bool operator==(const G4GeometryCell &k1, const G4GeometryCell &k2)
{
  if (k1.fVPhysiclaVolume != k2.fVPhysiclaVolume) return false;
  if (k1.fRepNum != k2.fRepNum) return false;
  return true;
}

G4bool operator!=(const G4GeometryCell &k1, const G4GeometryCell &k2)
{
  if (k1.fVPhysiclaVolume != k2.fVPhysiclaVolume) return true;
  if (k1.fRepNum != k2.fRepNum) return true;
  return false;  
}
