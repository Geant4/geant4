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
// $Id: G4GeometryCell.cc,v 1.3 2002-09-02 15:22:32 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeometryCell.cc
//
// ----------------------------------------------------------------------

#include "G4GeometryCell.hh"
#include "G4VPhysicalVolume.hh"

G4GeometryCell::G4GeometryCell(const G4VPhysicalVolume &aVolume,
                                 G4int RepNum)
 : fVPhysiclaVolume(&aVolume),
   fRepNum(RepNum)
{}

G4GeometryCell::~G4GeometryCell()
{}

G4GeometryCell::G4GeometryCell(const G4GeometryCell &rhs){
  *this = rhs;
}

G4GeometryCell &G4GeometryCell::operator=(const G4GeometryCell &rhs){
  if (this == &rhs) {
    return *this;
  }
  fVPhysiclaVolume = rhs.fVPhysiclaVolume; // this is treated 
                                           // as identifyer
  fRepNum = rhs.fRepNum;
  return *this;
}


G4bool G4GeometryCellComp::operator() (const G4GeometryCell &k1,
                              const G4GeometryCell &k2) const
{
  if (&(k1.GetPhysicalVolume()) != &(k2.GetPhysicalVolume())) {
    return  &(k1.GetPhysicalVolume()) < &(k2.GetPhysicalVolume());
  } else {
    return k1.GetReplicaNumber() < k2.GetReplicaNumber();
  }
}

G4bool operator==(const G4GeometryCell &k1, const G4GeometryCell &k2)
{
  if (&(k1.GetPhysicalVolume()) != &(k2.GetPhysicalVolume())) return false;
  if (k1.GetReplicaNumber() != k2.GetReplicaNumber()) return false;
  return true;
}

G4bool operator!=(const G4GeometryCell &k1, const G4GeometryCell &k2)
{
  if (&(k1.GetPhysicalVolume()) != &(k2.GetPhysicalVolume())) return true;
  if (k1.GetReplicaNumber() != k2.GetReplicaNumber()) return true;
  return false;  
}
