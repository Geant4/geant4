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
// $Id: G4GeometryCellStepStream.cc,v 1.1 2002-10-22 13:18:47 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeometryCellStepStream.cc
//
// ----------------------------------------------------------------------

#include "G4GeometryCellStepStream.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeometryCellStep.hh"

G4std::ostream& operator<<(G4std::ostream &out, const G4GeometryCell &tk)
{
  out << "Volume name = " << tk.GetPhysicalVolume().GetName() << ", ";
  out << "Replica number = " << tk.GetReplicaNumber();
  return out;
}

G4std::ostream& operator<<(G4std::ostream &out, const G4GeometryCellStep &ps)
{
  out << "PreGeometryCell : " <<  ps.GetPreGeometryCell() << " ";
  out << "PostGeometryCell: " <<  ps.GetPostGeometryCell() << " ";
  out << "CrossBoundary   : " <<  ps.GetCrossBoundary() << "\n";
  return out;
}
