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
// $Id: G4PStepStream.cc,v 1.3 2002-04-09 16:23:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4PStepStream.cc
//
// ----------------------------------------------------------------------

#include "G4PStepStream.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PStep.hh"

G4std::ostream& operator<<(G4std::ostream &out, const G4PTouchableKey &tk)
{
  out << "Volume name = " << tk.fVPhysiclaVolume->GetName() << ", ";
  out << "Replica number = " << tk.fRepNum;
  return out;
}

G4std::ostream& operator<<(G4std::ostream &out, const G4PStep &ps)
{
  out << "PreTouchableKey : " <<  ps.fPreTouchableKey << " ";
  out << "PostTouchableKey: " <<  ps.fPostTouchableKey << " ";
  out << "CrossBoundary   : " <<  ps.fCrossBoundary << "\n";
  return out;
}
