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
// $Id: G4ImportanceFinder.cc,v 1.7 2002-10-14 12:36:03 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ImportanceFinder.cc
//
// ----------------------------------------------------------------------

#include "g4std/strstream"
#include "G4VParallelStepper.hh"
#include "G4ImportanceFinder.hh"
#include "G4VIStore.hh"
#include "G4PStepStream.hh"

G4ImportanceFinder::G4ImportanceFinder(const G4VIStore &aIStore)
  : fIStore(aIStore)
{}

G4ImportanceFinder::~G4ImportanceFinder()
{}

G4double
G4ImportanceFinder::GetImportance(const G4GeometryCell &gCell) const
{  
  G4double  imp = fIStore.GetImportance(gCell);
  // importances < 0 are not allowed
  if (imp < 0) {
    G4std::ostrstream os;
    os << "imp < 0: GeometryCell = " << gCell  << '\0';
    Error(os.str());
  }
  return imp;
}

void G4ImportanceFinder::Error(const G4String &m) const
{
   G4std::G4cout << "ERROR - G4ImportanceFinder::" << m << G4endl;
   G4std::G4Exception("Program aborted.");
}
