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
// $Id: G4ImportanceFinder.cc,v 1.5 2002-08-29 15:30:51 dressel Exp $
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
G4ImportanceFinder::GetIPre_over_IPost(const G4GeometryCell &prekey,
		                       const G4GeometryCell &postkey) 
  const
{  
  G4double ratio = 1;
  // if either the pro or the post gCell is not known
  // the ratio of pre over post importance is set to 1
  // so no splitting od RR is done
  if ( fIStore.IsKnown(prekey) && fIStore.IsKnown(postkey) ) {
    G4double  ipre = fIStore.GetImportance(prekey);
    G4double ipost = fIStore.GetImportance(postkey);
    // importances < 0 are not allowed
    if (ipre < 0 || ipost < 0 ) {
      G4std::ostrstream os;
      os << "ipre < 0 || ipost < 0, preGeometryCell = " << prekey 
	 << ", postGeometryCell = " << postkey << '\0';
      Error(os.str());
    }
    // importances == 0 mean don't do any biaisng here
    else if (ipre == 0 || ipost == 0) {
      ratio = 1;
    }
    else {
      ratio =  ipre/ipost;
    }
  }
  else {
    ratio = 1;
  }

  return ratio;
}

void G4ImportanceFinder::Error(const G4String &m) const
{
   G4cout << "ERROR - G4ImportanceFinder::" << m << G4endl;
   G4Exception("Program aborted.");
}
