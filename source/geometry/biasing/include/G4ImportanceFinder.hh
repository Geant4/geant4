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
// $Id: G4ImportanceFinder.hh,v 1.6 2002-10-14 12:36:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ImportanceFinder
//
// Class description:
//
// This class is used internally by importance sampling in a parallel
// geometry. It derives the ratio of importances of the pre over the post
// G4GeometryCell. 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ImportanceFinder_hh
#define G4ImportanceFinder_hh G4ImportanceFinder_hh

#include "globals.hh"

class G4VParallelStepper;
class G4VIStore;
class G4GeometryCell;

class G4ImportanceFinder
{

public:  // with description

  explicit G4ImportanceFinder(const G4VIStore &aIStore);
    // initialise fIStore

  ~G4ImportanceFinder();
    // destruction

  G4double GetImportance(const G4GeometryCell &gCell) const;
  // get the importance to a G4GeometryCell

private:

  void Error(const G4String &m) const;

private:

  const G4VIStore &fIStore;
};

#endif
