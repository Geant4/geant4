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
// $Id: G4GeometryCellComp.hh,v 1.1 2002-10-22 13:18:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4GeometryCellComp
//
// Class description:
//
// A class needed for comparing G4GeometryCell e.g. in stl containers.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4GeometryCellComp_hh
#define G4GeometryCellComp_hh G4GeometryCellComp_hh 

#include "globals.hh"

class G4GeometryCell;

class G4GeometryCellComp
{

public:  // without description
  G4GeometryCellComp();

  G4bool operator() (const G4GeometryCell &g1,
                     const G4GeometryCell &g2) const;
    // returns true if g1 < g2
};


#endif
