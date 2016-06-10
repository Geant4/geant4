//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4GeometryCellComp.hh 66356 2012-12-18 09:02:32Z gcosmo $
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
