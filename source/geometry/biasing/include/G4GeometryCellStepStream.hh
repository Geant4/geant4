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
// $Id: G4GeometryCellStepStream.hh,v 1.1 2002-10-22 13:18:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Declarations
//
// Declaration description:
//
// declarations of streams for G4GeometryCell and G4GeometryCellStep

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4GeometryCellStepStream_hh
#define G4GeometryCellStepStream_hh G4GeometryCellStepStream_hh

#include "G4GeometryCell.hh"
#include "G4GeometryCellStep.hh"

G4std::ostream& operator<<(G4std::ostream &out, const G4GeometryCell &tk);
G4std::ostream& operator<<(G4std::ostream &out, const G4GeometryCellStep &ps);

#endif
