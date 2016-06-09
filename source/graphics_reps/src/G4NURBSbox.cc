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
// $Id$
//
// 
// Olivier Crumeyrolle  12 September 1996

// Box builder implementation (KidBox)
// OC 060996

#include "G4NURBSbox.hh"

G4NURBSbox::G4NURBSbox(G4double DX, G4double DY, G4double DZ)
  : G4NURBS  ( 2, 2,  // linear along U and V
               4, 5 ) // line with two 90 degrees folds along U
                      // rectangle along V (3 folds)
{
  // let's it Generate regular knots vector
  // (note we are calling the second constructor)

  t_indCtrlPt i = 0;

  CP(mpCtrlPts[i++], DX, DY, DZ, 1 ); 
  CP(mpCtrlPts[i++], DX, DY, DZ, 1 ); 
  CP(mpCtrlPts[i++], DX, DY,-DZ, 1 ); 
  CP(mpCtrlPts[i++], DX, DY,-DZ, 1 ); 

  CP(mpCtrlPts[i++], DX, DY, DZ, 1 ); 
  CP(mpCtrlPts[i++],-DX, DY, DZ, 1 ); 
  CP(mpCtrlPts[i++],-DX, DY,-DZ, 1 ); 
  CP(mpCtrlPts[i++], DX, DY,-DZ, 1 ); 

  CP(mpCtrlPts[i++], DX,-DY, DZ, 1 ); 
  CP(mpCtrlPts[i++],-DX,-DY, DZ, 1 ); 
  CP(mpCtrlPts[i++],-DX,-DY,-DZ, 1 ); 
  CP(mpCtrlPts[i++], DX,-DY,-DZ, 1 ); 

  CP(mpCtrlPts[i++], DX,-DY, DZ, 1 ); 
  CP(mpCtrlPts[i++], DX,-DY, DZ, 1 ); 
  CP(mpCtrlPts[i++], DX,-DY,-DZ, 1 ); 
  CP(mpCtrlPts[i++], DX,-DY,-DZ, 1 ); 

  CP(mpCtrlPts[i++], DX, DY, DZ, 1 ); 
  CP(mpCtrlPts[i++], DX, DY, DZ, 1 ); 
  CP(mpCtrlPts[i++], DX, DY,-DZ, 1 ); 
  CP(mpCtrlPts[i++], DX, DY,-DZ, 1 ); 
}

const char*  G4NURBSbox::Whoami() const
{
  return "Box (as a folded piece)";
}
