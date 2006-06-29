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
// $Id: G4NURBStube.cc,v 1.8 2006-06-29 19:06:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Tube builder implementation
// OC 090796

#include "G4NURBStube.hh"

G4NURBStube::G4NURBStube(G4double r, G4double R, G4double DZ)
  : G4NURBS( 2, 3,  // linear along U, quadratic along V
             5, 9,  // rectangle along U, circle along V
             Regular,     // the knot vector along U will be generated
             RegularRep ) // knot vector for the circle also
{
  // define control points

  const G4double sr2o2 = std::sqrt(2.)/2.;

  CP(mpCtrlPts[ 0] ,  r, 0,  DZ, 1 );
  CP(mpCtrlPts[ 1] ,  R, 0,  DZ, 1 );
  CP(mpCtrlPts[ 2] ,  R, 0, -DZ, 1 );
  CP(mpCtrlPts[ 3] ,  r, 0, -DZ, 1 );
  CP(mpCtrlPts[ 4] ,  r, 0,  DZ, 1 );

  CP(mpCtrlPts[ 5] ,  r, r,  DZ, 1 , sr2o2);
  CP(mpCtrlPts[ 6] ,  R, R,  DZ, 1 , sr2o2);
  CP(mpCtrlPts[ 7] ,  R, R, -DZ, 1 , sr2o2);
  CP(mpCtrlPts[ 8] ,  r, r, -DZ, 1 , sr2o2);
  CP(mpCtrlPts[ 9] ,  r, r,  DZ, 1 , sr2o2);

  CP(mpCtrlPts[10] ,  0, r,  DZ, 1 );
  CP(mpCtrlPts[11] ,  0, R,  DZ, 1 );
  CP(mpCtrlPts[12] ,  0, R, -DZ, 1 );
  CP(mpCtrlPts[13] ,  0, r, -DZ, 1 );
  CP(mpCtrlPts[14] ,  0, r,  DZ, 1 );

  CP(mpCtrlPts[15] , -r, r,  DZ, 1 , sr2o2);
  CP(mpCtrlPts[16] , -R, R,  DZ, 1 , sr2o2);
  CP(mpCtrlPts[17] , -R, R, -DZ, 1 , sr2o2);
  CP(mpCtrlPts[18] , -r, r, -DZ, 1 , sr2o2);
  CP(mpCtrlPts[19] , -r, r,  DZ, 1 , sr2o2);

  CP(mpCtrlPts[20] , -r, 0,  DZ, 1 );
  CP(mpCtrlPts[21] , -R, 0,  DZ, 1 );
  CP(mpCtrlPts[22] , -R, 0, -DZ, 1 );
  CP(mpCtrlPts[23] , -r, 0, -DZ, 1 );
  CP(mpCtrlPts[24] , -r, 0,  DZ, 1 );

  CP(mpCtrlPts[25] , -r,-r,  DZ, 1 , sr2o2);
  CP(mpCtrlPts[26] , -R,-R,  DZ, 1 , sr2o2);
  CP(mpCtrlPts[27] , -R,-R, -DZ, 1 , sr2o2);
  CP(mpCtrlPts[28] , -r,-r, -DZ, 1 , sr2o2);
  CP(mpCtrlPts[29] , -r,-R,  DZ, 1 , sr2o2);

  CP(mpCtrlPts[30] ,  0,-r,  DZ, 1 );
  CP(mpCtrlPts[31] ,  0,-R,  DZ, 1 );
  CP(mpCtrlPts[32] ,  0,-R, -DZ, 1 );
  CP(mpCtrlPts[33] ,  0,-r, -DZ, 1 );
  CP(mpCtrlPts[34] ,  0,-r,  DZ, 1 );

  CP(mpCtrlPts[35] ,  r,-r,  DZ, 1 , sr2o2);
  CP(mpCtrlPts[36] ,  R,-R,  DZ, 1 , sr2o2);
  CP(mpCtrlPts[37] ,  R,-R, -DZ, 1 , sr2o2);
  CP(mpCtrlPts[38] ,  r,-r, -DZ, 1 , sr2o2);
  CP(mpCtrlPts[39] ,  r,-r,  DZ, 1 , sr2o2);

  CP(mpCtrlPts[40] ,  r, 0,  DZ, 1 );
  CP(mpCtrlPts[41] ,  R, 0,  DZ, 1 );
  CP(mpCtrlPts[42] ,  R, 0, -DZ, 1 );
  CP(mpCtrlPts[43] ,  r, 0, -DZ, 1 );
  CP(mpCtrlPts[44] ,  r, 0,  DZ, 1 );
}

const char* G4NURBStube::Whoami() const
{
  return "Tube";
}
