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

// Cylinder builder implementation
// OC 090796

#include "G4NURBScylinder.hh"

// the cylinder constructor use the first G4NURBS constructor
// look in G4NURBStube if you want to see how to use the second one.

G4NURBScylinder::G4NURBScylinder(G4double R, G4double DZ)
  : G4NURBS  ( 2, 3,  // linear along U, quadratic along V
               4, 9,  // half rectangle along U, circle along V
               (new t_CtrlPt [ 4 * 9 ]),  // the array for CtrlPts
               0,     // the knot vector along U will be generated
               (new t_Knot [ 3 + 9 ]) )   // knot vector for the circle
{
  // define the V knot vector
  m[V].pKnots[ 0] = 0;
  m[V].pKnots[ 1] = 0;
  m[V].pKnots[ 2] = 0;
  m[V].pKnots[ 3] = 0.25;
  m[V].pKnots[ 4] = 0.25;
  m[V].pKnots[ 5] = 0.5;
  m[V].pKnots[ 6] = 0.5;
  m[V].pKnots[ 7] = 0.75;
  m[V].pKnots[ 8] = 0.75;
  m[V].pKnots[ 9] = 1;
  m[V].pKnots[10] = 1;
  m[V].pKnots[11] = 1;

  // define control points

  const G4double sr2o2 = std::sqrt(2.)/2. ;

  CP(mpCtrlPts[ 0] ,  0,  0,  DZ, 1 );
  CP(mpCtrlPts[ 1] ,  R, 0,  DZ, 1 );
  CP(mpCtrlPts[ 2] ,  R, 0, -DZ, 1 );
  CP(mpCtrlPts[ 3] ,  0,  0, -DZ, 1 );

  CP(mpCtrlPts[ 4] ,  0,  0,  DZ, 1 );
  CP(mpCtrlPts[ 5] ,  R, R, DZ, 1 , sr2o2);
  CP(mpCtrlPts[ 6] ,  R, R,-DZ, 1 , sr2o2);
  CP(mpCtrlPts[ 7] ,  0,  0, -DZ, 1 );

  CP(mpCtrlPts[ 8] ,  0,  0,  DZ, 1 );
  CP(mpCtrlPts[ 9] ,  0, R,  DZ, 1 );
  CP(mpCtrlPts[10] ,  0, R, -DZ, 1 );
  CP(mpCtrlPts[11] ,  0,  0, -DZ, 1 );

  CP(mpCtrlPts[12] ,  0,  0,  DZ, 1 );
  CP(mpCtrlPts[13] , -R, R, DZ, 1 , sr2o2);
  CP(mpCtrlPts[14] , -R, R,-DZ, 1 , sr2o2);
  CP(mpCtrlPts[15] ,  0,  0, -DZ, 1 );

  CP(mpCtrlPts[16] ,  0,  0,  DZ, 1 );
  CP(mpCtrlPts[17] , -R, 0,  DZ, 1 );
  CP(mpCtrlPts[18] , -R, 0, -DZ, 1 );
  CP(mpCtrlPts[19] ,  0,  0, -DZ, 1 );

  CP(mpCtrlPts[20] ,  0,  0,  DZ, 1 );
  CP(mpCtrlPts[21] , -R,-R, DZ, 1 , sr2o2);
  CP(mpCtrlPts[22] , -R,-R,-DZ, 1 , sr2o2);
  CP(mpCtrlPts[23] ,  0,  0, -DZ, 1 );

  CP(mpCtrlPts[24] ,  0,  0,  DZ, 1 );
  CP(mpCtrlPts[25] ,  0, -R, DZ, 1 );
  CP(mpCtrlPts[26] ,  0, -R,-DZ, 1 );
  CP(mpCtrlPts[27] ,  0,  0, -DZ, 1 );

  CP(mpCtrlPts[28] ,  0,  0,  DZ, 1 );
  CP(mpCtrlPts[29] ,  R,-R, DZ, 1 , sr2o2);
  CP(mpCtrlPts[30] ,  R,-R,-DZ, 1 , sr2o2);
  CP(mpCtrlPts[31] ,  0,  0, -DZ, 1 );

  CP(mpCtrlPts[32] ,  0,  0,  DZ, 1 );
  CP(mpCtrlPts[33] ,  R, 0,  DZ, 1 );
  CP(mpCtrlPts[34] ,  R, 0, -DZ, 1 );
  CP(mpCtrlPts[35] ,  0,  0, -DZ, 1 );
}

const char* G4NURBScylinder::Whoami() const
{
  return "Cylinder";
}
