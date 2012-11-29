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

// Tubesector builder implementation
// OC 290896

#include <sstream>

#include "G4NURBStubesector.hh"
#include "G4PhysicalConstants.hh"

G4NURBStubesector::G4NURBStubesector(G4double r, G4double R,
                                     G4double DZ, G4double PHI1, G4double PHI2)
  : G4NURBS( 2, 3,  // linear along U, quadratic along V
             5, DecideNbrCtrlPts(PHI1, PHI2),  
                    // rectangle along U,  required stuff along V
                    // we must use a static function which
                    // take the two angles because the
                    // mother constructor is initialised
                    // before everything
             Regular,     // the knot vector along U will be generated
             RegularRep ) // circular like knot vector also
{  
  // check angles
  G4double deltaPHI = PHI2-PHI1;
  while (deltaPHI <= 0) { PHI2 += twopi; deltaPHI += twopi; };

  G4int f = (int)floor(deltaPHI / (halfpi));  //number of pi/2 arcs

  const G4double mr = (r+R)/2;

  const G4double cp1 = std::cos(PHI1);
  const G4double sp1 = std::sin(PHI1);
  const G4double cp2 = std::cos(PHI2);
  const G4double sp2 = std::sin(PHI2);


  // define control points
  CP(mpCtrlPts[ 0] ,  cp1*mr, sp1*mr,  0, 1 );
  CP(mpCtrlPts[ 1] ,  cp1*mr, sp1*mr,  0, 1 );
  CP(mpCtrlPts[ 2] ,  cp1*mr, sp1*mr,  0, 1 );
  CP(mpCtrlPts[ 3] ,  cp1*mr, sp1*mr,  0, 1 );
  CP(mpCtrlPts[ 4] ,  cp1*mr, sp1*mr,  0, 1 );

  CP(mpCtrlPts[ 5] ,  cp1*mr, sp1*mr,  0, 1 );
  CP(mpCtrlPts[ 6] ,  cp1*mr, sp1*mr,  0, 1 );
  CP(mpCtrlPts[ 7] ,  cp1*mr, sp1*mr,  0, 1 );
  CP(mpCtrlPts[ 8] ,  cp1*mr, sp1*mr,  0, 1 );
  CP(mpCtrlPts[ 9] ,  cp1*mr, sp1*mr,  0, 1 );

  CP(mpCtrlPts[10] ,  cp1*r, sp1*r,  DZ, 1 );
  CP(mpCtrlPts[11] ,  cp1*R, sp1*R,  DZ, 1 );
  CP(mpCtrlPts[12] ,  cp1*R, sp1*R, -DZ, 1 );
  CP(mpCtrlPts[13] ,  cp1*r, sp1*r, -DZ, 1 );
  CP(mpCtrlPts[14] ,  cp1*r, sp1*r,  DZ, 1 );

  t_indCtrlPt  i = 15;
  G4double  srcAngle = PHI1;
  G4double  deltaAngleo2;

  G4double destAngle = halfpi + PHI1;

  for(; f > 0; f--)
  {
    // the first arc CP is already Done

    deltaAngleo2 = (destAngle - srcAngle) / 2;
    const G4double csa = std::cos(srcAngle);
    const G4double ssa = std::sin(srcAngle);
    const G4double tdao2 = std::tan(deltaAngleo2); 

    // to calculate the intermediate CP :
    // rotate by srcAngle the (1, tdao2) point
    const t_Coord x = csa - ssa*tdao2;
    const t_Coord y = ssa + csa*tdao2;

    // weight of the CP
    const G4Float weight = (std::cos(deltaAngleo2));

    // initialization. postfix ++ because i initialized to 15
    CP(mpCtrlPts[i++], x*r, y*r,  DZ, 1, weight);
    CP(mpCtrlPts[i++], x*R, y*R,  DZ, 1, weight);
    CP(mpCtrlPts[i++], x*R, y*R, -DZ, 1, weight);
    CP(mpCtrlPts[i++], x*r, y*r, -DZ, 1, weight);
    CP(mpCtrlPts[i++], x*r, y*r,  DZ, 1, weight);

    // end CP (which is the first CP of the next arc)
    const G4double cda = std::cos(destAngle);
    const G4double sda = std::sin(destAngle);
    CP(mpCtrlPts[i++], cda*r, sda*r,  DZ, 1);
    CP(mpCtrlPts[i++], cda*R, sda*R,  DZ, 1);
    CP(mpCtrlPts[i++], cda*R, sda*R, -DZ, 1);
    CP(mpCtrlPts[i++], cda*r, sda*r, -DZ, 1);
    CP(mpCtrlPts[i++], cda*r, sda*r,  DZ, 1);

    // prepare next arc
    srcAngle = destAngle;
    destAngle += halfpi;
  }

  // f == 0, final Arc
  // could be handled in the loops

  destAngle = PHI2;
  deltaAngleo2 = (destAngle - srcAngle) / 2;
  const G4double csa = std::cos(srcAngle);
  const G4double ssa = std::sin(srcAngle);
  const G4double tdao2 = std::tan(deltaAngleo2); 

  // to calculate the intermediate CP :
  // rotate by srcAngle the (1, tdao2) point
  const t_Coord x = csa - ssa*tdao2;
  const t_Coord y = ssa + csa*tdao2;

  // weight of the CP
  const G4Float weight = (std::cos(deltaAngleo2));

  // initialization.
  CP(mpCtrlPts[i++], x*r, y*r,  DZ, 1, weight);
  CP(mpCtrlPts[i++], x*R, y*R,  DZ, 1, weight);
  CP(mpCtrlPts[i++], x*R, y*R, -DZ, 1, weight);
  CP(mpCtrlPts[i++], x*r, y*r, -DZ, 1, weight);
  CP(mpCtrlPts[i++], x*r, y*r,  DZ, 1, weight);

  // end CP
  const G4double cda = std::cos(destAngle);
  const G4double sda = std::sin(destAngle);
  CP(mpCtrlPts[i++], cda*r, sda*r,  DZ, 1);
  CP(mpCtrlPts[i++], cda*R, sda*R,  DZ, 1);
  CP(mpCtrlPts[i++], cda*R, sda*R, -DZ, 1);
  CP(mpCtrlPts[i++], cda*r, sda*r, -DZ, 1);
  CP(mpCtrlPts[i++], cda*r, sda*r,  DZ, 1);

  if (i != (mtotnbrCtrlPts - 10) ) 
  { 
    G4cerr << "\nERROR: G4NURBStubesector::G4NURBStubesector: wrong index,"
           << i << " instead of " << (mtotnbrCtrlPts - 10)
           << "\n\tThe tubesector won't be correct."
           << G4endl;
  }

  CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
  CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
  CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
  CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
  CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);

  CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
  CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
  CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
  CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);
  CP(mpCtrlPts[i++] ,  cp2*mr, sp2*mr, 0, 1);

  // possible to put a DZ DZ -DZ -DZ DZ column to scratch
  // to a line instead of a point

  // creating the nurbs identity
  std::ostringstream  tmpstr;
  tmpstr << "Tubs" << " \tPHI1=" << PHI1 << " ; PHI2=" << PHI2;
  mpwhoami = new char [tmpstr.str().length() + 1];
  mpwhoami = std::strcpy(mpwhoami, tmpstr.str().c_str());
}

const char* G4NURBStubesector::Whoami() const
{
  return mpwhoami;
}

G4NURBStubesector::~G4NURBStubesector()
{
  if (mpwhoami) { delete [] mpwhoami; mpwhoami = 0; }
}

G4NURBStubesector::t_inddCtrlPt
G4NURBStubesector::DecideNbrCtrlPts(G4double PHI1, G4double PHI2)
{
  // check angles
  G4double deltaPHI = PHI2-PHI1;
  while (deltaPHI <= 0) { PHI2 += twopi; deltaPHI += twopi; }
  G4double k = deltaPHI / (halfpi);

  //    G4cerr << " k " << k << G4endl;
  //    G4cerr << " fk " << std::floor(k) << G4endl;
  //    G4cerr <<  " ifk " << ((int)(std::floor(k))) << G4endl;
  //    G4cerr << " n " << (2*((int)(std::floor(k))) + 7) << G4endl;

  return ( 2*((int)(std::floor(k))) + 7 );     
}
