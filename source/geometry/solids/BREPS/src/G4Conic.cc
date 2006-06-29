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
// $Id: G4Conic.cc,v 1.8 2006-06-29 18:41:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4conic.cc
//
// ----------------------------------------------------------------------

#include "G4Conic.hh"

G4Conic::G4Conic () : pShift(0)
{
}

G4Conic::~G4Conic()
{
}

G4Conic::G4Conic(const G4Conic& right)
  : G4Curve(), position(right.position), pShift(right.pShift)
{
  bBox      = right.bBox;
  start     = right.start;
  end       = right.end;
  pStart    = right.pStart;
  pEnd      = right.pEnd;
  pRange    = right.pRange;
  bounded   = right.bounded;
  sameSense = right.sameSense;
}

G4Conic& G4Conic::operator=(const G4Conic& right)
{
  if (&right == this) return *this;

  pShift   = right.pShift;
  position = right.position;
  bBox      = right.bBox;
  start     = right.start;
  end       = right.end;
  pStart    = right.pStart;
  pEnd      = right.pEnd;
  pRange    = right.pRange;
  bounded   = right.bounded;
  sameSense = right.sameSense;

  return *this;
}

/*
void
G4ConicalCurve::ProjectCurve(const G4Plane& Pl1, const G4Plane& Pl2)
{  
  // Curve start
  Project(ProjStart, Start, Pl1, Pl2);
  // Curve end
  Project(ProjEnd, End, Pl1, Pl2);
  // Placement
  Position.ProjectPlacement(Pl1,Pl2);
}

G4int
G4ConicalCurve::HitPartOfCurve(G4double Angle, G4double Solution,
                               const G4Point2d& ProjHit)
{
    	// Check if Solution1 is part of the curve i.e. in the "pie"
	G4double TmpSol1 = Solution - ProjHit.X();
	G4Point2d ArcHit1(TmpSol1, ProjHit.Y());
	G4double Cross1 = CrossProduct( ProjStart, ArcHit1);
	G4double Cross2 = CrossProduct( ArcHit1  , ProjEnd);
	if( (Angle<=0 && Cross1<=0 && Cross2 <=0) ||
	    (Angle> 0 && !(Cross1>=0 && Cross2 >=0)) )
	    // Solution1 is on the curve
	  return 1;
	return 0;
}
*/
