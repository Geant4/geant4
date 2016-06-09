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
// $Id: G4Conic.hh,v 1.7 2006-06-29 18:38:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Conic
//
// Class description:
// 
// Definition of a generic conical curve.

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __CONIC_H
#define __CONIC_H 

#include "G4Curve.hh"
#include "G4Axis2Placement3D.hh"

class G4Conic : public G4Curve
{

public:  // with description
 
  G4Conic();
  virtual ~G4Conic();
    // Constructor & destructor.

  G4Conic(const G4Conic& right);
  G4Conic& operator=(const G4Conic& right);
    // Copy contructor and assignment operator.

  inline const G4Axis2Placement3D* GetPosition() const;
  inline G4double GetPShift() const;
  inline void SetPShift(G4double pShift0); 
    // Get/Set geometric data.

public:  // without description

  //inline G4Placement GetPosition() {return Position;}
  //virtual const char *Name(){return "G4ConicalCurve";}

protected:

  //void ProjectCurve(const G4Plane&, const G4Plane&);
  //G4int HitPartOfCurve(G4double, G4double, const G4Point2d&);
  //G4Placement Position;

  G4Axis2Placement3D position;
    // Geometric data.

  G4double pShift;
    // pShift must be added/subtracted from the parameter.
    // No STEP I/O if not 0!!!
    // Set by Project members.

};

#include "G4Conic.icc"

#endif
