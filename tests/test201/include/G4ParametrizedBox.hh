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
// $Id: G4ParametrizedBox.hh,v 1.3 2001-07-11 10:10:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Original class written by Hans-Peter Wellisch.

#include "G4VPVParameterisation.hh"
 class G4ParametrizedBox: public G4VPVParameterisation
 {
  virtual void ComputeTransformation(const G4int n,
                                     G4VPhysicalVolume* pRep) const
  {
    pRep->SetTranslation(G4ThreeVector(0,(n-1)*15*m,0));
  }

  virtual void ComputeDimensions(G4Box &pBox,
                                 const G4int n,
                                 const G4VPhysicalVolume* pRep) const
  {
    pBox.SetXHalfLength(10*m*n);
    pBox.SetYHalfLength(5*m*n);
    pBox.SetZHalfLength(5*m*n);
  }
 };
