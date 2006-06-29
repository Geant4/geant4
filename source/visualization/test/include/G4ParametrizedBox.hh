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
// $Id: G4ParametrizedBox.hh,v 1.4 2006-06-29 21:33:46 gunter Exp $
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
