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
// $Id: GammaRayTelTrackerROGeometry.hh,v 1.3 2001-11-29 11:19:17 griccard Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelTrackerROGeometry  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#ifndef GammaRayTelTrackerROGeometry_h
#define GammaRayTelTrackerROGeometry_h 1

#include "G4VReadOutGeometry.hh"

class GammaRayTelDetectorConstruction;

class GammaRayTelTrackerROGeometry : public G4VReadOutGeometry
{
public:
  GammaRayTelTrackerROGeometry();
  GammaRayTelTrackerROGeometry(G4String);
  ~GammaRayTelTrackerROGeometry();

private:
  G4VPhysicalVolume* Build();
  GammaRayTelDetectorConstruction* GammaRayTelDetector;
  //pointer to the geometry
};

#endif


