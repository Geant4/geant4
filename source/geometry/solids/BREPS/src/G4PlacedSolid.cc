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
// $Id: G4PlacedSolid.cc,v 1.6 2007-05-11 13:49:32 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4PlacedSolid.cc
//
// ----------------------------------------------------------------------

#include "G4PlacedSolid.hh"
#include "G4GeometryTolerance.hh"

G4PlacedSolid::G4PlacedSolid()
{
  solid =0;
  solidRotation =0;
  solidTranslation =0;
}


G4PlacedSolid::G4PlacedSolid(G4BREPSolid* s, G4Axis2Placement3D* p)
{
  solid =s;
  if(p)
  {
    G4double x,y,z;
    G4Point3D srfpoint = p->GetLocation();
    
    x = srfpoint.x();
    y = srfpoint.y();
    z = srfpoint.z();
    solidTranslation = new G4ThreeVector(x,y,z);
    
    G4Vector3D tmpvec = p->GetAxis();
    x = tmpvec.x();
    y = tmpvec.y();
    z = tmpvec.z();
    G4ThreeVector x_axis(x,y,z);

    G4double kCarTolerance = G4GeometryTolerance::GetInstance()
                             ->GetSurfaceTolerance();
    if( (x<kCarTolerance)&&
	(y<kCarTolerance)&&
	(z<kCarTolerance)   )
      solidRotation=0;
    else
    {
      tmpvec = p->GetRefDirection();
      x = tmpvec.x();
      y = tmpvec.y();
      z = tmpvec.z();

      G4ThreeVector y_axis(x,y,z);
      solidRotation = new G4RotationMatrix();
      solidRotation->rotateAxes(x_axis, y_axis, x_axis.cross(y_axis));
    }
  }
  else
  {
    solidTranslation=0;
    solidRotation=0;
  }
}


G4PlacedSolid::~G4PlacedSolid()
{
  //delete solid;
  if(solidRotation)
    delete solidRotation;
  
  if(solidTranslation)
    delete solidTranslation;
}
