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
// $Id: G4PointRat.cc,v 1.6 2006-06-29 18:42:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4PointRat.cc
//
// ----------------------------------------------------------------------

#include "G4PointRat.hh"

G4PointRat::G4PointRat()
 : pt3d(), s(1)
{
}

G4PointRat::G4PointRat(const G4Point3D& tmp)
 : pt3d(tmp), s(1)
{
}

G4PointRat::~G4PointRat()
{
}

G4PointRat& G4PointRat::operator=(const G4PointRat& a)
{
    pt3d.setX(a.x());
    pt3d.setY(a.y());
    pt3d.setZ(a.z());
    s=a.w();
    
    return *this;
}

G4PointRat& G4PointRat::operator=(const G4Point3D& a)
{
    pt3d = a;
    s=1;
    
    return *this;
}
