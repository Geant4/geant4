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
// $Id: G4PointRat.cc,v 1.5 2001-07-11 09:59:46 gunter Exp $
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
