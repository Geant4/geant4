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
// $Id: G4GRSSolid.cc,v 1.4 2001-07-11 10:00:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4GRSSolid Implementation

#include "G4GRSSolid.hh"

G4GRSSolid::~G4GRSSolid()
{
  delete frot;			// safe if null
}

G4GRSSolid::G4GRSSolid(const G4GRSSolid& right)
{
  if ((&right) && (&right != this))
  {
    fsolid = right.fsolid;
    ftlate = right.ftlate;
    if (frot)
    {
      delete frot;
      frot = 0;
    }
    if (right.frot)
      frot = new G4RotationMatrix(*(right.frot));
  }
}

G4GRSSolid& G4GRSSolid::operator=(const G4GRSSolid& right)
{
  if (&right == this) return *this;
  if (&right)
  {
    fsolid = right.fsolid;
    ftlate = right.ftlate;
    if (frot)
    {
      delete frot;
      frot = 0;
    }
    if (right.frot)
      frot = new G4RotationMatrix(*(right.frot));
  }
  return *this;
}
