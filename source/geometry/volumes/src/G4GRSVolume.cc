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
// $Id: G4GRSVolume.cc,v 1.5 2002-08-06 08:23:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4GRSVolume Implementation
//
// ********************************************************************

#include "G4GRSVolume.hh"

G4GRSVolume::~G4GRSVolume()
{
  delete frot;			// safe if null
}

G4GRSVolume::G4GRSVolume(const G4GRSVolume& right)
{
  if ((&right) && (&right != this))
  {
    fvol = right.fvol;
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

G4GRSVolume& G4GRSVolume::operator=(const G4GRSVolume& right)
{
  if (&right == this) return *this;
  if (&right)
  {
    fvol = right.fvol;
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
