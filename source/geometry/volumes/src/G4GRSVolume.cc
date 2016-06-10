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
// $Id: G4GRSVolume.cc 87615 2014-12-12 15:23:11Z gcosmo $
//
// 
// class G4GRSVolume Implementation
//
// ----------------------------------------------------------------------

#include "G4GRSVolume.hh"

G4GRSVolume::~G4GRSVolume()
{
  delete frot;             // safe if null
}

G4GRSVolume::G4GRSVolume(const G4GRSVolume& right)
  : G4VTouchable(), fvol(0)
{
  if (frot) { delete frot; }
  frot = 0;
  fvol = right.fvol;
  ftlate = right.ftlate;
  if (right.frot)
  {
    frot = new G4RotationMatrix(*(right.frot));
  }
}

G4GRSVolume& G4GRSVolume::operator=(const G4GRSVolume& right)
{
  if (&right == this)  { return *this; }

  fvol = right.fvol;
  ftlate = right.ftlate;
  if (frot)
  {
    delete frot;
    frot = 0;
  }
  if (right.frot)
  {
    frot = new G4RotationMatrix(*(right.frot));
  }

  return *this;
}
