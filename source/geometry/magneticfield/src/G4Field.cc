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
// First implementation class for G4Field 
//  J. Apostolakis,  4 Nov 2011 - to add fGravityActive data member
// -------------------------------------------------------------------

#include "G4Field.hh"

G4Field::G4Field( G4bool gravityOn):
  fGravityActive( gravityOn )
{
}
 
G4Field::~G4Field()
{
}

G4Field& G4Field::operator = (const G4Field &p)
{
   if (&p == this) return *this;
   fGravityActive= p.fGravityActive;
   return *this;
}

G4Field::G4Field (const G4Field &p)
  : fGravityActive(p.fGravityActive)
{
}

G4Field* G4Field::Clone() const
{
    G4ExceptionDescription msg;
    msg << "Derived class does not implement cloning,\n"
        << "but Clone method called.\n"
        << "Cannot continue;";
    G4Exception("G4Field::Clone", "GeomField004", FatalException,msg );
    return NULL;
}
// ------------------------------------------------------------------------
