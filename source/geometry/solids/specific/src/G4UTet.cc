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
// * This  code  implementation is the  intellectual property  of the *
// * Vanderbilt University Free Electron Laser Center                 *
// * Vanderbilt University, Nashville, TN, USA                        *
// * Development supported by:                                        *
// * United States MFEL program  under grant FA9550-04-1-0045         *
// * and NASA under contract number NNG04CT05P                        *
// * Written by Marcus H. Mendenhall and Robert A. Weller.            *
// *                                                                  *
// * Contributed to the Geant4 Core, January, 2005.                   *
// *                                                                  *
// ********************************************************************
//
// $Id:$
//
// 
// Implementation for G4UTet wrapper class
// --------------------------------------------------------------------

#include "G4Tet.hh"
#include "G4UTet.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor - create a tetrahedron
// This class is implemented separately from general polyhedra,
// because the simplex geometry can be computed very quickly,
// which may become important in situations imported from mesh generators,
// in which a very large number of G4Tets are created.
// A Tet has all of its geometrical information precomputed
//
G4UTet::G4UTet(const G4String& pName,
                     G4ThreeVector anchor,
                     G4ThreeVector p2,
                     G4ThreeVector p3,
                     G4ThreeVector p4, G4bool* degeneracyFlag)
  : G4USolid(pName, new UTet(pName,
                             UVector3(anchor.x(),anchor.y(),anchor.z()),
                             UVector3(p2.x(), p2.y(), p2.z()),
                             UVector3(p3.x(), p3.y(), p3.z()),
                             UVector3(p4.x(), p4.y(), p4.z()),
                             degeneracyFlag))
{
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UTet::G4UTet( __void__& a )
  : G4USolid(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UTet::~G4UTet()
{
}

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UTet::G4UTet(const G4UTet& rhs)
  : G4USolid(rhs)
{
}


///////////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UTet& G4UTet::operator = (const G4UTet& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4USolid::operator=(rhs);

   return *this;
}
