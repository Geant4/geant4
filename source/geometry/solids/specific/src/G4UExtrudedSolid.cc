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
// $Id:$
//
// 
// Implementation of G4UExtrudedSolid wrapper class
// --------------------------------------------------------------------

#include "G4ExtrudedSolid.hh"
#include "G4UExtrudedSolid.hh"
#include "G4Polyhedron.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructors
//
G4UExtrudedSolid::G4UExtrudedSolid(const G4String&          name,
                                   std::vector<G4TwoVector> polygon,
                                   std::vector<ZSection>    zsections)
  : G4USolid(name, new UExtrudedSolid())
{
  GetShape()->SetName(name);
  std::vector<UVector2> pvec;
  for (unsigned int i=0; i<polygon.size(); ++i)
  {
    pvec.push_back(UVector2(polygon[i].x(), polygon[i].y()));
  }
  std::vector<UExtrudedSolid::ZSection> svec;
  for (unsigned int i=0; i<zsections.size(); ++i)
  {
    ZSection sec = zsections[i];
    svec.push_back(UExtrudedSolid::ZSection(sec.fZ,
                   UVector2(sec.fOffset.x(), sec.fOffset.y()), sec.fScale));
  }
  GetShape()->Initialise(pvec, svec);
}


G4UExtrudedSolid::G4UExtrudedSolid(const G4String&          name,
                                   std::vector<G4TwoVector> polygon,
                                   G4double                 halfZ,
                                   G4TwoVector off1, G4double scale1,
                                   G4TwoVector off2, G4double scale2)
  : G4USolid(name, new UExtrudedSolid())
{ 
  GetShape()->SetName(name);
  std::vector<UVector2> pvec;
  for (unsigned int i=0; i<polygon.size(); ++i)
  {
    pvec.push_back(UVector2(polygon[i].x(), polygon[i].y()));
  }
  GetShape()->Initialise(pvec, halfZ, UVector2(off1.x(), off1.y()), scale1,
                                      UVector2(off2.x(), off2.y()), scale2);
}

////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UExtrudedSolid::G4UExtrudedSolid(__void__& a)
  : G4USolid(a)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UExtrudedSolid::~G4UExtrudedSolid()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UExtrudedSolid::G4UExtrudedSolid(const G4UExtrudedSolid &source)
  : G4USolid(source)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UExtrudedSolid&
G4UExtrudedSolid::operator=(const G4UExtrudedSolid &source)
{
  if (this == &source) return *this;
  
  G4USolid::operator=( source );
  
  return *this;
}
