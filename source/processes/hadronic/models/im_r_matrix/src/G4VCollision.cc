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
// $Id: G4VCollision.cc,v 1.1 2003-10-07 12:37:40 hpw Exp $ //

#include "globals.hh"
#include "G4ios.hh"
#include "G4VCollision.hh"
#include "G4VCrossSectionSource.hh"
#include "G4KineticTrack.hh"

G4VCollision::G4VCollision()
{ }

G4VCollision::~G4VCollision()
{
}


G4bool G4VCollision::operator==(const G4VCollision &right) const
{
  return (this == (G4VCollision *) &right);
}


G4bool G4VCollision::operator!=(const G4VCollision &right) const
{
  return (this != (G4VCollision *) &right);
}


G4double G4VCollision::CrossSection(const G4KineticTrack& aTrk1, 
				    const G4KineticTrack& aTrk2) const
{
  G4double sigma = 0.;
  
  const G4VCrossSectionSource* xSource = GetCrossSectionSource();

  if (xSource != 0)
    {
      // There is a cross section for this Collision
      sigma = xSource->CrossSection(aTrk1,aTrk2);
    }
  return sigma;
}


void G4VCollision::Print() const
{
  G4String name = GetName();

  G4cout << "---- " << name << "---- Cross section" << G4endl;

   const G4VCrossSectionSource* xSource = GetCrossSectionSource();
 if (xSource) xSource->Print();

  G4int nComponents = 0;
  const G4CollisionVector* components = GetComponents();
  if (components)
    {
      nComponents = components->size();
    }
  G4cout << "---- " << name << "---- has " << nComponents << " components" <<G4endl;
  G4int i = 0;
  G4CollisionVector::const_iterator iter;
  for (iter = components->begin(); iter != components->end(); ++iter)
    {
      G4cout << "---- " << name << " ---- Component " << i << G4endl;
      ((*iter)())->Print();
      i++;
    }
  
//  G4int i;
//  for (i=0; i<nComponents; i++)
//    {
//      G4VCollision* component = (*components)[i];
//    G4cout << "---- " << name << " ---- Component " << i << G4endl;
//      component->Print();
//    }
}


void G4VCollision::Print(const G4KineticTrack& trk1, 
			 const G4KineticTrack& trk2) const
{
  G4String name = GetName();
  
  if (IsInCharge(trk1,trk2))
    {    
      G4cout << "---- " << name << "is in charge ---- " << G4endl;
    }
  else
    {    
      G4cout << "---- " << name << "is not in charge ---- " << G4endl;
    }
  
  G4cout << "---- " << name << "---- Cross section" << G4endl;
  
  const G4VCrossSectionSource* xSource = GetCrossSectionSource();
  if (xSource) xSource->Print();
  G4cout << "Cross section = " << CrossSection(trk1,trk2) << G4endl;
  
  G4int nComponents = 0;
  const G4CollisionVector* components = GetComponents();
  if (components)
    {
      nComponents = components->size();
    }
  G4cout << "---- " << name << "has " << nComponents << " components" <<G4endl;
  
G4int i = 0;
  G4CollisionVector::const_iterator iter;
  for (iter = components->begin(); iter != components->end(); ++iter)
    {
      G4cout << "Component " << i << G4endl;
      ((*iter)())->Print();
      i++;
    }
}


