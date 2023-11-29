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

  if (xSource != nullptr)
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
    nComponents = (G4int)components->size();
  }
  G4cout << "---- " << name << "---- has " << nComponents << " components" <<G4endl;
  G4int i = 0;
  if (components)
  {
    for (auto iter = components->cbegin(); iter != components->cend(); ++iter)
    {
      G4cout << "---- " << name << " ---- Component " << i << G4endl;
      ((*iter))->Print();
      ++i;
    }
  }
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
    nComponents = (G4int)components->size();
  }
  G4cout << "---- " << name << "has " << nComponents << " components" <<G4endl;

  G4int i = 0;
  if (components)
  {
    for (auto iter = components->cbegin(); iter != components->cend(); ++iter)
    {
      G4cout << "Component " << i << G4endl;
      ((*iter))->Print();
      ++i;
    }
  }
}

G4VCollision::G4VCollision(void*, void*, void *, void *, void *, void *, void *)
{
}

void G4VCollision::establish_G4MT_TLS_G4VCollision()
{
}
