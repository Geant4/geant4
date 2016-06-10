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
// $Id: G4FermiConfiguration.hh 85841 2014-11-05 15:35:06Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
// 23.04.2011 V.Ivanchenko: make this class to be a simple container and no physics

#ifndef G4FermiConfiguration_h
#define G4FermiConfiguration_h 1

#include "globals.hh"
#include "G4VFermiFragment.hh"
#include "G4Fragment.hh"

class G4FermiConfiguration 
{
public:

  G4FermiConfiguration(const std::vector<const G4VFermiFragment*>&);

  ~G4FermiConfiguration();

  inline G4int GetA() const;
  inline G4int GetZ() const;
  inline G4double GetMass() const;
  
  inline const std::vector<const G4VFermiFragment*>* GetFragmentList() const;

private:

  inline G4FermiConfiguration(const G4FermiConfiguration &);
  inline const G4FermiConfiguration & operator=(const G4FermiConfiguration &);
  inline G4bool operator==(const G4FermiConfiguration &) const;
  inline G4bool operator!=(const G4FermiConfiguration &) const;
  
  G4int    totalZ;
  G4int    totalA;

  G4double totalMass;

  std::vector<const G4VFermiFragment*> Configuration;

};

inline G4int G4FermiConfiguration::GetA() const
{
  return totalA;
}

inline G4int G4FermiConfiguration::GetZ() const
{
  return totalZ;
}

inline G4double G4FermiConfiguration::GetMass() const
{
  return totalMass;
}

inline const std::vector<const G4VFermiFragment*>* 
G4FermiConfiguration::GetFragmentList() const
{
  return &Configuration;
}

#endif


