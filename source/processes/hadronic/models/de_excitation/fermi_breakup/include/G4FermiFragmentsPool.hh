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
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#ifndef G4FermiFragmentsPool_hh 
#define G4FermiFragmentsPool_hh

#include "G4VFermiFragment.hh"
#include "G4StableFermiFragment.hh"
#include "G4B9FermiFragment.hh"
#include "G4Be8FermiFragment.hh"
#include "G4He5FermiFragment.hh"
#include "G4Li5FermiFragment.hh"

#include <map>
#include <functional>

class G4FermiFragmentsPool
{
public:
  G4FermiFragmentsPool();
  ~G4FermiFragmentsPool();


  G4int Count(const std::pair<G4int,G4int>& az) 
  {
    return GetMap().count(az);
  }
  
  std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment*,
		       std::less<const std::pair<G4int,G4int> > >::iterator 
  LowerBound(const std::pair<G4int,G4int> & az) 
  {
    return GetMap().lower_bound(az);
  }

  std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment*,
			 std::less<const std::pair<G4int,G4int> > >::iterator 
  UpperBound(const std::pair<G4int,G4int> & az) 
  {
    return GetMap().upper_bound(az);
  }

  
private:
  G4FermiFragmentsPool(const G4FermiFragmentsPool&);
  const G4FermiFragmentsPool & operator=(const G4FermiFragmentsPool&);
  G4bool operator==(const G4FermiFragmentsPool&) const;
  G4bool operator!=(const G4FermiFragmentsPool&) const;

  std::multimap<const std::pair<G4int,G4int>, const G4VFermiFragment* , 
                std::less<const std::pair<G4int,G4int> > >  &
  GetMap();


private:

  static G4bool MapIsEmpty;
  
};
#endif

