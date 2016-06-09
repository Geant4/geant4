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

