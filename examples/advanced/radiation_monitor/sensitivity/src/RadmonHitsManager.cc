//
// File name:     RadmonHitsManager.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonHitsManager.cc,v 1.1 2005-11-24 02:31:47 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonHitsManager.hh"
#include "RadmonHit.hh"

                                                RadmonHitsManager :: RadmonHitsManager()
:
 idCounter(0),
 allocator(new Allocator)
{
}


 
 
 
G4int                                           RadmonHitsManager :: ReserveIdByLabel(const G4String & label)
{
 const IdLabels::const_iterator i(idLabels.find(label));
 
 if (i!=idLabels.end())
  return i->second;
  
 G4int id(idCounter);
 idCounter++;
 
 idLabels[label]=id;
 return id;
}

RadmonHitsManager *                             RadmonHitsManager :: instance(0);
