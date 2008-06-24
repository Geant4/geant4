#ifndef _G4STEPREAD_INCLUDED_
#define _G4STEPREAD_INCLUDED_

#include "G4TessellatedSolid.hh"
#include "G4QuadrangularFacet.hh"
#include "G4TriangularFacet.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"

#include <fstream>

class G4STEPRead {
   G4TessellatedSolid* ReadTessellatedSolid(const std::string&);
   G4VFacet* ReadFacet(const std::string&);
public:
   void ReadGeom(const G4String&);
};

#endif
