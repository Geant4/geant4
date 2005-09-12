//
// File name:     RadmonDetectorEnvironmentLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorEnvironmentLayout.cc,v 1.1 2005-09-12 17:13:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorEnvironmentLayout.hh"
#include "RadmonDetectorDumpStyle.hh"
#include <iomanip>



void                                            RadmonDetectorEnvironmentLayout :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 size_t width(RADMONDETECTORDUMPWIDTH-indent.length());

 if (!enabled)
  out << indent << std::setw(width) << "Disabled\n";
 else
  out << indent << std::setw(width) << "Type" << " = \""  << environmentType << "\"\n";
}  
