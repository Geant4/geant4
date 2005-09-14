//
// File name:     RadmonDetectorEnvironmentLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorEnvironmentLayout.cc,v 1.2 2005-09-14 12:28:31 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorEnvironmentLayout.hh"
#include "RadmonDetectorDumpStyle.hh"
#include <iomanip>



void                                            RadmonDetectorEnvironmentLayout :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 if (!enabled)
  out << indent << "Environment disabled\n";
 else
 {
  G4int width(RADMONDETECTORDUMP_INDENT_WIDTH-indent.length());
  if (width<0)
   width=0;
  
  out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Type" << " = \""  << environmentType << "\"\n"
      << indent << "Attributes:\n";

  G4String indent2(indent);
  indent2.prepend("  ");

  DumpAttributesLayout(out, indent2);
 }
}  
