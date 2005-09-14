//
// File name:     RadmonDetectorLayerLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerLayout.cc,v 1.2 2005-09-14 12:28:31 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorLayerLayout.hh"
#include "RadmonDetectorDumpStyle.hh"
#include "G4UnitsTable.hh"

#include <iomanip>



void                                            RadmonDetectorLayerLayout :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 G4int width(RADMONDETECTORDUMP_INDENT_WIDTH-indent.length());
 if (width<0)
  width=0;

 out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Label";     out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << layerLabel << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Type";      out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << layerType << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Thickness"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = "   << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << G4BestUnit(layerThickness, "Length") << '\n'
     << indent << "Attributes:\n";

 G4String indent2(indent);
 indent2.prepend("  ");

 DumpAttributesLayout(out, indent2);
}  
