//
// File name:     RadmonDetectorLayerLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerLayout.cc,v 1.1 2005-09-12 17:13:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorLayerLayout.hh"
#include "RadmonDetectorDumpStyle.hh"
#include <iomanip>



void                                            RadmonDetectorLayerLayout :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 size_t width(RADMONDETECTORDUMPWIDTH-indent.length());

 out << indent << std::setw(width) << "Label" << " = \"" << layerLabel << "\"\n"
     << indent << std::setw(width) << "Type" << " = \""  << layerType << "\"\n"
     << indent << std::setw(width) << "Thickness" << " = " << std::setprecision(2) << layerThickness/mm << " mm\n"
     << indent << "Attributes:\n";

 G4String indent2(indent);
 indent2.prepend("  ");

 DumpAttributesLayout(out, indent2);
}  
