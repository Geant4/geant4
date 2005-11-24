//
// File name:     RadmonDataAnalysisLayout.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDataAnalysisLayout.cc,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDataAnalysisLayout.hh"
#include "RadmonDumpStyle.hh"
#include "G4UnitsTable.hh"

#include <iomanip>



void                                            RadmonDataAnalysisLayout :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 G4int width(RADMONDUMP_INDENT_WIDTH-indent.length());
 if (width<0)
  width=0;

 out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Label"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << dataAnalysisLabel << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Type";  out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << dataAnalysisType << "\"\n"
     << indent << "Attributes:\n";

 G4String indent2(indent);
 indent2.prepend("  ");

 DumpAttributesLayout(out, indent2);
}  
