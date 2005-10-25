//
// File name:     RadmonGeneratorSourceAlgorithmLayout.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorSourceAlgorithmLayout.cc,v 1.1 2005-10-25 16:36:43 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonGeneratorSourceAlgorithmLayout.hh"
#include "RadmonDumpStyle.hh"

#include <iomanip>



                                                RadmonGeneratorSourceAlgorithmLayout :: RadmonGeneratorSourceAlgorithmLayout(const RadmonGeneratorSourceAlgorithmLayout & copy)
:
 RadmonLayoutEntityWithAttributes(copy),
 algorithmLabel(copy.algorithmLabel),
 algorithmType(copy.algorithmType)
{
}
 
 
 
 
 
RadmonGeneratorSourceAlgorithmLayout &          RadmonGeneratorSourceAlgorithmLayout :: operator=(const RadmonGeneratorSourceAlgorithmLayout & copy)
{
 RadmonLayoutEntityWithAttributes::operator=(copy);
 algorithmLabel=copy.algorithmLabel;
 algorithmType=copy.algorithmType;
 
 return (*this);
}





void                                            RadmonGeneratorSourceAlgorithmLayout :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 G4int width(RADMONDUMP_INDENT_WIDTH-indent.length());
 if (width<0)
  width=0;

 out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Label";     out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << algorithmLabel << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Type";      out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << algorithmType << "\"\n"
     << indent << "Attributes:\n";

 G4String indent2(indent);
 indent2.prepend("  ");

 DumpAttributesLayout(out, indent2);
}
