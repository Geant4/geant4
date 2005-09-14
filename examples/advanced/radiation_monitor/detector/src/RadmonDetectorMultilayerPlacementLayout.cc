//
// File name:     RadmonDetectorMultilayerPlacementLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementLayout.cc,v 1.2 2005-09-14 12:28:31 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorMultilayerPlacementLayout.hh"
#include "RadmonDetectorDumpStyle.hh"
#include "G4UnitsTable.hh"

#include <iomanip>



G4ThreeVector                                   RadmonDetectorMultilayerPlacementLayout :: GetRelativePosition(const RadmonDetectorMultilayerPlacementLayout & reference) const
{
 G4ThreeVector v(multilayerPosition);
 v-=reference.multilayerPosition;

 return reference.multilayerRotation.inverse()*v;
}



G4RotationMatrix                                RadmonDetectorMultilayerPlacementLayout :: GetRelativeRotation(const RadmonDetectorMultilayerPlacementLayout & reference) const
{
 return reference.multilayerRotation.inverse()*multilayerRotation; 
}





void                                            RadmonDetectorMultilayerPlacementLayout :: SetRelativePosition(const RadmonDetectorMultilayerPlacementLayout & reference, const G4ThreeVector & position)
{
 multilayerPosition=reference.multilayerRotation*position;
 multilayerPosition+=reference.multilayerPosition;
}



void                                            RadmonDetectorMultilayerPlacementLayout :: SetRelativeRotation(const RadmonDetectorMultilayerPlacementLayout & reference, const G4RotationMatrix & rotation)
{
 multilayerRotation=reference.multilayerRotation;
 multilayerRotation*=rotation;
}





void                                            RadmonDetectorMultilayerPlacementLayout :: DumpLayout(std::ostream & out, const G4String &indent) const
{
 G4int width(RADMONDETECTORDUMP_INDENT_WIDTH-indent.length());
 if (width<0)
  width=0;
  
 G4int width2(width+4+indent.length());
  
 out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Label";            out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << placementLabel << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Multilayer label"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << multilayerLabel << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Position";         out.setf(std::ostream::right, std::ostream::adjustfield); out << " = ("  << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << G4BestUnit(multilayerPosition.getX(), "Length") << ", " 
                                                                                                                                                                                                 << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << G4BestUnit(multilayerPosition.getY(), "Length") << ", "
                                                                                                                                                                                                 << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << G4BestUnit(multilayerPosition.getZ(), "Length") << ")\n"
               << std::setw(width2); out.setf(std::ostream::right, std::ostream::adjustfield); out <<                                                                                       "|"  << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << multilayerRotation.xx() << ", " << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << multilayerRotation.xy() << ", " << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << multilayerRotation.xz() << "|\n"
     << indent << std::setw(width);  out.setf(std::ostream::left, std::ostream::adjustfield);  out << "Rotation";       out.setf(std::ostream::right, std::ostream::adjustfield); out << " = |"  << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << multilayerRotation.yx() << ", " << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << multilayerRotation.yy() << ", " << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << multilayerRotation.yz() << "|\n"
               << std::setw(width2); out.setf(std::ostream::right, std::ostream::adjustfield); out <<                                                                                       "|"  << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << multilayerRotation.zx() << ", " << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << multilayerRotation.zy() << ", " << std::setprecision(RADMONDETECTORDUMP_DOUBLE_PRECISION) << std::setw(RADMONDETECTORDUMP_DOUBLE_WIDTH) << multilayerRotation.zz() << "|\n";
}
