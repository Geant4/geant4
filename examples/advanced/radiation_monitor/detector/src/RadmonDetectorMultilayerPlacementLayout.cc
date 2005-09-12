//
// File name:     RadmonDetectorMultilayerPlacementLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementLayout.cc,v 1.1 2005-09-12 17:13:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorMultilayerPlacementLayout.hh"
#include "RadmonDetectorDumpStyle.hh"
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
 size_t width(RADMONDETECTORDUMPWIDTH-indent.length());
 out << indent << std::setw(width) << "Label" << " = \"" << placementLabel << "\"\n"
     << indent << std::setw(width) << "Multilayer label" << " = \"" << multilayerLabel << "\"\n"
     << indent << std::setw(width) << "Position" << " = (" << multilayerPosition.getX()/mm << " mm, " << multilayerPosition.getY()/mm << " mm, " << multilayerPosition.getZ()/mm << " mm)\n"
     << indent << std::setw(width) << "Rotation" << " = |" << multilayerRotation.xx() << ", " << multilayerRotation.xy() << ", " << multilayerRotation.xz() << "|\n"
     << indent << std::setw(width+4) << std::setiosflags(std::ostream::right) << " = |" << multilayerRotation.yx() << ", " << multilayerRotation.yy() << ", " << multilayerRotation.yz() << "|\n"
     << indent << std::setw(width+4) << std::setiosflags(std::ostream::right) << " = |" << multilayerRotation.zx() << ", " << multilayerRotation.zy() << ", " << multilayerRotation.zz() << "|\n";
}
