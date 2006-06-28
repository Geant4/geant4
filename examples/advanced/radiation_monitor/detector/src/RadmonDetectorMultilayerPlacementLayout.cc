//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// File name:     RadmonDetectorMultilayerPlacementLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementLayout.cc,v 1.4 2006-06-28 13:50:44 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorMultilayerPlacementLayout.hh"
#include "RadmonDumpStyle.hh"
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
 G4int width(RADMONDUMP_INDENT_WIDTH-indent.length());
 if (width<0)
  width=0;
  
 G4int width2(width+4+indent.length());
  
 out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Label";            out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << placementLabel << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Multilayer label"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << multilayerLabel << "\"\n"
     << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Position";         out.setf(std::ostream::right, std::ostream::adjustfield); out << " = ("  << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << G4BestUnit(multilayerPosition.getX(), "Length") << ", " 
                                                                                                                                                                                                 << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << G4BestUnit(multilayerPosition.getY(), "Length") << ", "
                                                                                                                                                                                                 << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << G4BestUnit(multilayerPosition.getZ(), "Length") << ")\n"
               << std::setw(width2); out.setf(std::ostream::right, std::ostream::adjustfield); out <<                                                                                       "|"  << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << multilayerRotation.xx() << ", " << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << multilayerRotation.xy() << ", " << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << multilayerRotation.xz() << "|\n"
     << indent << std::setw(width);  out.setf(std::ostream::left, std::ostream::adjustfield);  out << "Rotation";       out.setf(std::ostream::right, std::ostream::adjustfield); out << " = |"  << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << multilayerRotation.yx() << ", " << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << multilayerRotation.yy() << ", " << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << multilayerRotation.yz() << "|\n"
               << std::setw(width2); out.setf(std::ostream::right, std::ostream::adjustfield); out <<                                                                                       "|"  << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << multilayerRotation.zx() << ", " << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << multilayerRotation.zy() << ", " << std::setprecision(RADMONDUMP_DOUBLE_PRECISION) << std::setw(RADMONDUMP_DOUBLE_WIDTH) << multilayerRotation.zz() << "|\n";
}
