//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonDetectorMultilayerPlacementLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMultilayerPlacementLayout.cc,v 1.5 2006/06/29 16:14:05 gunter Exp $
// Tag:           $Name: geant4-08-01 $
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
