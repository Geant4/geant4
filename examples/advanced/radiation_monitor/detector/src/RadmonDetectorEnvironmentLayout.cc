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
// File name:     RadmonDetectorEnvironmentLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorEnvironmentLayout.cc,v 1.4 2006-06-28 13:50:18 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorEnvironmentLayout.hh"
#include "RadmonDumpStyle.hh"
#include <iomanip>



void                                            RadmonDetectorEnvironmentLayout :: DumpLayout(std::ostream & out, const G4String & indent) const
{
 if (!enabled)
  out << indent << "Environment disabled\n";
 else
 {
  G4int width(RADMONDUMP_INDENT_WIDTH-indent.length());
  if (width<0)
   width=0;
  
  out << indent << std::setw(width); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Type" << " = \""  << environmentType << "\"\n"
      << indent << "Attributes:\n";

  G4String indent2(indent);
  indent2.prepend("  ");

  DumpAttributesLayout(out, indent2);
 }
}  
