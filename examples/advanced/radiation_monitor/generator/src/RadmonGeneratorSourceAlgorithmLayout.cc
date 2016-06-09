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
// File name:     RadmonGeneratorSourceAlgorithmLayout.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorSourceAlgorithmLayout.cc,v 1.3 2006/06/29 16:16:29 gunter Exp $
// Tag:           $Name: geant4-09-02 $
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
