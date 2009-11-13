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
// $Id: G4BREPSolidOfLinearExtrusion.cc,v 1.1 2009-11-13 14:29:52 gcamelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolidOfLinearExtrusion.cc
//
// ----------------------------------------------------------------------

#include "G4BREPSolidOfLinearExtrusion.hh"
#include "G4SurfaceOfLinearExtrusion.hh"

G4BREPSolidOfLinearExtrusion::G4BREPSolidOfLinearExtrusion(
                                const G4String& pName,
                                G4Curve* pCurve,
                                G4double pDz)
  : G4BREPSolid(pName)
{
  G4Curve* curve = pCurve;
  G4double dz = pDz;
  curve->GetEntityType();
  dz +=0.;
//  Initialize();
}

G4BREPSolidOfLinearExtrusion::G4BREPSolidOfLinearExtrusion( __void__& a )
  : G4BREPSolid(a)
{
}

G4BREPSolidOfLinearExtrusion::~G4BREPSolidOfLinearExtrusion()
{
}

// Streams solid contents to output stream.
std::ostream& G4BREPSolidOfLinearExtrusion::StreamInfo(std::ostream& os) const
{
  G4BREPSolid::StreamInfo( os )
//  << "\n curve:        " << constructorParams.curve
  << "\n-----------------------------------------------------------\n";

  return os;
}
