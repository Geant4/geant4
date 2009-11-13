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
// $Id: G4BREPSolidOfRevolution.hh,v 1.1 2009-11-13 14:29:52 gcamelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidOfRevolution
//
// Class description:
//
//  Definition of a generic BREP solid of revolution:
//  the 2d curve defined on x,y plane is rotated around z-axis in 3D
//  space to generate a surface of revolution closed with two planes
//
//  G4BREPSolidOfRevolution(const G4String& pName,
//                          const G4Curve& pCurve)
//
// Author: G.Camellini.
// Revisions by:
// ----------------------------------------------------------------------
#ifndef __G4BREPSolidOfRevolution
#define __G4BREPSolidOfRevolution

#include "G4BREPSolid.hh"

class G4BREPSolidOfRevolution : public G4BREPSolid
{

 public: // with description

  G4BREPSolidOfRevolution(const G4String& pName,
                          G4Curve* pCurve);
    // Constructor defining the solid through surface of revolution and
    // two planes.

  ~G4BREPSolidOfRevolution();
    // Empty destructor.

  virtual std::ostream& StreamInfo(std::ostream& os) const;
    // Streams solid contents to output stream.

 public:  // without description

  G4BREPSolidOfRevolution(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

 private:

    G4BREPSolidOfRevolution(const G4BREPSolidOfRevolution&);
    G4BREPSolidOfRevolution& operator=(const G4BREPSolidOfRevolution&);
      // Private copy constructor and assignment operator.

//  struct G4BREPSolidOfRevolutionParams
//  {
//  } fConstructorParams;
};

#endif
