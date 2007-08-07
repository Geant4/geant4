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
// $Id: pyG4IntersctionSolid.cc,v 1.1 2007-08-07 03:35:53 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4IntersectionSolid.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4IntersectionSolid.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4IntersectionSolid()
{
  class_<G4IntersectionSolid,G4IntersectionSolid*,
    bases<G4BooleanSolid>,boost::noncopyable>
    ("G4IntersectionSolid", "intersection solid class", no_init)
    // ---
    .def(init<const G4String&, G4VSolid*, G4VSolid*>())
    .def(init<const G4String&, G4VSolid*, G4VSolid*,
         G4RotationMatrix*, const G4ThreeVector&>())
    .def(init<const G4String&, G4VSolid*, G4VSolid*,
         const G4Transform3D&>())   
    ;
}

