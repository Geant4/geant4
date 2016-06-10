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
// $Id: pyG4Polycone.cc 81291 2014-05-26 09:31:19Z gcosmo $
// ====================================================================
//   pyG4Polycone.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Polycone.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4Polycone {

// create solid methods

G4Polycone* f1_CreatePolycone(const G4String& name, G4double phiStart,
                              G4double phiTotal, G4int numZPlanes,
                              const std::vector<G4double>& zPlane,
                              const std::vector<G4double>& rInner,
                              const std::vector<G4double>& rOuter)
{
  G4double zlist[numZPlanes];
  G4double r0list[numZPlanes];
  G4double r1list[numZPlanes];

  for (G4int i=0; i< numZPlanes; i++) {
    zlist[i]= zPlane[i];
    r0list[i]= rInner[i];
    r1list[i]= rOuter[i];
  }

  return new G4Polycone(name, phiStart, phiTotal, numZPlanes,
                        zlist, r0list, r1list);
}


G4Polycone* f2_CreatePolycone(const G4String& name, G4double phiStart,
                              G4double phiTotal, G4int numRZ,
                              const std::vector<G4double>& r,
                              const std::vector<G4double>& z)
{
  G4double zlist[numRZ];
  G4double rlist[numRZ];

  for (G4int i=0; i< numRZ; i++) {
    rlist[i]= r[i];
    zlist[i]= z[i];
  }

  return new G4Polycone(name, phiStart, phiTotal, numRZ,
                        rlist, zlist);

}

}

using namespace pyG4Polycone;


// ====================================================================
// module definition
// ====================================================================
void export_G4Polycone()
{
  class_<G4Polycone, G4Polycone*, bases<G4VSolid> >
    ("G4Polycone", "Polycone solid class", no_init)
    // ---
    .def("GetStartPhi",    &G4Polycone::GetStartPhi)
    .def("GetEndPhi",      &G4Polycone::GetEndPhi)
    .def("IsOpen",         &G4Polycone::IsOpen)
    .def("GetNumRZCorner", &G4Polycone::GetNumRZCorner)

    // operators
    .def(self_ns::str(self))
    ;

  // Create solid
  def("CreatePolycone", f1_CreatePolycone,
      return_value_policy<manage_new_object>());
  def("CreatePolycone", f2_CreatePolycone,
      return_value_policy<manage_new_object>());

}

