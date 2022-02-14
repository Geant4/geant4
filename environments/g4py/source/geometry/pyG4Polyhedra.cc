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
// ====================================================================
//   pyG4Polyhedra.cc
//
//                                         2007 Q
// ====================================================================
#include <boost/python.hpp>
#include <memory>
#include "G4Polyhedra.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4Polyhedra {

// create solid methods

G4Polyhedra* f1_CreatePolyhedra(const G4String& name,
                                G4double phiStart, G4double phiTotal,
                                G4int numSide, G4int numZPlanes,
                                const std::vector<G4double>& zPlane,
                                const std::vector<G4double>& rInner,
                                const std::vector<G4double>& rOuter)
{
  std::unique_ptr<G4double[]> zlist(new G4double[numZPlanes]);
  std::unique_ptr<G4double[]> r0list(new G4double[numZPlanes]);
  std::unique_ptr<G4double[]> r1list(new G4double[numZPlanes]);

  for (G4int i=0; i< numZPlanes; i++) {
    zlist[i]= zPlane[i];
    r0list[i]= rInner[i];
    r1list[i]= rOuter[i];
  }

  return new G4Polyhedra(name, phiStart, phiTotal, numSide, numZPlanes,
                        zlist.get(), r0list.get(), r1list.get());
}


G4Polyhedra* f2_CreatePolyhedra(const G4String& name,
                                G4double phiStart, G4double phiTotal,
                                G4int numSide, G4int numRZ,
                                const std::vector<G4double>& r,
                                const std::vector<G4double>& z)
{
  std::unique_ptr<G4double[]> zlist(new G4double[numRZ]);
  std::unique_ptr<G4double[]> rlist(new G4double[numRZ]);

  for (G4int i=0; i< numRZ; i++) {
    zlist[i]= z[i];
    rlist[i]= r[i];
  }

  return new G4Polyhedra(name, phiStart, phiTotal, numSide, numRZ,
                         rlist.get(), zlist.get());

}

}

using namespace pyG4Polyhedra;

// ====================================================================
// module definition
// ====================================================================
void export_G4Polyhedra()
{
  class_<G4Polyhedra, G4Polyhedra*, bases<G4VSolid> >
    ("G4Polyhedra", "Polyhedra solid class", no_init)
    // ---
    .def("GetStartPhi",    &G4Polyhedra::GetStartPhi)
    .def("GetEndPhi",      &G4Polyhedra::GetEndPhi)
    .def("GetNumSide",     &G4Polyhedra::GetNumSide)
    .def("GetNumRZCorner", &G4Polyhedra::GetNumRZCorner)
    .def("IsOpen",         &G4Polyhedra::IsOpen)
    .def("IsGeneric",      &G4Polyhedra::IsGeneric)

    // operators
    .def(self_ns::str(self))
    ;

  // Create solid
  def("CreatePolyhedra", f1_CreatePolyhedra,
      return_value_policy<manage_new_object>());
  def("CreatePolyhedra", f2_CreatePolyhedra,
      return_value_policy<manage_new_object>());
}

