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
// $Id: pyG4ParticleTable.cc 86749 2014-11-17 15:03:05Z gcosmo $
// ====================================================================
//   pyG4ParticleTable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4ParticleTable.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4ParticleTable {

// contains...
G4bool(G4ParticleTable::*f1_contains)(const G4ParticleDefinition*) const
  = &G4ParticleTable::contains;

G4bool(G4ParticleTable::*f2_contains)(const G4String&) const
  = &G4ParticleTable::contains;

// FindParticle...
G4ParticleDefinition*(G4ParticleTable::*f1_FindParticle)(G4int)
  = &G4ParticleTable::FindParticle;

G4ParticleDefinition*(G4ParticleTable::*f2_FindParticle)(const G4String&)
  = &G4ParticleTable::FindParticle;

G4ParticleDefinition*(G4ParticleTable::*f3_FindParticle)(
  const G4ParticleDefinition*)= &G4ParticleTable::FindParticle;

// FindAntiParticle...
G4ParticleDefinition*(G4ParticleTable::*f1_FindAntiParticle)(G4int)
  = &G4ParticleTable::FindAntiParticle;

G4ParticleDefinition*(G4ParticleTable::*f2_FindAntiParticle)(const G4String&)
  = &G4ParticleTable::FindAntiParticle;

G4ParticleDefinition*(G4ParticleTable::*f3_FindAntiParticle)(
  const G4ParticleDefinition*)= &G4ParticleTable::FindAntiParticle;

// DumpTable
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_DumpTable, DumpTable, 0, 1)


// --------------------------------------------------------------------
// GetParticleList (returning python list)

list GetParticleList(G4ParticleTable* particleTable)
{
  list particleList;
  G4ParticleTable::G4PTblDicIterator*
    theParticleIterator= particleTable-> GetIterator();
  theParticleIterator-> reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle= theParticleIterator-> value();
    particleList.append(&particle);
  }

  return particleList;
}

}

using namespace pyG4ParticleTable;

// ====================================================================
// module definition
// ====================================================================
void export_G4ParticleTable()
{
  class_<G4ParticleTable, G4ParticleTable*, boost::noncopyable>
    ("G4ParticleTable", "particle table", no_init)
    // ---
    .def("GetParticleTable",  &G4ParticleTable::GetParticleTable,
         return_value_policy<reference_existing_object>())
    .staticmethod("GetParticleTable")
    .def("contains",          f1_contains)
    .def("contains",          f2_contains)
    .def("entries",           &G4ParticleTable::entries)
    .def("size",              &G4ParticleTable::size)
    // ---
    .def("GetParticle",       &G4ParticleTable::GetParticle,
         return_value_policy<reference_existing_object>())
    .def("GetParticleName",   &G4ParticleTable::GetParticleName,
	 return_value_policy<return_by_value>())
    .def("FindParticle",      f1_FindParticle,
         return_value_policy<reference_existing_object>())
    .def("FindParticle",      f2_FindParticle,
         return_value_policy<reference_existing_object>())
    .def("FindParticle",      f3_FindParticle,
         return_value_policy<reference_existing_object>())
    .def("FindAntiParticle",  f1_FindAntiParticle,
         return_value_policy<reference_existing_object>())
    .def("FindAntiParticle",  f2_FindAntiParticle,
         return_value_policy<reference_existing_object>())
    .def("FindAntiParticle",  f3_FindAntiParticle,
         return_value_policy<reference_existing_object>())
    .def("DumpTable",         &G4ParticleTable::DumpTable, f_DumpTable())
    //.def("GetIonTable",     &G4ParticleTable::GetIonTable,
    //...)
    //.def("GetShortLivedTable", &G4ParticleTable::GetShortLivedTable,
    //...)
    .def("SetVerboseLevel",   &G4ParticleTable::SetVerboseLevel)
    .def("GetVerboseLevel",   &G4ParticleTable::GetVerboseLevel)
    .def("SetReadiness",      &G4ParticleTable::SetReadiness)
    .def("GetReadiness",      &G4ParticleTable::GetReadiness)
    // ---
    // additionals
    .def("GetParticleList",   GetParticleList)
    ;
}
