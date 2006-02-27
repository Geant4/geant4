// $Id: pyG4ParticleTable.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ParticleTable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ParticleTable.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4ParticleTable {

// contains...
G4bool(G4ParticleTable::*f1_contains)(const G4ParticleDefinition*)
  = &G4ParticleTable::contains;

G4bool(G4ParticleTable::*f2_contains)(const G4String&)
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

// FindIon
G4ParticleDefinition*(G4ParticleTable::*f1_FindIon)(G4int, G4int, G4double)
  = &G4ParticleTable::FindIon;

G4ParticleDefinition*(G4ParticleTable::*f2_FindIon)(
G4int, G4int, G4int, G4int)= &G4ParticleTable::FindIon;

// DumpTable
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_DumpTable, DumpTable, 0, 1);

};

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
    .def("FindIon",           f1_FindIon,
         return_value_policy<reference_existing_object>())
    .def("FindIon",           f2_FindIon,
         return_value_policy<reference_existing_object>())
    .def("GetIon",            &G4ParticleTable::GetIon,
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
    ;
}
