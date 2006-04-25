// $Id: pymodG4particles.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4particles.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4ParticleDefinition();
void export_G4DynamicParticle();
void export_G4ParticleTable();
void export_G4DecayTable();
void export_G4PrimaryParticle();
void export_G4PrimaryVertex();
void export_PyG4ParticleList();


BOOST_PYTHON_MODULE(G4particles)
{
  export_G4ParticleDefinition();
  export_G4DynamicParticle();
  export_G4ParticleTable();
  export_G4DecayTable();
  export_G4PrimaryParticle();
  export_G4PrimaryVertex();
  export_PyG4ParticleList();
}

