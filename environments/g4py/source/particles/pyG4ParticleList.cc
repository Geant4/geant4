// $Id: pyG4ParticleList.cc,v 1.1 2006-04-25 08:13:51 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ParticleList.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ParticleTable.hh"

using namespace boost::python;

// ====================================================================
// internal class
// ====================================================================

class PyG4ParticleList {
public:
  typedef std::vector<G4ParticleDefinition*> ParticleList;
  typedef ParticleList::iterator p_iterator;

  static ParticleList particleTableCache;

  p_iterator p_begin() {
    G4ParticleTable* particleTable= G4ParticleTable::GetParticleTable();
    if(particleTableCache.size() != particleTable-> size() ) {
      particleTableCache.clear();
      G4ParticleTable::G4PTblDicIterator* 
	theParticleIterator= particleTable-> GetIterator();
      theParticleIterator-> reset();
      while( (*theParticleIterator)() ){
	G4ParticleDefinition* particle= theParticleIterator-> value();
	particleTableCache.push_back(particle);
      }
    }
    return particleTableCache.begin();
  }

  p_iterator p_end() {
    G4ParticleTable* particleTable= G4ParticleTable::GetParticleTable();
    if(particleTableCache.size() != particleTable-> size() ) {
      particleTableCache.clear();
      G4ParticleTable::G4PTblDicIterator* 
	theParticleIterator= particleTable-> GetIterator();
      theParticleIterator-> reset();
      while( (*theParticleIterator)() ){
	G4ParticleDefinition* particle= theParticleIterator-> value();
	particleTableCache.push_back(particle);
      }
    }
    return particleTableCache.end();
  }
};

PyG4ParticleList::ParticleList PyG4ParticleList::particleTableCache;


// ====================================================================
// module definition
// ====================================================================
void export_PyG4ParticleList()
{
  class_<PyG4ParticleList>("PyG4ParticleList", "particle list")
    .def("__iter__",  iterator<PyG4ParticleList::ParticleList>())
    .add_property("particles", range(&PyG4ParticleList::p_begin, 
				     &PyG4ParticleList::p_end))
    ;
}

