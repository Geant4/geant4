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
// $Id: pyG4ParticleList.cc 66892 2013-01-17 10:57:59Z gunter $
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

