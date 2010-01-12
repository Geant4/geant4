#ifndef G4INUCL_ELEMENTARY_PARTICLE_HH
#define G4INUCL_ELEMENTARY_PARTICLE_HH
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
// $Id: G4InuclElementaryParticle.hh,v 1.17 2010-01-12 06:27:15 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $

#include "G4InuclParticle.hh"
#include "globals.hh"
#include "G4Allocator.hh"

class G4ParticleDefinition;

class G4InuclElementaryParticle : public G4InuclParticle {

//                     known particle types:
//      1 - proton          11 - k+         111 - quasideuteron PP
//      2 - neutron         13 - k-         112 - quasideuteron PN
//      3 - pi+             15 - k0         122 - quasideuteron NN
//      5 - pi-             17 - k0bar
//      7 - pi 0            21 - lambda 
//     10 - photon          23 - sigma+
//                          25 - sigma0
//                          27 - sigma-
//                          29 - xi0
//                          31 - xi-

public:
  G4InuclElementaryParticle() 
    : G4InuclParticle("InuclElemPart"), generation(0) {}

  explicit G4InuclElementaryParticle(G4int type) 
    : G4InuclParticle("InuclElemPart", makeDefinition(type)), generation(0) {}

  G4InuclElementaryParticle(const G4CascadeMomentum& mom,
			    G4int type, G4int model=0) 
    : G4InuclParticle("InuclElemPart", makeDefinition(type), mom),
      generation(0) {
    setModel(model);
  }

  G4InuclElementaryParticle(G4double ekin, G4int type) 
    : G4InuclParticle("InuclElemPart", makeDefinition(type), ekin),
      generation(0) {}

  /******
  //  new/delete operators are overridden to use G4Allocator
  inline void *operator new(size_t);
  inline void operator delete(void *inuclElemPart);
  ******/

  // Copy and assignment constructors for use with std::vector<>
  G4InuclElementaryParticle(const G4InuclElementaryParticle& right)
    : G4InuclParticle(right), generation(right.generation) {}

  G4InuclElementaryParticle& operator=(const G4InuclElementaryParticle& right);

  void setType(G4int ityp);
  G4int type() const;

  G4bool photon() const { return (type() == 10); }

  G4bool nucleon() const { return (type() <= 2); }

  G4bool baryon() const { 
    return (getDefinition()->GetBaryonNumber() != 0);
  }

  G4bool pion() const { 
    return (type() == 3 || type() == 5 || type() == 7); 
  }

  G4bool quasi_deutron() const { return (type() > 100); }

  G4bool valid() const { return type()>0; }

  virtual void printParticle() const {
    G4InuclParticle::printParticle();

    G4cout << " Particle: type " << type() << " mass " << getMass()
	   << " ekin " << getKineticEnergy() << G4endl; 
  }

  void setGeneration(G4int gen) { generation = gen; }

  G4int getGeneration() { return generation; }

  static G4double getStrangeness(G4int type);
  static G4double getParticleMass(G4int type);

protected:
  // Convert internal type code to standard GEANT4 pointer
  static G4ParticleDefinition* makeDefinition(G4int ityp);

private: 
  G4int generation;
};        

/******
//  new/delete operators are overloaded to use G4Allocator

#ifdef G4INUCL_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<G4InuclElementaryParticle> anInuclElemPartAllocator;
#else
  extern G4DLLIMPORT G4Allocator<G4InuclElementaryParticle> anInuclElemPartAllocator;
#endif

inline void *G4InuclElementaryParticle::operator new(size_t) {
  void* temp = anInuclElemPartAllocator.MallocSingle();
  G4cout << "G4InuclElemPart::new returning @ " << temp << G4endl;
  return temp;
}

inline void G4InuclElementaryParticle::operator delete(void *inuclElemPart) {
  G4cout << "G4InuclElemPart::delete @ " << inuclElemPart << G4endl;
  anInuclElemPartAllocator.FreeSingle((G4InuclElementaryParticle*)inuclElemPart);
}
******/

#endif // G4INUCL_ELEMENTARY_PARTICLE_HH 
