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
// $Id: G4CollisionComposite.hh 82960 2014-07-21 09:19:17Z gcosmo $
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      For information related to this code contact:
//
//      File name:     G4CollisionComposite
//
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4CollisionComposite_h
#define G4CollisionComposite_h

#include "globals.hh"
#include "G4HadronicException.hh"
#include "G4VCollision.hh"
#include "G4CollisionVector.hh"
#include "G4KineticTrackVector.hh"
#include "G4CrossSectionBuffer.hh"
#include "G4Pair.hh"
#include "G4ParticleTable.hh"
#include "G4Threading.hh"

class G4KineticTrack;
class G4VCrossSectionSource;

class G4CollisionComposite : public G4VCollision
{
public:

  G4CollisionComposite();
  virtual ~G4CollisionComposite();

  virtual G4double CrossSection(const G4KineticTrack& trk1, 
				const G4KineticTrack& trk2) const;

  virtual G4KineticTrackVector* FinalState(const G4KineticTrack& trk1, 
					      const G4KineticTrack& trk2) const;

  virtual G4bool IsInCharge(const G4KineticTrack& trk1, 
			    const G4KineticTrack& trk2) const;
  void AddComponent(G4VCollision * aC) {components.push_back(aC);}

public:
  virtual const G4VCrossSectionSource* GetCrossSectionSource() const { return 0; }
  virtual const G4VAngularDistribution* GetAngularDistribution() const { return 0; }

  virtual const G4CollisionVector* GetComponents() const  { return &components;}
  struct Register
  {
    template <class T> void operator()(T*, G4CollisionComposite * aC)
    {
      aC->AddComponent(new T());
    }
  };
  struct Resolve
  {
//     template <class t1, int t2, int t3, int t4, int t5> 
    template <class T> 
    void operator()(T * , G4CollisionComposite * aC)
    {
      G4ParticleDefinition * p2, *p3, *p4, *p5;
      G4int pdg = 0;
      pdg = T::i1;
      p2=G4ParticleTable::GetParticleTable()->FindParticle(pdg);
      pdg = T::i2;
      p3=G4ParticleTable::GetParticleTable()->FindParticle(pdg);
      pdg = T::i3;
      p4=G4ParticleTable::GetParticleTable()->FindParticle(pdg);
      pdg = T::i4;
      p5=G4ParticleTable::GetParticleTable()->FindParticle(pdg);
      if(p2->GetPDGCharge()+p3->GetPDGCharge() != p4->GetPDGCharge()+p5->GetPDGCharge())
      {
        G4cerr << "charge-unbalance in collision composite"<<G4endl;
      }
      aC->AddComponent(new typename T::it(p2, p3, p4, p5));  
    }
  };

private:  

  G4CollisionComposite(const G4CollisionComposite &right);

  const G4CollisionComposite& operator=(const G4CollisionComposite &right);
  void BufferCrossSection(const G4ParticleDefinition * aP, const G4ParticleDefinition * bP);
  G4double BufferedCrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;

private:

  G4CollisionVector components;
  std::vector<G4CrossSectionBuffer> theBuffer;
  
  static const G4int nPoints;
  static const G4double theT[];
  
  mutable G4Mutex bufferMutex; //Protects concurrent access to theBuffer in MT
};

#endif
