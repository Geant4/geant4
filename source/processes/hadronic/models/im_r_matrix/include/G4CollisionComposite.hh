//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4CollisionComposite.hh,v 1.5 2003-12-12 15:38:21 hpw Exp $
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

protected:
  virtual const G4VCrossSectionSource* GetCrossSectionSource() const { return 0; }
  virtual const G4VAngularDistribution* GetAngularDistribution() const { return 0; }

  void AddComponent(G4VCollision * aC) {components.push_back(aC);}
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
    template <class t1, int t2, int t3, int t4, int t5> 
    void operator()(INT4(t1,t2,t3,t4,t5) * , G4CollisionComposite * aC)
    {
      G4ParticleDefinition * p2, *p3, *p4, *p5;
      p2=G4ParticleTable::GetParticleTable()->FindParticle(t2);
      p3=G4ParticleTable::GetParticleTable()->FindParticle(t3);
      p4=G4ParticleTable::GetParticleTable()->FindParticle(t4);
      p5=G4ParticleTable::GetParticleTable()->FindParticle(t5);
      if(p2->GetPDGCharge()+p3->GetPDGCharge() != p4->GetPDGCharge()+p5->GetPDGCharge())
      {
        G4cerr << "charge-unbalance in collision composite"<<G4endl;
      }
      aC->AddComponent(new t1(p2, p3, p4, p5));  
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
  static G4double theT[];
  
};

#endif
