// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PrimaryVertex.hh,v 1.2 1999-11-05 04:16:17 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//


#ifndef G4PrimaryVertex_h
#define G4PrimaryVertex_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4PrimaryParticle.hh"

// class description:
//
//  This is the class which represents a primary vertex. The ofject of this
// class is set to G4Event objct by G4VPrimaryGenerator concrete class.
// This class object has one or more G4PrimaryParticle objects as primary
// particles.

class G4PrimaryVertex 
{
  public:
      inline void *operator new(size_t);
      inline void operator delete(void *aStackedTrack);

      G4PrimaryVertex();
      G4PrimaryVertex(G4double x0,G4double y0,G4double z0,G4double t0);
      G4PrimaryVertex(G4ThreeVector xyz0,G4double t0);
      ~G4PrimaryVertex();

      const G4PrimaryVertex & operator=(const G4PrimaryVertex &right);
      int operator==(const G4PrimaryVertex &right) const;
      int operator!=(const G4PrimaryVertex &right) const;

      void Print() const;

  private:
      G4double X0;
      G4double Y0;
      G4double Z0;
      G4double T0;
      G4PrimaryParticle * theParticle;
      G4PrimaryParticle * theTail;
      G4PrimaryVertex* nextVertex;
      G4int numberOfParticle;

  public:
      inline G4ThreeVector GetPosition() const
      { return G4ThreeVector(X0,Y0,Z0); }
      inline G4double GetX0() const
      { return X0; }
      inline G4double GetY0() const
      { return Y0; }
      inline G4double GetZ0() const
      { return Z0; }
      inline G4double GetT0() const
      { return T0; }
      inline G4int GetNumberOfParticle() const
      { return numberOfParticle; }
      inline void SetPrimary(G4PrimaryParticle * pp)
      { 
        if(theParticle == NULL)
        { 
          theParticle = pp;
          theTail = pp;
        }
        else
        { 
          theTail->SetNext(pp); 
          theTail = pp;
        }
        numberOfParticle++;
      }
      inline G4PrimaryParticle* GetPrimary(G4int i=0) const
      { 
        if( i == 0 )
        { return theParticle; }
        else if( i > 0 && i < numberOfParticle )
        {
          G4PrimaryParticle* particle = theParticle;
          for( int j=0; j<i; j++ )
          { 
            if( particle == NULL ) return NULL;
            particle = particle->GetNext();
          }
          return particle;
        }
        else
        { return NULL; }
      }
      inline void SetNext(G4PrimaryVertex* nv)
      { 
        if(nextVertex == NULL)
        { nextVertex = nv; }
        else
        { nextVertex->SetNext(nv); }
      }
      inline G4PrimaryVertex* GetNext() const
      { return nextVertex; }
};

extern G4Allocator<G4PrimaryVertex> aPrimaryVertexAllocator;

inline void * G4PrimaryVertex::operator new(size_t)
{
  void * aPrimaryVertex;
  aPrimaryVertex = (void *) aPrimaryVertexAllocator.MallocSingle();
  return aPrimaryVertex;
}

inline void G4PrimaryVertex::operator delete(void * aPrimaryVertex)
{
  aPrimaryVertexAllocator.FreeSingle((G4PrimaryVertex *) aPrimaryVertex);
}


#endif

