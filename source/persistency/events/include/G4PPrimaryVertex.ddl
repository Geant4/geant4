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
// $Id: G4PPrimaryVertex.ddl,v 1.10.2.1 2001/06/28 19:11:26 gunter Exp $
// GEANT4 tag $Name:  $
//

// Class Description:
//      This is a class which represent a persistent primary vertex
//    in Geant4.  It is constructed by a G4PEvent object with a pointer
//    of transient G4PrimaryVertex object.
//


#ifndef G4PPrimaryVertex_h
#define G4PPrimaryVertex_h 1

#include "G4Pglobals.hh"
#include "G4ThreeVector.hh"
#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

class G4PrimaryVertex;
class G4PPrimaryParticle;
#pragma ooclassref G4PPrimaryParticle "G4PPrimaryParticle_ref.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PPrimaryVertex 
: public HepPersObj
{
  public: // with description
      G4PPrimaryVertex();
      G4PPrimaryVertex(const G4PrimaryVertex* vertex);
      ~G4PPrimaryVertex();
      // Constructors and destructor.

      G4PrimaryVertex* MakeTransientObject();
      // Creates a transient primary vertex object and returns the pointer.

  private:
      G4Pdouble X0;
      G4Pdouble Y0;
      G4Pdouble Z0;
      G4Pdouble T0;
      d_Ref<G4PPrimaryParticle> theParticle;
      d_Ref<G4PPrimaryParticle> theTail;
      d_Ref<G4PPrimaryVertex> nextVertex;
      G4Pint numberOfParticle;

  public:  // with description
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
      // Get methods to read the vertex values.

      HepRef(G4PPrimaryParticle) GetPrimary(G4int i) const;
      inline void SetNext(HepRef(G4PPrimaryVertex) nv)
      {
        if(nextVertex == 0)
        { nextVertex = nv; }
        else
        { nextVertex->SetNext(nv); }
      }
      inline HepRef(G4PPrimaryVertex) GetNext() const
      { return nextVertex; }
      // Set and Get methods to navigate thru the vertex chain.

};

#endif

