// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsLnVector.hh,v 1.6 2001-03-09 03:39:26 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//--------------------------------------------------------------------
//      GEANT 4 class header file
//
//  G4PhysicsLnVector.hh
//
//  Class description:
//
//    A physics vector which has values of energy-loss, cross-section, 
//    and other physics values of a particle in matter in a given 
//    range of the energy, momentum, etc. The scale of energy/momentum
//    bins is natural logarithmic.
//
//  History:
//    27 Apr. 1999, M.G. Pia: Created, copying from G4PhysicsLogVector 
//    11 Nov. 2000, H.Kurashige : Use STL vector for dataVector and binVector
//
//--------------------------------------------------------------------

#ifndef G4PhysicsLnVector_h
#define G4PhysicsLnVector_h 1

#include "globals.hh"
#include "G4PhysicsVector.hh"

class G4PhysicsLnVector : public G4PhysicsVector  
{
  public:
    // Constructors
    G4PhysicsLnVector();
    G4PhysicsLnVector(size_t theNbin);

  public: // with description

    // Constructor
    G4PhysicsLnVector(G4double theEmin, G4double theEmax, size_t theNbin);
         // Because of logarithmic scale, note that 'theEmin' has to be 
         // greater than zero. No protection exists against this error.

    // Destructor
    ~G4PhysicsLnVector();

    virtual G4bool Retrieve(G4std::ifstream& fIn, G4bool ascii);

  protected:

    size_t FindBinLocation(G4double theEnergy) const;
         // Find bin# in which theEnergy belongs - pure virtual function

  private:

    G4double dBin;          // Bin width - useful only for fixed binning
    G4double baseBin;       // Set this in constructor for performance

};


inline 
 size_t G4PhysicsLnVector::FindBinLocation(G4double theEnergy) const 
{
 
  // For G4PhysicsLnVector, FindBinLocation is implemented using
  // a simple arithmetic calculation.
  //
  // Because this is a virtual function, it is accessed through a
  // pointer to the G4PhyiscsVector object for most usages. In this
  // case, 'inline' will not be invoked. However, there is a possibility 
  // that the user access to the G4PhysicsLnVector object directly and 
  // not through pointers or references. In this case, the 'inline' will
  // be invoked. (See R.B.Murray, "C++ Strategies and Tactics", Chap.6.6)

  return size_t( log(theEnergy)/dBin - baseBin );
}

#endif
