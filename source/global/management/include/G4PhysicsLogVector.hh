// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsLogVector.hh,v 1.2 1999-11-16 17:40:41 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//--------------------------------------------------------------------
//      GEANT 4 class header file
//
//  G4PhysicsLogVector.hh
//
//  Class description:
//
//    A physics vector which has values of energy-loss, cross-section, 
//    and other physics values of a particle in matter in a given 
//    range of the energy, momentum, etc. The scale of energy/momentum
//    bins is in logarithmic.

//  History:
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    03 Mar. 1996, K.Amako : Implemented the 1st version
//    27 Apr. 1996, K.Amako : Cache mechanism added
//    01 Jul. 1996, K.Amako : Hidden bin from the user introduced
//    26 Sep. 1996, K.Amako : Constructor with only 'bin size' added.
//
//--------------------------------------------------------------------

#ifndef G4PhysicsLogVector_h
#define G4PhysicsLogVector_h 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4PhysicsVector.hh"


class G4PhysicsLogVector : public G4PhysicsVector  
{
  public:

    // Constructors
    G4PhysicsLogVector();
    G4PhysicsLogVector(size_t theNbin);
    G4PhysicsLogVector(G4double theEmin, G4double theEmax, size_t theNbin);
         // Because of logarithmic scale, note that 'theEmin' has to be 
         // greater than zero. No protection exists against this error.

    // Destructor
    ~G4PhysicsLogVector();

  protected:

    size_t FindBinLocation(G4double theEnergy) const;
         // Find bin# in which theEnergy belongs - pure virtual function

  private:

    G4double dBin;          // Bin width - useful only for fixed binning
    G4double baseBin;       // Set this in constructor for performance

};


inline size_t G4PhysicsLogVector::FindBinLocation(G4double theEnergy) const
{
  // For G4PhysicsLogVector, FindBinLocation is implemented using
  // a simple arithmetic calculation.
  //
  // Because this is a virtual function, it is accessed through a
  // pointer to the G4PhyiscsVector object for most usages. In this
  // case, 'inline' will not be invoked. However, there is a possibility 
  // that the user access to the G4PhysicsLogVector object directly and 
  // not through pointers or references. In this case, the 'inline' will
  // be invoked. (See R.B.Murray, "C++ Strategies and Tactics", Chap.6.6)

  return size_t( log10(theEnergy)/dBin - baseBin );
}

#endif
