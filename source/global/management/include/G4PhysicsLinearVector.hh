// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsLinearVector.hh,v 1.1 1999-01-07 16:09:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//--------------------------------------------------------------------
//      GEANT 4 class header file
//
//  G4PhysicsLinearVector.hh
//
//  Description:
//    A physics vector which has values of energy-loss, cross-section, 
//    and other physics values of a particle in matter in a given 
//    range of the energy, momentum, etc. The scale of energy/momentum
//    bins is in linear.
//
//  History:
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    03 Mar. 1996, K.Amako : Implemented the 1st version
//    01 Jul. 1996, K.Amako : Cache mechanism and hidden bin from the 
//                            user introduced.
//    26 Sep. 1996, K.Amako : Constructor with only 'bin size' added.
//
//--------------------------------------------------------------------

#ifndef G4PhysicsLinearVector_h
#define G4PhysicsLinearVector_h 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4PhysicsVector.hh"


class G4PhysicsLinearVector : public G4PhysicsVector  
{
  public:

    // Constructors
    G4PhysicsLinearVector();
    G4PhysicsLinearVector(size_t theNbin);
    G4PhysicsLinearVector(G4double theEmin, G4double theEmax, size_t theNbin);

    // Destructor
    ~G4PhysicsLinearVector();

  protected:

    size_t FindBinLocation(G4double theEnergy) const;
         // Find bin# in which theEnergy belongs - pure virtual function

  private:

    G4double dBin;          // Bin width - useful only for fixed binning
    G4double baseBin;       // Set this in constructor to gain performance
};


inline size_t G4PhysicsLinearVector::FindBinLocation(G4double theEnergy) const {

  // For G4PhysicsLinearVector, FindBinLocation is implemented using
  // a simple arithmetic calculation.
  //
  // Because this is a virtual function, it is accessed through a
  // pointer to the G4PhyiscsVector object for most usages. In this
  // case, 'inline' will not be invoked. However, there is a possibility 
  // that the user access to the G4PhysicsLinearVector object directly and 
  // not through pointers or references. In this case, the 'inline' will
  // be invoked. (See R.B.Murray, "C++ Strategies and Tactics", Chap.6.6)

  return size_t( theEnergy/dBin - baseBin ); 
}

#endif














