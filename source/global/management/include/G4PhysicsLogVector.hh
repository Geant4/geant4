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
// $Id: G4PhysicsLogVector.hh,v 1.15 2010-05-28 05:13:43 kurasige Exp $
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
//    26 Sep. 1996, K.Amako : Constructor with only 'bin size' added
//    11 Nov. 2000, H.Kurashige : Use STL vector for dataVector and binVector
//
//--------------------------------------------------------------------

#ifndef G4PhysicsLogVector_h
#define G4PhysicsLogVector_h 1

#include "globals.hh"
#include "G4PhysicsVector.hh"

class G4PhysicsLogVector : public G4PhysicsVector  
{
  public:

    G4PhysicsLogVector();
    explicit G4PhysicsLogVector(size_t theNbin);
      // Constructors


  public:  // with description

    G4PhysicsLogVector(G4double theEmin, G4double theEmax, size_t theNbin);
       // Because of logarithmic scale, note that 'theEmin' has to be 
       // greater than zero. No protection exists against this error.

    ~G4PhysicsLogVector();
      // Destructor

    G4PhysicsLogVector(const G4PhysicsLogVector&);
    G4PhysicsLogVector& operator=(const G4PhysicsLogVector&);
      // Copy constructor and assignment operator.
      
    G4bool Retrieve(std::ifstream& fIn, G4bool ascii);
      // To retrieve persistent data from file stream.
  
  protected:

    size_t FindBinLocation(G4double theEnergy) const;
      // Find bin# in which theEnergy belongs - pure virtual function

  private:

    G4double dBin;          // Bin width - useful only for fixed binning
    G4double baseBin;       // Set this in constructor for performance

};


inline 
 size_t G4PhysicsLogVector::FindBinLocation(G4double theEnergy) const
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

  return size_t( std::log10(theEnergy)/dBin - baseBin );
}

#endif
