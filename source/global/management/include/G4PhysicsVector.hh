// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsVector.hh,v 1.6 2001-01-09 11:27:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//---------------------------------------------------------------
//      GEANT 4 class header file
//
//  G4PhysicsVector.hh
//
//  Class description:
//
//    A physics vector which has values of energy-loss, cross-section, 
//    and other physics values of a particle in matter in a given 
//    range of the energy, momentum, etc.
//    This class serves as the base class for a vector having various 
//    energy scale, for example like 'log', 'linear', 'free', etc.

//  History:
//    02 Dec. 1995, G.Cosmo : Structure created based on object model
//    03 Mar. 1996, K.Amako : Implemented the 1st version
//    27 Apr. 1996, K.Amako : Cache mechanism added
//    01 Jul. 1996, K.Amako : Now GetValue not virtual. 
//    21 Sep. 1996, K.Amako : Added [] and () operators. 
//    11 Nov. 2000, H.Kurashige : use STL vector for dataVector and binVector
//
//---------------------------------------------------------------

#ifndef G4PhysicsVector_h
#define G4PhysicsVector_h 1

#include "g4std/vector"
#include "globals.hh"
#include "G4PhysicsTable.hh"

class G4PhysicsVector 
{
  public:

    // Constructor and destructor
    G4PhysicsVector();
    virtual ~G4PhysicsVector();

    // Public functions
    inline G4double GetValue(G4double theEnergy, G4bool& isOutRange);
         // Get the cross-section/energy-loss value corresponding to the
         // given energy. An appropriate interpolation is used to calculate
         // the value. 
         // [Note] isOutRange is not used anymore. This argument is kept
         //        for the compatibility reason.
    // Public operators
    G4int operator==(const G4PhysicsVector &right) const ;
    G4int operator!=(const G4PhysicsVector &right) const ;
    inline G4double operator[](const size_t binNumber) const ;
         // Returns simply the value in the bin specified by 'binNumber'
         // of the dataVector. The boundary check will be Done. If you
         // don't want this check, use the operator ().
    inline G4double operator()(const size_t binNumber) const ;
         // Returns simply the value in the bin specified by 'binNumber'
         // of the dataVector. The boundary check will not be Done. If 
         // you want this check, use the operator [].

    // Public functions
    inline void PutValue(size_t binNumber, G4double theValue);
         // Put 'theValue' into the bin specified by 'binNumber'.
         // Take note that the 'binNumber' starts from '0'.
         // To fill the vector, you have beforehand to Construct a vector
         // by the constructor with Emin, Emax, Nbin. 'theValue' should
         // be the crosssection/energyloss value corresponding to the low 
         // edge energy of the bin specified by 'binNumber'. You can get
         // the low edge energy value of a bin by GetLowEdgeEnergy().
    virtual G4double GetLowEdgeEnergy(size_t binNumber) const;
         // Get the energy value at the low edge of the specified bin.
         // Take note that the 'binNumber' starts from '0'.
         // This value is defined when a physics vector is constructed
         // by a constructor of a derived class. Use this function
         // when you fill physis vector by PutValue().
    inline size_t GetVectorLength() const;
         // Get the toal length (bin number) of the vector. 
    inline G4bool IsFilledVectorExist() const;
         // Is non-empty physics vector already exist?

    inline void LinkPhysicsTable(G4PhysicsTable& theTable);
         // Link the given G4PhysicsTable to the current G4PhysicsVector.
    inline G4bool IsLinkedTableExist() const;
         // Has this physics vector an extended physics table?
    inline const G4PhysicsTable* GetNextTable() const;
         // Returns the pointer to a physics table created for elements 
         // or isotopes (when the cross-sesctions or energy-losses 
         // depend explicitly on them).

    inline void PutComment(const G4String& theComment);
         // Put a comment to the G4PhysicsVector. This may help to check
         // whether your are accessing to the one you want. 
    inline G4String GetComment() const;
         // Retrieve the comment of the G4PhysicsVector.

  protected:

    typedef G4std::vector<G4double> G4PVDataVector;

    G4double edgeMin;           // Lower edge value of the lowest bin
    G4double edgeMax;           // Lower edge value of the highest bin
    size_t numberOfBin;

    G4double lastEnergy;        // Cache the last input value
    G4double lastValue;         // Cache the last output value   
    size_t lastBin;             // Cache the last bin location

    G4PVDataVector dataVector;    // Vector to keep the crossection/energyloss
    G4PVDataVector binVector;     // Vector to keep the low edge value of bin

    G4PhysicsTable* ptrNextTable;  
                                // Link to the connected physics table

    inline G4double LinearInterpolation(G4double theEnergy, size_t theLocBin);
         // Linear interpolation function
    virtual size_t FindBinLocation(G4double theEnergy) const=0;
         // Find the bin# in which theEnergy belongs - pure virtual function

  private:

    G4PhysicsVector(const G4PhysicsVector&);
    G4PhysicsVector& operator=(const G4PhysicsVector&);
         // Private copy constructor and assignment operator.

  private:

    G4String comment;
};

#include "G4PhysicsVector.icc"

#endif
