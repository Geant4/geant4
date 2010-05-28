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
// $Id: G4PhysicsVector.hh,v 1.31 2010-05-28 05:13:43 kurasige Exp $
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
//    01 Jul. 1996, K.Amako : Now GetValue not virtual
//    21 Sep. 1996, K.Amako : Added [] and () operators
//    11 Nov. 2000, H.Kurashige : Use STL vector for dataVector and binVector
//    09 Mar. 2001, H.Kurashige : Added G4PhysicsVectorType & Store/Retrieve()
//    02 Apr. 2008, A.Bagulya : Added SplineInterpolation() and SetSpline()
//    11 May  2009, V.Ivanchenko : Added ComputeSecondDerivatives
//    19 Jun. 2009, V.Ivanchenko : Removed hidden bin 
//    22 Dec. 2009  H.Kurashige  : Use pointers to G4PVDataVector
//    04 May. 2010  H.Kurashige  : Use G4PhysicsVectorCache
//    28 May  2010  H.Kurashige  : Stop using  pointers to G4PVDataVector
//---------------------------------------------------------------

#ifndef G4PhysicsVector_h
#define G4PhysicsVector_h 1

#include <vector>
#include "globals.hh"
#include "G4ios.hh"
#include <iostream>
#include <fstream>

#include "G4PhysicsVectorCache.hh"
#include "G4PhysicsVectorType.hh"

typedef std::vector<G4double> G4PVDataVector;

class G4PhysicsVector 
{
  public:  

    G4PhysicsVector(G4bool spline = false);
         // constructor  
         // This class is an abstract class with pure virtual method of
         // virtual size_t FindBinLocation(G4double theEnergy) const
         // So, default constructor is not supposed to be invoked explicitly

    G4PhysicsVector(const G4PhysicsVector&);
    G4PhysicsVector& operator=(const G4PhysicsVector&);
         // Copy constructor and assignment operator.

  public:  // with description

    virtual ~G4PhysicsVector();
         // destructor

    G4double Value(G4double theEnergy);
         // Get the cross-section/energy-loss value corresponding to the
         // given energy. An appropriate interpolation is used to calculate
         // the value. 

    inline G4double GetValue(G4double theEnergy, G4bool& isOutRange);
         // Obolete method to get value, isOutRange is not used anymore. 
         // This method is kept for the compatibility reason.

    G4int operator==(const G4PhysicsVector &right) const ;
    G4int operator!=(const G4PhysicsVector &right) const ;

    inline G4double operator[](const size_t binNumber) const ;
         // Returns simply the value in the bin specified by 'binNumber'
         // of the dataVector. The boundary check will not be done. 

    inline G4double operator()(const size_t binNumber) const ;
         // Returns simply the value in the bin specified by 'binNumber'
         // of the dataVector. The boundary check will not be Done. 

    inline void PutValue(size_t index, G4double theValue);
         // Put 'theValue' into the bin specified by 'binNumber'.
         // Take note that the 'index' starts from '0'.
         // To fill the vector, you have beforehand to construct a vector
         // by the constructor with Emin, Emax, Nbin. 'theValue' should
         // be the crosssection/energyloss value corresponding to the  
         // energy of the index. You can get this energy by the next method
         // or by the old method GetLowEdgeEnergy().

    void ScaleVector(G4double factorE, G4double factorV);
         // Scale all values of the vector and second derivatives
         // by factorV, energies by vectorE. This method may be applied 
         // for example after Retrieve a vector from an external file to 
         // convert values into Geant4 units

    inline G4double Energy(size_t index) const;
         // Returns simply the value in the energy specified by 'index'
         // of the energy vector. The boundary check will not be done. 
         // Use this function when you fill physis vector by PutValue().

    virtual G4double GetLowEdgeEnergy(size_t binNumber) const;
         // Obsolete method
         // Get the energy value at the low edge of the specified bin.
         // Take note that the 'binNumber' starts from '0'.
         // This value should be defined before the call.
         // The boundary check will not be done.

    inline size_t GetVectorLength() const;
         // Get the toal length (bin number) of the vector. 

    void FillSecondDerivatives();
        // Initialise second derivatives for spline keeping 
        // 3d derivative continues - default algorithm

    void ComputeSecDerivatives();
         // Initialise second derivatives for spline using algorithm 
         // which garantee only 1st derivative continues 
         // Warning: this method should be called when the vector 
         // is already filled

    void ComputeSecondDerivatives(G4double firstPointDerivative, 
                                  G4double endPointDerivative);
         // Initialise second derivatives for spline using 
         // user defined 1st derivatives at edge points
         // Warning: this method should be called when the vector 
         // is already filled

    inline G4bool IsFilledVectorExist() const;
         // Is non-empty physics vector already exist?

    inline void PutComment(const G4String& theComment);
         // Put a comment to the G4PhysicsVector. This may help to check
         // whether your are accessing to the one you want. 

    inline const G4String& GetComment() const;
         // Retrieve the comment of the G4PhysicsVector.

    inline G4PhysicsVectorType GetType() const;
         // Get physics vector type
  
    inline void SetSpline(G4bool);
         // Activate/deactivate Spline interpolation.

    virtual G4bool Store(std::ofstream& fOut, G4bool ascii=false);
    virtual G4bool Retrieve(std::ifstream& fIn, G4bool ascii=false);
         // To store/retrieve persistent data to/from file streams.

    friend std::ostream& operator<<(std::ostream&, const G4PhysicsVector&);

    
    G4double GetLastEnergy() const;
    G4double GetLastValue() const;
    size_t GetLastBin() const;
         // Get cache values 

  protected:

    virtual size_t FindBinLocation(G4double theEnergy) const=0;
         // Find the bin# in which theEnergy belongs - pure virtual function

    void DeleteData();
    void CopyData(const G4PhysicsVector& vec);
         // Internal methods for allowing copy of objects

  protected:

    G4PhysicsVectorType type;   // The type of PhysicsVector (enumerator)

    G4double edgeMin;           // Energy of first point
    G4double edgeMax;           // Energy of the last point

    size_t numberOfNodes;

    G4PhysicsVectorCache*  cache;

    G4PVDataVector  dataVector;    // Vector to keep the crossection/energyloss
    G4PVDataVector  binVector;     // Vector to keep energy
    G4PVDataVector  secDerivative; // Vector to keep second derivatives 

  private:

    G4bool SplinePossible();

    inline G4double LinearInterpolation(G4int lastBin);
         // Linear interpolation function
    inline G4double SplineInterpolation(G4int lastBin);
         // Spline interpolation function

    inline void Interpolation(G4int lastBin);

    G4String   comment;
    G4bool     useSpline;
};

#include "G4PhysicsVector.icc"

#endif
