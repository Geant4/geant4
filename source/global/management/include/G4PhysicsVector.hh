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
// $Id: G4PhysicsVector.hh 98864 2016-08-15 11:53:26Z gcosmo $
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
//    16 Aug. 2011  H.Kurashige  : Add dBin, baseBin and verboseLevel
//    02 Oct. 2013  V.Ivanchenko : FindBinLocation method become inlined;
//                                 instead of G4Pow G4Log is used
//---------------------------------------------------------------

#ifndef G4PhysicsVector_h
#define G4PhysicsVector_h 1

#include <iostream>
#include <fstream>
#include <vector>

#include "globals.hh"
#include "G4ios.hh"
#include "G4PhysicsVectorType.hh"
#include "G4Log.hh"

typedef std::vector<G4double> G4PVDataVector;

class G4PhysicsVector 
{
  public:// with description

    explicit G4PhysicsVector(G4bool spline = false);
         // default constructor - vector will be filled via Retrieve() method 

    G4PhysicsVector(const G4PhysicsVector&);
    G4PhysicsVector& operator=(const G4PhysicsVector&);
         // Copy constructor and assignment operator.

    virtual ~G4PhysicsVector();

    G4double Value(G4double theEnergy, size_t& lastidx) const; 
         // Get the cross-section/energy-loss value corresponding to the
         // given energy. An appropriate interpolation is used to calculate
         // the value. Consumer code got changed index and may reuse it
         // for the next call to save CPU for bin location. 

    inline G4double Value(G4double theEnergy) const; 
         // Get the cross-section/energy-loss value corresponding to the
         // given energy. An appropriate interpolation is used to calculate
         // the value. This method is kept for backward compatibility reason,
         // it should be used instead of the previous method if bin location 
         // cannot be kept thread safe

    inline G4double GetValue(G4double theEnergy, G4bool& isOutRange) const;
         // Obsolete method to get value, isOutRange is not used anymore. 
         // This method is kept for the compatibility reason.

    G4int operator==(const G4PhysicsVector &right) const ;
    G4int operator!=(const G4PhysicsVector &right) const ;

    inline G4double operator[](const size_t index) const ;
         // Returns the value for the specified index of the dataVector
         // The boundary check will not be done. 

    inline G4double operator()(const size_t index) const ;
         // Returns the value for the specified index of the dataVector
         // The boundary check will not be done. 

    inline void PutValue(size_t index, G4double theValue);
         // Put 'theValue' into the dataVector specified by 'index'.
         // Take note that the 'index' starts from '0'.
         // To fill the vector, you have beforehand to construct a vector
         // by the constructor with Emin, Emax, Nbin. 'theValue' should
         // be the crosssection/energyloss value corresponding to the  
         // energy of the index. 

    virtual void ScaleVector(G4double factorE, G4double factorV);
         // Scale all values of the vector and second derivatives
         // by factorV, energies by vectorE. This method may be applied 
         // for example after Retrieve a vector from an external file to 
         // convert values into Geant4 units

    inline G4double Energy(size_t index) const;
         // Returns the value in the energy specified by 'index'
         // of the energy vector. The boundary check will not be done. 
         // Use this function when compute cross section or dEdx 
         // before filling the vector by PutValue(..).

    inline G4double GetMaxEnergy() const;
         // Returns the energy of the last point of the vector

    G4double GetLowEdgeEnergy(size_t binNumber) const;
         // Obsolete method
         // Get the energy value at the low edge of the specified bin.
         // Take note that the 'binNumber' starts from '0'.
         // The boundary check will not be done.

    inline size_t GetVectorLength() const;
         // Get the total length of the vector. 

    inline size_t FindBin(G4double energy, size_t idx) const;
         // find low edge index of a bin for given energy
         // min value 0, max value VectorLength-1
         // idx is suggested bin number from user code

    void FillSecondDerivatives();
         // Initialise second derivatives for spline keeping 
         // 3d derivative continues - default algorithm
         // Warning: this method should be called when the vector 
         // is already filled

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

    G4double FindLinearEnergy(G4double rand) const;
         // Find energy using linear interpolation for vector
         // filled by cumulative probability function 
         // value of rand should be between 0 and 1

    inline G4bool IsFilledVectorExist() const;
         // Is non-empty physics vector already exist?

    inline G4PhysicsVectorType GetType() const;
         // Get physics vector type
  
    inline void SetSpline(G4bool);
         // Activate/deactivate Spline interpolation.

    G4bool Store(std::ofstream& fOut, G4bool ascii=false) const;
    virtual G4bool Retrieve(std::ifstream& fIn, G4bool ascii=false);
         // To store/retrieve persistent data to/from file streams.

    friend std::ostream& operator<<(std::ostream&, const G4PhysicsVector&);
    void DumpValues(G4double unitE=1.0, G4double unitV=1.0) const;
         // print vector

    inline void SetVerboseLevel(G4int value);

  protected:

    void DeleteData();
    void CopyData(const G4PhysicsVector& vec);
         // Internal methods for allowing copy of objects

    void PrintPutValueError(size_t index);

  protected:

    G4PhysicsVectorType type;   // The type of PhysicsVector (enumerator)

    G4double edgeMin;           // Energy of first point
    G4double edgeMax;           // Energy of the last point

    size_t numberOfNodes;

    G4PVDataVector  dataVector;    // Vector to keep the crossection/energyloss
    G4PVDataVector  binVector;     // Vector to keep energy
    G4PVDataVector  secDerivative; // Vector to keep second derivatives 

  private:

    G4bool SplinePossible();

    inline G4double LinearInterpolation(size_t idx, G4double energy) const;
         // Linear interpolation function
    inline G4double SplineInterpolation(size_t idx, G4double energy) const;
         // Spline interpolation function

    inline G4double Interpolation(size_t idx, G4double energy) const;

    inline size_t FindBinLocation(G4double theEnergy) const;
         // find low edge index of a bin for given energy
         // min value 0, max value VectorLength-1

    G4bool     useSpline;

  protected:

    G4double dBin;          // Bin width - useful only for fixed binning
    G4double baseBin;       // Set this in constructor for performance

    G4int verboseLevel;
};

#include "G4PhysicsVector.icc"

#endif
