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
// $Id$
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
//---------------------------------------------------------------

#ifndef G4PhysicsVector_h
#define G4PhysicsVector_h 1

#include "globals.hh"
#include "G4ios.hh"

#include <iostream>
#include <fstream>
#include <vector>

#include "G4Allocator.hh"
#include "G4PhysicsVectorCache.hh"
#include "G4PhysicsVectorType.hh"

typedef std::vector<G4double> G4PVDataVector;

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//The class PhysicsVectorPrivateSubclass is introduced to
//encapsulate the fields associated to the class G4PhysicsVector that
//may not be read-only.
#ifndef PHYSICSVECTORPRIVATESUBCLASS_HH
#define PHYSICSVECTORPRIVATESUBCLASS_HH

class PhysicsVectorPrivateSubclass
{
public:
  G4PhysicsVectorCache*  cache;
  void initialize() {
    cache = new G4PhysicsVectorCache();
    //    cache->lastEnergy = 0;
    //    cache->lastValue = 0;
    //    cache->lastBin = 0;
  };
};
#endif

//01.25.2009 Xin Dong: Phase II change for Geant4 multithreading.
//The class G4LogicalVolumeSubInstanceManager is introduced to
//encapsulate the methods used by both the master thread and
//worker threads to allocate memory space for the fields encapsulated
//by the class PhysicsVectorPrivateSubclass. When each thread                  
//initializes the value for these fields, it refers to them using a macro           
//definition defined below. For every G4PhysicsVectorCache instance, there is            
//a corresponding PhysicsVectorPrivateSubclass instance. All                   
//PhysicsVectorPrivateSubclass instances are organized by the                  
//class G4PhysicsVectorSubInstanceManager as an array. The field "                       
//int g4physicsVectorInstanceID" is added to the class G4PhysicsVectorCache.            
//The value of this field in each G4PhysicsVectorCache instance is the subscript         
//of the corresponding PhysicsVectorPrivateSubclass instance. In order         
//to use the class G4PhysicsVectorSubInstanceManager, we add a static member in          
//the class G4PhysicsVectorCache as follows: "                                           
//static G4PhysicsVectorSubInstanceManager g4physicsVectorSubInstanceManager".                
//Both the master and worker threads change the length of the array because         
//threads do not share all physics vectors. For the master thread, the array        
//for PhysicsVectorPrivateSubclass instances grows dynamically along with      
//G4PhysicsVector instances are created. For each worker thread, it copies the      
//array of PhysicsVectorPrivateSubclass instances from the master thread       
//first. For some physics vectors, worker threads share them and each thread        
//just uses the array copied to hold thread private data. However, each thread      
//will still create some physics vectors which extend the array to hold thread      
//private data although these physics vectors are not shared. It makes some         
//elements in the thread private array useless and consumes more memory space.      
//However, we determine to share almost all large physics vectors. Even if          
//we ignore some small physics vectors, the waste due to replication of these       
//small physics vectors is neglectable.
#ifndef G4PHYSICSVECTORSUBINSTANCEMANAGER_HH
#define G4PHYSICSVECTORSUBINSTANCEMANAGER_HH

#include "G4MTTransitoryPhysicsVector.hh"
typedef G4MTPrivatePhysicsVectorCounter<PhysicsVectorPrivateSubclass>  G4PhysicsVectorSubInstanceManager;

//01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.
//These macros changes the references to fields that are now encapsulated
//in the class PhysicsVectorPrivateSubclass.
#define cacheG4MTThreadPrivate ((g4physicsVectorSubInstanceManager.offset[g4physicsVectorInstanceID]).cache)

#endif

class G4PhysicsVector 
{
  public:  // with description

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.              
    //This new field is used as instance ID.                                        
    int g4physicsVectorInstanceID;

    //01.25.2009 Xin Dong: Phase II change for Geant4 multi-threading.              
    //This new field helps to use the class G4PhysicsVectorSubInstanceManager            
    //introduced above.                                                             
    static G4PhysicsVectorSubInstanceManager g4physicsVectorSubInstanceManager;

    G4PhysicsVector(G4bool spline = false);
         // constructor  
         // This class is an abstract class with pure virtual method of
         // virtual size_t FindBinLocation(G4double theEnergy) const
         // So, default constructor is not supposed to be invoked explicitly

    G4PhysicsVector(const G4PhysicsVector&);
    G4PhysicsVector& operator=(const G4PhysicsVector&);
         // Copy constructor and assignment operator.

    virtual ~G4PhysicsVector();
         // destructor

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    inline G4double Value(G4double theEnergy);
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

    virtual void ScaleVector(G4double factorE, G4double factorV);
         // Scale all values of the vector and second derivatives
         // by factorV, energies by vectorE. This method may be applied 
         // for example after Retrieve a vector from an external file to 
         // convert values into Geant4 units

    inline G4double Energy(size_t index) const;
         // Returns simply the value in the energy specified by 'index'
         // of the energy vector. The boundary check will not be done. 
         // Use this function when you fill physis vector by PutValue().

    inline G4double GetMaxEnergy() const;
         // Returns the energy of last point

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

    inline G4PhysicsVectorType GetType() const;
         // Get physics vector type
  
    inline void SetSpline(G4bool);
         // Activate/deactivate Spline interpolation.

    virtual G4bool Store(std::ofstream& fOut, G4bool ascii=false);
    virtual G4bool Retrieve(std::ifstream& fIn, G4bool ascii=false);
         // To store/retrieve persistent data to/from file streams.

    friend std::ostream& operator<<(std::ostream&, const G4PhysicsVector&);

    
    inline G4double GetLastEnergy() const;
    inline G4double GetLastValue() const;
    inline size_t GetLastBin() const;
         // Get cache values 

    inline void SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel(G4int);
         // Set/Get Verbose level

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

    //Change to thread private
    //    G4PhysicsVectorCache*  cacheG4MTThreadPrivate;

    G4PVDataVector  dataVector;    // Vector to keep the crossection/energyloss
    G4PVDataVector  binVector;     // Vector to keep energy
    G4PVDataVector  secDerivative; // Vector to keep second derivatives 

  private:

    void ComputeValue(G4double theEnergy);

    G4bool SplinePossible();

    inline G4double LinearInterpolation(G4int lastBin);
         // Linear interpolation function
    inline G4double SplineInterpolation(G4int lastBin);
         // Spline interpolation function

    inline void Interpolation(G4int lastBin);

    G4bool     useSpline;

  protected:
    G4double dBin;          // Bin width - useful only for fixed binning
    G4double baseBin;       // Set this in constructor for performance

    G4int verboseLevel;
};

#include "G4PhysicsVector.icc"

#endif
