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
// $Id:$
//
// 
//---------------------------------------------------------------
//      GEANT 4 class header file
//
//  G4Physics2DVector.hh 
//
//  Class description:
//
//  A 2-dimentional vector with linear interpolation.

//  Author:        Vladimir Ivanchenko
//
//  Creation date: 25.09.2011
//
//  Modified:
//  16.05.2013 V.Ivanchenko removed cache; changed signature of 
//             several methods; all run time methods become const;
//             the class become read only in run time
//---------------------------------------------------------------

#ifndef G4Physics2DVector_h
#define G4Physics2DVector_h 1

#include <iostream>
#include <fstream>
#include <vector>

#include "globals.hh"
#include "G4ios.hh"
#include "G4PhysicsVectorType.hh"

typedef std::vector<G4double> G4PV2DDataVector;

class G4Physics2DVector 
{
public:  // with description

  G4Physics2DVector();
  // Vector will be filled via Retrieve method

  explicit G4Physics2DVector(size_t nx, size_t ny);
  // Vector will be filled via Put methods

  G4Physics2DVector(const G4Physics2DVector&);
  G4Physics2DVector& operator=(const G4Physics2DVector&);
  // Copy constructor and assignment operator.

  ~G4Physics2DVector();
  // destructor
 
  G4double Value(G4double x, G4double y, 
		 size_t& lastidx, size_t& lastidy) const;
  G4double Value(G4double x, G4double y) const;
  // Main method to interpolate 2D vector
  // Consumer class should provide initial values 
  // lastidx and lastidy,  

  inline void PutX(size_t idx, G4double value);
  inline void PutY(size_t idy, G4double value);
  inline void PutValue(size_t idx, size_t idy, G4double value);
  void PutVectors(const std::vector<G4double>& vecX,
		  const std::vector<G4double>& vecY);
  // Methods to fill vector 
  // Take note that the 'index' starts from '0'.

  void ScaleVector(G4double factor);
  // Scale all values of the vector by factor, 
  // This method may be applied 
  // for example after Retrieve a vector from an external file to 
  // convert values into Geant4 units

  G4double 
  FindLinearX(G4double rand, G4double y, size_t& lastidy) const;
  inline G4double FindLinearX(G4double rand, G4double y) const;
  // Find Y using linear interpolation for Y-vector
  // filled by cumulative probability function 
  // value of rand should be between 0 and 1

  inline G4double GetX(size_t index) const;
  inline G4double GetY(size_t index) const;
  inline G4double GetValue(size_t idx, size_t idy) const;
  // Returns simply the values of the vector by index
  // of the energy vector. The boundary check will not be done. 

  inline size_t FindBinLocationX(G4double x, size_t lastidx) const;
  inline size_t FindBinLocationY(G4double y, size_t lastidy) const;
  // Find the bin# in which theEnergy belongs
  // Starting from 0 

  inline size_t GetLengthX() const;
  inline size_t GetLengthY() const;
  // Get the lengths of the vector. 

  inline G4PhysicsVectorType GetType() const;
  // Get physics vector type

  inline void SetBicubicInterpolation(G4bool);
  // Activate/deactivate bicubic interpolation.

  void Store(std::ofstream& fOut) const;
  G4bool Retrieve(std::ifstream& fIn);
  // To store/retrieve persistent data to/from file streams.

  inline void SetVerboseLevel(G4int value);

protected:

  void PrepareVectors();

  void ClearVectors();

  void CopyData(const G4Physics2DVector& vec);

  G4double BicubicInterpolation(G4double x, G4double y,
				size_t idx, size_t idy) const;
  // Bicubic interpolation of 2D vector  

  size_t FindBinLocation(G4double z, const G4PV2DDataVector&) const;
  // Main method to local bin

  inline size_t FindBin(G4double z, const G4PV2DDataVector&, 
                        size_t idz, size_t idzmax) const;

private:

  G4double InterpolateLinearX(G4PV2DDataVector& v, G4double rand) const;

  inline G4double DerivativeX(size_t idx, size_t idy, G4double fac) const;
  inline G4double DerivativeY(size_t idx, size_t idy, G4double fac) const;
  inline G4double DerivativeXY(size_t idx, size_t idy, G4double fac) const;
  // computation of derivatives

  G4int operator==(const G4Physics2DVector &right) const = delete;
  G4int operator!=(const G4Physics2DVector &right) const = delete;

  G4PhysicsVectorType type;   // The type of PhysicsVector (enumerator)    

  size_t numberOfXNodes;
  size_t numberOfYNodes;

  G4PV2DDataVector  xVector;
  G4PV2DDataVector  yVector;
  std::vector<G4PV2DDataVector*> value;

  G4int  verboseLevel;
  G4bool useBicubic;
};

#include "G4Physics2DVector.icc"

#endif
