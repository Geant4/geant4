// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyUtilities.hh,v 1.7 2001-02-05 17:45:16 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ------------ G4LowEnergyPhotoElectric physics process ------
//                   by A.Forti 1999/06/28
//
// Class description:
// Class for the handling of the data libraries
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// ************************************************************

#ifndef G4LowEnergyUtilities_h
#define G4LowEnergyUtilities_h 1

//#include "G4FirstLevel.hh"
//#include "G4SecondLevel.hh"
#include "G4ThirdLevel.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsTable.hh"

typedef G4FirstLevel oneShellTable;
typedef G4SecondLevel oneAtomTable;
typedef G4ThirdLevel allAtomTable;

class G4LowEnergyUtilities
{

public:
  
  G4LowEnergyUtilities();

  ~G4LowEnergyUtilities();

  G4FirstLevel* BuildFirstLevelTables(const G4int, const G4int, const char*);
  G4SecondLevel* BuildSecondLevelTables(const G4int, const G4int, const char*);

  //inline functions
  G4int FindBinLocation(const G4double BinValue, const G4DataVector& arg);
  G4int FindBinLocation(const G4double, const G4PhysicsVector*);

  G4double DataLogInterpolation(const G4double, 
				const G4DataVector&,
				const G4DataVector&);

  G4double DataLogInterpolation(const G4double Argument, 
				const G4double AtomicNumber, 
				const G4PhysicsTable* Table);

  G4double DataSemiLogInterpolation(const G4double, 
				    const G4DataVector&, 
				    const G4DataVector&);

};

inline G4int G4LowEnergyUtilities::FindBinLocation(const G4double arg, const G4DataVector& vec){

  G4int numberOfBin = vec.size();
  G4int lowerBound = 0;
  G4int upperBound = numberOfBin-1;

  do {

    G4int midBin = (lowerBound + upperBound)/2;

    if( arg < vec[midBin] )

       upperBound = midBin-1;
    else

       lowerBound = midBin+1;

  } while (lowerBound <= upperBound); 

  return upperBound;
}


inline G4double G4LowEnergyUtilities::DataLogInterpolation(const G4double Argument, 
                                                               const G4DataVector& argVec, 
                                                               const G4DataVector& valVec){

  G4int theLoc = FindBinLocation(Argument, argVec); 

  if(theLoc == argVec.size()-1){
    return valVec[theLoc];
  }

  G4double val1 = valVec[theLoc], val2 = valVec[theLoc+1];
  G4double arg1 = argVec[theLoc], arg2 = argVec[theLoc+1];

  if(arg1 == 0.0) arg1 = 1e-17;	if(val1 == 0.0) val1 = 1e-17;
 
  G4double theVal = (log10(val1)*log10(arg2/Argument)
                     +log10(val2)*log10(Argument/arg1))/log10(arg2/arg1);
  
  theVal = pow(10,theVal);

  return theVal;
}

inline G4int G4LowEnergyUtilities::FindBinLocation(const G4double arg, 
						   const G4PhysicsVector* vec){

  if(!vec){

    G4Exception("G4LowEnergy: FindBinLocation: Vector Empty "
        "probably the program hasn't found data files or data files are empty");
  }

  G4int numberOfBin = vec->GetVectorLength();
  G4int lowerBound = 0;
  G4int upperBound = numberOfBin-1;
  do {

    G4int midBin = (lowerBound + upperBound)/2;

    if( arg < vec->GetLowEdgeEnergy(midBin) )

       upperBound = midBin-1;

    else

       lowerBound = midBin+1;

  } while (lowerBound <= upperBound); 

  return upperBound;
}

inline G4double G4LowEnergyUtilities::DataLogInterpolation(const G4double Argument, 
							       const G4double TableIndex, 
							       const G4PhysicsTable* Table){

  G4PhysicsVector* theVec = 0;
  theVec = (*Table)(TableIndex);

  G4int theLoc = FindBinLocation(Argument, theVec); 

  G4double val1 = (*theVec)(theLoc), val2 = (*theVec)(theLoc+1);
  G4double arg1 = theVec->GetLowEdgeEnergy(theLoc), arg2 = theVec->GetLowEdgeEnergy(theLoc+1);

  G4double theVal = (log10(val1)*log10(arg2/Argument)
		     +log10(val2)*log10(Argument/arg1))/log10(arg2/arg1);
  
  theVal = pow(10,theVal);
  return theVal;
}

inline G4double G4LowEnergyUtilities::DataSemiLogInterpolation(const G4double Argument, 
                                                                const G4DataVector& argVec, 
                                                                const G4DataVector& valVec){

  G4int theLoc = FindBinLocation(Argument, argVec); 

  if(theLoc == argVec.size()-1){
    return valVec[theLoc];
  }

  G4double val1 = valVec[theLoc], val2 = valVec[theLoc+1];
  G4double arg1 = argVec[theLoc], arg2 = argVec[theLoc+1];

  G4double theVal = (val1*log10(arg2/Argument)
                     +val2*log10(Argument/arg1))/log10(arg2/arg1);
  
  return theVal;
}

#endif










