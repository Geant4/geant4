// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyUtilities.hh,v 1.2 1999-12-15 14:51:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4LowEnergyPhotoElectric physics process ------
//                   by Michel Maire, April 1996
// ************************************************************
// 12-06-96, Added SelectRandomAtom() method and new data member
//           for cumulative total cross section, by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 17-09-96, Dynamic array PartialSumSigma
//           split ComputeBindingEnergy(), M.Maire
// 08-01-97, crossection table + meanfreepath table, M.Maire
// 13-03-97, adapted for the new physics scheme, M.Maire
// ------------------------------------------------------------

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
  G4int FindBinLocation(const G4double BinValue, const G4Data& arg);
  G4int FindBinLocation(const G4double, const G4PhysicsVector*);

  G4double DataLogInterpolation(const G4double, 
				const G4Data&,
				const G4Data&);

  G4double DataLogInterpolation(const G4double Argument, 
				const G4double AtomicNumber, 
				const G4PhysicsTable* Table);

  G4double DataSemiLogInterpolation(const G4double, 
				    const G4Data&, 
				    const G4Data&);

};

inline G4int G4LowEnergyUtilities::FindBinLocation(const G4double arg, const G4Data& vec){

  G4int numberOfBin = vec.length();
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
};


inline G4double G4LowEnergyUtilities::DataLogInterpolation(const G4double Argument, 
                                                               const G4Data& argVec, 
                                                               const G4Data& valVec){

  G4int theLoc = FindBinLocation(Argument, argVec); 

  if(theLoc == argVec.length()-1){
    return valVec[theLoc];
  }

  G4double val1 = valVec[theLoc], val2 = valVec[theLoc+1];
  G4double arg1 = argVec[theLoc], arg2 = argVec[theLoc+1];

  if(arg1 == 0.0) arg1 = 1e-17;	if(val1 == 0.0) val1 = 1e-17;
 
  G4double theVal = (log10(val1)*log10(arg2/Argument)
                     +log10(val2)*log10(Argument/arg1))/log10(arg2/arg1);
  
  theVal = pow(10,theVal);

  return theVal;
};

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
                                                                const G4Data& argVec, 
                                                                const G4Data& valVec){

  G4int theLoc = FindBinLocation(Argument, argVec); 

  if(theLoc == argVec.length()-1){
    return valVec[theLoc];
  }

  G4double val1 = valVec[theLoc], val2 = valVec[theLoc+1];
  G4double arg1 = argVec[theLoc], arg2 = argVec[theLoc+1];

  G4double theVal = (val1*log10(arg2/Argument)
                     +val2*log10(Argument/arg1))/log10(arg2/arg1);
  
  return theVal;
}

#endif










