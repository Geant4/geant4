// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyCompton.hh,v 1.1 1999-01-08 14:16:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file --- Copyright CERN 1995
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4LowEnergyCompton physics process ------
//                   by Michel Maire, April 1996
// ************************************************************
// 10-06-96, updated by M.Maire 
// 21-06-96, SetCuts implementation, M.Maire
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 17-02-97, New Physics scheme
// 25-02-97, GetMeanFreePath() now is public function
// 12-03-97, new physics scheme again
// ------------------------------------------------------------

#ifndef G4LowEnergyCompton_h
#define G4LowEnergyCompton_h 

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh" 
#include "G4ComptonScattering.hh"
#include "G4OrderedTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsFreeVector.hh" 
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Step.hh"

class G4LowEnergyCompton : public G4ComptonScattering{

private: 

  // hide assignment operator as private 
  G4LowEnergyCompton& operator=(const G4LowEnergyCompton &right);
  G4LowEnergyCompton(const G4LowEnergyCompton& );
 
public:
  
  G4LowEnergyCompton(const G4String& processName ="LowEnCompt");
  
  ~G4LowEnergyCompton();
  
  void BuildPhysicsTable(const G4ParticleDefinition& GammaType);
  void BuildMeanFreePathTable();
  
  G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize, G4ForceCondition* condition);

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  

protected:
  
  G4double DataLogInterpolation(G4double Argument, G4double AtomicNumber, G4PhysicsTable* Table);

  inline void ComputeMeanFreePath(){cout<<"ComputeMFP not available in this class"<<endl;}
  inline void ComputeMicroscopicCrossSection(){cout<<"ComputeMCS not available in this class"<<endl;}

  private:

  G4Element* SelectRandomAtom(const G4DynamicParticle*, G4Material*);
  G4int FindBinLocation(G4double,G4PhysicsVector*);

  G4PhysicsTable* theCrossSectionTable;
  G4PhysicsTable* theMeanFreePathTable;
  G4PhysicsTable* theScatteringFunctionTable;

  G4double LowestEnergyLimit; // low  energy limit of the crosssection data 
  G4double HighestEnergyLimit; // high energy limit of the crosssection data
  G4int NumbBinTable; // number of bins in the data tables

  G4double MeanFreePath; // actual Mean Free Path (current medium)
};

#endif
 
inline G4double 
G4LowEnergyCompton::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*){

  // returns the gamma mean free path in GEANT4 internal units
  const G4DynamicParticle* aDynamicGamma = aTrack.GetDynamicParticle();
  G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
  G4Material* aMaterial = aTrack.GetMaterial();
  
  G4bool isOutRange ;
  //  cout<<"************ the MFP ****************"<<endl;
  //  for(G4int U =0; U<theMeanFreePathTable->length(); U++){
    //    cout<<"ELEMENT "<<U+1<<":"<<endl;
  //
  //for(G4int f =0; f<(*theMeanFreePathTable)(U)->GetVectorLength() ; f++){
      
      //      cout<<"datvec["<<f<<"]: "<<(*(*theMeanFreePathTable)(U))(f)<<"    binvec["<<f<<"]: "<<(*theMeanFreePathTable)(U)->GetLowEdgeEnergy(f)<<endl;
  // }
  // }
    if(GammaEnergy > HighestEnergyLimit){
    MeanFreePath = DBL_MAX;
  }
  else {
    if (GammaEnergy < LowestEnergyLimit) GammaEnergy = LowestEnergyLimit;   
    MeanFreePath = (*theMeanFreePathTable)(aMaterial->GetIndex())->GetValue(GammaEnergy, isOutRange);
  }                                     
  
    //  cout<<"END OF LOW ENERGY COMPTON GETMEANFREEPATH!!!!!!!!!!!"<<endl;
    //  cout<<"----------------------------------------------------------"<<endl;
    //  cout<<"MeanFreePath: "<<MeanFreePath<<endl;

  return MeanFreePath;
}











