// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleWithCuts.hh,v 1.4 1999-10-28 23:24:13 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: first implementation, based on object model of
//      Hisaya Kurashige, 21 Oct 1996
//      Calculation of Range Table is based on 
//      implementeation for Muon by L.Urban, 10 May 1996
// ----------------------------------------------------------------
//      modified by Hisaya Kurashige, 04 Jan 1997
//      added verboseLevel by Hisaya Kurashige 13 Nov 1997
//      added ReCalcCuts() and ResetCuts  H.Kurashige 15 Nov.1996
//      BuildPhysicsTable() becomes dummy H.Kurashige 06 June 1998
//      added  GetEnergyThreshold  H.Kurashige 08 June 1998
//      change Lowest/HighestEnergy as static H.Kurashige 18 June 1998 
// ----------------------------------------------------------------
// Class Description
// "theCutInMaxInteractionLength", for charged particles, is
//    coincident with a cut in stopping range; for neutral
//    particles it corresponds to the distance allowed by
//    the total cross section.
// "theKineticEnergyCuts" is the vector of cuts in kinetic
//    energy (a value per material); it is allocated/computed
//    in the specific SetCuts method.  
// void SetCuts(G4double aCut):
//    Sets the cuts values relative to this particle type in stopping
//    range (or absorption length) and converts those cuts into energy cuts
//    for all the materials defined in the Material table.
//    [in protected method of CalcEnergyCuts(G4double aCut) ;
//    It also triggers the recomputation of the physics tables of the
// void ReCalcCuts():
//    Re calculate energy cut values with the previous cut value in range
// void ResetCuts():
//    Reset alll cut values in energy, though theCutInMaxInteractionLength
//    remains unchanged
// G4double GetCuts():
//    Returns value of the cut in interaction length for all the
//    processes of this particle type.
// const G4double* GetCutsInEnergy();
//    Returns vector of energy cuts (ordered per material)
//    corresponding to the current stopping range or absorption length
//

#ifndef G4ParticleWithCuts_h
#define G4ParticleWithCuts_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

#include "G4PhysicsTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"
class G4PhysicsLogVector;


class G4ParticleWithCuts : public G4ParticleDefinition
{
  public:
     G4ParticleWithCuts(const G4String&  aName,  
                G4double         mass,     
                G4double         width,
                G4double         charge,   
                G4int            iSpin,
                G4int            iParity,
                G4int            iConjugation,
                G4int            iIsospin,   
                G4int            iIsospinZ, 
                G4int            gParity,
                const G4String&  pType,
                G4int            lepton,
                G4int            baryon,
                G4int            encoding,
                G4bool           stable,
                G4double         lifetime,
                G4DecayTable     *decaytable,
		G4bool           resonance = false);
      virtual ~G4ParticleWithCuts();
   
  //--------------for SetCuts-------------------------------------------
  protected:
    G4double  theCutInMaxInteractionLength;
    G4double* theKineticEnergyCuts;

  public:  // With Description
   // virtual methods derived from G4ParticleDefinition
   virtual void          ResetCuts();
   // Reset alll cut values in energy 
   // but theCutInMaxInteractionLength remain unchanged
   virtual void          ReCalcCuts();
   // Set cut values in energy derived from theCutInMaxInteractionLength
   virtual void          SetCuts(G4double aCut);
   // Set cut values in energy derived from the cut range of aCut

   virtual G4double        	GetLengthCuts() const;
   virtual G4double* 	        GetEnergyCuts() const;
  
   virtual G4double      	GetEnergyThreshold(const G4Material* aMaterial) const;
   static  void          	SetEnergyRange(G4double, G4double);

 protected:
    virtual   void  CalcEnergyCuts(G4double aCut);

    // BuildPhysicsTable is defined as a dummy routine
    void  BuildPhysicsTable() {};

  //-------------- Loss Table ------------------------------------------
  // theLossTable is a collection of loss vectors for all elements.
  // Each loss vector has energy loss values (cross section values 
  // for neutral particles) which are calculated by  
  // ComputeLoss(G4double AtomicNumber,G4double KineticEnergy).

  protected:
    typedef G4PhysicsTable     G4LossTable;
    G4LossTable*               theLossTable;
    G4int                      NumberOfElements;

    typedef G4PhysicsLogVector G4LossVector;
    static G4double            LowestEnergy, HighestEnergy;
    G4int                      TotBin;

  protected:  
    virtual void BuildLossTable();
    virtual G4double ComputeLoss(G4double AtomicNumber,
                                 G4double KineticEnergy
                                ) const;
//-------------- Range Table ------------------------------------------
  protected:
    typedef G4PhysicsLogVector G4RangeVector;
    virtual void BuildRangeVector(const G4Material* aMaterial,
				  const G4LossTable* aLossTable,
				  G4double       maxEnergy,     
				  G4double       aMass,
                                  G4RangeVector* rangeVector);

  protected:  
    G4double ConvertCutToKineticEnergy(
                            G4RangeVector* theRangeVector
                          ) const;

    static G4double RangeLinSimpson(
			    const G4ElementVector* elementVector,
                            const G4double* atomicNumDensityVector,
                            const G4LossTable* aLossTable,
                            G4double aMass,   
			    G4double taulow, G4double tauhigh,
                            G4int nbin, G4int NumEl
                          );
   static G4double RangeLogSimpson(
			    const G4ElementVector* elementVector,
                            const G4double* atomicNumDensityVector,
                            const G4LossTable* aLossTable,
                            G4double aMass,   
			    G4double ltaulow, G4double ltauhigh,
                            G4int nbin, G4int NumEl
                          );
};

inline G4double	G4ParticleWithCuts::GetLengthCuts() const 
{
	return theCutInMaxInteractionLength;
}

inline G4double* 	G4ParticleWithCuts::GetEnergyCuts() const 
{
	return theKineticEnergyCuts;
}
  
inline void G4ParticleWithCuts::ResetCuts()
{
   if(theKineticEnergyCuts) delete [] theKineticEnergyCuts;
   theKineticEnergyCuts = 0;
}

inline void G4ParticleWithCuts::ReCalcCuts()
{
  if (theCutInMaxInteractionLength>0.0) {
    CalcEnergyCuts(theCutInMaxInteractionLength);
  } else {
    if (GetVerboseLevel()>0) {
      G4cout << "G4ParticleWithCuts::ReCalcCuts() :";
      G4cout << "theCutInMaxInteractionLength is not defined " << endl; 
    }
  }
}

inline void G4ParticleWithCuts::SetCuts(G4double aCut)
{
  CalcEnergyCuts(aCut);
}

inline 
 void G4ParticleWithCuts::SetEnergyRange(G4double lowedge, G4double highedge)
{
  LowestEnergy = lowedge;
  HighestEnergy = highedge;
}

#include "G4Material.hh"

inline 
 G4double  G4ParticleWithCuts::GetEnergyThreshold(const G4Material* aMaterial) const
{
   return theKineticEnergyCuts[aMaterial->GetIndex()]; 
}
#endif









