//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ParticleWithCuts.hh,v 1.13 2001-10-28 05:08:37 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//       first implementation, based on object model of Hisaya Kurashige, 
//                                                      21 Oct 1996
//       calculation of Range Table is based on implementeation for Muon 
//                                           by L.Urban, 10 May 1996
//       added  RestoreCuts  H.Kurashige 09 Mar. 2001
//       introduced material dependent range cuts   08 Sep. 2001
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
//    It also triggers the recomputation of the physics tables of the
// void SetRangeCut(G4double aCut, const G4Material*):
//    Set a cut value in range for the specified material and converts 
//    it into an energy cut.
// void ReCalcCuts():
//    Re calculate energy cut values with the previous cut value in range
// void ResetCuts():
//    Reset alll cut values in energy, though theCutInMaxInteractionLength
//    remains unchanged
// G4double GetCuts():
//    Returns value of the cut in interaction length for all the
//    processes of this particle type.
// const G4double* GetLengthCuts():
//    Returns an array of cuts in range (ordered per material)
// G4double  GetRangeThreshold(const G4Material* ) const:
//    Returns a range cut for a material   
// const G4double* GetEnergyCuts():
//    Returns an array of energy cuts (ordered per material)
// G4double GetEnergyThreshold(const G4Material* ) const:
//    Returns a energy cut for a material      
//

#ifndef G4ParticleWithCuts_h
#define G4ParticleWithCuts_h 1

#include "globals.hh"
#include "g4std/vector"
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
   G4double*  theCutInMaxInteractionLength;
   G4double*  theKineticEnergyCuts;

  public:  // With Description
   // G4ParticleWithCuts  
   // virtual methods derived from G4ParticleDefinition

  // Set Cuts methods
      virtual void              SetCuts(G4double aCut);
      // Set the range of aCut for all materials
      virtual void              SetRangeCut(G4double aCut, const G4Material*);
      // Set the cut range of aCut for a material
      virtual void              SetRangeCutVector(G4std::vector<G4double>&);
      // Set the vector of range cuts for all material

     // Get cuts methods
      virtual G4double*         GetLengthCuts() const;
      // Get an array of range cuts for all materials      
      virtual G4double          GetRangeThreshold(const G4Material* ) const;
      // Get a range cut for a material      
      virtual G4double*        GetEnergyCuts() const;
      // Get an array of energy cuts for all materials      
      virtual G4double         GetEnergyThreshold(const G4Material* ) const;
      // Get a energy cut for a material      

     // Other methods related with cuts
      virtual void              ResetCuts();
      // Reset alll cut values in energy 
      //  but theCutInMaxInteractionLength remain unchanged
      virtual void              ReCalcCuts();
      // Set cut values in energy derived from theCutInMaxInteractionLength

     //  set energy range  
      static void SetEnergyRange(G4double lowedge, G4double highedge) ;    

  public:  // With Description
   // This method concerning cut values is supposed to be used by
   // G4VUserPhysicsList to restore cutvalues witout calculation

   virtual void                  RestoreCuts(const G4double* cutInLength,
					     const G4double* cutInEnergy );
      
 protected:
    void    SetCutInMaxInteractionLength(G4double aCut);
    // Set a value of theCutInMaxInteractionLength for all materials
    void    SetCutInMaxInteractionLength(G4double aCut , G4int matrialIndex);
    // Set a value of theCutInMaxInteractionLength for a material
    void    SetCutInMaxInteractionLength(G4double aCut , 
                                         const G4Material* aMaterial);
    // Set a value of theCutInMaxInteractionLength for a material       
    void    SetEnergyCutValues(G4double energyCuts);
    // Set a energy cut value for all materials

    virtual   void  CalcEnergyCuts(const G4Material* material=0);
    // Calculate energy cut values by using range cuts

    // BuildPhysicsTable is defined as a dummy routine
    void  BuildPhysicsTable() {};

  protected:
  // cut values for proton is used for all heavy charged particles 
    static G4ParticleDefinition* theProton;
    G4bool  UseProtonCut();

    G4bool  CheckEnergyBinSetting() const;

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
				       G4RangeVector* theRangeVector,
                                       size_t         materialIndex
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

inline G4double*        G4ParticleWithCuts::GetLengthCuts() const 
{
        return theCutInMaxInteractionLength;
}

inline G4double*        G4ParticleWithCuts::GetEnergyCuts() const 
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
  if (theCutInMaxInteractionLength != 0) {
    CalcEnergyCuts();
  } else {
    if (GetVerboseLevel()>0) {
      G4cout << "G4ParticleWithCuts::ReCalcCuts() :";
      G4cout << "theCutInMaxInteractionLength is not defined " << G4endl; 
    }
  }
}

#include "G4Material.hh"

inline 
 G4double  G4ParticleWithCuts::GetEnergyThreshold(const G4Material* aMaterial) const
{
   return theKineticEnergyCuts[aMaterial->GetIndex()]; 
}

inline 
 G4double  G4ParticleWithCuts::GetRangeThreshold(const G4Material* aMaterial) const
{
   return theCutInMaxInteractionLength[aMaterial->GetIndex()]; 
}


#endif









