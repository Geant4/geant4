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
// $Id: G4ParticleWithCuts.hh,v 1.15 2002-12-16 11:15:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------
//      GEANT 4 class header file
//
//      History: 
//       first implementation, based on object model of Hisaya Kurashige, 
//                                                      21 Oct 1996
//       calculation of Range Table is based on implementeation for Muon 
//                                           by L.Urban, 10 May 1996
//       added  RestoreCuts  H.Kurashige 09 Mar. 2001
//       introduced material dependent range cuts   08 Sep. 2001
//       
//       restructuring for Cuts per Region  by Hisaya    07 Oct.2002 
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
class G4ProductionCutsTable;

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
      G4ProductionCutsTable* theCutsTable;

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
					     const G4double* cutInEnergy )
   {
    G4cerr << "WARNING !" << G4endl;
    G4cerr << " RestoreCuts is currently not supported." << G4endl;
   }
      
  protected: 
   G4int     GetParticleIndex() const;
  
  protected: 
    G4double* theCutInMaxInteractionLength;
    G4double* theKineticEnergyCuts;

};

#endif
