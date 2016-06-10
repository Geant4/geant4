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
// $Id: G4VRangeToEnergyConverter.hh 70745 2013-06-05 10:54:00Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//
// Class Description
//  This class is base class for Range to Energy Converters.
//  Cut in energy corresponding to given cut value in range
//  is calculated for a material by using Convert method
//  
// ------------------------------------------------------------
//   First Implementation          5 Oct. 2002  H.Kurahige
// ------------------------------------------------------------

#ifndef G4VRangeToEnergyConverter_h
#define G4VRangeToEnergyConverter_h 1

#include "globals.hh"
#include <cmath>
#include "G4ios.hh"
#include <vector>

#include "G4ParticleDefinition.hh"

#include "G4PhysicsTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"
class G4PhysicsLogVector;

class G4VRangeToEnergyConverter
{
  public: // with description
  //  constructor
  G4VRangeToEnergyConverter();

  //  copy constructor
  G4VRangeToEnergyConverter(const G4VRangeToEnergyConverter &right);

  G4VRangeToEnergyConverter & operator=(const G4VRangeToEnergyConverter &right);

  public:
  //  destructor
  virtual ~G4VRangeToEnergyConverter();

  // equal opperators
  G4int operator==(const G4VRangeToEnergyConverter &right) const;
  G4int operator!=(const G4VRangeToEnergyConverter &right) const;

  public: // with description 
  // calculate energy cut from given range cut for the material
  virtual G4double Convert(G4double rangeCut, const G4Material* material);

  //  set energy range for all particle type
  static void SetEnergyRange(G4double lowedge, G4double highedge);

  //  get energy range for all particle type
  static G4double GetLowEdgeEnergy();
  static G4double GetHighEdgeEnergy();

  //  get/set max cut energy for all particle type
  static G4double GetMaxEnergyCut();
  static void SetMaxEnergyCut(G4double value);
  
  // return pointer to the particle type which this converter takes care
  const G4ParticleDefinition* GetParticleType() const;

  // return the Loss Table
  const  G4PhysicsTable* GetLossTable() const;   
   //-------------- Loss Table ------------------------------------------
   // theLossTable is a collection of loss vectors for all elements.
   // Each loss vector has energy loss values (cross section values
   // for neutral particles) which are calculated by
   // ComputeLoss(G4double AtomicNumber,G4double KineticEnergy).
   // ComputeLoss method is pure virtual and should be provided for each 
   // particle type

  // reset Loss Table and Range Vectors
  virtual void Reset();
    
 protected:

    static G4double               LowestEnergy, HighestEnergy;
    static G4double               MaxEnergyCut; 
    G4double                      fMaxEnergyCut;
   
    const G4ParticleDefinition*   theParticle;
    typedef G4PhysicsTable        G4LossTable;
    G4LossTable*                  theLossTable;
    G4int                         NumberOfElements;
  
    typedef G4PhysicsLogVector    G4LossVector;
    const G4int                   TotBin;

  protected:// with description  
    virtual void BuildLossTable();

    virtual G4double ComputeLoss(G4double AtomicNumber,
                                 G4double KineticEnergy
			       ) =0;

  //-------------- Range Table ------------------------------------------
  protected:
    typedef G4PhysicsLogVector G4RangeVector;

    virtual void BuildRangeVector(const G4Material* aMaterial,
                                   G4RangeVector* rangeVector);

    std::vector< G4RangeVector* > fRangeVectorStore;   
      
  protected:
    G4double ConvertCutToKineticEnergy(
                                       G4RangeVector* theRangeVector,
				       G4double       theCutInLength, 
                                       size_t         materialIndex
                                      ) const;

  public: // with description  
      void  SetVerboseLevel(G4int value);
      G4int GetVerboseLevel() const;
      // controle flag for output message
      //  0: Silent
      //  1: Warning message
      //  2: More

 private:
   G4int verboseLevel;

};

inline 
 void G4VRangeToEnergyConverter::SetVerboseLevel(G4int value)
{
   verboseLevel = value;
}

inline 
 G4int G4VRangeToEnergyConverter::GetVerboseLevel() const
{
   return verboseLevel;
}


inline 
 const G4ParticleDefinition* G4VRangeToEnergyConverter::GetParticleType() const
{
   return theParticle;
}
#endif








