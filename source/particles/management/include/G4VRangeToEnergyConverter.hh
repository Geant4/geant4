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
// $Id: G4VRangeToEnergyConverter.hh,v 1.2 2002-12-16 11:15:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4ios.hh"
#include "g4std/vector"

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
  
  // return pointer to the particle type which this converter takes care
  const G4ParticleDefinition* GetParticleType() const;

  // return the Loss Table
  const  G4PhysicsTable* GetLossTable() const;   
   //-------------- Loss Table ------------------------------------------
   // theLossTable is a collection of loss vectors for all elements.
   // Each loss vector has energy loss values (cross section values
   // for neutral particles) which are calculated by
   // ComputeLoss(G4double AtomicNumber,G4double KineticEnergy).
 
 protected:
    static G4double               LowestEnergy, HighestEnergy;

    const G4ParticleDefinition*   theParticle;
    typedef G4PhysicsTable        G4LossTable;
    G4LossTable*                  theLossTable;
    G4int                         NumberOfElements;

    typedef G4PhysicsLogVector    G4LossVector;
    G4int                         TotBin;

  protected:// with description  
   virtual void BuildLossTable();

    virtual G4double ComputeLoss(G4double AtomicNumber,
                                 G4double KineticEnergy
                                ) const;

  //-------------- Range Table ------------------------------------------
  protected:
    typedef G4PhysicsLogVector G4RangeVector;
    virtual void BuildRangeVector(const G4Material* aMaterial,
                                  G4double       maxEnergy,
                                  G4double       aMass,
                                  G4RangeVector* rangeVector);

  protected:
    G4double ConvertCutToKineticEnergy(
                                       G4RangeVector* theRangeVector,
				       G4double       theCutInLength, 
                                       size_t         materialIndex
                                      ) const;

    G4double RangeLinSimpson(
                            const G4ElementVector* elementVector,
                            const G4double* atomicNumDensityVector,
			    G4double aMass,
                            G4double taulow, G4double tauhigh
                          ); 

    G4double RangeLogSimpson(
                            const G4ElementVector* elementVector,
                            const G4double* atomicNumDensityVector,
                            G4double aMass,
                            G4double ltaulow, G4double ltauhigh
                          );

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


#endif








