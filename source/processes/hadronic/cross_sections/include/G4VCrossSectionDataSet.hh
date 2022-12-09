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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4VCrossSectionDataSet
//
// Author  F.W. Jones, TRIUMF, 20-JAN-97
//
// Modifications:
// 23.01.2009 V.Ivanchenko move constructor and destructor to source
// 05.07.2010 V.Ivanchenko added name, min and max energy limit and
//            corresponding access methods
// 12.08.2011 G.Folger, V.Ivanchenko, T.Koi, D.Wright redesign the class
// 
//
// Class Description
// This is a base class for hadronic cross section data sets.  Users may 
// derive specialized cross section classes and register them with the
// appropriate process, or use provided data sets.
//
// Each cross section should have unique name
// Minimal and maximal energy for the cross section will be used in run
// time before IsApplicable method is called
// 
// Both the name and the energy interval will be used for documentation 
//
// Class Description - End

#ifndef G4VCrossSectionDataSet_h
#define G4VCrossSectionDataSet_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include <iostream>

class G4DynamicParticle;
class G4Isotope;
class G4Material;
class G4CrossSectionDataSetRegistry;

class G4VCrossSectionDataSet
{
public: //with description

  G4VCrossSectionDataSet(const G4String& nam = "");

  virtual ~G4VCrossSectionDataSet();

  //============== Is Applicable methods ===============================
  // The following three methods have default implementations returning
  // "false".  Derived classes should implement only needed methods.
 
  // Element-wise cross section
  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z, 
			     const G4Material* mat = nullptr);

  // Derived classes should implement this method if they provide isotope-wise 
  // cross sections.  Default arguments G4Element and G4Material are needed to 
  // access low-energy neutron cross sections, but are not required for others.
  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,    
			 const G4Element* elm = nullptr,
			 const G4Material* mat = nullptr);

  //============== GetCrossSection methods ===============================

  // This is a generic method to access cross section per element
  // This method should not be overwritten in a derived class
  inline G4double GetCrossSection(const G4DynamicParticle*, const G4Element*,
                                  const G4Material* mat = nullptr);

  // This is a generic method to compute cross section per element
  // If the DataSet is not applicable the method returns zero
  // This method should not be overwritten in a derived class
  G4double ComputeCrossSection(const G4DynamicParticle*, 
                               const G4Element*,
                               const G4Material* mat = nullptr);

  // Implement element cross section, IsApplicable does not checked.
  // In the default implementation a sum of isotope cross sections is computed
  virtual
  G4double ComputeCrossSectionPerElement(G4double kinEnergy, G4double loge,
                                         const G4ParticleDefinition*, 
                                         const G4Element*,
                                         const G4Material* mat = nullptr);

  // Implement these methods for element-wise cross section 
  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z,
				  const G4Material* mat = nullptr);

  // Derived classes should implement these methods if they provide isotope-wise
  // cross sections. Extra arguments G4Isotope, G4Element, and G4Material are 
  // needed to access low-energy neutron cross sections, but not in other cases. 
  virtual
  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,  
			      const G4Isotope* iso = nullptr,
			      const G4Element* elm = nullptr,
			      const G4Material* mat = nullptr);

  virtual
  G4double ComputeIsoCrossSection(G4double kinEnergy, G4double loge,
                                  const G4ParticleDefinition*, G4int Z, G4int A,  
			          const G4Isotope* iso = nullptr,
			          const G4Element* elm = nullptr,
			          const G4Material* mat = nullptr);

  //=====================================================================

  // Implement this method if needed
  // This method is called for element-wise cross section
  // Default implementation assumes equal cross sections for all isotopes 
  virtual const G4Isotope* SelectIsotope(const G4Element*, G4double kinEnergy, 
                                         G4double logE);

  // Implement this method if needed
  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  // Implement this method if needed
  // Default implementation will provide a dump of the cross section 
  // in logarithmic scale in the interval of applicability 
  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&);

  virtual void CrossSectionDescription(std::ostream&) const;

  virtual void SetVerboseLevel(G4int value);

  inline G4double GetMinKinEnergy() const;

  inline void SetMinKinEnergy(G4double value);

  inline G4double GetMaxKinEnergy() const;

  inline void SetMaxKinEnergy(G4double value);

  inline bool ForAllAtomsAndEnergies() const;

  inline void SetForAllAtomsAndEnergies(G4bool val);

  inline const G4String& GetName() const;

  inline void SetName(const G4String& nam);

  G4VCrossSectionDataSet & operator=
  (const G4VCrossSectionDataSet &right) = delete;
  G4VCrossSectionDataSet(const G4VCrossSectionDataSet&) = delete;

protected:

  G4int verboseLevel;

  G4String name;

private:

  G4CrossSectionDataSetRegistry* registry;

  G4double minKinEnergy;
  G4double maxKinEnergy;

  G4bool isForAllAtomsAndEnergies;

};

inline G4double 
G4VCrossSectionDataSet::GetCrossSection(const G4DynamicParticle* dp, 
                                        const G4Element* elm,
                                        const G4Material* mat)
{
  return ComputeCrossSection(dp, elm, mat);
}

inline void G4VCrossSectionDataSet::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline void G4VCrossSectionDataSet::SetMinKinEnergy(G4double value)
{
  minKinEnergy = value;
}

inline G4double G4VCrossSectionDataSet::GetMinKinEnergy() const
{
  return minKinEnergy;
}

inline void G4VCrossSectionDataSet::SetMaxKinEnergy(G4double value)
{
  maxKinEnergy = value;
}

inline G4double G4VCrossSectionDataSet::GetMaxKinEnergy() const
{
  return maxKinEnergy;
}

inline const G4String& G4VCrossSectionDataSet::GetName() const
{
  return name;
}

inline bool G4VCrossSectionDataSet::ForAllAtomsAndEnergies() const
{
  return isForAllAtomsAndEnergies;
}

inline void G4VCrossSectionDataSet::SetForAllAtomsAndEnergies(G4bool val)
{
  isForAllAtomsAndEnergies = val;
}

inline void G4VCrossSectionDataSet::SetName(const G4String& nam)
{
  name = nam;
}

#endif
