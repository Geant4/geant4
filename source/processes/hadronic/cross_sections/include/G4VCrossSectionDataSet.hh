// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VCrossSectionDataSet.hh,v 1.4 2000-12-14 11:25:57 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics abstract class: G4VCrossSectionData -- header file
// F.W. Jones, TRIUMF, 20-JAN-97
//
// Class Description
// This class serves as base class for cross-section data sets in geant4
// hadronic physics. Users can derive their specialized classes, and register with
// the system, or use provided data sets.
// Class Description - End

#ifndef G4VCrossSectionDataSet_h
#define G4VCrossSectionDataSet_h 1

#include "G4DynamicParticle.hh"
#include "G4Element.hh"


class G4VCrossSectionDataSet
{
public:

   G4VCrossSectionDataSet() :
      verboseLevel(0)
   {
   }

   virtual ~G4VCrossSectionDataSet()
   {
   }

public: //with description

   // the following methods need to be implemented for a new data-set.
   virtual
   G4bool IsApplicable(const G4DynamicParticle*, const G4Element*) = 0;

   virtual
   G4double GetCrossSection(const G4DynamicParticle*, const G4Element*) = 0;

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&) = 0;

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) = 0;

public: // Without Description

   void SetVerboseLevel(G4int value)
   {
      verboseLevel = value;
   }

   G4int GetVerboseLevel(G4int value)
   {
      return verboseLevel;
   }

protected:

   G4int verboseLevel;
};

#endif
