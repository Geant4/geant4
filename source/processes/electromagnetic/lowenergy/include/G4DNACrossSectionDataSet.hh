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
// $Id: G4DNACrossSectionDataSet.hh,v 1.4 2007/10/15 08:31:49 pia Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// Author: Riccardo Capra <capra@ge.infn.it>
//
// History:
// -----------
// 30 Jun 2005  RC         Created
// 14 Oct 2007  MGP        Removed inheritance from concrete class G4ShellEMDataSet
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Data set for an electromagnetic physics process
// A strategy pattern is used to encapsulate algorithms for data interpolation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef   G4DNACROSSSECTIONDATASET_HH
 #define  G4DNACROSSSECTIONDATASET_HH 1

 #include "G4ShellEMDataSet.hh"

 class G4DNACrossSectionDataSet : public G4VEMDataSet
 { 
  public:
   G4DNACrossSectionDataSet(G4VDataSetAlgorithm* argAlgorithm, 
			    G4double argUnitEnergies=MeV, 
			    G4double argUnitData=barn);

   virtual ~G4DNACrossSectionDataSet();

   virtual G4double FindValue(G4double argEnergy, G4int argComponentId=0) const;
  
   virtual void PrintData(void) const;

   virtual const G4VEMDataSet*  GetComponent(G4int argComponentId) const 
   { return components[argComponentId]; }

   virtual void AddComponent(G4VEMDataSet * argDataSet) 
   { components.push_back(argDataSet); }
   virtual size_t NumberOfComponents(void) const 
   { return components.size(); }

   virtual const G4DataVector& GetEnergies(G4int argComponentId) const 
   { return GetComponent(argComponentId)->GetEnergies(0); }

   virtual const G4DataVector& GetData(G4int argComponentId) const 
   { return GetComponent(argComponentId)->GetData(0); }

   virtual void SetEnergiesData(G4DataVector* argEnergies, G4DataVector* argData, G4int argComponentId);

   virtual G4bool LoadData(const G4String & argFileName);
   virtual G4bool SaveData(const G4String & argFileName) const;

   //   void CleanUpComponents();
   
  private:

   G4String FullFileName(const G4String & argFileName) const;

   // Hide copy constructor and assignment operator 
   G4DNACrossSectionDataSet();
   G4DNACrossSectionDataSet(const G4DNACrossSectionDataSet & copy);
   G4DNACrossSectionDataSet& operator=(const G4DNACrossSectionDataSet & right);

   std::vector<G4VEMDataSet*> components;          // Owned pointers

   G4int z;

   G4VDataSetAlgorithm* algorithm;           // Owned pointer 
  
   G4double unitEnergies;
   G4double unitData; 

   G4double GetUnitEnergies() const { return unitEnergies; }
   G4double GetUnitData() const { return unitData; }
   const G4VDataSetAlgorithm* GetAlgorithm() const { return algorithm; }
   
   void CleanUpComponents(void);


};
#endif /* G4DNACROSSSECTIONDATASET_HH */
