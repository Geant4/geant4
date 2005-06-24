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
// $Id: G4ShellEMDataSet.hh,v 1.6 2005-06-24 09:55:05 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Data set for an electromagnetic physics process
// A strategy pattern is used to encapsulate algorithms for data interpolation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef   G4SHELLEMDATASET_HH
 #define  G4SHELLEMDATASET_HH 1

 #include "globals.hh"
 #include "G4VEMDataSet.hh"
 #include <vector>

 class G4VDataSetAlgorithm;

 class G4ShellEMDataSet : public G4VEMDataSet 
 { 
  public:
                                                G4ShellEMDataSet(G4int argZ, G4VDataSetAlgorithm* argAlgorithm, G4double argUnitEnergies=MeV, G4double argUnitData=barn);
   virtual                                     ~G4ShellEMDataSet();
 
   virtual G4double                             FindValue(G4double argEnergy, G4int argComponentId=0) const;
  
   virtual void                                 PrintData(void) const;

   virtual const G4VEMDataSet *                 GetComponent(G4int argComponentId) const { return components[argComponentId]; }
   virtual void                                 AddComponent(G4VEMDataSet * argDataSet) { components.push_back(argDataSet); }
   virtual size_t                               NumberOfComponents(void) const { return components.size(); }

   virtual const G4DataVector &                 GetEnergies(G4int argComponentId) const { return GetComponent(argComponentId)->GetEnergies(0); }
   virtual const G4DataVector &                 GetData(G4int argComponentId) const { return GetComponent(argComponentId)->GetData(0); }
   virtual void                                 SetEnergiesData(G4DataVector * argEnergies, G4DataVector * argData, G4int argComponentId);

   virtual G4bool                               LoadData(const G4String & argFileName);
   virtual G4bool                               SaveData(const G4String & argFileName) const;
   
  private:
   void                                         CleanUpComponents(void);

   G4String                                     FullFileName(const G4String & argFileName) const;
  
   // Hide copy constructor and assignment operator 
                                                G4ShellEMDataSet();
                                                G4ShellEMDataSet(const G4ShellEMDataSet & copy);
   G4ShellEMDataSet &                           operator=(const G4ShellEMDataSet & right);

   std::vector<G4VEMDataSet *>                  components;          // Owned pointers

   G4int                                        z;

   const G4VDataSetAlgorithm *                  algorithm;           // Owned pointer 
  
   G4double                                     unitEnergies;
   G4double                                     unitData;
 };
#endif /* G4SHELLEMDATASET_HH */
