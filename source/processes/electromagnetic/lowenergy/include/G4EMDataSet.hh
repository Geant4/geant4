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
// $Id: G4EMDataSet.hh,v 1.7 2006/06/29 19:35:31 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
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

#ifndef   G4EMDATASET_HH
 #define  G4EMDATASET_HH 1

 #include "globals.hh"
 #include "G4VEMDataSet.hh"

 class G4VDataSetAlgorithm;

 class G4EMDataSet : public G4VEMDataSet
 {
  public:
                                                G4EMDataSet(G4int argZ, G4VDataSetAlgorithm * argAlgorithm, G4double argUnitEnergies=MeV, G4double argUnitData=barn);
                                                G4EMDataSet(G4int argZ, G4DataVector * argEnergies, G4DataVector * argData, G4VDataSetAlgorithm * argAlgorithm, G4double argUnitEnergues=MeV, G4double argUnitData=barn);
   virtual                                     ~G4EMDataSet();
 
   virtual G4double                             FindValue(G4double argEnergy, G4int argComponentId=0) const;
  
   virtual void                                 PrintData(void) const;

   virtual const G4VEMDataSet *                 GetComponent(G4int /* argComponentId */) const { return 0; }
   virtual void                                 AddComponent(G4VEMDataSet * /* argDataSet */) {}
   virtual size_t                               NumberOfComponents(void) const { return 0; }

   virtual const G4DataVector &                 GetEnergies(G4int /* argComponentId */) const { return *energies; }
   virtual const G4DataVector &                 GetData(G4int /* argComponentId */) const { return *data; }
   virtual void                                 SetEnergiesData(G4DataVector * argEnergies, G4DataVector * argData, G4int argComponentId);

   virtual G4bool                               LoadData(const G4String & argFileName);
   virtual G4bool                               SaveData(const G4String & argFileName) const;
   
  private:
   size_t                                       FindLowerBound(G4double argEnergy) const;
   
   G4String                                     FullFileName(const G4String & argFileName) const;

   // Hide copy constructor and assignment operator 
                                                G4EMDataSet();
                                                G4EMDataSet(const G4EMDataSet & copy);
   G4EMDataSet &                                operator=(const G4EMDataSet & right);

   G4int                                        z;

   G4DataVector *                               energies;            // Owned pointer
   G4DataVector *                               data;                // Owned pointer

   G4VDataSetAlgorithm *                        algorithm;           // Owned pointer 
  
   G4double                                     unitEnergies;
   G4double                                     unitData;
 };
#endif /* G4EMDATASET_HH */
