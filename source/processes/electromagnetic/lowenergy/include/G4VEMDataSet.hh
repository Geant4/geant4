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
// $Id: G4VEMDataSet.hh,v 1.7 2005-06-24 09:55:05 capra Exp $
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

#ifndef   G4VEMDATASET_HH
 #define  G4VEMDATASET_HH 1

 #include "globals.hh"
 #include "G4DataVector.hh"

 class G4VEMDataSet 
 { 
  public:
                                                G4VEMDataSet() { }
   virtual                                     ~G4VEMDataSet() { }
  
   virtual G4double                             FindValue(G4double argEnergy, G4int argComponentId=0) const = 0;
 
   virtual void                                 PrintData(void) const = 0;
  
   virtual const G4VEMDataSet *                 GetComponent(G4int argComponentId) const = 0;
   virtual void                                 AddComponent(G4VEMDataSet * argDataSet) = 0;
   virtual size_t                               NumberOfComponents(void) const = 0;
 
   virtual const G4DataVector &                 GetEnergies(G4int argComponentId) const = 0;
   virtual const G4DataVector &                 GetData(G4int argComponentId) const = 0;
   virtual void                                 SetEnergiesData(G4DataVector * argEnergies, G4DataVector * argData, G4int argComponent=0) = 0;
 
   virtual G4bool                               LoadData(const G4String & argFileName) = 0;
   virtual G4bool                               SaveData(const G4String & argFileName) const = 0;
   
  private:
   // Hide copy constructor and assignment operator 
                                                G4VEMDataSet(const G4VEMDataSet & copy);
   G4VEMDataSet &                               operator=(const G4VEMDataSet & right);
 };
#endif /* G4VEMDATASET_HH */
