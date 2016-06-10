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
// $Id: G4IDataSet.hh 70904 2013-06-07 10:34:25Z gcosmo $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created as G4VDataSet
// 31 Jul 2008   MGP        Revised and renamed to G4IDataSet
//
// -------------------------------------------------------------------

// Class description:
// Data set abstract interface
// Applies a Composite design pattern for data library management

// -------------------------------------------------------------------

#ifndef  G4IDATASET_HH
#define  G4IDATASET_HH 1

#include "globals.hh"
#include "G4DataVector.hh"

class G4IDataSet 
{ 
public:
  G4IDataSet() { }
  virtual ~G4IDataSet() { }
  
  virtual G4double FindValue(G4double x, G4int componentId = 0) const = 0;
 
  virtual void PrintData(void) const = 0;
  
  virtual const G4IDataSet* GetComponent(G4int componentId) const = 0;
  virtual void AddComponent(G4IDataSet* dataSet) = 0;
  virtual size_t NumberOfComponents(void) const = 0;
 
  virtual const G4DataVector& GetEnergies(G4int componentId) const = 0;
  virtual const G4DataVector& GetData(G4int componentId) const = 0;
  virtual void SetEnergiesData(G4DataVector* x, G4DataVector* data, G4int component=0) = 0;
 
  virtual G4bool LoadData(const G4String& fileName) = 0;
  virtual G4bool SaveData(const G4String& fileName) const = 0;

  virtual G4double RandomSelect(G4int componentId = 0) const = 0;
   
private:
  // Hide copy constructor and assignment operator 
  G4IDataSet(const G4IDataSet& copy);
  G4IDataSet& operator=(const G4IDataSet& right);
};
#endif /* G4IDATASET_HH */
