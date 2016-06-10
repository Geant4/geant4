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
// $Id: G4VEMDataSet.hh 66241 2012-12-13 18:34:42Z gunter $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP                 Created
//
// 15 Jul 2009   N.A.Karakatsanis    New methods added for loading logarithmic data
//                                   to enhance computing performance of interpolation
//
// 24 Dec 2010   V.Ivanchenko moved constructor and destructor to the source
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Data set for an electromagnetic physics process
// A strategy pattern is used to encapsulate algorithms for data interpolation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef  G4VEMDATASET_HH
#define  G4VEMDATASET_HH 1

#include "globals.hh"
#include "G4DataVector.hh"

class G4VEMDataSet 
{ 
public:

  G4VEMDataSet();
  virtual ~G4VEMDataSet();
  
  virtual G4double FindValue(G4double x, G4int componentId = 0) const = 0;
 
  virtual void PrintData(void) const = 0;
  
  virtual const G4VEMDataSet* GetComponent(G4int componentId) const = 0;
  virtual void AddComponent(G4VEMDataSet* dataSet) = 0;
  virtual size_t NumberOfComponents(void) const = 0;
 
  virtual const G4DataVector& GetEnergies(G4int componentId) const = 0;
  virtual const G4DataVector& GetData(G4int componentId) const = 0;
  virtual const G4DataVector& GetLogEnergies(G4int componentId) const = 0;
  virtual const G4DataVector& GetLogData(G4int componentId) const = 0;

  virtual void SetEnergiesData(G4DataVector* x, G4DataVector* data, G4int component=0) = 0;
  virtual void SetLogEnergiesData(G4DataVector* x,
                                  G4DataVector* data,
                                  G4DataVector* Log_x, 
                                  G4DataVector* Log_data,
                                  G4int component=0) = 0;
 
  virtual G4bool LoadData(const G4String& fileName) = 0;
  virtual G4bool LoadNonLogData(const G4String& fileName) = 0;

  virtual G4bool SaveData(const G4String& fileName) const = 0;

  virtual G4double RandomSelect(G4int componentId = 0) const = 0;
   
private:
  // Hide copy constructor and assignment operator 
  G4VEMDataSet(const G4VEMDataSet& copy);
  G4VEMDataSet& operator=(const G4VEMDataSet& right);
};
#endif /* G4VEMDATASET_HH */
