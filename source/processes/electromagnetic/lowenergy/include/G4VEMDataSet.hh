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
// $Id: G4VEMDataSet.hh,v 1.2 2001-08-31 15:43:09 vnivanch Exp $
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

#ifndef G4VEMDATASET_HH
#define G4VEMDATASET_HH 1

#include "globals.hh"
#include "G4DataVector.hh"

class G4VDataSetAlgorithm;

class G4VEMDataSet 
{ 
public:

  G4VEMDataSet() { }

  virtual ~G4VEMDataSet() { }
 
  virtual G4double FindValue(G4double e, G4int id = 0) const = 0;

  virtual void PrintData() const = 0;

  virtual const G4VEMDataSet* GetComponent(G4int i) const { return 0; }

  virtual void AddComponent(G4VEMDataSet* dataSet) { }

  virtual size_t NumberOfComponents() const { return 0; }

  virtual const G4DataVector& GetEnergies(G4int i) const = 0;
  virtual const G4DataVector& GetData(G4int i) const = 0;
  
protected:

private:

  // hide assignment operator 
  G4VEMDataSet(const  G4VEMDataSet&);
  G4VEMDataSet & operator=(const G4VEMDataSet &right);

};
 
#endif
 










