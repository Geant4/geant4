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
// $Id: 
// GEANT4 tag $Name: 
//
// 
//---------------------------------------------------------------
//
//  G4FastSimulationVector.hh
//
//  Description:
//    Extends the STL vector to replace RW.
//
//  History:
//    May 00: Verderi && MoraDeFreitas - First Implementation.
//
//---------------------------------------------------------------

#ifndef G4FastSimulationVector_h
#define G4FastSimulationVector_h 1

#include <vector>
#include "G4Types.hh"

template<class T>
class G4FastSimulationVector : public std::vector<T*>
{
  typedef std::vector<T*> std_pvector;
  typedef typename std_pvector::iterator iterator;
  typedef typename std_pvector::const_iterator const_iterator;

public:

  G4FastSimulationVector(){};
  //  G4FastSimulationVector(const G4FastSimulationVector<T>&){};

  virtual ~G4FastSimulationVector(){};
  T* remove (const T*);	
  T* removeAt (G4int);
  void clearAndDestroy ();
};

#include "G4FastSimulationVector.icc"

#endif
