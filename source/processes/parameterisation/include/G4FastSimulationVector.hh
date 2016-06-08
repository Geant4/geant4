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

#include "g4std/vector"

template<class T>
class G4FastSimulationVector : public G4std::vector<T*>
{
  typedef G4std::vector<T*> std_pvector;
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
