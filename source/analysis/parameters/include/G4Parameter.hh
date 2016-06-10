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
// $Id$

// Class template for parameters handled by Geant4 analysis
//
// Author: Ivana Hrivnacova, 07/09/2015  (ivana@ipno.in2p3.fr)

#ifndef G4Parameter_h
#define G4Parameter_h 1

#include "G4VParameter.hh"
#include "globals.hh"

template <typename T>
class G4Parameter : public G4VParameter
{
  public:
    G4Parameter(const G4String& name, T initValue, 
                G4MergeMode mergeMode = G4MergeMode::kAddition);
    G4Parameter(const G4Parameter& rhs);
    G4Parameter() = delete;
    virtual ~G4Parameter();

    // operators
    G4Parameter<T>& operator= (const G4Parameter<T>& rhs);
    G4Parameter<T>& operator+=(const G4Parameter<T>& rhs);
    G4Parameter<T>& operator*=(const G4Parameter<T>& rhs);

    G4Parameter<T>& operator= (const T& rhs);
    G4Parameter<T>& operator+=(const T& rhs);
    G4Parameter<T>& operator*=(const T& rhs);

    // methods
    virtual void Merge(const G4VParameter& other) final;
    virtual void Reset() final;

    // get methods
    T  GetValue() const;

  private:
    // data members
    T  fValue;
    T  fInitValue;
 };

// inline functions

#include "G4Parameter.icc"

#endif
