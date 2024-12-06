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

// Value type independent base class for accumulables handled by Geant4 analysis
//
// Author: Ivana Hrivnacova, IJCLab IN2P3/CNRS, 04/09/2015

#ifndef G4VAccumulable_h
#define G4VAccumulable_h 1

#include "G4AccType.hh"
#include "G4MergeMode.hh"

#include "globals.hh"

class G4VAccumulable
{
  public:
    G4VAccumulable(G4MergeMode mergeMode = G4MergeMode::kAddition);
    G4VAccumulable(const G4String& name,
                   G4MergeMode mergeMode = G4MergeMode::kAddition);
    G4VAccumulable(const G4VAccumulable& rhs) = default;
    G4VAccumulable(G4VAccumulable&& rhs) = default;
    virtual ~G4VAccumulable() = default;

    // Operators
    G4VAccumulable& operator=(const G4VAccumulable& rhs) = default;
    G4VAccumulable& operator=(G4VAccumulable&& rhs) = default;

    // Methods
    virtual void Merge(const G4VAccumulable& other) = 0;
    virtual void Reset() = 0;
    virtual void Print(G4PrintOptions options = G4PrintOptions()) const;
    virtual void SetMergeMode(G4MergeMode value);

    // Set methods
    void SetName(const G4String& name);
    void SetId(G4int id);

    // Get methods
    G4String  GetName() const;
    G4MergeMode GetMergeMode() const;
    G4int GetId() const;

    virtual G4AccType GetType() const { return G4AccType::kUser; }

  protected:
    void PrintBase(G4PrintOptions options) const;

    G4String  fName;
    G4MergeMode  fMergeMode = G4MergeMode::kAddition;
    G4int fId = G4Accumulables::kInvalidId;
 };

#include "G4VAccumulable.icc"

#endif

