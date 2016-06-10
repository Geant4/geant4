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
// $Id: G4CsvNtupleManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Value type independent base class for parameters handled by Geant4 analysis
//
// Author: Ivana Hrivnacova, 04/09/2015  (ivana@ipno.in2p3.fr)

#ifndef G4VParameter_h
#define G4VParameter_h 1

#include "tools/cids" 

#include "G4MergeMode.hh"
#include "globals.hh"


class G4VParameter
{
  public:
    G4VParameter(const G4String& name,
                 G4MergeMode mergeMode = G4MergeMode::kAddition);
    G4VParameter(const G4VParameter& rhs);
    G4VParameter() = delete;
    virtual ~G4VParameter() {}

    // operators
    G4VParameter& operator=(const G4VParameter& rhs);

    // methods
    virtual void Merge(const G4VParameter& other) = 0;
    virtual void Reset() = 0;

    // get methods
    G4String    GetName() const;
    G4MergeMode GetMergeMode() const;

  protected:
    G4String    fName;
    G4MergeMode fMergeMode;
 };

#include "G4VParameter.icc"

#endif

