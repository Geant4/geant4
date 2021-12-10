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

// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file RootIO.hh
/// \brief Definition of the CaTS::RootIO class

#ifdef WITH_ROOT
#pragma once

#include "RtypesCore.h"
#include <G4Types.hh>
#include "G4ThreadLocalSingleton.hh"

class TFile;
class TTree;
class TBranch;
class Event;

class RootIO {
    friend class G4ThreadLocalSingleton< RootIO >;
public:
    virtual ~RootIO()= default;
    static RootIO* GetInstance();
    void Write(Event*);
    void Close();
    void Merge();
protected:
    RootIO();

private:
    TFile* fFile{ nullptr };
    G4int fNevents{ 0 };
    TTree* ftree{ nullptr };
    TBranch* fevtbranch{ nullptr };
    Long64_t fnb{ 0 };
    G4bool fevtinitialized{ false };
    static G4ThreadLocal RootIO* fgInstance;
};
#endif /* WITH_ROOT */
