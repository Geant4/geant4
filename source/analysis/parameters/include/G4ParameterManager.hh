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
// $Id: G4ParameterManager.hh 91012 2015-06-15 10:38:07Z ihrivnac $

// The common implementation of analysis manager classes.

// Author: Ivana Hrivnacova, 07/09/2015  (ivana@ipno.in2p3.fr)

#ifndef G4ParameterManager_h
#define G4ParameterManager_h 1

#include "G4Parameter.hh"
#include "G4AnalysisManagerState.hh"
#include "globals.hh"

#include <map>
#include <vector>

class G4AnalysisManagerState;
class G4VParameter;

class G4ParameterManager
{
  public:
    G4ParameterManager(G4bool isMaster = true);
    virtual ~G4ParameterManager();
     
    // static methods
    static G4ParameterManager* Instance();
   
    // Methods

    // Requirements for parameter type:
    // Must have defined operators: +, ==, +=, *=
    template <typename T> 
    void CreateParameter(const G4String& name, T value, 
                         G4MergeMode mergeMode = G4MergeMode::kAddition);

    template <typename T> 
    void RegisterParameter(G4Parameter<T>& parameter);

    template <typename T> 
    G4Parameter<T>*  GetParameter(const G4String& name, G4bool warn = true) const;

    // Handling user defined parameters
    void RegisterParameter(G4VParameter* parameter);
    G4VParameter*  GetParameter(const G4String& name, G4bool warn = true) const;


    // Iterators not (yet) provided

    // Methods for merging
    //G4bool MergeImpl(tools::histo::hmpi* hmpi);
    void Merge();
    void Reset();

  private:
    // static data members
    static G4ParameterManager* fgMasterInstance;
    static G4ThreadLocal G4ParameterManager* fgInstance;    

    // data members
    const G4AnalysisManagerState      fState;
    std::map<G4String, G4VParameter*> fMap;
    std::vector<G4VParameter*>        fParametersToDelete;
 };

#include "G4ParameterManager.icc"

#endif

