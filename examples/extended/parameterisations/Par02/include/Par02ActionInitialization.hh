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
/// \file Par02ActionInitialization.hh
/// \brief Definition of the Par02ActionInitialization class

#ifndef PAR02_ACTION_INITIALIZATION_H
#define PAR02_ACTION_INITIALIZATION_H

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

/// Actions initialization.
///
/// Initialization of all user defined actions as well as mandatory
/// particle generator. 
/// Based on G4 examples/extended/parametrisations/Par01/include/
///                Par01ActionInitialization.hh .
/// @author Anna Zaborowska

class Par02ActionInitialization : public G4VUserActionInitialization {
  public:
    
    /// A default constructor: sets the output file name fFileName 
    /// (passed to Par02RunAction) to DefaultOutput.root as well as
    /// flag fSmear (passed to Par02EventAction) to true indicating
    /// that smearing will be performed.
    Par02ActionInitialization();
    
    /// A constructor: sets the flag fSmear (passed to Par02EventAction)
    ///                to true indicating that smearing will be performed.
    /// @param aOutName The output file name passed to Par02RunAction.
    Par02ActionInitialization( const G4String aOutName );
    
    /// A constructor.
    /// @param aOutName The output file name passed to Par02RunAction.
    /// @param aSmear The flag indicating if smearing should be done,
    ///               passed to Par02EventAction.
    Par02ActionInitialization( const G4String aOutName, const G4String aSmear );

    virtual ~Par02ActionInitialization();

    virtual void BuildForMaster() const;
    
    /// A method where all the user actions are created. 
    /// One of them, Par02PrimaryGeneratorAction is a mandatory class.
    virtual void Build() const;

  private:
    /// An output file name. Passed in Build() to the Par02RunAction.
    G4String fFileName;

    /// A flag indicating if smearing should be performed. 
    /// Passed in Build() to the Par02EventAction.    
    G4bool fSmear;
};

#endif

