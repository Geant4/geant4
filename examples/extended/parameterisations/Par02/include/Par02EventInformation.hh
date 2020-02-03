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
/// \file Par02EventInformation.hh
/// \brief Definition of the Par02EventInformation class

#ifndef PAR02_EVENT_INFORMATION_H
#define PAR02_EVENT_INFORMATION_H

#include "G4VUserEventInformation.hh"
#include "globals.hh"

/// Event information.
///
/// Describes the information that can be associated with a G4Event class object.
/// @author Anna Zaborowska

class Par02EventInformation : public G4VUserEventInformation {
  public:
    
    /// A default constructor. Sets flag fDoSmearing to true.
    Par02EventInformation();
    
    /// A constructor.
    /// @param aSmear The flag indicating if smearing should be done.
    Par02EventInformation( G4bool aSmear );

    virtual ~Par02EventInformation();
    
    /// Prints event information.
    virtual void Print() const;
    
    /// Sets the flag indicating if smearing should be done.
    /// @param aSmear A boolean flag.
    void SetDoSmearing( G4bool aSmear );
    
    /// Gets the flag indicating if smearing should be done.
    G4bool GetDoSmearing();

  private:
    
    /// A flag indicating if smearing should be performed. 
    /// It is read by implementations of G4VFastSimulationModel.
    G4bool fDoSmearing;
};

#endif

