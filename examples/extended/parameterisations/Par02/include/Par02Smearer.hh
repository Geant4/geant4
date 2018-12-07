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
/// \file Par02Smearer.hh
/// \brief Definition of the Par02Smearer class

#ifndef PAR02_SMEARER_H
#define PAR02_SMEARER_H

#include "Par02Output.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandGauss.h"

/// Smearing of the particle momentum or energy.
///
/// A singleton class used to smear (alter) the particle momentum (for tracking
/// detectors) and energy (for calorimeters). In case the resolution is given,
/// the momentum (energy) is smeared with Gaussian distribution.
/// @author Anna Zaborowska

class Par02Smearer {
  public:

    /// Allows the access to the unique Par02Smearer class object.
    /// @return A pointer to the Par02Smearer class.
    static Par02Smearer* Instance();
    
    /// Smears the momentum with a given resolution.
    /// @param aTrack A track to smear.
    /// @param aResolution A resolution. Gaussian smearing is done with a 
    ///                    given resolution as a standard deviation.
    G4ThreeVector SmearMomentum( const G4Track* aTrack, G4double aResolution = -1 );
    
    /// Smears the energy deposit with a given resolution.
    /// @param aTrack A track to smear.
    /// @param aResolution A resolution. Gaussian smearing is done with a
    ///                    given resolution as a standard deviation.
    G4double SmearEnergy( const G4Track* aTrack, G4double aResolution = -1 );
    
    /// First possible type of smearing. Smears the momentum with a given resolution.
    /// @param aTrackOriginal A track to smear.
    /// @param aResolution A resolution taken as a standard deviation of a
    ///                    Gaussian distribution.
    G4ThreeVector SmearGaussian( const G4Track* aTrackOriginal, G4double aResolution );
    
    /// Returns a random number from a Gaussian distribution.
    /// @param aMean The mean of the Gaussian distribution.
    /// @param aStandardDeviation The standard deviation of a Gaussian distribution.
    G4double Gauss( G4double aMean, G4double aStandardDeviation );

  protected:
    
    /// A default constructor.
    Par02Smearer();

    ~Par02Smearer();

  private:
    
    /// A pointer to Par02Smearer object.
    static Par02Smearer* fPar02Smearer;
    
    /// CLHEP random engine.
    CLHEP::HepRandomEngine* fRandomEngine;
    
    /// CLHEP random engine used in gaussian smearing.
    CLHEP::RandGauss* fRandomGauss;
};

#endif

