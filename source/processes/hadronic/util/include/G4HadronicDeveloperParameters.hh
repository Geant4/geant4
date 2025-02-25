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
#ifndef G4HadronicDeveloperParameters_h
#define G4HadronicDeveloperParameters_h

#include "globals.hh"

#include<map>
#include<string>
#include<cfloat>


class G4HadronicDeveloperParameters
{
   public:
      static G4HadronicDeveloperParameters& GetInstance();

   //protected:
   private:
      G4HadronicDeveloperParameters();
      G4HadronicDeveloperParameters( const G4HadronicDeveloperParameters& );
      G4HadronicDeveloperParameters &operator=( const G4HadronicDeveloperParameters& );

   public:
      G4bool Set( const G4String& name , const G4bool );
      G4bool Set( const G4String& name , const G4int );
      G4bool Set( const G4String& name , const G4double );
      G4bool GetDefault( const G4String& name , G4bool& value );
      G4bool GetDefault( const G4String& name , G4int& value );
      G4bool GetDefault( const G4String& name , G4double& value );
      G4bool Get( const G4String& name , G4bool& value );
      G4bool Get( const G4String& name , G4int& value );
      G4bool Get( const G4String& name , G4double& value );
      G4bool DeveloperGet( const G4String& name , G4bool& value );
      G4bool DeveloperGet( const G4String& name , G4int& value );
      G4bool DeveloperGet( const G4String& name , G4double& value );
      void Dump( const G4String& name );

   //protected:
   public:
      G4bool SetDefault( const G4String& name , const G4bool value );
      G4bool SetDefault( const G4String& name , const G4int value , G4int lower_limit = -INT_MAX , G4int upper_limit = INT_MAX );
      G4bool SetDefault( const G4String& name , const G4double value , G4double lower_limit = -DBL_MAX , G4double upper_limit = DBL_MAX );

   private:
      G4bool get( const G4String& name , G4bool& value , G4bool check_change = false );
      G4bool get( const G4String& name , G4int& value , G4bool check_change = false );
      G4bool get( const G4String& name , G4double& value , G4bool check_change= false );

      std::map<G4String,G4bool> b_values;
      std::map<G4String,const G4bool> b_defaults;

      std::map<G4String,G4int> i_values;
      std::map<G4String,const G4int> i_defaults;
      std::map<G4String,std::pair<const G4int,const G4int>> i_limits;

      std::map<G4String,G4double> values;
      std::map<G4String,const G4double> defaults;
      std::map<G4String,std::pair<const G4double,const G4double>> limits;

      G4bool check_value_within_limits( std::pair<const G4double,const G4double>& , G4double );
      G4bool check_value_within_limits( std::pair<const G4int,const G4int>& , G4int );

      void issue_no_param( const G4String& name );
      void issue_has_changed( const G4String& name );
      void issue_non_eligible_value( const G4String& name );
      void issue_is_already_defined( const G4String& name );
      void issue_is_modified( const G4String& name );
};

#endif
