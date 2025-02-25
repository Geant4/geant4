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
#include"G4HadronicDeveloperParameters.hh"

G4HadronicDeveloperParameters& G4HadronicDeveloperParameters::GetInstance(){
    static  G4HadronicDeveloperParameters instance;
    return instance;
}

G4HadronicDeveloperParameters::G4HadronicDeveloperParameters(){
}

G4HadronicDeveloperParameters::G4HadronicDeveloperParameters( const G4HadronicDeveloperParameters& ){ 
}

G4bool G4HadronicDeveloperParameters::SetDefault( const G4String& name , const G4bool value ){
   G4bool status = false;
   const std::map< G4String , const G4bool >::const_iterator it = b_defaults.find( name );
   if ( it == b_defaults.cend() ) {
      status = true;
      b_defaults.insert( std::pair<G4String, const G4bool>( name , value ) );
      b_values.insert( std::pair<G4String, G4bool>( name , value ) );
   } else {
            /*error*/
      issue_is_already_defined( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::SetDefault( const G4String& name , const G4int value , G4int lower_limit , G4int upper_limit ){
   G4bool status = false;
   const std::map< G4String , const G4int >::const_iterator it = i_defaults.find( name );
   if ( it == i_defaults.cend() ) {
      status = true;
      i_defaults.insert( std::pair<G4String, const G4int>( name , value ) );
      i_values.insert( std::pair<G4String, G4int>( name , value ) );
      i_limits.insert( std::pair< G4String,std::pair< const G4int , const G4int> >( name, std::pair< const G4int , const G4int> ( lower_limit , upper_limit ) ) );
   } else {
            /*error*/
      issue_is_already_defined( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::SetDefault( const G4String& name , const G4double value , G4double lower_limit , G4double upper_limit ){
   G4bool status = false;
   const std::map< G4String , const G4double >::const_iterator it = defaults.find( name );
   if ( it == defaults.cend() ) {
      status = true;
      defaults.insert( std::pair<G4String, const G4double>( name , value ) );
      values.insert( std::pair<G4String, G4double>( name , value ) );
      limits.insert( std::pair< G4String,std::pair< const G4double , const G4double> >( name, std::pair< const G4double , const G4double> ( lower_limit , upper_limit ) ) );
   } else {
            /*error*/
      issue_is_already_defined( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::Set( const G4String& name , const G4bool value ){
   G4bool status = false;
   const std::map<G4String,G4bool>::iterator it = b_values.find( name );
   if ( it != b_values.cend() ) {
      if ( it->second == b_defaults.find(name)->second ) {
         status = true;
         it->second = value;
      } else {
         /*change more than once*/
         issue_has_changed( name );
      }
   } else {
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::Set( const G4String& name , const G4int value ){
   G4bool status = false;
   const std::map<G4String,G4int>::iterator it = i_values.find( name );
   if ( it != i_values.cend() ) {
      if ( it->second == i_defaults.find(name)->second ) {
         if ( check_value_within_limits( i_limits.find(name)->second , value ) ) {
            /*value is OK*/
            status = true;
            it->second = value;
         } else {
            /*Value is outside of valid range*/
            issue_non_eligible_value( name );
         }
      } else {
         /*change more than once*/
         issue_has_changed( name );
      }
   } else {
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::Set( const G4String& name , const G4double value ){
   G4bool status = false;
   const std::map<G4String,G4double>::iterator it = values.find( name );
   if ( it != values.cend() ) {
      if ( it->second == defaults.find(name)->second ) {
         if ( check_value_within_limits( limits.find(name)->second , value ) ) {
            /*value is OK*/
            status = true;
            it->second = value;
         } else {
            /*Value is outside of valid range*/
            issue_non_eligible_value( name );
         }
      } else {
         /*change more than once*/
         issue_has_changed( name );
      }
   } else {
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::GetDefault( const G4String& name , G4bool& value ) {
   G4bool status = false;
   const std::map<G4String,const G4bool>::const_iterator it = b_defaults.find( name );
   if ( it != b_defaults.cend() ) {
      status = true;
      value = it->second;
   } else { 
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::GetDefault( const G4String& name , G4int& value ) {
   G4bool status = false;
   const std::map<G4String,const G4int>::const_iterator it = i_defaults.find( name );
   if ( it != i_defaults.cend() ) {
      status = true;
      value = it->second;
   } else { 
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::GetDefault( const G4String& name , G4double& value ) {
   G4bool status = false;
   const std::map<G4String,const G4double>::const_iterator it = defaults.find( name );
   if ( it != defaults.cend() ) {
      status = true;
      value = it->second;
   } else { 
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::Get( const G4String& name , G4double& value ) {
   return get( name , value );
}

G4bool G4HadronicDeveloperParameters::DeveloperGet( const G4String& name , G4double& value ) {
   return get( name , value , true );
}

G4bool G4HadronicDeveloperParameters::get( const G4String& name , G4bool& value , G4bool check_change ) {
   G4bool status = false;
   const std::map<G4String,G4bool>::const_iterator it = b_values.find( name );
   if ( it != b_values.cend() ) {
      status = true;
      value = it->second;
      if ( check_change && value != b_defaults.find(name)->second ) {
         //Parameter "name" has changed from default
         issue_is_modified( name );
      }
   } else { 
      //Parameter of "name" does not exist
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::Get( const G4String& name , G4int& value ) {
   return get( name , value );
}

G4bool G4HadronicDeveloperParameters::DeveloperGet( const G4String& name , G4int& value ) {
   return get( name , value , true );
}

G4bool G4HadronicDeveloperParameters::get( const G4String& name , G4int& value , G4bool check_change ) {
   G4bool status = false;
   const std::map<G4String,G4int>::const_iterator it = i_values.find( name );
   if ( it != i_values.cend() ) {
      status = true;
      value = it->second;
      if ( check_change && value != i_defaults.find(name)->second ) {
         //Parameter "name" has changed from default
         issue_is_modified( name );
      }
   } else { 
      //Parameter of "name" does not exist
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::Get( const G4String& name , G4bool& value ) {
   return get( name , value );
}

G4bool G4HadronicDeveloperParameters::DeveloperGet( const G4String& name , G4bool& value ) {
   return get( name , value , true );
}

G4bool G4HadronicDeveloperParameters::get( const G4String& name , G4double& value , G4bool check_change ) {
   G4bool status = false;
   const std::map<G4String,G4double>::const_iterator it = values.find( name );
   if ( it != values.cend() ) {
      status = true;
      value = it->second;
      if ( check_change && value != defaults.find(name)->second ) {
         /*Parameter "name" has changed from default*/
         issue_is_modified( name );
      }
   } else { 
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
   return status;
}

void G4HadronicDeveloperParameters::Dump( const G4String& name ) {
   //const std::map<G4String,G4double>::const_iterator it = values.find( name );
   if ( b_values.find( name ) != b_values.cend() ) {
      G4cout << "G4HadronicDeveloperParameters: "
      << "name = " << name 
      << ", default value = " << b_defaults.find( name )->second
      << ", current value = " << b_values.find( name )->second
      << "." << G4endl;
   } else if ( i_values.find( name ) != i_values.cend() ) {
      G4cout << "G4HadronicDeveloperParameters: "
      << "name = " << name 
      << ", default value = " << i_defaults.find( name )->second
      << ", lower limit = " << i_limits.find( name )->second.first
      << ", upper limit = " << i_limits.find( name )->second.second
      << ", current value = " << i_values.find( name )->second
      << "." << G4endl;
   } else if ( values.find( name ) != values.cend() ) {
      G4cout << "G4HadronicDeveloperParameters: "
      << "name = " << name 
      << ", default value = " << defaults.find( name )->second
      << ", lower limit = " << limits.find( name )->second.first
      << ", upper limit = " << limits.find( name )->second.second
      << ", current value = " << values.find( name )->second
      << "." << G4endl;
   } else {
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
}

G4bool G4HadronicDeveloperParameters::check_value_within_limits( std::pair<const G4double,const G4double>& alimits , const G4double value ){
   if ( alimits.first <= value && value <= alimits.second ) {
      return true;
   } else {
      return false;
   }
}

G4bool G4HadronicDeveloperParameters::check_value_within_limits( std::pair<const G4int,const G4int>& alimits , const G4int value ){
   if ( alimits.first <= value && value <= alimits.second ) {
      return true;
   } else {
      return false;
   }
}

void G4HadronicDeveloperParameters::issue_no_param( const G4String& name ){
   G4String text("Parameter ");
   text += name;  
   text += " does not exist.";
   G4Exception( "G4HadronicDeveloperParameters" , "HadDevPara_001", FatalException , text );
}

void G4HadronicDeveloperParameters::issue_has_changed( const G4String& name ) {
   G4String text("Parameter ");
   text += name;  
   text += " has already been changed once.";
   G4Exception( "G4HadronicDeveloperParameters" , "HadDevPara_002", FatalException , text );
}
void G4HadronicDeveloperParameters::issue_non_eligible_value( const G4String& name ) {
   G4String text("The value of the parameter ");
   text += name;  
   text += " is outside the allowable range.";
   G4Exception( "G4HadronicDeveloperParameters" , "HadDevPara_003", FatalException , text );
}
void G4HadronicDeveloperParameters::issue_is_already_defined( const G4String& name ) {
   G4String text("Parameter ");
   text += name;  
   text += " is already defined.";
   G4Exception( "G4HadronicDeveloperParameters" , "HadDevPara_004", FatalException , text );
}
void G4HadronicDeveloperParameters::issue_is_modified( const G4String& name ) {
   G4String text("Parameter ");
   text += name;  
   text += " has changed from default value.";
   G4Exception( "G4HadronicDeveloperParameters" , "HadDevPara_005", JustWarning , text );
}
