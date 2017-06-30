#include"G4HadronicDeveloperParameters.hh"

G4HadronicDeveloperParameters& G4HadronicDeveloperParameters::GetInstance(){
    static  G4HadronicDeveloperParameters instance;
    return instance;
}

G4HadronicDeveloperParameters::G4HadronicDeveloperParameters(){
}

G4HadronicDeveloperParameters::G4HadronicDeveloperParameters( const G4HadronicDeveloperParameters& ){ 
}

G4bool G4HadronicDeveloperParameters::SetDefault( const std::string name , const G4bool value ){
   G4bool status = false;
   const std::map< std::string , const G4bool >::iterator it = b_defaults.find( name );
   if ( it == b_defaults.end() ) {
      status = true;
      b_defaults.insert( std::pair<std::string, const G4bool>( name , value ) );
      b_values.insert( std::pair<std::string, G4bool>( name , value ) );
   } else {
            /*error*/
      issue_is_already_defined( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::SetDefault( const std::string name , const G4int value , G4int lower_limit , G4int upper_limit ){
   G4bool status = false;
   const std::map< std::string , const G4int >::iterator it = i_defaults.find( name );
   if ( it == i_defaults.end() ) {
      status = true;
      i_defaults.insert( std::pair<std::string, const G4int>( name , value ) );
      i_values.insert( std::pair<std::string, G4int>( name , value ) );
      i_limits.insert( std::pair< std::string,std::pair< const G4int , const G4int> >( name, std::pair< const G4int , const G4int> ( lower_limit , upper_limit ) ) );
   } else {
            /*error*/
      issue_is_already_defined( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::SetDefault( const std::string name , const G4double value , G4double lower_limit , G4double upper_limit ){
   G4bool status = false;
   const std::map< std::string , const G4double >::iterator it = defaults.find( name );
   if ( it == defaults.end() ) {
      status = true;
      defaults.insert( std::pair<std::string, const G4double>( name , value ) );
      values.insert( std::pair<std::string, G4double>( name , value ) );
      limits.insert( std::pair< std::string,std::pair< const G4double , const G4double> >( name, std::pair< const G4double , const G4double> ( lower_limit , upper_limit ) ) );
   } else {
            /*error*/
      issue_is_already_defined( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::Set( const std::string name , const G4bool value ){
   G4bool status = false;
   const std::map<std::string,G4bool>::iterator it = b_values.find( name );
   if ( it != b_values.end() ) {
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

G4bool G4HadronicDeveloperParameters::Set( const std::string name , const G4int value ){
   G4bool status = false;
   const std::map<std::string,G4int>::iterator it = i_values.find( name );
   if ( it != i_values.end() ) {
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

G4bool G4HadronicDeveloperParameters::Set( const std::string name , const G4double value ){
   G4bool status = false;
   const std::map<std::string,G4double>::iterator it = values.find( name );
   if ( it != values.end() ) {
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

G4bool G4HadronicDeveloperParameters::GetDefault( const std::string name , G4bool& value ) {
   G4bool status = false;
   const std::map<std::string,const G4bool>::iterator it = b_defaults.find( name );
   if ( it != b_defaults.end() ) {
      status = true;
      value = it->second;
   } else { 
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::GetDefault( const std::string name , G4int& value ) {
   G4bool status = false;
   const std::map<std::string,const G4int>::iterator it = i_defaults.find( name );
   if ( it != i_defaults.end() ) {
      status = true;
      value = it->second;
   } else { 
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::GetDefault( const std::string name , G4double& value ) {
   G4bool status = false;
   const std::map<std::string,const G4double>::iterator it = defaults.find( name );
   if ( it != defaults.end() ) {
      status = true;
      value = it->second;
   } else { 
      /*Parameter of "name" does not exist*/
      issue_no_param( name );
   }
   return status;
}

G4bool G4HadronicDeveloperParameters::Get( const std::string name , G4double& value ) {
   return get( name , value );
}

G4bool G4HadronicDeveloperParameters::DeveloperGet( const std::string name , G4double& value ) {
   return get( name , value , true );
}

G4bool G4HadronicDeveloperParameters::get( const std::string name , G4bool& value , G4bool check_change ) {
   G4bool status = false;
   const std::map<std::string,G4bool>::iterator it = b_values.find( name );
   if ( it != b_values.end() ) {
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

G4bool G4HadronicDeveloperParameters::Get( const std::string name , G4int& value ) {
   return get( name , value );
}

G4bool G4HadronicDeveloperParameters::DeveloperGet( const std::string name , G4int& value ) {
   return get( name , value , true );
}

G4bool G4HadronicDeveloperParameters::get( const std::string name , G4int& value , G4bool check_change ) {
   G4bool status = false;
   const std::map<std::string,G4int>::iterator it = i_values.find( name );
   if ( it != i_values.end() ) {
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

G4bool G4HadronicDeveloperParameters::Get( const std::string name , G4bool& value ) {
   return get( name , value );
}

G4bool G4HadronicDeveloperParameters::DeveloperGet( const std::string name , G4bool& value ) {
   return get( name , value , true );
}

G4bool G4HadronicDeveloperParameters::get( const std::string name , G4double& value , G4bool check_change ) {
   G4bool status = false;
   const std::map<std::string,G4double>::iterator it = values.find( name );
   if ( it != values.end() ) {
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

void G4HadronicDeveloperParameters::Dump( const std::string name ) {
   //const std::map<std::string,G4double>::iterator it = values.find( name );
   if ( b_values.find( name ) != b_values.end() ) {
      G4cout << "G4HadronicDeveloperParameters: "
      << "name = " << name 
      << ", default value = " << b_defaults.find( name )->second
      << ", current value = " << b_values.find( name )->second
      << "." << G4endl;
   } else if ( i_values.find( name ) != i_values.end() ) {
      G4cout << "G4HadronicDeveloperParameters: "
      << "name = " << name 
      << ", default value = " << i_defaults.find( name )->second
      << ", lower limit = " << i_limits.find( name )->second.first
      << ", upper limit = " << i_limits.find( name )->second.second
      << ", current value = " << i_values.find( name )->second
      << "." << G4endl;
   } else if ( values.find( name ) != values.end() ) {
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

void G4HadronicDeveloperParameters::issue_no_param( const std::string& name ){
   std::string text("Parameter ");
   text += name;  
   text += " does not exist.";
   G4Exception( "G4HadronicDeveloperParameters" , "HadDevPara_001", FatalException , text.c_str() );
}

void G4HadronicDeveloperParameters::issue_has_changed( const std::string& name ) {
   std::string text("Parameter ");
   text += name;  
   text += " has already been changed once.";
   G4Exception( "G4HadronicDeveloperParameters" , "HadDevPara_002", FatalException , text.c_str() );
}
void G4HadronicDeveloperParameters::issue_non_eligible_value( const std::string& name ) {
   std::string text("The value of the parameter ");
   text += name;  
   text += " is outside the allowable range.";
   G4Exception( "G4HadronicDeveloperParameters" , "HadDevPara_003", FatalException , text.c_str() );
}
void G4HadronicDeveloperParameters::issue_is_already_defined( const std::string& name ) {
   std::string text("Parameter ");
   text += name;  
   text += " is already defined.";
   G4Exception( "G4HadronicDeveloperParameters" , "HadDevPara_004", FatalException , text.c_str() );
}
void G4HadronicDeveloperParameters::issue_is_modified( const std::string& name ) {
   std::string text("Parameter ");
   text += name;  
   text += " has changed from default value.";
   G4Exception( "G4HadronicDeveloperParameters" , "HadDevPara_005", JustWarning , text.c_str() );
}
