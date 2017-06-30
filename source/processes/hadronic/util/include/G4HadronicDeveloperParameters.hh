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
      G4bool Set( const std::string name , const G4bool );
      G4bool Set( const std::string name , const G4int );
      G4bool Set( const std::string name , const G4double );
      G4bool GetDefault( const std::string name , G4bool& value );
      G4bool GetDefault( const std::string name , G4int& value );
      G4bool GetDefault( const std::string name , G4double& value );
      G4bool Get( const std::string name , G4bool& value );
      G4bool Get( const std::string name , G4int& value );
      G4bool Get( const std::string name , G4double& value );
      G4bool DeveloperGet( const std::string name , G4bool& value );
      G4bool DeveloperGet( const std::string name , G4int& value );
      G4bool DeveloperGet( const std::string name , G4double& value );
      void Dump( const std::string name );

   //protected:
   public:
      G4bool SetDefault( const std::string name , const G4bool value );
      G4bool SetDefault( const std::string name , const G4int value , G4int lower_limit = -INT_MAX , G4int upper_limit = INT_MAX );
      G4bool SetDefault( const std::string name , const G4double value , G4double lower_limit = -DBL_MAX , G4double upper_limit = DBL_MAX );

   private:
      G4bool get( const std::string name , G4bool& value , G4bool check_change = false );
      G4bool get( const std::string name , G4int& value , G4bool check_change = false );
      G4bool get( const std::string name , G4double& value , G4bool check_change= false );

      std::map<std::string,G4bool> b_values;
      std::map<std::string,const G4bool> b_defaults;

      std::map<std::string,G4int> i_values;
      std::map<std::string,const G4int> i_defaults;
      std::map<std::string,std::pair<const G4int,const G4int>> i_limits;

      std::map<std::string,G4double> values;
      std::map<std::string,const G4double> defaults;
      std::map<std::string,std::pair<const G4double,const G4double>> limits;

      G4bool check_value_within_limits( std::pair<const G4double,const G4double>& , G4double );
      G4bool check_value_within_limits( std::pair<const G4int,const G4int>& , G4int );

      void issue_no_param( const std::string& name );
      void issue_has_changed( const std::string& name );
      void issue_non_eligible_value( const std::string& name );
      void issue_is_already_defined( const std::string& name );
      void issue_is_modified( const std::string& name );
};

#endif
