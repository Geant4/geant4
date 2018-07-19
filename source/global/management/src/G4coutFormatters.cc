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
// $Id: G4coutFormatters.cc 103582 2017-04-18 17:24:45Z adotti $
//
// --------------------------------------------------------------------
//
// G4coutFormatters.cc
//
// Author: A.Dotti (SLAC), April 2017
// --------------------------------------------------------------------

#include "G4coutFormatters.hh"

namespace G4coutFormatters
{
  // Internal functions and utilites used to setup default formatters
  namespace
  {
    // Split a single string in an array of strings
    String_V split(const G4String& input, char separator='\n')
    {
      String_V output;
      G4String::size_type prev_pos=0, pos=0;
      while( (pos=input.find(separator,pos))!=G4String::npos)
      {
          G4String substr( input.substr(prev_pos,pos-prev_pos)) ;
          output.push_back(substr);
          prev_pos = ++pos;
      }
      // output.push_back( input.substr(prev_pos,pos-prev_pos));
      return output;
    }

    // Return a syslog style message with input message, type identifies
    // the type of the message
    G4bool transform( G4String& input , const G4String& type)
    {
      std::time_t result = std::time(nullptr);
      std::ostringstream newm;
#if __GNUC__ >= 5
      newm << std::put_time(std::localtime(&result),"%d/%b/%Y:%H:%M:%S %z");
#else
      std::tm* time_ = std::localtime(&result);
      newm << time_->tm_mday << "/" << time_->tm_mon << "/" << time_->tm_year;
      newm << ":" << time_->tm_hour << ":"<<time_->tm_min<<":"<<time_->tm_sec;
#endif
      newm<<" "<<type<<" [";
      G4String delimiter = "";
      for (const auto& el : split(input) )
      {
          if ( !el.empty() )
          {
            newm << delimiter << el ;
            delimiter = "\\n";
          }
      }
      newm<<" ]"<<G4endl;
      input = newm.str();
      return true;
    }

    // Style used in master thread
    G4String masterStyle = "";

    // Modify output to look like syslog messages:
    // DATE TIME **LOG|ERROR** [ "multi","line","message"]
    SetupStyle_f SysLogStyle = [](G4coutDestination* dest)->G4int
    {
      if ( dest != nullptr )
      {
          dest->AddCoutTransformer(std::bind(&transform,std::placeholders::_1,
                                             "INFO"));
          dest->AddCerrTransformer(std::bind(&transform,std::placeholders::_1,
                                             "ERROR"));
      }
      return 0;
    };

    // Bring back destination to original state
    SetupStyle_f DefaultStyle = [](G4coutDestination* dest)->G4int
    {
      if ( dest != nullptr )
      {
          dest->ResetTransformers();
      }
      return 0;
    };

    std::unordered_map<std::string,SetupStyle_f> transformers =
    {
        {ID::SYSLOG,SysLogStyle},
        {ID::DEFAULT,DefaultStyle}
    };
  }

  void SetMasterStyle(const G4String& news )
  {
    masterStyle = news;
  }

  G4String GetMasterStyle()
  {
    return masterStyle;
  }

  void SetupStyleGlobally(const G4String& news)
  {
    static G4coutDestination ss;
    G4coutbuf.SetDestination(&ss);
    G4cerrbuf.SetDestination(&ss);
    G4coutFormatters::HandleStyle(&ss,news);
    G4coutFormatters::SetMasterStyle(news);
  }

  String_V Names()
  {
    String_V result;
    for ( const auto& el : transformers )
    {
        result.push_back(el.first);
    }
    return result;
  }

  G4int HandleStyle( G4coutDestination* dest , const G4String& style)
  {
    const auto& handler = transformers.find(style);
    return ( handler != transformers.end() ) ? (handler->second)(dest) : -1;
  }

  void RegisterNewStyle(const G4String& name , SetupStyle_f& fmt)
  {
    if ( transformers.find(name) != transformers.end() )
    {
        G4ExceptionDescription msg;
        msg << "Format Style with name " << name
            << " already exists. Replacing existing.";
        G4Exception("G4coutFormatters::RegisterNewStyle()",
                    "FORMATTER001", JustWarning, msg);
    }
    // transformers.insert(std::make_pair(name,fmt));
    transformers[name]=fmt;
  }
}
