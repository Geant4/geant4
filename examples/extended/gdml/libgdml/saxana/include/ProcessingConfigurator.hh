//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ProcessingConfigurator.hh,v 1.2 2002-06-03 12:09:32 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef PROCESSING_CONFIGURATOR_H
#define PROCESSING_CONFIGURATOR_H 1

#include <string>

class ProcessingConfigurator
{
public:
  typedef enum { eBatch, eProgressive, eSearch } SAXMode;
  
public:
  ProcessingConfigurator()
  : fMode( eBatch ), fDataFileURI( " " ), fSetupName( "" ), fType( " " ), fVersion( "" )
  {
  }
  
  ProcessingConfigurator( const ProcessingConfigurator& right )
  {
    fMode        = right.Mode();
    fDataFileURI = right.URI();
    fSetupName   = right.SetupName();
    fType        = right.Type();
    fVersion     = right.SetupVersion();
  }
  
  ~ProcessingConfigurator()
  {
  }

  SAXMode Mode() const
	{
	  return fMode;
	}
	
  const std::string& URI() const
  {
    return fDataFileURI;
  }
  
  const std::string& SetupName() const
  {
    return fSetupName;
  }
  
  const std::string& Type() const
  {
    return fType;
  }
  
  const std::string& SetupVersion() const
  {
    return fVersion;
  }
	
	void SetMode( SAXMode mode )
	{
	  fMode = mode;
	}
	
	void SetURI( const std::string& uri )
	{
	  fDataFileURI = uri;
	}
	
	void SetSetupName( const std::string& name )
	{
	  fSetupName = name;
	}
	
	void SetType( const std::string& type )
	{
	  fType = type;
	}
	
	void SetSetupVersion( const std::string& version )
	{
	  fVersion = version;
	}
	
private:
  SAXMode     fMode;
  std::string fDataFileURI;
  std::string fSetupName;
  std::string fType;
  std::string fVersion;
};

#endif // PROCESSING_CONFIGURATOR_H

