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

