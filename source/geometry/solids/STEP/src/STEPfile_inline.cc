

//



//
// $Id: STEPfile_inline.cc,v 1.2 1999-05-21 20:20:51 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

/*
* NIST STEP Core Class Library
* cleditor/STEPfile.inline.cc
* May 1995
* Peter Carr
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 

#include <STEPfile.h>

#ifdef __GNUG__
extern "C" { 
double  ceil(double);
}
#endif

#ifdef __SUNCPLUSPLUS__
extern "C" { 
double  ceil(double);
}
#endif

#ifdef __OBJECTCENTER__
extern "C" { 
double  ceil(double);
}
#endif

extern void HeaderSchemaInit (Registry & reg);

//To Be inline functions

//constructor & destructor

STEPfile::STEPfile(Registry& r, InstMgr& i, const char *filename)

#ifdef __O3DB__
: _reg(&r), _instances(&i),
#else
: _reg(r), _instances(i), 
#endif
  _headerId(0), _maxErrorCount(5000), 
  _fileName (0), _entsNotCreated(0), _entsInvalid(0), 
  _entsIncomplete(0), _entsWarning(0), 
  _errorCount (0), _warningCount (0),
  _fileIdIncr (0)

{ 
    SetFileType(VERSION_CURRENT);
    SetFileIdIncrement(); 
    _currentDir = new DirObj("");
//    _headerRegistry = new Registry(&s_HeaderSchemaInit);
    _headerRegistry = new Registry(&HeaderSchemaInit);
    _headerInstances = new InstMgr;
    if (filename) ReadExchangeFile(filename);
}

STEPfile::~STEPfile() 
{
    delete _currentDir;

//  remove everything from the Registry before deleting it
    _headerRegistry -> DeleteContents ();
    delete _headerRegistry;

//DAS    delete _headerRegistryOld;
    _headerInstances -> ClearInstances ();
    delete _headerInstances;
}

int
STEPfile::SetFileType(FileTypeCode ft) 
{
    FileType (ft);

    switch (_fileType) 
      {
	case (VERSION_OLD):
	  ENTITY_NAME_DELIM = '@';
	  FILE_DELIM = "STEP;";
	  END_FILE_DELIM = "ENDSTEP;";
/*DAS
	  if (!_headerRegistryOld) 
	      _headerRegistryOld = 
			new Registry(& s_Header_Section_Schema_N279Init);
*/

	  break;
	case (VERSION_UNKNOWN):
	case (VERSION_CURRENT):
	  ENTITY_NAME_DELIM = '#';
	  FILE_DELIM = "ISO-10303-21;";
	  END_FILE_DELIM = "END-ISO-10303-21;";
	  break;
	case (WORKING_SESSION):
	  ENTITY_NAME_DELIM = '#';
	  FILE_DELIM = "STEP_WORKING_SESSION;";
	  END_FILE_DELIM = "END-STEP_WORKING_SESSION;";
	  break;

	default:
	  // some kind of error
	  G4cerr << "Internal error:  " << __FILE__ <<  __LINE__ 
	       << "\n" << _POC_ "\n";
	  return 0;
      }
    return 1;
}

 
/******************************************************/
const char*
STEPfile::TruncFileName(const char* filename) const
{
    char* tmp = strrchr(filename,'/');
    if (tmp) return tmp++;
    else return filename;
    
}


/******************************************************/
Severity
STEPfile::ReadExchangeFile(const char* filename)
{
    _error.ClearErrorMsg();
    _errorCount = 0;
    istream* in = OpenInputFile(filename);
    if (_error.severity() < SEVERITY_WARNING) 
      {	  
	CloseInputFile(in);
	return _error.severity();  
      }
    
    instances ().ClearInstances ();
    if (_headerInstances)
	_headerInstances->ClearInstances ();
    _headerId = 5;
    Severity rval = AppendFile (in);
    CloseInputFile(in);
    return rval;
}

Severity 
STEPfile::AppendExchangeFile (const char* filename)
{
    _error.ClearErrorMsg();
    _errorCount = 0;
    istream* in = OpenInputFile(filename);
    if (_error.severity() < SEVERITY_WARNING) 
      {	      
	CloseInputFile(in);
	return _error.severity();  
      }
    Severity rval = AppendFile (in);
    CloseInputFile(in);
    return rval;
}

/******************************************************/
Severity
STEPfile::ReadWorkingFile(const char* filename) 
{
    _error.ClearErrorMsg();
    _errorCount = 0;
    istream* in = OpenInputFile(filename);
    if (_error.severity() < SEVERITY_WARNING) 
      {	  
	CloseInputFile(in);
	return _error.severity();  
      }

    instances ().ClearInstances();
    _headerInstances->ClearInstances ();
    SetFileType(WORKING_SESSION);

    Severity rval = AppendFile (in);
    SetFileType();
    CloseInputFile(in);
    return rval;
}


Severity
STEPfile::AppendWorkingFile(const char* filename)
{
    _error.ClearErrorMsg();
    _errorCount = 0;
    istream* in = OpenInputFile(filename);
    if (_error.severity() < SEVERITY_WARNING) 
      {	      
	CloseInputFile(in);
	return _error.severity();  
      }
    SetFileType(WORKING_SESSION);
    Severity rval = AppendFile (in);
    SetFileType();
    CloseInputFile(in);
    return rval;
}



/******************************************************/
istream*
STEPfile::OpenInputFile (const char* filename)
{
    //  if there's no filename to use, fail
    if (! (strcmp (filename, "") || strcmp (FileName (), "")) ) 
      {
	  _error.AppendToUserMsg("Unable to open file for input. No current file name.\n");
	  _error.GreaterSeverity(SEVERITY_INPUT_ERROR);
	  return (0);
      }		
    else  {
	if (!SetFileName (filename)) 
	  {
	      char msg[BUFSIZ];
	      sprintf(msg,"Unable to find file for input: \'%s\'. File not read.\n",filename);
	      _error.AppendToUserMsg(msg);
	      _error.GreaterSeverity(SEVERITY_INPUT_ERROR);
	      return (0);
	  }
    }
    //  istream* in = new istream(FileName(), io_readonly, a_useonly);
    // port 29-Mar-1994 kcm
    istream* in = new ifstream(FileName());
    // default for ostream is readonly and protections are set to 644 
//    if ( !in || !(in -> readable ()) )
    if ( !in || !(in -> good ()) )
      {
	      char msg[BUFSIZ];
	      sprintf(msg,"Unable to open file for input: \'%s\'. File not read.\n",filename);
	      _error.AppendToUserMsg(msg);
	      _error.GreaterSeverity(SEVERITY_INPUT_ERROR);
	      return (0);
      }
    return in;
}

/******************************************************/
void
STEPfile::CloseInputFile(istream* in)
{
    delete in;
}

    
/******************************************************/

/*
void
STEPfile::ReadWhiteSpace (istream& in)
{

  char c = ' ';
  while ((c == ' ') || (c == '\n') || (c == '\t'))  {
    in >> c; 
  }
  in.putback (c);
}
*/

/***************************
***************************/
ofstream*
STEPfile::OpenOutputFile(const char* filename)
{
    if (!filename) 
      {
	  if (!FileName())
	    { 
		_error.AppendToUserMsg("No current file name.\n");
		_error.GreaterSeverity(SEVERITY_INPUT_ERROR);
	    }
      }
    else 
      {
	  if (!SetFileName (filename)) 
	    {
		char msg[BUFSIZ];
		sprintf(msg,"can't find file: %s\nFile not written.\n",filename);
		_error.AppendToUserMsg(msg);
		_error.GreaterSeverity(SEVERITY_INPUT_ERROR);
	    }	      
      }

    if (_currentDir->FileExists(TruncFileName(FileName())))
	MakeBackupFile();
//    ostream* out  = new ostream(FileName(), io_writeonly, a_create);  
    // - port 29-Mar-1994 kcm
    ofstream* out  = new ofstream(FileName());  
    // default for ostream is writeonly and protections are set to 644 
    if (!out) 
      {
	  _error.AppendToUserMsg("unable to open file for output\n");
	  _error.GreaterSeverity(SEVERITY_INPUT_ERROR);
      }
    return out;
}

void
STEPfile::CloseOutputFile(ostream* out) 
{
    delete out;
}



/******************************************************/
int 
STEPfile::IncrementFileId (int fileid) 
{ 
    return (fileid + FileIdIncr()); 
}


void 
STEPfile::SetFileIdIncrement()
{

  if (instances ().MaxFileId() < 0) _fileIdIncr = 0;
  else
    {
      double instMaxFileId = (instances ().MaxFileId()+ 99)/1000;
      _fileIdIncr = (int)((ceil(instMaxFileId) + 1) * 1000);
    }
}
