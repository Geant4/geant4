

//



//
// $Id: STEPfile.h,v 1.3 1999-12-15 14:50:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef _STEPFILE_H
#define	_STEPFILE_H

/*
* NIST STEP Core Class Library
* cleditor/STEPfile.h
* May 1995
* Peter Carr
* K. C. Morris
* David Sauder

* Development of this software was funded by the United States Government,
* and is not subject to copyright.
*/

/*  */ 
#ifdef __O3DB__
#include <OpenOODB.h>
#endif
#include "globals.hh"
//#include <math.h>


//typedef unsigned boolean;

#include <instmgr.h>
#include <Registry.h>
#include "g4std/fstream"
#include <dirobj.h>
#include <errordesc.h>
#include <time.h>

#include <read_func.h>

//error reporting level
#define READ_COMPLETE    10
#define READ_INCOMPLETE  20

enum  FileTypeCode {
    TYPE_UNKNOWN	= -2,
    VERSION_OLD		= -1,
    VERSION_UNKNOWN	=  0,
    VERSION_CURRENT	=  1,
    WORKING_SESSION	=  2,
    OLD_WORKING_SESSION =  3
  };

class STEPfile
{
  protected:
    //data members

#ifdef __O3DB__
    InstMgr *  _instances;
    Registry * _reg;

    InstMgr & instances ()  { return *_instances; }
    Registry & reg () { return *_reg; }
#else
    InstMgr&  _instances;
    Registry& _reg;

    InstMgr & instances ()  { return _instances; }
    Registry & reg () { return _reg; }
#endif
    int _fileIdIncr;   //Increment value to be added to FileId Numbers on input

//header information
    InstMgr*  _headerInstances;
    Registry *_headerRegistry;
//Registry *_headerUserDefined;
//    Registry *_headerRegistryOld;
    
    int _headerId;     //STEPfile_id given to STEPentity from header section

//file information
    DirObj* _currentDir;
    char* _fileName;

//error information
    ErrorDescriptor _error;

    // new errors
    int _entsNotCreated; // num entities not created in first pass
    int _entsInvalid;    // num entities that had invalid attr values
    int _entsIncomplete; // num entities that had missing attr values
			 // (includes entities that had invalid values
			 // for required attrs)
    int _entsWarning; // num entities that may have had problems 
		      // with attrs - reported as an attr user msg

    // old errors
    int _errorCount;
    int _warningCount;

    int _maxErrorCount;

  protected:
    
//file type information
    FileTypeCode _fileType;
    char ENTITY_NAME_DELIM;
    char* FILE_DELIM;    
    char* END_FILE_DELIM;
    
//public member functions
  public:

//public access to member variables
//header information
    InstMgr* HeaderInstances() { return _headerInstances; }

//file information
    const char* FileName() const { return _fileName; }
    const char* SetFileName (const char* name = "");
    const char* TruncFileName (const char* name) const;

//error information
    ErrorDescriptor& Error() /* const */  { return _error;        }
    int ErrorCount() const  { return _errorCount;   }
    int WarningCount() const { return _warningCount; }
    Severity AppendEntityErrorMsg (ErrorDescriptor *e);	
    
//version information
    FileTypeCode FileType() const   { return _fileType; }
    void FileType (FileTypeCode ft) { _fileType = ft; }	
    int SetFileType (FileTypeCode ft = VERSION_CURRENT);
    
//Reading and Writing 
    Severity ReadExchangeFile (const char* filename =0);
    Severity AppendExchangeFile (const char* filename =0); 

    Severity ReadWorkingFile (const char* filename =0);
    Severity AppendWorkingFile (const char* filename =0);

    Severity AppendFile (G4std::istream* in) ;

    Severity WriteExchangeFile (G4std::ostream& out, int validate =1);
    Severity WriteExchangeFile (const char* filename =0, int validate =1);

    Severity WriteWorkingFile (G4std::ostream& out);
    Severity WriteWorkingFile (const char* filename =0);

    stateEnum EntityWfState(char c);
    
    void Renumber ();

//constructors
    STEPfile (Registry& r, InstMgr& i, const char *filename = (const char*)0);
    ~STEPfile();

 protected:    
//member functions
    
//called by ReadExchangeFile
    G4std::istream* OpenInputFile (const char* filename = "");
    void CloseInputFile(G4std::istream* in);
    
    Severity ReadHeader(G4std::istream& in);

    InstMgr* HeaderConvertToNew(InstMgr& oldinst);
    Severity HeaderVerifyInstances(InstMgr* im);
    void HeaderMergeInstances(InstMgr* im);
    STEPentity* HeaderDefaultFileName();	
    STEPentity* HeaderDefaultFileDescription();	
    STEPentity* HeaderDefaultFileSchema();	
  
    int HeaderId (int Increment ); // =1);
    int HeaderId (const char* nm); //  ="\0");  // DELETED both defaults !!
                                                //   J.Apostolakis Oct 30, 98
    int HeaderIdOld (const char* nm ="\0");

    int ReadData1 (G4std::istream& in); // first pass to create instances
    int ReadData2 (G4std::istream& in); // second pass to read instances

// obsolete
    int ReadWorkingData1 (G4std::istream& in);
    int ReadWorkingData2 (G4std::istream& in);

    void ReadRestOfFile(G4std::istream& in);

	// create instance - used by ReadData1()
    STEPentity *  CreateInstance(G4std::istream& in, G4std::ostream& out);
	// create complex instance - used by CreateInstance()
    STEPentity * CreateSubSuperInstance(G4std::istream& in, int fileid);

	// read the instance - used by ReadData2()
    STEPentity * ReadInstance (G4std::istream& in, G4std::ostream& out, SCLstring &cmtStr);

  //  reading scopes are still incomplete
  //  these functions are stubs
    Severity CreateScopeInstances(G4std::istream& in, STEPentityH ** scopelist);
    Severity ReadScopeInstances(G4std::istream& in);
//    Severity ReadSubSuperInstance(G4std::istream& in);

    int FindDataSection (G4std::istream& in);
    int FindHeaderSection (G4std::istream& in);

// writing working session files
    void WriteWorkingData(G4std::ostream& out);

//called by WriteExchangeFile
    G4std::ofstream* OpenOutputFile(const char* filename =0);
    void CloseOutputFile(G4std::ostream* out);

    void WriteHeader (G4std::ostream& out);
    void WriteHeaderInstance (STEPentity *obj, G4std::ostream& out);
    void WriteHeaderInstanceFileName (G4std::ostream& out);
    void WriteHeaderInstanceFileDescription (G4std::ostream& out);
    void WriteHeaderInstanceFileSchema (G4std::ostream& out);

    void WriteData (G4std::ostream& out);
    
    int IncrementFileId (int fileid);
    int FileIdIncr() { return _fileIdIncr; }
    void SetFileIdIncrement ();
    void MakeBackupFile();

//    void ReadWhiteSpace(G4std::istream& in);
};

//inline functions
#endif


