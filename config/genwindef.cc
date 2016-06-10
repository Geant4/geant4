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
//$Id: genconf.cpp,v 1.35 2008/10/15 21:51:24 marcocle Exp $	//

#ifdef _WIN32
  // Disable a warning in Boost program_options headers:
  // inconsistent linkage in program_options/variables_map.hpp
  #pragma warning ( disable : 4273 )
  #define popen _popen
  #define pclose _pclose 
  #define fileno _fileno 
  #include <stdlib.h>
#endif

// Include files----------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

// LibSymbolinfo----------------------------------------------------------------
#if !defined(AFX_LIBSYMBOLINFO_H__1A7003B4_BA53_11D1_AE46_1CFB51000000__INCLUDED_)
#define AFX_LIBSYMBOLINFO_H__1A7003B4_BA53_11D1_AE46_1CFB51000000__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include <string>
#include <iostream>

#include <stdio.h>
#include <assert.h>
#include <windows.h>

class CLibSymbolInfo  
{
public:
  CLibSymbolInfo();
  virtual ~CLibSymbolInfo();
  BOOL DumpSymbols(LPTSTR lpszLibPathName, std::ostream& pFile);
  std::string GetLastError() const;

protected:
  std::string m_strResultsString;
  std::string m_strErrorMsg;

  BOOL Dump(LPTSTR lpszLibPathName, std::ostream& pFile); 
  BOOL IsRegularLibSymbol( PSTR pszSymbolName );
  BOOL IsFiltedSymbol( std::string& pszSymbolName );
  DWORD ConvertBigEndian(DWORD bigEndian);
};

enum errMMF {   errMMF_NoError, errMMF_FileOpen,
                errMMF_FileMapping, errMMF_MapView };

class MEMORY_MAPPED_FILE
{
  public:   
  MEMORY_MAPPED_FILE( PSTR pszFileName );
  ~MEMORY_MAPPED_FILE(void);
    
  PVOID   GetBase( void ){ return m_pMemoryMappedFileBase; }
  DWORD   GetFileSize( void ){ return m_cbFile; }
  BOOL    IsValid( void ) { return errMMF_NoError == m_errCode; } 
  errMMF  GetErrorType(){ return m_errCode; }

  private:

  HANDLE      m_hFile;
  HANDLE      m_hFileMapping; // Handle of memory mapped file
  PVOID       m_pMemoryMappedFileBase;
  DWORD       m_cbFile;
  errMMF      m_errCode;  
};

typedef MEMORY_MAPPED_FILE* PMEMORY_MAPPED_FILE;

#endif // !defined(AFX_LIBSYMBOLINFO_H__1A7003B4_BA53_11D1_AE46_1CFB51000000__INCLUDED_)

using namespace std;

#define MakePtr( cast, ptr, addValue ) (cast)( (DWORD)(ptr) + (DWORD)(addValue))

/////////////////////////////////////////////////////////////////////////////
// CLibSymbolInfo

CLibSymbolInfo::CLibSymbolInfo()
{
}

CLibSymbolInfo::~CLibSymbolInfo()
{
}

//=============================================================================
// Function: DumpSymbols
//
// Parameters:
//      LPTSTR lpszLibPathName    - The library file path name
//      CStdioFile* pFile         - Address of the file in which to dump the symbols
//
// Description:
//
// Given a library file path name, the function dumps the symbol info into the file
// pointed to by pFile.
//=============================================================================
BOOL CLibSymbolInfo::DumpSymbols(LPTSTR lpszLibPathName, ostream& pFile)
{
  if(lpszLibPathName == NULL || pFile.bad() ) {
    assert(lpszLibPathName != NULL);
    assert(pFile.good());
    m_strErrorMsg.assign("NULL <lpszLibPathName> or Invalid file handle.");
    return FALSE;
  }

  if(!Dump(lpszLibPathName, pFile))  return FALSE;
  return TRUE;
}

//=============================================================================
// Function: Dump
//
// Parameters:
//      As mentioned above.
//
// Description:
//
// Depending on the value specified in <m_bDump>/<m_bGrep> the routine dumps/greps
// the symbo info.
//=============================================================================
BOOL CLibSymbolInfo::Dump(LPTSTR lpszLibPathName, ostream& pFile) 
{
  string sBuff;
  MEMORY_MAPPED_FILE libFile(lpszLibPathName);
        
  // Ensure that the file mapping worked
  if( FALSE == libFile.IsValid() ) {
        m_strErrorMsg = "Unable to access file ";
    m_strErrorMsg+= lpszLibPathName;
    return FALSE;
  }
  // All COFF libraries start with the string "!<arch>\n".  Verify that this
  // string is at the beginning of the mapped file

  PSTR pArchiveStartString = (PSTR)libFile.GetBase();

  if ( 0 != strncmp( pArchiveStartString, IMAGE_ARCHIVE_START,
                       IMAGE_ARCHIVE_START_SIZE )  )  {
      m_strErrorMsg.assign("Not a valid COFF LIB file.");
        return FALSE;
  }
    
  // Point to the first archive member.  This entry contains the LIB symbols,
  // and immediately follows the archive start string ("!<arch>\n")
  PIMAGE_ARCHIVE_MEMBER_HEADER pMbrHdr;
  pMbrHdr = MakePtr( PIMAGE_ARCHIVE_MEMBER_HEADER, pArchiveStartString,
                     IMAGE_ARCHIVE_START_SIZE );

  // First DWORD after this member header is a symbol count
  PDWORD pcbSymbols = (PDWORD)(pMbrHdr + 1);  // Pointer math!

  // The symbol count is stored in big endian format, so adjust as
  // appropriate for the target architecture
  DWORD cSymbols = ConvertBigEndian( *pcbSymbols );

  // Following the symbol count is an array of offsets to archive members
  // (essentially, embedded .OBJ files)
  PDWORD pMemberOffsets = pcbSymbols + 1;     // Pointer math!

  // Following the array of member offsets is an array of offsets to symbol
  // names.
  PSTR pszSymbolName = MakePtr( PSTR, pMemberOffsets, 4 * cSymbols );

  //
  // Loop through every symbol in the first archive member
  //      
  for ( unsigned i = 0; i < cSymbols; i++ )
  {
    DWORD offset;

    // The offsets to the archive member that contains the symbol is stored
    // in big endian format, so convert it appropriately.        
    offset = ConvertBigEndian( *pMemberOffsets );

    // Call DisplayLibInfoForSymbol, which figures out what kind of symbol
    // it is.  The "IsRegularLibSymbol" filters out symbols that are
    // internal to the linking process
    if ( IsRegularLibSymbol( pszSymbolName ) ) {
      string symbol(pszSymbolName);
      if (IsFiltedSymbol(symbol) ) {
        pFile << symbol << endl;
	  }
    }            
    // Advanced to the next symbol offset and name.  The MemberOffsets
    // array has fixed length entries, while the symbol names are
    // sequential null-terminated strings
    pMemberOffsets++;
    pszSymbolName += strlen(pszSymbolName) + 1;
  }
  return TRUE;
}

//=============================================================================
// Filters out symbols that are internal to the linking process, and which
// the programmer never explicitly sees.
//=============================================================================
BOOL CLibSymbolInfo::IsRegularLibSymbol( PSTR pszSymbolName )
{
  if ( 0 == strncmp( pszSymbolName, "__IMPORT_DESCRIPTOR_", 20 ) )
      return FALSE;

  if ( 0 == strncmp( pszSymbolName, "__NULL_IMPORT_DESCRIPTOR", 24 ) )
      return FALSE;
      
  if ( strstr( pszSymbolName, "_NULL_THUNK_DATA" ) )
      return FALSE;
        
  return TRUE;
}
//=============================================================================
// Filters out symbols that are not needed....
//=============================================================================
BOOL CLibSymbolInfo::IsFiltedSymbol( string& symbolName )
{ 
  // Filter problematic symbols for Win64  
  if ( symbolName.substr(0,3) == "_CT" ) return FALSE;
  if ( symbolName.substr(0,3) == "_TI" ) return FALSE;
  // Filter other symbols
  if ( symbolName.substr(0,2) == "__" ) 
    return FALSE;
  if ( symbolName.substr(0,3) == "??_" && symbolName[3] != '0') // Keep 'operator/='  [??_0]
    return FALSE;
  if( symbolName[0] == '_') {
    symbolName.erase(0, 1);  // C functions ...
  }
  // Filter the internal Boost symbols
  if (symbolName.find ("detail@boost") != string::npos )
        return FALSE;
  if (symbolName.find ("details@boost") != string::npos ) 
        return FALSE;
  return TRUE;
}

//=============================================================================
// Converts from big endian to little endian numbers.
//=============================================================================
DWORD CLibSymbolInfo::ConvertBigEndian(DWORD bigEndian)
{
  DWORD temp = 0;

  temp |= bigEndian >> 24;
  temp |= ((bigEndian & 0x00FF0000) >> 8);
  temp |= ((bigEndian & 0x0000FF00) << 8);
  temp |= ((bigEndian & 0x000000FF) << 24);

  return temp;
}

string CLibSymbolInfo::GetLastError() const
{
  return m_strErrorMsg;
}


MEMORY_MAPPED_FILE::MEMORY_MAPPED_FILE( PSTR pszFileName ) {

   //
   // Given a filename, the constructor opens a file handle, creates a file
   // mapping, and maps the entire file into memory.
   //
   m_hFile = INVALID_HANDLE_VALUE;
   m_hFileMapping = 0;
   m_pMemoryMappedFileBase = 0;
   m_cbFile = 0;
   m_errCode = errMMF_FileOpen;    // Initial error code: not found
    // First get a file handle      
   m_hFile = CreateFile(pszFileName, GENERIC_READ, FILE_SHARE_READ, NULL,
                       OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, (HANDLE)0);
                   
   if ( m_hFile == INVALID_HANDLE_VALUE )
   {
       m_errCode = errMMF_FileOpen;
       return;
   }
    m_cbFile = ::GetFileSize( m_hFile, 0 );
    // Now, create a file mapping   
   m_hFileMapping = CreateFileMapping(m_hFile,NULL, PAGE_READONLY, 0, 0,NULL);
   if ( m_hFileMapping == 0 )
   {
       // Oops.  Something went wrong.  Clean up.
       CloseHandle(m_hFile);
       m_hFile = INVALID_HANDLE_VALUE;
       m_errCode = errMMF_FileMapping;
       return;
   }
    m_pMemoryMappedFileBase = (PCHAR)MapViewOfFile( m_hFileMapping,
                                                   FILE_MAP_READ, 0, 0, 0);
    if ( m_pMemoryMappedFileBase == 0 )
   {
       // Oops.  Something went wrong.  Clean up.
       CloseHandle(m_hFileMapping);
       m_hFileMapping = 0;
       CloseHandle(m_hFile);
       m_hFile = INVALID_HANDLE_VALUE;
       m_errCode = errMMF_MapView;
       return;
   }
    m_errCode = errMMF_NoError;
}

MEMORY_MAPPED_FILE::~MEMORY_MAPPED_FILE(void)
{
    // Clean up everything that was created by the constructor
    if ( m_pMemoryMappedFileBase )
        UnmapViewOfFile( m_pMemoryMappedFileBase );

    if ( m_hFileMapping )
        CloseHandle( m_hFileMapping );

    if ( m_hFile != INVALID_HANDLE_VALUE )
        CloseHandle( m_hFile ); 

    m_errCode = errMMF_FileOpen;
}

namespace windef {
  void usage(){
    cerr << "Usage: genwindef [-l <dllname>] [-o <output-file> | exports.def]  <obj or lib filenames>" << endl;
    exit(1);
  }
}


//--- Command main program-----------------------------------------------------
int main ( int argc, char** argv )
//-----------------------------------------------------------------------------
{
  string outfile("exports.def");
  string library("UnknownLib");
  string objfiles;
  bool debug(false);

  int arg;
  if (argc < 3) windef::usage();
  arg = 1;
  while (argv[arg][0] == '-') {
    if (strcmp(argv[arg], "--") == 0) {
      windef::usage();
    } 
    else if (strcmp(argv[arg], "-l") == 0) {
      arg++; 
      if (arg == argc) windef::usage();
      library = argv[arg];
    } 
    else if (strcmp(argv[arg], "-o") == 0) {
      arg++; 
      if (arg == argc) windef::usage();
      outfile = argv[arg];
    } 
    arg++;
  }
  if (arg == argc) windef::usage();
  for (arg; arg < argc; arg++) {
     objfiles += argv[arg];
     if( arg+1 < argc) objfiles += " ";
  }

  CLibSymbolInfo libsymbols;
  ofstream out(outfile.c_str());
  if(out.fail()) {
    cerr << "windef: Error opening file " << outfile << endl;
    return 1;
  }
  out << "LIBRARY " << library << endl;
  out << "EXPORTS" << endl;

  libsymbols.DumpSymbols(const_cast<char*>(objfiles.c_str()), out);

  out.close();


  return 0;
}
