//==========================================
// Matt Pietrek
// Microsoft Systems Journal, April 1998
// Program: LibDump
// FILE: LibSymbolInfo.CPP
//==========================================

// LibSymbolInfo.cpp: implementation of the CLibSymbolInfo class.
//
//////////////////////////////////////////////////////////////////////


#include "LibSymbolInfo.h"

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
