// LibSymbolInfo.h: interface for the CLibSymbolInfo class.
//
//////////////////////////////////////////////////////////////////////

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