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
