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
// $Id: StatusCode.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef STATUS_CODE_H
#define STATUS_CODE_H 1

class StatusCode{
public:
  typedef enum
  {
    eError = 0,
    eOk = 1,
    eWarning,
    eParserInitFailure
  } EInfo;

public:
  StatusCode( EInfo info = StatusCode::eOk )
  : fInfo( info )
  {
  }  
  ~StatusCode()
  {
  }
  
  bool IsOk()
  {
    return( (fInfo == eOk) ? true : false );
  } 
  bool IsFailure()
  {
    return( (fInfo != eOk) ? true : false );
  }
  StatusCode::EInfo Info()
  {
    return fInfo;
  }
  
private:
  EInfo fInfo;
};

#endif // STATUS_CODE_H

