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
// $Id: StateProcessMap.hh,v 1.2 2002-06-03 12:09:33 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------
//
#ifndef STATE_PROCESS_MAP_H
#define STATE_PROCESS_MAP_H 1

#include <string>
#include <map>

#include "RCObjectHandle.hh"
#include "SAXStateProcess.hh"

class StateProcessMap
{
public:
  typedef RCObjectHandle<SAXStateProcess>                 Process;
  typedef std::map<std::string,Process>                   ProcessMap;

public:
  StateProcessMap();
  ~StateProcessMap();
  
  void Initialize();
  void Reset();
  
  void AddProcess( const std::string& tag, Process process );
  void AddProcess( char*        tag, Process process );
  bool Check( const std::string& tag );
  bool Check( char*        tag );
  Process GetProcess( const std::string& tag );
  Process GetProcess( char*        tag );

private:
  ProcessMap fMap;
};

#endif // STATE_PROCESS_MAP_H

