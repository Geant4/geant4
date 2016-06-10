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
//
// $Id: G4MTcoutDestination.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 
// ----------------------------------------------------------------------
// G4MTcoutDestination
//

#include "G4MTcoutDestination.hh"
#include "G4strstreambuf.hh"
#include "G4AutoLock.hh"

G4MTcoutDestination::G4MTcoutDestination(const G4int& threadId,
                              std::ostream& co, std::ostream&  ce)
: finalcout(co), finalcerr(ce), id(threadId), useBuffer(false),
  threadCoutToFile(false), threadCerrToFile(false),
  ignoreCout(false), ignoreInit(true)
{
  G4coutbuf.SetDestination(this);
  G4cerrbuf.SetDestination(this);
  prefix = "G4WT";
}

G4MTcoutDestination::~G4MTcoutDestination()
{
  if( useBuffer ) DumpBuffer();
  if( threadCoutToFile ) CloseCoutFile();
  if( threadCerrToFile ) CloseCerrFile();
}

namespace  { G4Mutex coutm = G4MUTEX_INITIALIZER; }
#include "G4StateManager.hh"

G4int G4MTcoutDestination::ReceiveG4cout(const G4String& msg)
{
  if( threadCoutToFile )
  { coutFile<<msg<<std::flush; }
  else if( useBuffer )
  { cout_buffer<<msg; }
  else if( !ignoreCout )
  {
    if(!ignoreInit || 
       G4StateManager::GetStateManager()->GetCurrentState() != G4State_Idle )
    {
      G4AutoLock l(&coutm);
        finalcout<<prefix;
        if ( id!=G4Threading::GENERICTHREAD_ID ) finalcout<<id;
        finalcout<<" > "<<msg<<std::flush;
    }
  }
  //forward message to master G4coutDestination if set
  if ( masterG4coutDestination &&  !ignoreCout &&
       ( !ignoreInit || G4StateManager::GetStateManager()->GetCurrentState() != G4State_Idle )
    ){
        G4AutoLock l(&coutm);
        std::stringstream ss;
      ss<<prefix;
      if ( id!=G4Threading::GENERICTHREAD_ID ) ss<<id;
      ss<<" > "<<msg;
      masterG4coutDestination->ReceiveG4cout(ss.str());
    }
  return 0;
}

G4int G4MTcoutDestination::ReceiveG4cerr(const G4String& msg)
{
  if( threadCerrToFile )
  { cerrFile<<msg<<std::flush; }
  if( useBuffer )
  { cerr_buffer<<msg; }
  else
    {   G4AutoLock l(&coutm);
        finalcerr<<prefix;
        if ( id!=G4Threading::GENERICTHREAD_ID ) finalcerr<<id;
        finalcerr<<" > "<<msg<<std::flush;
    }
  //forward message to master G4coutDestination if set
    if ( masterG4coutDestination &&  !ignoreCout &&
        ( !ignoreInit || G4StateManager::GetStateManager()->GetCurrentState() != G4State_Idle )
    ){
        G4AutoLock l(&coutm);
        std::stringstream ss;
        ss<<prefix;
        if ( id!=G4Threading::GENERICTHREAD_ID ) ss<<id;
        ss<<" > "<<msg;
        masterG4coutDestination->ReceiveG4cerr(ss.str());
    }
  return 0;
}

void G4MTcoutDestination::SetCoutFileName(const G4String& fileN, G4bool ifAppend)
{
  if( threadCoutToFile ) CloseCoutFile();
  if( fileN == "**Screen**" ) return;
  if( ! coutFile.is_open() )
  {
    std::ios::openmode mode = std::ios::out;
    if ( ifAppend ) mode |= std::ios::app;
    coutFile.open(fileN,mode);
  }
  threadCoutToFile = true;
}

void G4MTcoutDestination::SetCerrFileName(const G4String& fileN, G4bool ifAppend)
{
  if( threadCerrToFile ) CloseCerrFile();
  if( fileN == "**Screen**" ) return;
  if( ! cerrFile.is_open() )
  {
    std::ios::openmode mode = std::ios::out;
    if ( ifAppend ) mode |= std::ios::app;
    cerrFile.open(fileN,mode);
  }
  threadCerrToFile = true;
}

void G4MTcoutDestination::EnableBuffering(G4bool flag)
{
  if(useBuffer && !flag) DumpBuffer();
  useBuffer = flag;
}

void G4MTcoutDestination::SetPrefixString(const G4String& wd)
{ prefix = wd; }

void G4MTcoutDestination::SetIgnoreCout(G4int tid)
{
 if(tid<0)
 { ignoreCout = false; }
 else
 { ignoreCout = (tid!=id); }
}

void G4MTcoutDestination::CloseCoutFile()
{
  if( coutFile.is_open() ) coutFile.close();
  threadCoutToFile = false;
}

void G4MTcoutDestination::CloseCerrFile()
{
  if( cerrFile.is_open() ) cerrFile.close(); 
  threadCerrToFile = false;
}

void G4MTcoutDestination::DumpBuffer()
{
  G4AutoLock l(&coutm);
  finalcout<<"====================="<<std::endl;
  finalcout<<"cout buffer for worker with ID:"<<id<<std::endl;
  finalcout<<cout_buffer.str()<<std::endl;
  finalcerr<<"====================="<<std::endl;
  finalcerr<<"cerr buffer for worker with ID:"<<id<<std::endl;
  finalcerr<<cerr_buffer.str()<<std::endl;
  finalcerr<<"====================="<<std::endl;
}

