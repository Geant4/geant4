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
// $Id: G4MTcoutDestination.cc 66241 2012-12-13 18:34:42Z gunter $
// 
// --------------------------------------------------------------------
//
// G4MTcoutDestination.cc
//
// --------------------------------------------------------------------

#include <sstream>
#include <assert.h>

#include "G4MTcoutDestination.hh"
#include "G4LockcoutDestination.hh"
#include "G4MasterForwardcoutDestination.hh"
#include "G4FilecoutDestination.hh"
#include "G4BuffercoutDestination.hh"
#include "G4strstreambuf.hh"
#include "G4AutoLock.hh"

namespace
{
  G4String empty = "";
}

G4MTcoutDestination::G4MTcoutDestination(const G4int& threadId)
  : ref_defaultOut(nullptr), ref_masterOut(nullptr),
    masterDestinationFlag(true),masterDestinationFmtFlag(true),
    id(threadId), useBuffer(false), ignoreCout(false), ignoreInit(true),
    prefix("G4WT")
{
  // TODO: Move these two out of here and in the caller
  G4coutbuf.SetDestination(this);
  G4cerrbuf.SetDestination(this);

  stateMgr=G4StateManager::GetStateManager();
  SetDefaultOutput(masterDestinationFlag,masterDestinationFmtFlag);
}

void G4MTcoutDestination::SetDefaultOutput( G4bool addmasterDestination ,
                                            G4bool formatAlsoMaster )
{
  masterDestinationFlag = addmasterDestination;
  masterDestinationFmtFlag = formatAlsoMaster;
  // Formatter: add prefix to each thread
  const auto f = [this](G4String& msg)->G4bool {
    std::ostringstream str;
    str<<prefix;
    if ( id!=G4Threading::GENERICTHREAD_ID ) str<<id;
    str<<" > "<<msg;
    msg = str.str();
    return true;
  };
  // Block cout if not in correct state
  const auto filter_out = [this](G4String&)->G4bool {
    if (this->ignoreCout ||
        ( this->ignoreInit &&
          this->stateMgr->GetCurrentState() == G4State_Idle ) )
      { return false; }
    return true;
  };

  // Default behavior, add a destination that uses cout and uses a mutex
  auto output = G4coutDestinationUPtr( new G4LockcoutDestination );
  ref_defaultOut = output.get();
  output->AddCoutTransformer(filter_out);
  output->AddCoutTransformer(f);
  output->AddCerrTransformer(f);
  push_back( std::move(output) );
  if ( addmasterDestination )
  {
      AddMasterOutput(formatAlsoMaster);
  }
}

void G4MTcoutDestination::AddMasterOutput(G4bool formatAlsoMaster )
{
  // Add a destination, that forwards the message to the master thread
  auto forwarder = G4coutDestinationUPtr( new G4MasterForwardcoutDestination );
  ref_masterOut = forwarder.get();
  const auto filter_out = [this](G4String&)->G4bool {
    if (this->ignoreCout ||
        ( this->ignoreInit &&
          this->stateMgr->GetCurrentState() == G4State_Idle ) )
      { return false; }
    return true;
  };
  forwarder->AddCoutTransformer(filter_out);
  if ( formatAlsoMaster )
  {
      // Formatter: add prefix to each thread
      const auto f = [this](G4String& msg)->G4bool {
        std::ostringstream str;
        str<<prefix;
        if ( id!=G4Threading::GENERICTHREAD_ID ) str<<id;
        str<<" > "<<msg;
        msg = str.str();
        return true;
      };
      forwarder->AddCoutTransformer(f);
      forwarder->AddCerrTransformer(f);
  }
  push_back( std::move(forwarder ) );

}

G4MTcoutDestination::~G4MTcoutDestination()
{
  if ( useBuffer ) DumpBuffer();
}

void G4MTcoutDestination::Reset()
{
  clear();
  SetDefaultOutput(masterDestinationFlag,masterDestinationFmtFlag);
}

void G4MTcoutDestination::HandleFileCout(G4String fileN, G4bool ifAppend,
                                         G4bool suppressDefault)
{
  // Logic: we create a file destination. We want this to get only the G4cout
  // stream and should discard everything in G4cerr.
  // First we create the destination with the appropriate open mode

  std::ios_base::openmode mode = (ifAppend ? std::ios_base::app
                                           : std::ios_base::trunc);
  auto output = G4coutDestinationUPtr( new G4FilecoutDestination(fileN,mode));

  // This reacts only to G4cout, so let's make a filter that removes everything
  // from G4cerr
  output->AddCerrTransformer( [](G4String&) { return false;} );
  push_back(std::move(output));
  // Silence G4cout from default formatter
  if ( suppressDefault )
  {
     ref_defaultOut->AddCoutTransformer( [](G4String&) { return false; } );
     if ( ref_masterOut )
       ref_masterOut->AddCoutTransformer( [](G4String&) { return false; } );
  }
}

void G4MTcoutDestination::HandleFileCerr(G4String fileN, G4bool ifAppend,
                                         G4bool suppressDefault)
{
  // See HandleFileCout for explanation, switching cout with cerr

  std::ios_base::openmode mode = (ifAppend ? std::ios_base::app
                                           : std::ios_base::trunc);
  auto output = G4coutDestinationUPtr( new G4FilecoutDestination(fileN,mode));
  output->AddCoutTransformer( [](G4String&) { return false;} );
  push_back(std::move(output));
  if ( suppressDefault )
  {
     ref_defaultOut->AddCerrTransformer( [](G4String&) { return false; } );
     if ( ref_masterOut )
       ref_masterOut->AddCerrTransformer( [](G4String&) { return false; } );
  }
}

void G4MTcoutDestination::SetCoutFileName(const G4String& fileN,
                                                G4bool ifAppend)
{
  // First let's go back to the default
  Reset();
  if ( fileN != "**Screen**" )
  {
      HandleFileCout(fileN,ifAppend,true);
  }
}

void G4MTcoutDestination::EnableBuffering(G4bool flag)
{
  // I was using buffered output and now I want to turn it off, dump current
  // buffer content and reset output
  if ( useBuffer && !flag )
  {
      DumpBuffer();
      Reset();
  }
  else if ( useBuffer && flag ) { /* do nothing: already using */ }
  else if ( !useBuffer  && !flag ) { /* do nothing: not using */ }
  else if ( !useBuffer && flag )
  {
     // Remove everything, in this case also removing the forward to the master
     // thread, we want everything to be dumpled to a file
     clear();
     const size_t infiniteSize = 0;
     push_back(G4coutDestinationUPtr(new G4BuffercoutDestination(infiniteSize)));
  }
  else { assert(false); } // Should never happen
  useBuffer = flag;
}

void G4MTcoutDestination::AddCoutFileName(const G4String& fileN,
                                                G4bool ifAppend)
{
  // This is like the equivalent SetCoutFileName, but in this case we do not
  // remove or silence what is already exisiting
  HandleFileCout(fileN,ifAppend,false);
}

void G4MTcoutDestination::SetCerrFileName(const G4String& fileN,
                                                G4bool ifAppend)
{
  // See SetCoutFileName for explanation
  Reset();
  if ( fileN != "**Screen**")
  {
      HandleFileCerr(fileN,ifAppend,true);
  }
}

void G4MTcoutDestination::AddCerrFileName(const G4String& fileN,
                                                G4bool ifAppend)
{
  HandleFileCerr(fileN,ifAppend,false);
}

void G4MTcoutDestination::SetIgnoreCout(G4int tid)
{
 if (tid<0) { ignoreCout = false; }
 else   { ignoreCout = (tid!=id); }
}

namespace
{
  G4Mutex coutm = G4MUTEX_INITIALIZER;
}

void G4MTcoutDestination::DumpBuffer()
{
  G4AutoLock l(&coutm);
  std::ostringstream msg;
  msg << "=======================\n";
  msg << "cout buffer(s) for worker with ID:" << id << std::endl;
  G4coutDestination::ReceiveG4cout(msg.str());
  G4bool sep = false;
  std::for_each( begin() , end(),
      [this,&sep](G4coutDestinationUPtr& el) {
        auto cout = dynamic_cast<G4BuffercoutDestination*>(el.get());
        if ( cout != nullptr ) {
            cout->FlushG4cout();
            if ( sep ) { G4coutDestination::ReceiveG4cout("==========\n"); }
            else { sep = true; }
        }
  } );
  sep = false;
  msg.str("");
  msg.clear();
  msg << "=======================\n";
  msg << "cerr buffer(s) for worker with ID:" << id
      << " (goes to std error)" << std::endl;
  G4coutDestination::ReceiveG4cout(msg.str());
  std::for_each( begin() , end(),
      [this,&sep](G4coutDestinationUPtr& el) {
        auto cout = dynamic_cast<G4BuffercoutDestination*>(el.get());
        if ( cout != nullptr ) {
            cout->FlushG4cerr();
            if (sep ) { G4coutDestination::ReceiveG4cout("==========\n"); }
            else { sep = true; }
        }
  } );
  G4coutDestination::ReceiveG4cout("=======================\n");
}
