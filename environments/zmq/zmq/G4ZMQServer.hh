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
// $Id: G4RunManager.hh 95232 2016-02-01 14:31:22Z gcosmo $
//
//
#ifndef G4ZMQ_SERVER_H_
#define G4ZMQ_SERVER_H_

#include "G4VBasicShell.hh"

class G4UItcsh;

class G4ZMQServer : public  G4VBasicShell {
public:
  G4ZMQServer();
   ~G4ZMQServer();

  void SetEndpoint(const G4String& endpoint);
  G4String GetEndpoint() const;

  void SetDebug(G4bool flag);

  virtual G4UIsession* SessionStart();
  virtual void PauseSessionStart(const G4String& message);

  virtual G4int ReceiveG4cout(const G4String& coutString);
  virtual G4int ReceiveG4cerr(const G4String& cerrString);

private:
  G4bool qdebug_;
  G4String endpoint_;
  G4UItcsh* shell_;

  G4String GetCommand(const G4String& input);

  virtual void ExecuteCommand(const G4String& command);
  virtual G4bool GetHelpChoice(G4int& );
  virtual void ExitHelp() const;

};

// ==========================================================================
inline void G4ZMQServer::SetEndpoint(const G4String& endpoint)
{
  endpoint_ = endpoint;
}

inline G4String G4ZMQServer::GetEndpoint() const
{
  return endpoint_;
}

inline void G4ZMQServer::SetDebug(G4bool flag)
{
  qdebug_ = flag;
}

#endif
