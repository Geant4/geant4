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
#include <cstring>
#include <sstream>
#include <zmq.hpp>
#include "G4UItcsh.hh"
#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"
#include "G4ZMQServer.hh"

// --------------------------------------------------------------------------
namespace {

G4UImanager* ui_manager = nullptr;
G4bool qexit = false;
std::stringstream cout_stream;
std::string black_str = "\033[30m";
std::string command_list = "";

void ThrowException(const std::string& message)
{
  std::stringstream ss;
  ss << "[ERROR] " << message << std::endl;
  throw std::runtime_error(ss.str());
}

// --------------------------------------------------------------------------
void GetCommandTree(G4UIcommandTree* ctree)
{
  command_list += (ctree-> GetPathName() + " ");

  auto n_cmd = ctree-> GetCommandEntry();
  for ( auto icmd = 1; icmd <= n_cmd; icmd++ ) {
    auto cmd_path = ctree-> GetCommand(icmd)-> GetCommandPath();
    command_list += (cmd_path + " ");
  }

  auto n_tree = ctree-> GetTreeEntry();
  for ( auto itr = 1; itr <= n_tree ; itr++ ) {
    G4UIcommandTree* atree = ctree-> GetTree(itr);
    ::GetCommandTree(atree);
  }
}

} // end of namespace

// --------------------------------------------------------------------------
G4ZMQServer::G4ZMQServer()
{
  endpoint_ = "tcp://127.0.0.1:5555";
  qdebug_ = false;
  shell_= new G4UItcsh();
  shell_-> SetLsColor(BLUE, RED);

  ::ui_manager = G4UImanager::GetUIpointer();
  ::ui_manager-> SetSession(this);
  ::ui_manager-> SetCoutDestination(this);

  ::qexit = false;
}

// --------------------------------------------------------------------------
G4ZMQServer::~G4ZMQServer()
{
  delete shell_;
}

// --------------------------------------------------------------------------
G4UIsession* G4ZMQServer::SessionStart()
{
  zmq::context_t context(1);
  zmq::socket_t socket( context, ZMQ_REP );
  socket.bind(endpoint_);

  enum { kBufferSize = 4096 };
  char buffer[kBufferSize];

  while ( ! ::qexit ) {
    if ( qdebug_ ) {
      std::cout << "@@ Waiting..." << std::endl;
    }

    // waiting command
    zmq::message_t request;
    G4bool qok = socket.recv(&request);
    if ( qok ==  false ) ::ThrowException("G4ZMQSever: socket recv error");
    auto end_pos = request.size();
    if ( end_pos >= kBufferSize  ) end_pos = kBufferSize - 1;
    std::memcpy(buffer, request.data(), end_pos);
    buffer[end_pos] = '\0';
    std::string cmd_str = buffer;

    if ( qdebug_ ) {
      std::cout << "@@ Recv=" << cmd_str << "<<" << std::endl;
    }

    // store output & send back response
    ::cout_stream.str("");

    if ( cmd_str == "@@ping" ) {
      G4cout << "pong" << G4endl;

    } else if ( cmd_str == "@@debug") {
      qdebug_ = true;
      G4cout << "G4ZMQ debug activated" << G4endl;

    } else if ( cmd_str == "@@nodebug") {
      qdebug_ = false;
      G4cout << "G4ZMQ debug deactivated" << G4endl;

    } else if ( cmd_str == "@@get_command_tree" ) {
      auto cwd_name = GetCurrentWorkingDirectory();
      auto cwd_tree = FindDirectory(cwd_name.c_str());
      ::command_list = "";
      ::GetCommandTree(cwd_tree);
      G4cout << ::command_list << std::flush;

    } else if ( cmd_str == "@@get_fullcommand_tree" ) {
      auto root = ::ui_manager-> GetTree();
      ::command_list = "";
      ::GetCommandTree(root);
      G4cout << ::command_list << std::flush;

    } else if ( cmd_str == "help" ) {
      G4cout << "help <command>" << G4endl;

    } else {
      G4String new_command = GetCommand(cmd_str);
      if ( qdebug_ ) {
        std::cout << ::black_str << "@@ Cmd="
                  << new_command << "<<" << std::endl;
      }
      ExecuteCommand(new_command);
    }

    std::string reply = ::cout_stream.str();
    size_t cout_size = reply.size();
    zmq::message_t message(cout_size);
    std::strncpy((char*)message.data(), reply.c_str(), cout_size);
    qok = socket.send(message);
    if ( qok ==  false ) ::ThrowException("G4ZMQServer: socket send error");
  }

  return nullptr;
}

// --------------------------------------------------------------------------
void G4ZMQServer::PauseSessionStart(const G4String& msg)
{
}

// --------------------------------------------------------------------------
G4int G4ZMQServer::ReceiveG4cout(const G4String& coutString)
{
  if ( qdebug_ ) {
    std::cout << coutString << std::flush;
  }

  ::cout_stream << coutString << std::flush;

  return 0;
}

// --------------------------------------------------------------------------
G4int G4ZMQServer::ReceiveG4cerr(const G4String& cerrString)
{
  if ( qdebug_ ) {
    std::cerr << cerrString << std::flush;
  }

  ::cout_stream << cerrString << std::flush;

  return 0;
}

// --------------------------------------------------------------------------
G4String G4ZMQServer::GetCommand(const G4String& input)
{
  const std::string nullstr = "";
  G4String cmdstr = input;

  G4String cstr = cmdstr.strip(G4String::leading);
  if ( cstr.length() == 0 ) {
    cmdstr = nullstr;

  // define built-in shell commands...
  } else if ( cstr(0) == '#' ) {
    G4cout << cstr << G4endl;
    cmdstr = nullstr;

  } else if ( cstr == "ls" || cstr.substr(0,3) == "ls " ) {
    ListDirectory(cstr);
    cmdstr = nullstr;

  } else if ( cstr == "lc" || cstr.substr(0,3) == "lc " ) {
    shell_-> ListCommand(cstr.remove(0,2));
    cmdstr = nullstr;

  } else if (cstr == "pwd" ) {
    G4cout << "Current Command Directory : "
           << GetCurrentWorkingDirectory() << G4endl;
    cmdstr = nullstr;

  } else if ( cstr == "cwd" ) {
    shell_-> ShowCurrentDirectory();
    cmdstr = nullstr;

  } else if (cstr == "cd" || cstr.substr(0,3) == "cd " ) {
    ChangeDirectoryCommand(cstr);
    shell_-> SetCurrentDirectory(GetCurrentWorkingDirectory());
    cmdstr = nullstr;

  } else if ( cstr == "help" || cstr.substr(0,5) == "help " ) {
    TerminalHelp(cstr);
    cmdstr = nullstr;

  } else if ( cstr(0) == '?' ) {
    ShowCurrent(cstr);
    cmdstr = nullstr;

  } else if ( cstr == "history" ) {
    auto nh= ::ui_manager-> GetNumberOfHistory();
    for (auto i = 0; i < nh; i++) {
      G4cout << i << ": " << ::ui_manager->GetPreviousCommand(i) << G4endl;
    }
    cmdstr = nullstr;

  } else if ( cstr == "exit" ) {
    ::qexit = true;
    cmdstr = nullstr;
  }

  return ModifyToFullPathCommand(cmdstr);
}

// --------------------------------------------------------------------------
void G4ZMQServer::ExecuteCommand(const G4String& command)
{
  auto rc = ::ui_manager-> ApplyCommand(command);
  auto pcode = rc % 100;
  auto status = rc - pcode;

  G4UIcommand* cmd = nullptr;
  if( status != fCommandSucceeded ) cmd = FindCommand(command);

  switch ( status ) {
  case fCommandSucceeded:
    break;
  case fCommandNotFound:
    G4cerr << "command <" << ::ui_manager-> SolveAlias(command)
           << "> not found" << G4endl;
    break;
  case fIllegalApplicationState:
    G4cerr << "illegal application state -- command refused" << G4endl;
    break;
  case fParameterOutOfRange:
    G4cerr << "Parameter is out of range" << G4endl;
   break;
  case fParameterOutOfCandidates:
    G4cerr << "Parameter is out of candidate list (index "
           << pcode << ")" << G4endl;
    G4cerr << "Candidates : "
           << cmd-> GetParameter(pcode)-> GetParameterCandidates()
           << G4endl;
    break;
  case fParameterUnreadable:
    G4cerr << "Parameter is wrong type and/or is not omittable (index "
           << pcode << ")" << G4endl;
    break;
  case fAliasNotFound:
    break;
  default:
    G4cerr << "command refused (" << status << ")" << G4endl;
    break;
 }
}

// --------------------------------------------------------------------------
G4bool G4ZMQServer::GetHelpChoice(G4int&)
{
  return true;
}

// --------------------------------------------------------------------------
void G4ZMQServer::ExitHelp() const
{
}
