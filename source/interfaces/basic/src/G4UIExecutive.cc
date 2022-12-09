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

#include "G4UIExecutive.hh"
#include "G4UIsession.hh"
#include "G4UImanager.hh"

#if defined(G4UI_BUILD_QT_SESSION)
#include "G4UIQt.hh"
#include "G4Qt.hh"
#endif

#if defined(G4UI_BUILD_XM_SESSION)
#include "G4UIXm.hh"
#endif

#if defined(G4UI_BUILD_WIN32_SESSION)
#include "G4UIWin32.hh"
#endif

#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIcsh.hh"
#include "G4TiMemory.hh"

// --------------------------------------------------------------------------
// build flags as variables

#if defined(G4UI_BUILD_QT_SESSION)
static const G4bool qt_build = true;
#else
static const G4bool qt_build = false;
#endif

#if defined(G4UI_BUILD_XM_SESSION)
static const G4bool xm_build = true;
#else
static const G4bool xm_build = false;
#endif

#if defined(G4UI_BUILD_WIN32_SESSION)
static const G4bool win32_build = true;
#else
static const G4bool win32_build = false;
#endif

#ifndef WIN32
static const G4bool tcsh_build = true;
#else
static const G4bool tcsh_build = false;
#endif

#define DISCARD_PARAMETER(p) (void)p

// --------------------------------------------------------------------------
G4UIExecutive::G4UIExecutive(G4int argc, char** argv, const G4String& type)
  : selected(kNone), session(NULL), shell(NULL), isGUI(false), verbose(true)
{
  if ( verbose ) {
    G4cout << "Available UI session types: [ ";
    if ( qt_build ) G4cout << "Qt, ";
    if ( xm_build ) G4cout << "Xm, ";
    if ( win32_build) G4cout << "Win32, ";
    if (tcsh_build ) G4cout << "tcsh, ";
    G4cout << "csh ]" << G4endl;
  }

  // selecting session type...
  // 1st priority : in case argumant specified
  G4String stype = G4StrUtil::to_lower_copy(type); // session type is case-insensitive.
  if (type != "") SelectSessionByArg(stype);

  // 2nd priority : refer environment variables (as backword compatibility)
  if ( selected == kNone ) SelectSessionByEnv();

  // 3rd priority : refer $HOME/.g4session
  if ( selected == kNone ) {
    G4String appinput = argv[0];
    G4String appname = "";
    size_t islash = appinput.find_last_of("/\\");
    if (islash == G4String::npos)
      appname = appinput;
    else
      appname = appinput.substr(islash+1, appinput.size()-islash-1);

    SelectSessionByFile(appname);
  }

  // 4th, best guess of session type
  if ( selected == kNone) SelectSessionByBestGuess();

  // instantiate a session...
  switch ( selected ) {
  case kQt:
#if defined(G4UI_BUILD_QT_SESSION)
    session = new G4UIQt(argc, argv);
    isGUI = true;
#endif
    break;
  case kXm:
#if defined(G4UI_BUILD_XM_SESSION)
    session = new G4UIXm(argc, argv);
    isGUI = true;
#endif
    break;
  case kWin32:
#if defined(G4UI_BUILD_WIN32_SESSION)
    DISCARD_PARAMETER(argc);
    DISCARD_PARAMETER(argv);
    session = new G4UIWin32();
    isGUI = true;
#endif
    break;
 case kTcsh:
#if !(defined(WIN32) || defined(__MINGW32__))
    DISCARD_PARAMETER(argc);
    DISCARD_PARAMETER(argv);
    shell = new G4UItcsh;
    session = new G4UIterminal(shell);
#endif
    break;
  case kCsh:
    DISCARD_PARAMETER(argc);
    DISCARD_PARAMETER(argv);
    shell = new G4UIcsh;
    session = new G4UIterminal(shell);
  default:
    break;
  }

  // fallback (csh)
  if ( session == NULL ) {
    G4Exception("G4UIExecutive::G4UIExecutive()",
                "UI0002",
                JustWarning,
                "Specified session type is not build in your system,\n"
                "or no session type is specified.\n"
                "A fallback session type is used.");

    selected = kCsh;
    DISCARD_PARAMETER(argc);
    DISCARD_PARAMETER(argv);
    shell = new G4UIcsh;
    session = new G4UIterminal(shell);
  }

  TIMEMORY_INIT(argc, argv);
}

// --------------------------------------------------------------------------
G4UIExecutive::~G4UIExecutive()
{
  delete session;
}

// --------------------------------------------------------------------------
void G4UIExecutive::SelectSessionByArg(const G4String& stype)
{
  if ( qt_build && stype == "qt" ) selected = kQt;
  else if ( xm_build && stype == "xm" ) selected = kXm;
  else if ( win32_build && stype == "win32" ) selected = kWin32;
  else if ( tcsh_build && stype == "tcsh" ) selected = kTcsh;
  else if ( stype == "csh" ) selected = kCsh;
}

// --------------------------------------------------------------------------
void G4UIExecutive::SelectSessionByEnv()
{
  if ( qt_build && std::getenv("G4UI_USE_QT") ) selected = kQt;
  else if ( xm_build && std::getenv("G4UI_USE_XM") ) selected = kXm;
  else if ( win32_build && std::getenv("G4UI_USE_WIN32") ) selected = kWin32;
  else if ( tcsh_build && std::getenv("G4UI_USE_TCSH") ) selected = kTcsh;
}

// --------------------------------------------------------------------------
void G4UIExecutive::SelectSessionByFile(const G4String& appname)
{
  const char* path = std::getenv("HOME");
  if( path == NULL ) return;
  G4String homedir = path;

#ifndef WIN32
  G4String fname= homedir + "/.g4session";
#else
  G4String fname= homedir + "\\.g4session";
#endif

  std::ifstream fsession;
  enum { BUFSIZE= 1024 }; char linebuf[BUFSIZE];

  fsession.open(fname, std::ios::in);

  G4String default_session = "";
  G4int iline = 1;
  sessionMap.clear();
  while( fsession.good() ) {
    if( fsession.eof()) break;
    fsession.getline(linebuf, BUFSIZE);
    G4String aline = G4StrUtil::strip_copy(linebuf);
    if ( aline[0] == '#' ) continue;
    if ( aline == "" ) continue;
    if ( iline == 1 )
      default_session = aline;
    else {
      size_t idx = aline.find_first_of(" ");
      if ( idx == G4String::npos ) break;
      G4String aname = aline.substr(0, idx);
      idx = aline.find_first_not_of(" ", idx);
      if (idx == G4String::npos ) break;
      G4String sname = aline.substr(idx, aline.size()-idx);
      sessionMap[aname] = sname;
    }
    iline++;
  }
  fsession.close();

  G4String stype = "";
  std::map<G4String, G4String>::iterator it = sessionMap.find(appname);
  if ( it != sessionMap.end() ) stype = sessionMap[appname];
  else stype = default_session;
  G4StrUtil::to_lower(stype);

  // select session...
  if ( qt_build && stype == "qt" ) selected = kQt;
  else if ( xm_build && stype == "xm" ) selected = kXm;
  else if ( win32_build && stype == "win32" ) selected = kWin32;
  else if ( tcsh_build && stype == "tcsh" ) selected = kTcsh;
  else if ( stype == "csh" ) selected = kCsh;
}

// --------------------------------------------------------------------------
void G4UIExecutive::SelectSessionByBestGuess()
{
  if ( qt_build ) selected = kQt;
  else if ( win32_build ) selected = kWin32;
  else if ( tcsh_build ) selected = kTcsh;
  else if ( xm_build ) selected = kXm;
}

// --------------------------------------------------------------------------
void G4UIExecutive::SetPrompt(const G4String& prompt)
{
  if(shell) shell-> SetPrompt(prompt);
}

// --------------------------------------------------------------------------
void G4UIExecutive::SetLsColor(TermColorIndex dirColor,
                               TermColorIndex cmdColor)
{
  if(shell) shell-> SetLsColor(dirColor, cmdColor);
}

// --------------------------------------------------------------------------
void G4UIExecutive::SessionStart()
{
  session-> SessionStart();
}
