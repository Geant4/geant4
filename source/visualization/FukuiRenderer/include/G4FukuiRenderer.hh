// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FukuiRenderer.hh,v 1.5 1999-12-15 14:54:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Satoshi TANAKA
// Fukui Renderer factory.

//=================//
#if defined (G4VIS_BUILD_DAWN_DRIVER) || defined (G4VIS_USE_DAWN)
//=================//

#ifndef G4FUKUI_RENDERER_HH
#define G4FUKUI_RENDERER_HH

#include "G4VGraphicsSystem.hh"
#include "G4FRClientServer.hh"

	//----- prototype
class G4VSceneHandler   ;

	//---------------------------------//
	//----- class G4FukuiRenderer -----// 
	//---------------------------------//
class G4FukuiRenderer: public G4VGraphicsSystem {

public:
  G4FukuiRenderer ();
  virtual ~G4FukuiRenderer ();
  G4VSceneHandler* CreateSceneHandler (const G4String& name = "");
  G4VViewer*  CreateViewer  (G4VSceneHandler&, const G4String& name = "");

  G4FRClientServer& GetPrimDest() { return fPrimDest ; }

	//----- inter-process communication
  enum { IP_UNIX, IP_INET };
  enum { FR_MAX_PORT_INCR = 10 };
  G4bool IsUnixDomain() const { return (fIPMode == G4FukuiRenderer::IP_UNIX); }
  G4bool IsInetDomain() const { return (fIPMode == G4FukuiRenderer::IP_INET); }
  void                  UseInetDomainAuto();
  void                  UseInetDomain();
  void                  UseBSDUnixDomainAuto();
  void                  ConnectPort( int max_port_incr = FR_MAX_PORT_INCR);

  G4bool	        IsGUIMode   (void) { return flag_use_gui   ; }
  G4bool	        IsConnected (void) { return flag_connected ; }

public:
	//----- inter-process communication
  G4FRClientServer  fPrimDest ;

private:
  int		fIPMode   ;
  G4bool	flag_use_gui ;
  G4int		flag_connected ;
};

#endif
#endif //G4VIS_BUILD_DAWN_DRIVER
