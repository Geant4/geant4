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
// $Id: G4FukuiRenderer.hh 66373 2012-12-18 09:41:34Z gcosmo $
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
