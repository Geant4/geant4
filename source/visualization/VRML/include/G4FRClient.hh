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
// $Id: G4FRClient.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// G4FRClient.hh
// Yasuhide Sawada and Satoshi Tanaka

#ifndef WIN32

#ifdef  G4VIS_BUILD_VRML_DRIVER

#ifndef G4_FR_CLIENT_HH
#define G4_FR_CLIENT_HH

#include "globals.hh"

class FRClient;

class G4FRClient {
	// G4FRClient can SEND only! 
public:
	G4FRClient();
	virtual ~G4FRClient();

	G4bool connect(const char *hostname, G4int port);
	void close();

	G4int getPort() const;

	G4bool isConnected() const { return connected; }
	G4bool is_open    () const { return connected; }

	G4FRClient& operator << (G4int);
	G4FRClient& operator << (G4double);
	G4FRClient& operator << (const char *);
	G4FRClient& operator << (G4FRClient& (*)(G4FRClient&));
private:
	FRClient *fFRClient;
	G4bool connected;
	G4int fPort;
};

//manipulator
//G4FRClient& endl(G4FRClient&);

#endif //G4_FR_CLIENT_HH
#endif //G4VIS_BUILD_VRML_DRIVER
#endif //WIN32
