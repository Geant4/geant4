// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FRClient.hh,v 1.2 1999-01-09 16:27:29 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4FRClient.hh
// Yasuhide Sawada and Satoshi Tanaka

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
G4FRClient& endl(G4FRClient&);

#endif //G4_FR_CLIENT_HH
#endif //G4VIS_BUILD_VRML_DRIVER
