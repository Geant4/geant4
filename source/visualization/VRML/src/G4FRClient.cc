// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FRClient.cc,v 1.3 1999-12-15 14:54:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4FRClient.cc
// Yasuhide Sawada & Satoshi Tanaka

//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


#include <stdio.h>
#include "G4FRClient.hh"
#include "FRClient.h"

G4FRClient::G4FRClient()
{
	fFRClient = NULL;
	fPort = -1;
	connected = false;
}

G4FRClient::~G4FRClient()
{
	if (connected)
		this->close();
}

G4bool G4FRClient::connect(const char *hostname, G4int port)
{
	if (connected)
		return false;

	delete fFRClient;
	fFRClient = new FRClient();

	fPort = port;
	connected = (fFRClient->connect(hostname, port) < 0) ? false : true ;

	return connected;
}

void G4FRClient::close()
{
	delete fFRClient;
	fFRClient = NULL;
	connected = false;
}

G4int G4FRClient::getPort() const
{
	return fPort;
}

G4FRClient& G4FRClient::operator << (G4int val)
{
	char buf[64];
	sprintf(buf, "%d", val);
	fFRClient->send(buf);
	return *this;
}

G4FRClient& G4FRClient::operator << (G4double val)
{
	char buf[64];
	sprintf(buf, "%g", val);
	fFRClient->send(buf);
	return *this;
}

G4FRClient& G4FRClient::operator << (const char *pval)
{
	fFRClient->send(pval);
	return *this;
}

G4FRClient& G4FRClient::operator << (G4FRClient& (*func)(G4FRClient&))
{
	return func(*this);
}

//manipulator
G4FRClient& endl(G4FRClient& c)
{
	return c << "\n";
}

#endif
