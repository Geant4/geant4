//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4FRClient.cc,v 1.7 2002-06-23 03:31:50 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4FRClient.cc
// Yasuhide Sawada & Satoshi Tanaka


//=================//
#ifndef WIN32
//=================//


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

////////////////////////////////////////
////manipulator
//G4FRClient& endl(G4FRClient& c)
//{
//	return c << "\n";
//}
///////////////////////////////////////

#endif //G4VIS_BUILD_VRML_DRIVER

#endif //WIN32

