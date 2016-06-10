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
// $Id: G4FRClient.cc 66373 2012-12-18 09:41:34Z gcosmo $
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

