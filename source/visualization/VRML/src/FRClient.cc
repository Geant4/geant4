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
// $Id: FRClient.cc 66870 2013-01-14 23:38:59Z adotti $
//
// FRClient.cc
// FukuiRenderer Client
// Yasuhide Sawada & Satoshi Tanaka


//=================//
#ifndef WIN32
//=================//


//=================//
#ifdef G4VIS_BUILD_VRML_DRIVER
//=================//


#include "G4ios.hh"
#include <sys/time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <fcntl.h>

#include <unistd.h>
#include <string.h>
#include <stdio.h>

#include "G4VisManager.hh"
#include "FRClient.h"

FRClient::FRClient()
{
	fd = -1;
	create();
}

FRClient::~FRClient()
{
	close();
}

int FRClient::create()
{
	/* stream socket */
	fd = socket(AF_INET, SOCK_STREAM, 0);
	if (fd < 0)
		fputs("error: socket.\n", stderr);
	return fd;
}



int FRClient::connect(const char *hostname, int port_)
{
	// local variables
	struct sockaddr_in sa;
	struct hostent *hp;

	// set port ( sa.sin_family,  sa.sin_port )
	port = port_;  // Store port number to data member
	memset( (char *)&sa, '\0', sizeof(sa)) ;
	sa.sin_family = AF_INET;
	sa.sin_port = htons(port);

	// set server host ( sa.sin_addr )
	if (hostname == NULL) {
		hostname = "localhost"; 
			// reset arg 
	}
	hp = gethostbyname(hostname) ;
	if ( !hp ) {
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
		G4cout << "ERROR: gethostbyname() failed" << G4endl;
	  return -1; 
	}

	memcpy( (char * )&sa.sin_addr, (char * )hp->h_addr, hp->h_length );

	// make connection to server
	if (::connect(fd, (struct sockaddr *)&sa, sizeof(sa)) == -1) {
		fputs("error: connect\n", stderr);
		return -1;
	}

	// return file descripter (data member)
	return fd;
}



int FRClient::send(const char *sendbuf)
{
	int len = strlen(sendbuf);

	if (::send(fd, sendbuf, len, 0) < 0) {
		fputs("error: Send()\n", stderr);
		len = -1;
	}
	return len;
}

int FRClient::receive(char *recvbuf)
{
	int len;

	memset(recvbuf, '\0', FRSendLength + 1);
	len = ::recv(fd, recvbuf, FRSendLength, 0);
	if(len < 0) {
		fputs("error: Receive()\n", stderr);
		len = -1;
	}
	return len;
}

int FRClient::close()
{
	/*
		shutdown :argument '2' means shutdown both send and receive.
	*/
	if (::shutdown(fd, 2) < 0) {
		fputs("error: shutdown\n", stderr);
		return -1;
	}
	::close(fd);
	return 0;
}

#endif //G4VIS_BUILD_VRML_DRIVER

#endif // WIN32
