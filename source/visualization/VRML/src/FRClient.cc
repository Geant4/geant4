// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FRClient.cc,v 1.1 1999-01-07 16:15:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// FRClient.cc
// FukuiRenderer Client
// Yasuhide Sawada & Satoshi Tanaka

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
		G4cerr << "ERROR: gethostbyname() failed" << endl;
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

#endif
