// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FRClient.h,v 1.2 1999-12-15 14:54:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// FRClient.h
// FukuiRenderer Client
// Yasuhide Sawada & Satoshi Tanaka

#ifdef  G4VIS_BUILD_VRML_DRIVER

#ifndef FR_CLIENT_H
#define FR_CLIENT_H

#define FRSendLength 256

class FRClient {
public:
	FRClient();
	virtual ~FRClient();

	int create();
	int connect(const char *hostname, int port);

	int send(const char *sendbuf);
	int receive(char *recvbuf);

	int close();
protected:
	int fd;
	int port;
};

#endif //FR_CLIENT_H
#endif //G4VIS_BUILD_VRML_DRIVER
