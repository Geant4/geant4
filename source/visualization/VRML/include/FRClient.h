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
// $Id: FRClient.h,v 1.3 2001-07-11 10:09:09 gunter Exp $
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
