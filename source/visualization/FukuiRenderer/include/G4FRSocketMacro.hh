// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FRSocketMacro.hh,v 1.1 1999-01-07 16:14:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
///////////////////////////////
///// G4FRSocketMacro.hh  /////
///////////////////////////////
//----- MACRO for portability -----//

//=================//
#ifdef G4VIS_BUILD_DAWN_DRIVER
//=================//

#if !defined FR_SOCKET_MACRO_H
#define FR_SOCKET_MACRO_H


	//----- gethostname
#if defined FR_SOCKET_IRIX_SOLARIS
 #include <sys/systeminfo.h>
 #define  GET_HOSTNAME( hostname, length )  sysinfo( SI_HOSTNAME, hostname, length ) 
#else 
 #define  GET_HOSTNAME( hostname, length )  gethostname( hostname, length ) 
#endif			

#endif
#endif //G4VIS_BUILD_DAWN_DRIVER
