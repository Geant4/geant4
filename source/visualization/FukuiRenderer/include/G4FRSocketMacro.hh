// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FRSocketMacro.hh,v 1.3 1999-12-15 14:54:00 gunter Exp $
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
