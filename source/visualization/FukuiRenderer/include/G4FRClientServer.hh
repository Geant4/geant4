// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FRClientServer.hh,v 1.1 1999-01-07 16:14:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Satoshi TANAKA, Wed Jul  3 14:13:52 JST 1996
////////////////////////////////
///// G4FRClientServer.h  /////
////////////////////////////////

//=================//
#if defined (G4VIS_BUILD_DAWN_DRIVER) || defined (G4VIS_USE_DAWN)
//=================//


#if !defined G4FR_CLIENT_SERVER_H
#define G4FR_CLIENT_SERVER_H

#include<sys/types.h>
#include<sys/socket.h>
#include<netinet/in.h>
#include<arpa/inet.h>
#include<netdb.h>
#include<sys/un.h>
#include<unistd.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"G4ios.hh"


//----- constants
const  char	FR_ENV_SERVER_HOST_NAME[]	= "G4DAWN_HOST_NAME"  ;
const  char	FR_ENV_NAMED_PIPE_CONNECTION[]	= "G4DAWN_NAMED_PIPE" ;

	//-----------------------------------//
	//-----  class G4FRClientServer -----//
	//-----------------------------------//

class G4FRClientServer { 
 public:
	enum { SEND_BUFMAX = 1024 , RECV_BUFMAX = 1024 };
	enum { SUN_PATH_MAX = 16 };

 protected:	
	const char		TERMINATOR ;
	const char		END_OF_LINE;
	char			SUN_PATH[ SUN_PATH_MAX ];
	int			PORT_NUMBER ;
	int			fSocketFd ;

	char		fReceivedMessage [ RECV_BUFMAX ];
	char		fSendingMessage  [ SEND_BUFMAX ];

 protected:

	void		Err( char* message ) { perror(message) ;}
	void		SetSendingMessage( const char* message ) 
				{ strcpy( fSendingMessage, message );}
	void		Send() ; // send command in fSendingMessage

 public:	
		//----- Server
	int		AcceptUnix(){ return 0;}  // made unfunctioned  
	int		AcceptINET(){ return 0 ;}  // made unfunctioned  

		//----- Client
	int		ConnectUnix(); 
	int		ConnectINET(); 

		//----- Common to server and client

		//---------- (1)
	G4FRClientServer (	char terminator = '.'           ,
				char end_line = '\n'              ) ;  
	virtual ~G4FRClientServer () {;}
	void		SetSunPath( const char* sun_path ) 
			{ strcpy     ( SUN_PATH, sun_path ); }
	void		SetPortNumber( int port_num ) 
			{ PORT_NUMBER = port_num ; }
	void		IncrementPortNumber( int incr = 1 ) 
			{ PORT_NUMBER += incr ; }
	const char*	GetSendingMessage() const
				{ return  fSendingMessage            ;}
	int		GetSendingMessageLength() const
				{ return  strlen(fSendingMessage)    ;}
	void		SetReceivedMessage( const char* message ) 
				{ strcpy( fReceivedMessage, message );}
	const char*	GetReceivedMessage() const
				{ return fReceivedMessage            ;}
	int		GetReceivedMessageLength() const    
				{ return strlen(fReceivedMessage)    ;}
	int		GetSofd() const { return fSocketFd ; }
	int		GetPortNumber () const { return PORT_NUMBER ; }
	void		ClearReceivedMessage () 
			{ memset(fReceivedMessage, '\0', RECV_BUFMAX) ;	}

	int		IsTerminator(char ch ) { return ( ch == TERMINATOR ); }
	char		GetTerminator()  const  { return TERMINATOR ; }
	int		IsEndOfLine (char ch ) { return ( ch == END_OF_LINE ); }
	char		GetEndOfLine() const { return END_OF_LINE ; }
	void		DisConnect();
	void		Clear();
	void		WaitSendBack( const char* command_string ) ;

		//---------- (2) send and receive
	virtual	void	Send     ( const char* message ) ; 
	virtual void	SendLine ( const char* message ) ; 
				// Add END_OF_LINE if the message does not
				// Terminate with it. And then send the message.
	virtual	void	Receive() ;
	virtual void	ReceiveLine() ;

}; // G4FRClientServer


#endif
#endif //G4VIS_BUILD_DAWN_DRIVER

