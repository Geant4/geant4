// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FRClientServer.cc,v 1.1 1999-01-07 16:14:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Satoshi TANAKA, Wed Jul  3 14:14:29 JST 1996
////////////////////////////////
///// G4FRClientServer.cc /////
////////////////////////////////


//=================//
#ifdef G4VIS_BUILD_DAWN_DRIVER
//=================//

#include "G4FRSocketConfig.hh"
#include "G4FRSocketMacro.hh"
#include "G4FRClientServer.hh"

// #define DEBUG_CLIENT_SERVER

#include<sys/param.h>

	//----- const 	
const  char	DEFAULT_SUN_PATH[]		= "FR_TMP3"        ;
const  int	DEFAULT_PORT_NUMBER		= 40701            ;
//const  char	FR_ENV_SERVER_HOST_NAME[]	= "G4DAWN_HOST_NAME" ; // moved to .hh
const  int	MAX_BINDING_TRIAL		= 10               ;
const  int	BINDING_TRIAL_INTERVAL		= 5                ;
const  int	MAX_CONNECT_TRIAL		= 10               ;


	//----- G4FRClientServer::G4FRClientServer ()
G4FRClientServer::G4FRClientServer ( char terminator , char end_line ) :  
	TERMINATOR ( terminator ) , 
	END_OF_LINE( end_line   ) ,
	fSocketFd  ( -1 )           
{ 
	SetSunPath    ( DEFAULT_SUN_PATH ) ; // for Unix domain
	SetPortNumber ( DEFAULT_PORT_NUMBER ) ;
	ClearReceivedMessage () ;
}


	//----- G4FRClientServer::ConnectUnix()
int G4FRClientServer::ConnectUnix()
{
		//----- local
	int			flag_connected = 0 ; 
	struct sockaddr_un	server_address ;

		//----- make socket
	fSocketFd = socket( AF_UNIX, SOCK_STREAM, 0 );
	if( fSocketFd < 0 ) { Err("G4FRClientServer::ConnectUnix(),socket"); }

		//----- set server address
	memset( (char *)&server_address, '\0', sizeof(server_address)) ;
	server_address.sun_family = AF_UNIX ;
	strcpy( server_address.sun_path, SUN_PATH );

		//----- connection
	int	connection_status = -1  ;
	int	num_trial         = 0 ;
	while( connection_status < 0 && num_trial <= MAX_CONNECT_TRIAL ) {
		num_trial++ ; 
		connection_status = connect( fSocketFd, (struct sockaddr * )(&server_address), sizeof( server_address ) ) ;
		if( connection_status  <0 ) 
		{
#if defined DEBUG_CLIENT_SERVER
			Err("G4FRClientServer::ConnectUnix(),connect => RETRY");
#endif
			flag_connected = 0 ;
		} else {
			flag_connected = 1 ;
			break ;
		}

		sleep(1);
		
	} // while(connection_status...)

		//----- return status of connection
	return flag_connected ;

} // G4FRClientServer::ConnectUnix()


	//-----	G4FRClientServer::Receive()
void	G4FRClientServer::Receive()
{
		//-----
	ClearReceivedMessage () ;
	if( recv( fSocketFd, fReceivedMessage, G4FRClientServer::RECV_BUFMAX , 0 ) < 0 ) 
	{
		Err("G4FRClientServer::Receive(), recv");
	}

#if defined DEBUG_CLIENT_SERVER
	G4cerr << ">>>>> receivedMessage = " << fReceivedMessage << endl;
#endif

}


	//-----	G4FRClientServer::ReceiveLine()
void	G4FRClientServer::ReceiveLine()
{
		//----- local
	char    buf[1];
	int index = 0 ;

		//----- receive a line (until newline)
	memset(fReceivedMessage, '\0', RECV_BUFMAX) ;	
	while( read( fSocketFd, buf, 1 ) == 1 ) {
		fReceivedMessage[index++] = buf[0];
		if( IsEndOfLine(buf[0]) ) { break ;}
	}
} // G4FRClientServer::ReceiveLine()


	//-----	G4FRClientServer::Send()
void	G4FRClientServer::Send()
{
	if( send( fSocketFd, fSendingMessage, strlen(fSendingMessage) , 0 ) < 0 ) 
	{
		Err("G4FRClientServer::Send(), send");
	}

#if defined DEBUG_CLIENT_SERVER
	G4cerr << "<<<<< SentMessage = " << fSendingMessage << endl;
#endif

} // G4FRClientServer::Send()


	//-----	G4FRClientServer::Send( message )
void	G4FRClientServer::Send( const char* message ) 
{
	this->SetSendingMessage( message )      ;
	this->Send();

} // G4FRClientServer::Send( message )


	//-----	G4FRClientServer::SendLine()
void	G4FRClientServer::SendLine( const char* message ) 
{
		//----- local
	int	smsg_length ;

		//----- set message to sending buf
	this->SetSendingMessage( message )      ;
	smsg_length = GetSendingMessageLength() ;

		//----- add newline if necessary
	if( !IsEndOfLine( fSendingMessage[ (smsg_length - 1)] ) ) {
		fSendingMessage[ smsg_length ]      = GetEndOfLine() ;
		fSendingMessage[ (smsg_length +1) ] = '\0' ;
		smsg_length = GetSendingMessageLength();
	}

		//----- send
	this->Send();

}// G4FRClientServer::SendLine()


	//----- G4FRClientServer::DisConnect()
void G4FRClientServer::DisConnect()
{
		//----- close connection
	if( shutdown(fSocketFd,2) < 0 ) { 
		Err("G4FRClientServer::DisConnect,shutdown");
	}
	close( fSocketFd );

	this->Clear();
}



	//----- G4FRClientServer::Clear()
void G4FRClientServer::Clear()
{
	unlink(SUN_PATH) ;
	fSocketFd = -1   ;
}


	//----- G4FRClientServer::ConnectINET()
int G4FRClientServer::ConnectINET()
{
		//----- local
	int			flag_connected = 0 ;
	sockaddr_in		server_address ;
	char			server_hostname[32] ;
	hostent*		server_host_p ;

		//----- make socket
	fSocketFd = socket( AF_INET, SOCK_STREAM, 0 );
	if( fSocketFd < 0 ) { 
#if defined DEBUG_CLIENT_SERVER
	  Err("G4FRClientServer::ConnectINET(),socket"); 
#endif
	}

		//----- get IP address of server from its name
	if( getenv( FR_ENV_SERVER_HOST_NAME ) != NULL ) 
	{
			//----- get server name
		strcpy( server_hostname, getenv( FR_ENV_SERVER_HOST_NAME ) );

			//----- get IP address of server from its name,
			//..... reading /etc/hosts
		server_host_p = gethostbyname( server_hostname ) ;
			
			//----- If the host specified by FR_ENV_SERVER_HOST_NAME
			//.....  is not written in /etc/hosts, 
			//...... server host is set to the same as the
			//...... client host
		if( !server_host_p ) {
#if defined DEBUG_CLIENT_SERVER
			Err("G4FRClientServer::ConnectINET(), gethostbyname");
#endif
				//----- get server name
			if( GET_HOSTNAME( server_hostname, sizeof( server_hostname ) ) < 0 )
			{ Err("G4FRClientServer::ConnectNET(), GET_HOSTNAME");}

				//----- get IP address of server from its name
			if( !(server_host_p = gethostbyname( server_hostname ) ) )
			{ Err("G4FRClientServer::ConnectINET(), gethostbyname");}
		}

	} else {
			//----- If FR_ENV_SERVER_HOST_NAME is not set, 
			//...... server host is set to the same as the
			//...... client host.

			//----- get server name
		if( GET_HOSTNAME( server_hostname, sizeof( server_hostname ) ) < 0 )
		{ Err("G4FRClientServer::ConnectNET(), GET_HOSTNAME");}

			//----- get IP address of server from its name
		if( !(server_host_p = gethostbyname( server_hostname ) ) )
		{ Err("G4FRClientServer::ConnectINET(), gethostbyname");}
	}



// #if defined DEBUG_CLIENT_SERVER
	G4cerr << "***** Trying connection to  " << server_hostname << endl;
// #endif 
	

		//----- connection and binding 
	memset( (char *)&server_address, '\0', sizeof(server_address)) ; 
				// clear server_address
	server_address.sin_family = AF_INET ;
	server_address.sin_port   = htons( PORT_NUMBER );
	memcpy( (char *)(&server_address.sin_addr ), 
		(char *)( server_host_p->h_addr   ), 
		server_host_p->h_length             ); 

	int	connection_status = -1  ;
	int	num_trial         = 0 ;
	while( connection_status < 0 && num_trial <= MAX_CONNECT_TRIAL ) {
		num_trial++ ; 
		connection_status = connect( fSocketFd, (struct sockaddr * )(&server_address), sizeof( server_address ) ) ;
		if( connection_status  <0 ) 
		{
#if defined DEBUG_CLIENT_SERVER
			Err("G4FRClientServer::ConnectINET(),connect => RETRY");
#endif 
			flag_connected  = 0 ;
		} else {
			flag_connected  = 1 ;
			break ;
		}

		sleep(1);

	} // while(connection_status...)

		//----- return status of connection
	return flag_connected ;

} // G4FRClientServer::ConnectINET()


	//----- G4FRClientServer::WaitSendBack()
void	G4FRClientServer::WaitSendBack( const char* command_string ) 
{
		//----- wait for sending back
	while(1) { 
		this->ReceiveLine();

		if( !strncmp(	this->GetReceivedMessage(), \
				command_string                , \
				(strlen(command_string))       )   )
		{
			break;
		} else {
			sleep(2);
		}

	} // while

		//----- clear buffer to receive message
	this->ClearReceivedMessage();	

} // G4FRClientServer::WaitSendBack()

#endif // G4VIS_BUILD_DAWN_DRIVER
