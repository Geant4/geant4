// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FukuiRenderer.cc,v 1.1 1999-01-07 16:14:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Satoshi TANAKA
// FukuiRenderer factory.

//=================//
#ifdef G4VIS_BUILD_DAWN_DRIVER
//=================//

//#define DEBUG_FR_SYSTEM


#include "G4FukuiRenderer.hh"

#define __G_ANSI_C__

//#include "G4VisFeaturesOfFukuiRenderer.hh"
#include "G4VScene.hh"
#include "G4FukuiRendererScene.hh"
#include "G4FukuiRendererView.hh"
#include "G4FRClientServer.hh"
#include "G4FRConst.hh"
#include "G4FRFeatures.hh"

	//----- constant
const int FR_DEFAULT_CONNECTION_TIME = 5  ;
const char FR_ENV_CONNECTION_TIME [] = "G4DAWN_CONNECTION_TIME";
const char FR_ENV_DAWN_GUI_ALWAYS [] = "G4DAWN_GUI_ALWAYS";

	//----- G4FukuiRenderer, constructor
G4FukuiRenderer::G4FukuiRenderer ():
  G4VGraphicsSystem ("FukuiRenderer",
		     "DAWN",
		     FR_DAWNFILE_FEATURES,
		     G4VGraphicsSystem::threeD),
  fPrimDest() ,
  fIPMode( G4FukuiRenderer::IP_UNIX ),
  flag_use_gui     (false) ,
  flag_connected   (0) 
{
		//----- sorting of objects
	if( ( getenv( FR_ENV_DAWN_GUI_ALWAYS ) != NULL        ) && \
	    ( strcmp( getenv( FR_ENV_DAWN_GUI_ALWAYS ),"0"  ) )      ) 
	{
		flag_use_gui =	true ;
	} 
}


	//----- G4FukuiRenderer, destructor
G4FukuiRenderer::~G4FukuiRenderer () 
{
	if( flag_connected ) {
			//----- Terminate FukuiRenderer
		fPrimDest.SendLine( FR_TERMINATE_DAWND ) ; 

			//----- disconnection
		fPrimDest.DisConnect(); 
		flag_connected = 0 ;
	}
}

	//-----  G4FukuiRenderer::CreateScene () 
G4VScene* G4FukuiRenderer::CreateScene (const G4String& name) 
{
	G4VScene* p = new G4FukuiRendererScene (*this, name);

	if(!flag_connected) { delete p ;  p = NULL ; }

	G4cout	<< G4FukuiRendererScene::GetSceneCount ()
		<< ' ' << fName << " scenes extanct." << endl;

	return p;
}

	//-----  G4FukuiRenderer::CreateView () 
G4VView* G4FukuiRenderer::CreateView (G4VScene& scene, const G4String& name) 
{
	if(!flag_connected) return NULL;
       	G4VView* pView = 
	  new G4FukuiRendererView ((G4FukuiRendererScene&) scene, name);
	return pView;
}

	//----- G4FukuiRenderer::UseInetDomainAuto()
void G4FukuiRenderer::UseInetDomainAuto()
{
	int		pid ;

#if defined DEBUG_FR_SYSTEM
	G4cerr << "***** Unix Inet Domain (AUTO mode)" << endl;
#endif
	fIPMode = G4FukuiRenderer::IP_UNIX ;

	if( ( pid = fork() ) == 0 ) { // child
		if ( execlp ( "dawn", "dawn", "-G" , (char *)0 ) < 0 ) 
		{
			perror("dawn") ;
		}
	} else { // parent 

			//----- Set waiting time to ensure connection
		int connection_time = FR_DEFAULT_CONNECTION_TIME  ;
		if( getenv( FR_ENV_CONNECTION_TIME ) != NULL  ) {
			int sscanf_status = \
				sscanf( getenv( FR_ENV_CONNECTION_TIME ), "%d", &connection_time) ;
			if ( sscanf_status == EOF ) {  
			  connection_time = FR_DEFAULT_CONNECTION_TIME  ;
			}
		}


			//----- Wait for starting up of FukuiRenderer
		sleep(connection_time); 

			//----- establish connection
		this->ConnectPort();
	}

	if(!flag_connected) { 
	  G4cerr << "***** ERROR: Connection failed" << endl; 
	}
	else { 
	  G4cerr << "***** GEANT4 is connected to FukuiRenderer DAWN ";
	  G4cerr << "in the same host" << endl; 
	}

} //  G4FukuiRenderer::UseInetDomainAuto()


	//----- G4FukuiRenderer::UseInetDomain()
void
G4FukuiRenderer::UseInetDomain()
{
#if defined DEBUG_FR_SYSTEM
	G4cerr << "***** INET Domain " << endl;
#endif
	fIPMode = G4FukuiRenderer::IP_INET ;

	this->ConnectPort();

	if(!flag_connected) {
	  G4cerr << "***** ERROR: Connection failed" << endl; 
	}
	else { 
	  G4cerr << "GEANT4 is connected to FukuiRenderer DAWN via socket" ; 
	  G4cerr << endl; 
	}

} // G4FukuiRenderer::UseInetDomain()

void
G4FukuiRenderer::ConnectPort( int max_port_incr )
{
	//----- establish connection
  int connect_trial = 0 ;
  while(1) {
    if ( ++connect_trial > max_port_incr ) {
	this->flag_connected = 0 ;
	G4cerr << "***** INET Connection failed."                << endl;
	G4cerr << "      Maybe, you have not invoked DAWN yet,"  << endl;
	G4cerr << "      or too many DAWN's are already running" << endl;
	G4cerr << "      in the server host."                    << endl;
	fPrimDest.IncrementPortNumber( (- FR_MAX_PORT_INCR) );
	return ;
    } else if ( fPrimDest.ConnectINET() ) { 
	    // INET domain connection is established
	this->flag_connected = 1 ;
	G4cerr << "***** GEANT4 is connected to port  " ;
	G4cerr << fPrimDest.GetPortNumber() ; 
	G4cerr << "  of server" << endl;
	break ; 
    } else { 
	    // Connection failed. Try the next port.
      this->flag_connected = 0 ;
      fPrimDest.IncrementPortNumber();
      G4cerr << "***** GEANT4 incremented targeting port to " ;
      G4cerr << fPrimDest.GetPortNumber() << endl;

    } // if-elseif-else

  } // while (1) 
}


	//----- G4FukuiRenderer::UseBSDUnixDomainAuto()
void G4FukuiRenderer::UseBSDUnixDomainAuto()
{
	int     pid ;

#if defined DEBUG_FR_SYSTEM
	G4cerr << "***** UseBSDUnixDomainAuto " << endl;
#endif
	fIPMode = G4FukuiRenderer::IP_UNIX ; // Unix domain

	if( ( pid = fork() ) == 0 ) { // child
		if ( execlp ( "dawn", "dawn", "-g" , (char *)0 ) < 0 ) 
		{
			perror("dawn") ;
		}
	} else { // parent 

	  		//----- Sleep for a while to make sure that 
			//..... FukuiRenderer is ready
		int connection_time = FR_DEFAULT_CONNECTION_TIME  ;
		if( getenv( FR_ENV_CONNECTION_TIME ) != NULL  ) {
			int sscanf_status = \
				sscanf( getenv( FR_ENV_CONNECTION_TIME ), "%d", &connection_time) ;
			if ( sscanf_status == EOF ) {  connection_time = FR_DEFAULT_CONNECTION_TIME  ;}
		}
		sleep(connection_time); 

			//----- connect GEANT4 to FukuiRenderer
		this->flag_connected = fPrimDest.ConnectUnix();

			//----- display status
		if(!flag_connected) {
		  G4cerr << "***** ERROR: Connection failed" << endl; 
		} else { 
		  G4cerr << "*** GEANT4 is connected to FukuiRenderer DAWN ";
		  G4cerr <<  "in the same host" << endl; 
		}

	} // if--else

}// G4FukuiRenderer::UseBSDUnixDomainAuto()


#endif // G4VIS_BUILD_DAWN_DRIVER
