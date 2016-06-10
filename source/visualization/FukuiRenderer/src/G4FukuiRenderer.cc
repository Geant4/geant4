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
// $Id: G4FukuiRenderer.cc 66870 2013-01-14 23:38:59Z adotti $
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

#include "G4VisManager.hh"
//#include "G4VisFeaturesOfFukuiRenderer.hh"
#include "G4VSceneHandler.hh"
#include "G4FukuiRendererSceneHandler.hh"
#include "G4FukuiRendererViewer.hh"
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

	//-----  G4FukuiRenderer::CreateSceneHandler () 
G4VSceneHandler* G4FukuiRenderer::CreateSceneHandler (const G4String& name) 
{
	G4VSceneHandler* p = new G4FukuiRendererSceneHandler (*this, name);

	if(!flag_connected) { delete p ;  p = NULL ; }

	return p;
}

	//-----  G4FukuiRenderer::CreateViewer () 
G4VViewer* G4FukuiRenderer::CreateViewer (G4VSceneHandler& scene, const G4String& name) 
{
	if(!flag_connected) return NULL;
       	G4VViewer* pView = 
	  new G4FukuiRendererViewer ((G4FukuiRendererSceneHandler&) scene, name);
	return pView;
}

	//----- G4FukuiRenderer::UseInetDomainAuto()
void G4FukuiRenderer::UseInetDomainAuto()
{
	int		pid ;

#if defined DEBUG_FR_SYSTEM
	if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	  G4cout << "***** Unix Inet Domain (AUTO mode)" << G4endl;
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
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	    G4cout << "***** ERROR: Connection failed" << G4endl; 
	}
	else { 
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
	    G4cout << "***** GEANT4 is connected to FukuiRenderer DAWN ";
	    G4cout << "in the same host" << G4endl; 
	  }
	}

} //  G4FukuiRenderer::UseInetDomainAuto()


	//----- G4FukuiRenderer::UseInetDomain()
void
G4FukuiRenderer::UseInetDomain()
{
#if defined DEBUG_FR_SYSTEM
  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	G4cout << "***** INET Domain " << G4endl;
#endif
	fIPMode = G4FukuiRenderer::IP_INET ;

	this->ConnectPort();

	if(!flag_connected) {
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	    G4cout << "***** ERROR: Connection failed" << G4endl; 
	}
	else { 
	  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
	    G4cout << "GEANT4 is connected to FukuiRenderer DAWN via socket" ; 
	    G4cout << G4endl; 
	  }
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
	if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
	  G4cout << "***** INET Connection failed."                << G4endl;
	  G4cout << "      Maybe, you have not invoked DAWN yet,"  << G4endl;
	  G4cout << "      or too many DAWN's are already running" << G4endl;
	  G4cout << "      in the server host."                    << G4endl;
	}
	fPrimDest.IncrementPortNumber( (- FR_MAX_PORT_INCR) );
	return ;
    } else if ( fPrimDest.ConnectINET() ) { 
	    // INET domain connection is established
	this->flag_connected = 1 ;
	if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
	  G4cout << "***** GEANT4 is connected to port  " ;
	  G4cout << fPrimDest.GetPortNumber() ; 
	  G4cout << "  of server" << G4endl;
	}
	break ; 
    } else { 
	    // Connection failed. Try the next port.
      this->flag_connected = 0 ;
      fPrimDest.IncrementPortNumber();
      if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
	G4cout << "***** GEANT4 incremented targeting port to " ;
	G4cout << fPrimDest.GetPortNumber() << G4endl;
      }

    } // if-elseif-else

  } // while (1) 
}


	//----- G4FukuiRenderer::UseBSDUnixDomainAuto()
void G4FukuiRenderer::UseBSDUnixDomainAuto()
{
	int     pid ;

#if defined DEBUG_FR_SYSTEM
	if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
	  G4cout << "***** UseBSDUnixDomainAuto " << G4endl;
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
		  if (G4VisManager::GetVerbosity() >= G4VisManager::errors)
		    G4cout << "***** ERROR: Connection failed" << G4endl; 
		} else { 
		  if (G4VisManager::GetVerbosity() >= G4VisManager::errors) {
		    G4cout << "*** GEANT4 is connected to FukuiRenderer DAWN ";
		    G4cout <<  "in the same host" << G4endl; 
		  }
		}

	} // if--else

}// G4FukuiRenderer::UseBSDUnixDomainAuto()


#endif // G4VIS_BUILD_DAWN_DRIVER
