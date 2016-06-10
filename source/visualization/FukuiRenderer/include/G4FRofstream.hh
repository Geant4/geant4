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
// $Id: G4FRofstream.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
#include <fstream>

#if !defined G4_FR_OFSTREAM_HH
#define      G4_FR_OFSTREAM_HH


#include "globals.hh"

/////////////////////
//typedef int G4bool ;
//#define false 0 ;
//#define true  1 ;
////////////////////

class G4FRofstream {

 public:
	enum { SEND_BUFMAX = 1024 };

 public: 
	
	// constructors
	G4FRofstream ()                       { flag_file_open = false ; } 
	G4FRofstream ( const char* filename ) ;

	// destructor
	virtual ~G4FRofstream ();

	// open and close 
	void Open ( const char* filename );
	void Close() ;
	G4bool IsOpen() { return flag_file_open ;}

	// utilities
	void SendLine( const char* string ) ; // save string with new line

	// static functions
	static G4bool DoesFileExist( const char* filename ) ;

 protected:
	G4bool flag_file_open ;
	std::ofstream fout ;
} ;


inline  void G4FRofstream::Open ( const char* filename )
{ 
	if( !IsOpen() ) {
		fout.open( filename ) ; 
		flag_file_open = true ;
	}
}


inline  void G4FRofstream::Close ()
{ 
	if( IsOpen() ) {
		fout.close();
		flag_file_open = false ;
	}
}

inline  void    G4FRofstream::SendLine ( const char* message ) 
{
	if ( IsOpen() ) {
		fout << message << G4endl;
	}		
}


inline  G4bool  G4FRofstream::DoesFileExist ( const char* filename ) 
{
	G4bool status = false ;

	std::ifstream fout_tmp( filename ) ; 
	if( fout_tmp ) { status = true ; }
	fout_tmp.close();

	return status ;
}


inline 
G4FRofstream::G4FRofstream ( const char* filename ) 
{
	flag_file_open = false ; 
	Open( filename ); 
} 

inline 
G4FRofstream::~G4FRofstream () 
{
	Close() ;
} 


#endif
