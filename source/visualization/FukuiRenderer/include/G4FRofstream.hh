// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FRofstream.hh,v 1.2 1999-01-09 16:11:43 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include<fstream.h>

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
	ofstream fout ;
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
		fout << message << endl;
	}		
}


inline  G4bool  G4FRofstream::DoesFileExist ( const char* filename ) 
{
	G4bool status = false ;

	ifstream fout_tmp( filename ) ; 
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
