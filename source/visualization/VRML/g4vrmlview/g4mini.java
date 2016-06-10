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
// $Id: g4mini.java 66373 2012-12-18 09:41:34Z gcosmo $
//
import java.io.*;
import java.net.*;

public class g4mini
{
	public static void main( String[] args )
	{
	try{

		// port number 
		final int DEFALUT_PORT_NO = 40801 ;
		int   portNo =  DEFALUT_PORT_NO ;

		// argument checking
		if( args.length != 2 ) 
		{
			System.out.println( "Usage: java g4mini  src_file  server_hostname");
			return ;
		}
		String src         = args[0];
		String server      = args[1];

		// open connection
		Socket   socket    =  new Socket( server, portNo ) ;

		// get input stream from file
		BufferedReader br 
			= new BufferedReader ( new FileReader ( src ) ) ;

		// get output stream to socket
		BufferedWriter bw 
		= new BufferedWriter ( new OutputStreamWriter ( socket.getOutputStream() ) ) ;

		// file ==> socket
		String line ;
		while ( (line = br.readLine()) != null )
		{
			bw.write( line );   
			bw.newLine()    ;
			bw.flush  ()    ;
		}

		// close streams
		br.close();
		bw.close();
		socket.close() ; 
	 }

	catch( Exception e ) 
	{
		System.out.println( e.toString() );
	}	
	} // main
}
