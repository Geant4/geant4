// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: g4vrmlview.java,v 1.2 1999-12-15 14:54:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
import java.io.*;
import java.net.*;

//------------------//
// class g4vrmlview //
// ( main() )       // 
//------------------//
public class g4vrmlview
{
	public static void main( String[] args )
	{
	try{
		// CONST 
		final String VERSION             = "1.00"   ;
		final String DATE                = "August 19, 1997";

		final int    PORT_NO             = 40801    ;
		final String OUTPUT_FILE_HEAD    = "g4"     ;
		final String OUTPUT_FILE_EXT     = "wrl"    ;
		final int    MAX_TRIAL           = 10       ;

		// local
		int          portNo              = PORT_NO  ; 

		// argument checking
		if(  args.length != 1 && args.length != 2 )
		{
			System.out.println( "-------------------------------");
			System.out.println( " G4VRMLView version  " + VERSION );
			System.out.println( "     " + DATE                   );
			System.out.println( "-------------------------------");
			System.out.println( "Usage: java g4vrmlview browser_name [port_number]");
			System.out.println( "  Browser_name: netscape, vrweb, etc, or NONE");
			return ;
		}
		
		// VRML browser
		String       browser = new String ( args[0] ) ;

		// port number
		if( args.length == 2 ) 
		{
			portNo = Integer.parseInt( args[1] );
		}

		// make a server socket
		ServerSocket  ss = null ;
		for ( int i = 0 ; i < MAX_TRIAL ; i++ )
		{
			try 
			{
				ss = new ServerSocket( portNo );
				System.out.println( "Waiting for requests at port  " +portNo + " ...");
				break ;
			}
			catch ( Exception e ) 
			{
				portNo++ ; 
				if( i >= MAX_TRIAL ) 
				{ 
					System.out.println( "Sockets are not available.");
					return ;
				}
			}
		} // for 


		// open connection and invoke thread
		int nSpawn = 0 ;
		while( true )
		{
			Socket socket = ss.accept();	nSpawn++ ;

			System.out.println( "Connection accepted by thread " + nSpawn ); 

			( new g4vrmlviewThread( socket, OUTPUT_FILE_HEAD, OUTPUT_FILE_EXT , browser )).start() ;

		} // while

	}
	catch ( Exception e ) 
	{
		System.out.println( e.toString() );
	}	
	} // main()

} // g4vrmlview


//------------------------//
// class g4vrmlviewThread //
//------------------------//
class g4vrmlviewThread extends Thread
{
	private final String NONE = "NONE" ; // no browser
	private  Socket     m_socket     ;
	private  String     m_outputFile     ;		
	private  String     m_browser    ;

	// constuctor
	public g4vrmlviewThread(	Socket socket          , 
					String outputFileHead  , 
					String outputFileExt   , 
					String browser     ) 
	{
		m_socket     = socket ;
		SetOutputFileName ( outputFileHead  , outputFileExt );
		m_browser    = new String ( browser    );
	}

	private void SetOutputFileName(	String outputFileHead,
					String outputFileExt  )
	{
		// temporary file name
		String  outputFile_tmp 
		= new String (	outputFileHead + 
				"."            +
				outputFileExt   ) ; 

		// for modification of temporary filename
		int n = 1 ; 

		// make a non-existing filename
		while ( true )
		{
			File file = new File ( outputFile_tmp );

			if ( !file.exists() ) 
			{ 
				break ; 
			} else {

				outputFile_tmp 
				= new String (	outputFileHead + 
						"_"            +
						(n++)          +
						"."            +
						outputFileExt   ); 
			}

		} // while

		// set decided filename to data field
		m_outputFile = new String ( outputFile_tmp ); 

	} // g4vrmlviewThread::setOutputFileName()


	// run ()
	public void run ()
	{
	try{
		// get input stream from socket
		BufferedReader br 
		= new BufferedReader ( new InputStreamReader ( m_socket.getInputStream() ) ) ;

		// get output stream to file
		BufferedWriter bw 
		= new BufferedWriter ( new FileWriter ( m_outputFile ) ) ;

		// socket ==> file
		String line ;
		while ( (line = br.readLine()) != null )
		{
			bw.write( line );   
			bw.newLine()    ;
			bw.flush  ()    ;
		}
		System.out.println( "VRML data is saved to file " +m_outputFile );

		// close streams and socket
		br.close();
		bw.close();
		m_socket.close() ; 

		// invoke browser
		if( !m_browser.equals(NONE) )
		{
			try 
			{

				File file = new File ( m_outputFile );

				// visualize created VRML file with browser	
				if ( file.exists() ) 
				{
					String  outputFileAbs  = new String ( file.getAbsolutePath() ) ;
					String[] command = { m_browser, outputFileAbs };
					System.out.println( "Command: " + command[0] + " " + command[1] );

					Process  exec_process    = null ;
	
					Runtime runtime = Runtime.getRuntime();
					exec_process 	= runtime.exec( command );
					exec_process.waitFor() ;
				} else {
						System.out.println( "Error: Failed to open file" + m_outputFile );
				}
			} 
			catch ( Exception e ) 
			{  
				System.out.println( e.toString() );
			}

		} 
		else 
		{
			System.out.println( "No browser was invoked" );
		}

	}
	catch( Exception e ) 
	{
		System.out.println( e.toString() );
	}	
	} // run()

} // g4vrmlviewThread 
