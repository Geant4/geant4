// read exp.data
//
// FW data (forward = small angle with respect to projectile)
//
const int NSetsFW = 4;
//
// LA data (lateral = larger angle with respect to projectile)
//
const int NSetsLA = 9;

// general purpose counter
//
static int NSets = 0;

float** AngleBin = 0;
int*    NPoints = 0; 
float** XMin = 0;
float** XMax = 0;
float** Y = 0;
float** EY = 0;

void ReadHARPData( std::string beam, std::string target, std::string energy, 
                   std::string secondary, 
		   std::string region ) // should be either "FW" or "LA"
{

   std::string dirname = "./harp-exp-data/";
   
   std::string filename = beam + "_" + target + "_" + energy + "GeV_" + secondary + "_" + region + ".dat";
   
   std::string file = dirname + filename;
   
   ifstream infile;
   infile.open( file.c_str() );
      
   if (AngleBin) 
   {
      delete [] AngleBin[0];
      delete [] AngleBin[1];
      delete [] AngleBin;
   }    
   
   if ( NPoints ) delete [] NPoints;   
   
   int i = 0;
   
   if ( NSets > 0 )
   {
      if ( XMin )
      {
         for ( i=0; i<NSets; i++ )
         {
            delete [] XMin[i];
         }
         delete [] XMin;
      }   
      if ( XMax )
      {
         for ( i=0; i<NSets; i++ )
         {
            delete [] XMax[i];
         }
         delete [] XMax;
      }   
      if ( Y )
      {
         for ( i=0; i<NSets; i++ )
         {
            delete [] Y[i];
         }
         delete [] Y;
      }      
      if ( EY )
      {
         for ( i=0; i<NSets; i++ )
         {
            delete [] EY[i];
         }
         delete [] EY;
      }
      NSets = 0;
   }   

   if ( region == "FW" )
   {
      NSets = NSetsFW;
   }
   else if ( region == "LA" )
   {
      NSets = NSetsLA;
   }
   else
   {
      return;
   }
   
   AngleBin = new float*[2];
   AngleBin[0] = new float[NSets];
   AngleBin[1] = new float[NSets]; 

   NPoints = new int[NSets];

   XMin = new float*[NSets];
   XMax = new float*[NSets];
   Y    = new float*[NSets];
   EY   = new float*[NSets]; 

   for ( i=0; i<NSets; i++ )
   {
      infile >> AngleBin[0][i] >> AngleBin[1][i];
      infile >> NPoints[i];
      
      // std::cout << "Angle Bin: " << AngleBin[0][i] << " " << AngleBin[1][i] << std::endl;
      // std::cout << "NPoints: " << NPoints[i] << std::endl;
      
      XMin[i] = new float[NPoints[i]];
      XMax[i] = new float[NPoints[i]];
      Y[i]    = new float[NPoints[i]];
      EY[i]   = new float[NPoints[i]];
      for ( int j=0; j<NPoints[i]; j++ )
      {
         infile >> XMin[i][j] >> XMax[i][j] >> Y[i][j] >> EY[i][j];
	 Y[i][j] *= 1000.;
	 EY[i][j] *= 1000.;
      }
   }
   
   infile.close();

   return;

}
