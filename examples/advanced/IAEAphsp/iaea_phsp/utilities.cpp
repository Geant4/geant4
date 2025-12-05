/*
 * Copyright 2000-2003 Virginia Commonwealth University
 * -----------------------------------------------------------------------------
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 *-----------------------------------------------------------------------------
 *
 *   AUTHORS:
 *
 *   Jeffrey Vincent Siebers
 *   e-mail: jsiebers@vcu.edu
 *   Virginia Commonwealth University
 *   401 College Street, P.O.Box 980058
 *   Richmond, Viriginia 23298-0058
 *   Phone: +1-804-6287771
 *
 *
*/
/*  General Utilities for CPP programs
      File Created:
            18-December-1995: Combined ok_check.cpp and some open_file
      Modification History:
            01-Feb-1996: JVS: filename used in open_file has extension only in openfile
            09-Feb-1996: JVS: Add eprintf: outputs to screen and a buffer called pbuffer
                                                pbuffer is a global whose memory must be allocated
                                                this is useful for creating a "history" file
      17-June-1996: JVS: change latex_string so will work with win95/bc5.0
                         cannot have for(int i=0,j=0; ). j will not increment.
      05-Sept-1996: JVS: Add allocate_pbuffer and print_runtime_info
      23-April-1997: JVS: Add interpolate
      06-Jan-1998: JVS: Add array_read
      11-June-1998: jvs: add eprintf and view_errors
      22-Sept-1998: JVS: fix memory leak in eprintf
      07-Dec-1998: JVS: Add clean_name
      18-Feb-1999: JVS: eliminate atof in array_read because of failures
      25-Feb-1999: JVS: array_read will now read numbers that start with .
      26-Feb-1999: JVS: Add array_read for strings
      02-March-1999: JVS: Add global eprint_mode so can quite eprintf statements
      24-March-1999: JVS: Modify clean_name so names cannot have *'s in them
      July 20, 1999: JVS: eprintf modified to use fprintf(stdout), rather than printf
      Dec 3, 1999: JVS: Add check_byte_order
      Jan 11, 2000: JVS: Modify clean_name so names cannot have / in them
      Jun 16, 2000: JVS: Change open_file so will read in .extension properly when a . is in
                         the path name
      April 25, 2001: JVS: Add cp(SourceFile,DestinationFile)
      May 29, 2001: JVS: clean_name now removes & as well
      May 31, 2002: add reverse_short_byte_order
      Feb 18, 2004: JVS: Add writeBigEndianBinaryFile()
      Feb 10, 2005: JVS: Add writeLittleEndianBinaryFile() and writeBinaryFile()
      Feb 11, 2005: JVS: Add reverse_int_byte_order
      April 21, 2005: JVS: Add readBinaryDataFromFile
*/
#if (defined WIN32) || (defined WIN64)
#include <iostream>  // so that namespace std becomes defined
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cctype>
#include <ctime>

#if !(defined WIN32) && !(defined WIN64)
using namespace std;
#endif

#include "utilities.h"

/* ************************************************************************** */
int reverse_int_byte_order(int xold)
{
   int xnew;
   char *pn   = (char *) &xnew;
   char *po   = (char *) &xold;
   pn[0] = po[3];
   pn[1] = po[2];
   pn[2] = po[1];
   pn[3] = po[0];
   return(xnew);
}
/* *************************************************************************** */
float reverse_float_byte_order(float xold)
{
   float xnew;
   char *pn   = (char *) &xnew;
   char *po   = (char *) &xold;
   pn[0] = po[3];
   pn[1] = po[2];
   pn[2] = po[1];
   pn[3] = po[0];
   return(xnew);
}
/* **************************************************************************** */
short reverse_short_byte_order(short xold)
{
   short xnew;
   char *pn   = (char *) &xnew;
   char *po   = (char *) &xold;
   pn[0] = po[1];
   pn[1] = po[0];
   return(xnew);
}
/* **************************************************************************** */
int check_byte_order()
{
  /* Determine the byte order on this machine */
  float ftest=1.0f; /* assign a float to 1.0 */
  char *pf = (char *) &ftest;
  // printf("\n \t %x %x %x %x", pf[0],pf[1],pf[2],pf[3]);
  if(pf[0] == 0 && pf[3] != 0)
  {
    // printf("\n\n Byte order: INTEL / ALPHA,LINUX -> LITLE_ENDIAN \n");
    return(LITTLE_ENDIAN);
  }else if(pf[0] != 0 && pf[3] == 0)
  {
    // printf("\n\n Byte order: OTHER (SGI,SUN-SOLARIS) -> BIG_ENDIAN \n ");
    return(BIG_ENDIAN);
  }
  else
  {
        printf("\n\n ERROR: indeterminate byte order");
        printf("\n \t %x %x %x %x", pf[0],pf[1],pf[2],pf[3]);
        return(UNKNOWN_ENDIAN);
  }
}
/* ************************************************** */
void print_runtime_info(int argc, char *argv[])
{  // print file header stuff
        printf("\n Command Line: ");
        for(int i=0; i<argc; i++) printf(" %s", argv[i]);
        // printf("\n Program %s Revision %f",  Prog_Name,Revision);
        printf("\n \t Copyright XXXX MCV");
        time_t t;
        t = time(NULL);
        printf("\n Run on %s\n",  ctime(&t) );
}
/* ****************************************************************** */
void allocate_pbuffer()
{
   pbuffer = (char *)malloc(MAX_BUFFER_SIZE * sizeof(char) +2);
      if(pbuffer == NULL)
      {
                  printf("\n Error Allocating Memory Buffer");
                  exit(EXIT_FAILURE);
   }
}
/* ***************************************************************** */
/* *********************************************************************** */
int advance(char *istr, int *sval, int len)
{                           /* advances past white-space in file */
      while( !isspace(istr[*sval]) && (*sval < len) )
            *sval+=1; /* advance to space */
      while( isspace(istr[*sval]) && (*sval < len) )
            *sval+=1; /* advance to next thing */
      if(*sval > len) return(FAIL); /* return 0 when fails */
      return (OK);
}
/* ********************************************************************** */
int my_isascii( int c )
{
   return( !(c < 0 || c > 0177) );
}
/* ********************************************************************* */
int clean_name(char *name)
{
   int len = strlen(name);
   char *tname = (char *) calloc(len+1, sizeof(char));
   if(tname == NULL)
   {
      eprintf("\n ERROR: memory allocation error");
      return(FAIL);
   }
   strcpy(tname, name);
   if(clean_name(tname, name)!= OK)
   {
      eprintf("\n ERROR: cleaning Name");
      return(FAIL);
   }
   free(tname);
   return(OK);
}
/* ********************************************************************* */
int clean_name(char *tmp_path, char *opath)
{
   /* remove spaces, *'s, :'s &'s and commas from the name */
   int len = strlen(tmp_path);
   int o_index=0;
   for(int i=0; i<len; i++)
   {
      if( isspace(tmp_path[i] ) )
      {
          if( o_index &&              // add a _ if not first char
              opath[o_index-1] != '_' ) // and if previous char not a _
         opath[o_index++] = '_';
      }
      else
      if( my_isascii( tmp_path[i] ) &&
          tmp_path[i] != '&'   &&
          tmp_path[i] != ','   &&
          tmp_path[i] != '*'   &&
          tmp_path[i] != '/'   &&
          tmp_path[i] != ':' )
          opath[o_index++] = tmp_path[i];
   }
   opath[o_index] = '\0'; /* terminate the string */

   return(OK);
}
/* *********************************************************************** */
FILE *open_file(char *filename, const char*extension, const char *access)
{
  char string[MAX_STR_LEN];
  FILE *strm = NULL;

  if(filename[0]=='\0')
   {
      printf("\n INPUT FILENAME (%s) > ",access);
      //(MACG)-Wunused-result fgets(string,MAX_STR_LEN,stdin);
      char* dummy = fgets(string, MAX_STR_LEN, stdin);
      (void) dummy;  //(MACG) no -Wunused-variable
      sscanf(string,"%s",filename);
      printf(" FILE %s opened \n", filename);
   }
   int len=strlen(filename);

   if( len + strlen(extension) >= MAX_STR_LEN)
   {
      printf("\n ERROR: String Length  of %s.%s Exceeds Maximum",
              filename, extension);
      return(NULL);
   }

   // char *filename1 = new(char[len+strlen(extension)+1]);

   const int filenameLength = len+strlen(extension)+1;
   //(MACG)-Wvla char *filename1 = new(char[filenameLength]);
   char *filename1 = new char[filenameLength];

   strcpy(filename1,filename); // temp filename for appending extension

   /* check if file name has .extension    */
   /* if it does not, add .extension to it */
   int i=len-1;
   while(i > 0 && filename[i--] != '.');
   //   printf("\n Comparing %s to %s", extension, filename+i+1);
   if(strcmp(extension, filename+i+1)  )
      strcat(filename1,extension);
   if( (strm = fopen(filename1, access) ) == NULL )
   {
      printf("\n ERROR OPENING FILE %s (mode %s)", filename1,access);
   }
   //(MACG)-Wvla delete(filename1);
   delete[] filename1;
   return(strm);
}
/* *********************************************************************** */
int ok_check(void)                          /* GETS RESPONSE FROM USER      */
{                                           /* IF OK TO DO SOMETHING        */
   char reply[MAX_STR_LEN];                 /* RETURNS 1 ONLY IF REPLY Y    */
                                            /* OR y ELSE RETURNS 0          */
   //(MACG)-Wunused-result fgets(reply,MAX_STR_LEN,stdin);
   char* dummy = fgets(reply,MAX_STR_LEN,stdin);
   (void) dummy;  //(MACG) no -Wunused-variable
   if(  ( strncmp(reply,"Y",1)==0 )||
   ( strncmp(reply,"y",1)==0 ))
     return(1);
   return(0);
}
/* ***********************************************************************  */
int ok_checks(char *string)
{
   printf("\n %s", string);
   return(ok_check());
}
/* ********************************************************************** */
#include <stdarg.h> // for va function
int pprintf(char *fmt, ... )
{
  va_list  argptr;         /* Argument list pointer   */
  char str[MAX_STR_LEN];           /* Buffer to build sting into   */
  int cnt;            /* Result of SPRINTF for return */
  va_start( argptr, fmt );      /* Initialize va_ functions   */
  //(MACG) Apple SDK deprecated:
  //cnt = vsprintf( str, fmt, argptr );   /* prints string to buffer   */
  cnt = vsnprintf(str, MAX_STR_LEN, fmt, argptr);  /* prints string to buffer */
  if(str[0] == '\0') return(0);
  printf("%s", str);   /* Send  to screen */
  if(pbuffer != NULL && strlen(pbuffer) + strlen(str) < MAX_BUFFER_SIZE)
     strcat(pbuffer,str);
  else
     printf("\n ERROR: pbuffer is full");
  va_end( argptr );         /* Close va_ functions      */
  return( cnt );         /* Return the conversion count   */
}
/* *********************************************************************** */
/* eprintf: for buffering error reports, writes error messages to a buffer,
   and, also can echo them to the screen (if set at compile time)
   at first instance, allocates memory for the error buffer              */
static char *ebuffer = NULL;
int eprintf(const char *fmt, ... )
{
  va_list  argptr;         /* Argument list pointer   */
  char str[MAX_STR_LEN];           /* Buffer to build sting into   */
  int cnt;            /* Result of SPRINTF for return */
  va_start( argptr, fmt );      /* Initialize va_ functions   */
  //(MACG) Apple SDK deprecated:
  //cnt = vsprintf( str, fmt, argptr );   /* prints string to buffer   */
  cnt = vsnprintf(str, MAX_STR_LEN, fmt, argptr);  /* prints string to buffer */
  if(str[0] == '\0') return(0);

  if(eprintf_mode==ON)
     fprintf(stdout,"%s", str);   /* Send  to screen */

  // allocate memory for the error message
  int ilen = 0;
  if(ebuffer != NULL)
  {
     ilen+=strlen(ebuffer);
     ebuffer = (char *) realloc(ebuffer, (ilen+strlen(str)+1)*sizeof(char));
  }
  else
     ebuffer = (char *) calloc(ilen+strlen(str)+1,sizeof(char));

  if(ebuffer == NULL)
  {
     printf("\n ERROR: ebuffer cannot be allocated in eprintf");
  }
  else
     strcat(ebuffer,str);
  va_end( argptr );         /* Close va_ functions      */
  return( cnt );         /* Return the conversion count   */
}
int view_errors(void)
{
   printf("\n%s\n",ebuffer);
   return(OK);
}
/* ************************************************************************** */
int latex_string(char *string, char *nstring)
{
  // adds \\ in front of % so % will show up in the comment when printed
  // with LaTeX
   // must change all %'s to \% for latex output
   // also, must do the same for $, &, # _  { and }
   int len = strlen(string);
   int sval=0;
   int j;
   while(isspace(string[sval]) )sval++; // remove space from start of string
   while(isspace(string[len-1]))len--; // remove space from end to string

   j=0;
   for(int i=sval;i<len;i++)
   {
      if(string[i]=='%' ||
         string[i]=='$' ||
         string[i]=='&' ||
         string[i]=='#' ||
         string[i]=='_' ||
         string[i]=='{' ||
         string[i]=='}' )
      {
         nstring[j++]='\\';
      }
      else
      if(string[i]=='<' ||
         string[i]=='>' )
      {
         nstring[j++]='$';
      }
      nstring[j++] = string[i];
      if(string[i]=='<' ||
      string[i]=='>' )
      {
         nstring[j++]='$';
      }
   }
   nstring[j]='\0';
/*   printf("\n string: %s", string);
   printf("\n nstring: %d %s",j, nstring); */
   return(OK);
}
/* ************************************************************************** */
float interpolate(float xh, float xl, float xm, float yh, float yl)
{
  return(yh - (xh-xm)/(xh-xl)*(yh-yl));
}
/* *********************************************************************** */
// #define DEBUG_ARRAY
/* ********************************************************************** */
int array_read(char *in_string, float *array, int max_array)
{
   char delimeter_string[MAX_STR_LEN];
   //(MACG) Apple SDK deprecated:
   // sprintf(delimeter_string," ,\t"); /* spaces, commas, and tabs */
   int cnt =
     snprintf(delimeter_string,MAX_STR_LEN," ,\t"); /* spaces, commas, and tabs */
   (void) cnt; //(MACG) silent -Wunused-result and -Wunused-variable

   char *p; /* pointer to string read in */
   p = strtok(in_string,delimeter_string);
   int i=0;
   if(p!=NULL)
   {
      array[i++]=(float)atof(p); /* get the first value */
      // if( sscanf(p,"%f",&array[i]) == 1) i++; // sscanf rounds values....
      do{             /* get remaining values */
         p = strtok(NULL,delimeter_string);
         if(p!=NULL)
         {
         //array[i++] = atof(p);
            if( sscanf(p,"%f",&array[i]) == 1) i++;
            // printf("\n Got Value of %f", array[i-1]);
         }
      }while(p!=NULL && i < max_array);
   }

#ifdef DEBUG_ARRAY
   printf("\n atof %d", i);
   for(int j=0; j<i; j++) {
     //     array[j] = 0.0001*round(1000.0*array[j]);
      printf("\n i = %d, %f",j,array[j]);
   }
#endif
   return(i);

}
int array_read(FILE *istrm, float *array, int max_array)
{
   // reads in an array of floats from a single line of istrm
   // returns the number of elements read in

   char in_string[MAX_STR_LEN];

   if(fgets(in_string, MAX_STR_LEN, istrm) == NULL ) return(FAIL);
#ifdef DEBUG_ARRAY
   printf("\nInput String\n %s", in_string);
#endif

   int slen = strlen(in_string);

   int k=0;
   while(isspace(in_string[k]) && k < slen) k++;
   if(slen==0 || !(isdigit(in_string[k])
              ||    in_string[k] == '.'
              ||    in_string[k] == '+'
              ||    in_string[k] == '-'))
   {
      return(0); // skip blank and non-numerical lines
   }
   int nread = array_read(in_string,array,max_array);

   return(nread);  // return the number of elements read
}
/* ********************************************************************** */
int copy(char *SourceFile, char *DestinationFile)
{
  /* Copies sourceFile to destination file like unix cp command */
  FILE *sStream = fopen(SourceFile, "rb");
  if(sStream == NULL)
  {
     perror("\n ERROR: copy: ");
     printf("\n ERROR: copy: Opening Source File %s",SourceFile);return(FAIL);
  }
  FILE *dStream = fopen(DestinationFile,"wb");
  if(dStream == NULL)
  {
     perror("\n ERROR: copy:");
     printf("\n ERROR: copy: Opening Destination File %s",DestinationFile);return(FAIL);
  }
  char buffer[1000];
  int nRead;
  do{
     nRead = fread(buffer, sizeof(char), 1000, sStream);
     if( nRead )
        fwrite(buffer, sizeof(char), nRead, dStream);
  }while( !feof(sStream) && !ferror(dStream) && !ferror(sStream) );
  if(ferror(sStream) || ferror(dStream) )
  {
     perror("ERROR: Copy: ");
     printf("\n ERROR: source %s, destination %s", SourceFile, DestinationFile);
     return(FAIL);
  }
  fclose(sStream);
  fclose(dStream);
  return(OK);
}
/* ********************************************************************** */
/* ************************************************************************************ */
int readBinaryDataFromFile(FILE *iStream, int nItemsToRead, float **arrayToRead, int swab_flag)
{
  // Reads binary data to stream, swab's if requested (1=swab, 0=don't swab)
  // Swab if needed...Put swabbed results in different array so no need to "unswab" when done
  // Allocate memory to read array into
   float *inputArray;
   inputArray = (float *) calloc(nItemsToRead,sizeof(float));
   if(inputArray == NULL) {
      printf("\n ERROR: Allocating memory for inputArray in readBinaryDataFromFile");
      return(FAIL);
   }
   if(OK != readBinaryDataFromFile(iStream, nItemsToRead, inputArray, swab_flag)) {
     printf("\n ERROR: Reading binary data from file"); return(FAIL);
   }
   *arrayToRead = inputArray;
   return(OK);
}
/* ************************************************************************************ */
int readBinaryDataFromFile(FILE *iStream, int nItemsToRead, float *inputArray, int swab_flag)
{
  // Reads binary data to stream, swab's if requested (1=swab, 0=don't swab)
  // Swab if needed...Put swabbed results in different array so no need to "unswab" when done
  // Allocate memory to read array into
   // Read in the array....
   int nRead=fread(inputArray,sizeof(float),nItemsToRead, iStream);
   if(nRead != nItemsToRead) {
     eprintf("\n ERROR: Wrong number read from file (%d %d)\n",
           nRead, nItemsToRead); return(FAIL);
   }
   // Check if need to swab the data
   if(swab_flag) // swab if swab_flag != 0
   {
      for(int index=0; index<nItemsToRead;index++)
      {
         inputArray[index] = reverse_float_byte_order( inputArray[index] );
      }
   }
   return(OK);
}
/* ***************************************************************************************** */
int writeBinaryFile(char *binaryFileName, int nItemsToWrite, float *arrayToWrite, int swab_flag)
{
   FILE *outputStream= fopen(binaryFileName,"wb");
   if (outputStream == NULL) {
     eprintf("\n ERROR: Cannot open file %s for writing\n",binaryFileName); return(FAIL);
   }
   if(OK != writeBinaryDataToFile(outputStream, nItemsToWrite, arrayToWrite, swab_flag) )
   {
     eprintf("\n ERROR: Writing Binary File"); return(FAIL);
   }
   fclose(outputStream);
   return(OK);
}
/* ************************************************************************************ */
int writeBinaryDataToFile(FILE *outputStream, int nItemsToWrite, float *arrayToWrite, int swab_flag)
{
  // Writes binary data to stream, swab's if requested (1=swab, 0=don't swab)
   // Swab if needed...Put swabbed results in different array so no need to "unswab" when done
   float *swabbedArray;
   if(swab_flag) // swab if swab_flag != 0
   {
      //
      swabbedArray = (float *) calloc(nItemsToWrite,sizeof(float));
      if(swabbedArray == NULL) {
      eprintf("\n ERROR: Allocating memory for swabbedArray in writeBinaryFile");
        return(FAIL);
      }
      for(int index=0; index<nItemsToWrite;index++)
      {
         swabbedArray[index] = reverse_float_byte_order( arrayToWrite[index] );
      }
   } else {
     swabbedArray = arrayToWrite;
   }
   // Check that writing positive number of items
   if(nItemsToWrite < 0 )
   {
     eprintf("\n ERROR: writeBinaryDataToFile: nItemsToWrite= %d < 0", nItemsToWrite); return(FAIL);
   }
   // Write the dose distribution
   int nWrite=fwrite(swabbedArray,sizeof(float),nItemsToWrite, outputStream);
   if(nWrite != nItemsToWrite) {
     eprintf("\n ERROR: Wrong number written to file (%d %d)\n",
           nWrite, nItemsToWrite); return(FAIL);
   }
   // free swabbedArray if it was allocated here
   if(swab_flag) {
     free(swabbedArray);
   }
   return(OK);
}
/* ***************************************************************************************************** */

int writeBigEndianBinaryFile(char *binaryFileName,  int nItemsToWrite, float *arrayToWrite)
{
   // Pinnacle doses always written in BIG_ENDIAN format... Check if need to swab.....
   int swab_flag = 0;
   switch( (check_byte_order()) )
   {
      case BIG_ENDIAN:
         break;
      case LITTLE_ENDIAN:
         swab_flag=1;
        break;
      default:
         eprintf("\n ERROR: Indeterminate Byte Order\n");
         return(FAIL);
   }
   if(OK != writeBinaryFile(binaryFileName, nItemsToWrite, arrayToWrite, swab_flag) )
   {
     printf("\n ERROR: Writing dose file %s", binaryFileName); return(FAIL);
   }
   return(OK);
}
/* ***************************************************************************************************** */
int writeLittleEndianBinaryFile(char *binaryFileName,  int nItemsToWrite, float *arrayToWrite)
{
   // Pinnacle doses always written in BIG_ENDIAN format... Check if need to swab.....
   int swab_flag = 0;
   switch( (check_byte_order()) )
   {
      case BIG_ENDIAN:
      swab_flag=1;
         break;
      case LITTLE_ENDIAN:
        break;
      default:
         eprintf("\n ERROR: Indeterminate Byte Order\n");
         return(FAIL);
   }
   if(OK != writeBinaryFile(binaryFileName, nItemsToWrite, arrayToWrite, swab_flag) )
   {
     printf("\n ERROR: Writing dose file %s", binaryFileName); return(FAIL);
   }
   return(OK);
}
/* ***************************************************************************************************** */
char *strnset(char *s, int ch, size_t n)
{  /* mimic strnset command in dos/ os2/ win / ... */
   for(int i=0; i< (int) n; i++)
   {
     if(s[i] == STR_NULL ) return(s); // return when find null
      s[i] = ch;
   }
   return(s);
}
int get_string(FILE *fspec, char *string)
{
#ifdef DEBUG
  int rvalue=fget_c_string(string, MAX_STR_LEN, fspec);
  printf("\n fget_c_string returns %s", string);
  return(rvalue);
#else
  return(fget_c_string(string, MAX_STR_LEN, fspec));
#endif
}
#define REWIND_STREAM 100
/* ************************************************************************** */
int fget_c_string(char *string, int Max_Str_Len, FILE *fspec)
{
   /* gets a string from the input and removes comments from it */
   /* allows comments in standard "c" syntax,
           starting with / *
      ending with  * /  */
   /* also allows c++ type comments, // causes rest of line to be skipped */

   int check;
   char comment_start[4]="/*";  /* signals start of comment */
   char comment_stop[4]="*/";   /* signals end of comment   */
   int clen; /* length of string for start/stop*/
   int ilen; /* length of input string */
   char *istring; /* input string */
   int olen; /* location on output string */
   int icnt; /* location on string */

   //   int n_pass = 0; /* Number of passes through the file looking for a value */

   olen = 0;
   /* allocate memory for input string */
   istring = (char *)calloc(Max_Str_Len,sizeof(char));
   if(istring == NULL)
   {
      printf("\n ERROR: Allocating memory for input string if fget_c_string");
      return(FAIL);
   }
   strnset(string,'\0',Max_Str_Len); /* null entire output string */
   strnset(istring,'\0',Max_Str_Len); /* null entire input string */

#ifdef DEBUG
   printf ("\n --------------fget_c_string");
#endif
   clen = strlen(comment_start);
   /* read in the line, verify that it exists */
   do{
      /* read in a line from the file */
      while(fgets(istring, Max_Str_Len, fspec) == NULL) /* output warning if not a valid read */
      {
#ifdef DEBUG
        printf("\n istring: %s", istring);
#endif
#ifdef ALLOW_REWIND
   if(n_pass) /* if already gone through file once looking for value, quit */
#endif
   {
#ifdef DEBUG
           //  printf ("\n***End of Input File in get_string, closing");
#endif
           //  printf ("\nERROR: Reading File : End of File On Read ");
           // fclose(fspec);
           free(istring);
           return(FAIL);
        }
#ifdef ALLOW_REWIND
        n_pass++;      /* increment the number of times through the file */
        rewind(fspec); /* rewind to the beginning of the file */
        free(istring);
        return(REWIND_STREAM);
#endif
      }
#ifdef DEBUG
        printf("\n istring: %s", istring);
#endif
      ilen = strlen(istring); /* length of input string */
      istring[ilen]='\0'; /* null terminate the string */
      if(ilen < clen) /* not possible to have comment on the line */
      {               /* so output the string as is */
         strcpy(string,istring);
         olen = strlen(string);
      }
      else
      {
        /* strip comments out of input string */
        icnt=0;
        do{
          check = 1;
          if(icnt < ilen - clen)  /* make sure have enough characters for start of comment */
             check = strncmp(istring+icnt,comment_start,clen); /* check if start of comment */
          if(check == 0) /* comment found for standard c syntax */
          {
            /* find end of comment */
            icnt+=clen; /* advance past comment delimeter */
            clen=strlen(comment_stop); /* get length of end of comment delimiter */
            /* look for end of comment till end of string */
            do{
               check = 1;
               if(icnt < ilen - clen)  /* make sure have enough characters for end of comment */
                  check = strncmp(istring+icnt, comment_stop,clen);
               if(check != 0)  /* if not end of comment */
               {
                  icnt++; /* increment location on string */
                  if(icnt>ilen) /* if advance past end of string, get a new one */
                  {
                     if(fgets(istring, Max_Str_Len, fspec) == NULL) /* output warning if not a valid read */
                     {
                        printf ("\nERROR: Reading File, looking for end of comment %s",comment_stop);
                        // fclose(fspec);
                        free(istring);
                        return(FAIL);
                     }
                     ilen = strlen(istring); /* get length of this new string */
                     /* null terminate the string */
                     istring[ilen]='\0';
                     icnt = 0; /* reset the counter to the start of the string */
                  }
               }
               else
               {
                  icnt+=clen; /* advance past comment delimiter */
               }
            }while(check != 0); /* end of comment found */
          }  /* end if */
          else /* check if comment is in c++ format */
          {
             check = strncmp(istring+icnt, "//",2);
             if(check == 0) /* c++ style comment found */
             {  /* skip till end of string */
                icnt = ilen;
                string[olen++]='\n';
                string[olen]='\0';
             }
             else /* is a valid character for the string */
             {
                string[olen++] = istring[icnt++]; /* append value to the string */
                string[olen]='\0';
             }
          }
        }while(icnt < ilen          &&    /* do till end of string */
               olen < Max_Str_Len); /* and output string not too long */
        /* check for only carriage return (should have been caught above) */
        if(olen == 1 && string[0] == '\n') olen = 0;

      }  /* end else */
   }while(olen == 0); /* do till read in a string */
   if(olen == Max_Str_Len)
   {
      printf ("\nERROR: Input line too long");
      // fclose(fspec);
      free(istring);
      return(FAIL);
   }
   free(istring);
   return(OK);
}
