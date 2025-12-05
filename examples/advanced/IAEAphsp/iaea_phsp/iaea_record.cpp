/*
 * Copyright (C) 2006 International Atomic Energy Agency
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
 *   Roberto Capote Noy, PhD
 *   e-mail: R.CapoteNoy@iaea.org (rcapotenoy@yahoo.com)
 *   International Atomic Energy Agency
 *   Nuclear Data Section, P.O.Box 100
 *   Wagramerstrasse 5, Vienna A-1400, AUSTRIA
 *   Phone: +431-260021713; Fax: +431-26007
 *
 *   Iwan Kawrakow, PhD
 *   e-mail iwan@irs.phy.nrc.ca
 *   Ionizing Radiation Standards
 *   Institute for National Measurement Standards
 *   National Research Council of Canada Ottawa, ON, K1A 0R6 Canada
 *   Phone: +1-613-993 2197, ext.241; Fax: +1-613-952 9865
 *
 **********************************************************************************
 * For documentation
 * see http://www-nds.iaea.org/reports-new/indc-reports/indc-nds/indc-nds-0484.pdf
 **********************************************************************************/
//#define DEBUG // Comment to avoid printing for every particle write or read

#if (defined WIN32) || (defined WIN64)
#include <iostream>  // so that namespace std becomes defined
#endif
#include <math.h>
#include <cstdio>

#if !(defined WIN32) && !(defined WIN64)
using namespace std;
#endif

#include "iaea_record.h"

short iaea_record_type::initialize()
{
  if(p_file == NULL) {
     fprintf(stderr, "\n ERROR: Failed to open Phase Space file \n");
     return (FAIL);
  }

  // Defines i/o logic and variable quantities to be stored
  // If the value is zero, then corresponding quantity is fixed
  ix = 1;
  iy = 1;
  iz = 1;
  iu = 1;
  iv = 1;
  iw = 1;
  iweight = 1;
  // Defines a number of EXTRA long variables stored
  // (EGS need 2 for incremental history number and LATCH)
  iextralong = 1;
  if( iextralong >= NUM_EXTRA_LONG)
  {
     fprintf(stderr, "\n ERROR: Increase NUM_EXTRA_LONG number in iaea_record.h\n");
     return (FAIL);
  }
  // Defines a number of EXTRA float variables stored (EGS could need 1 for ZLAST)
  iextrafloat = 0;
  if( iextrafloat >= NUM_EXTRA_FLOAT)
  {
     fprintf(stderr, "\n ERROR: Increase NUM_EXTRA_FLOAT number in iaea_record.h\n");
     return (FAIL);
  }

  return (OK);
}

short iaea_record_type::write_particle()
{
  float floatArray[NUM_EXTRA_FLOAT+7];
  IAEA_I32 longArray[NUM_EXTRA_LONG];

  char ishort = (char) particle;
  if(w < 0) ishort = -ishort; // Sign of w is stored in particle type

  if( fwrite(&ishort, sizeof(char),   1, p_file) != 1)
  {
    fprintf(stderr, "\n ERROR: write_particle: Failed to write particle type\n");
    return (FAIL);;
  }

  int reclength = sizeof(char);

  if(IsNewHistory > 0) energy *= (-1); // New history is signaled by negative energy

  floatArray[0] = energy;

  int i = 0;

  if(ix > 0) floatArray[++i] = x;
  if(iy > 0) floatArray[++i] = y;
  if(iz > 0) floatArray[++i] = z;
  if(iu > 0) floatArray[++i] = u;
  if(iv > 0) floatArray[++i] = v;
  if(iweight > 0) floatArray[++i] = weight;

  int j;
  for(j=0;j<iextrafloat;j++) floatArray[++i] = extrafloat[j];

  reclength += (i+1)*sizeof(float);

  if( fwrite(floatArray, sizeof(float), (size_t)(i+1), p_file) != (size_t) (i+1))
  {
     fprintf(stderr, "\n ERROR: write_particle: Failed to write FLOAT phsp data\n");
     return (FAIL);
  }

  if(iextralong > 0)
  {
     for(j=0;j<iextralong;j++) longArray[j] = extralong[j];
     reclength += iextralong*sizeof(IAEA_I32);
     if( fwrite(longArray, sizeof(IAEA_I32), (size_t)iextralong, p_file) != (size_t)iextralong)
     {
        fprintf(stderr, "\n ERROR: write_particle: Failed to write LONG phsp data\n");
        return (FAIL);
     }
  }

  if(reclength == 0) return(FAIL);

  #ifdef DEBUG
  // charge defined
  int iaea_charge[MAX_NUM_PARTICLES]={0,-1,+1,0,+1};
  int charge = iaea_charge[particle - 1];

  printf("\n Wrote a particle with a record lenght %d",reclength);
  printf("\n Q %d E %f X %f Y %f Z %f \n\t u %f v %f w %f W %f Part %d \n",
  charge, energy, x, y, z, u, v, w, weight, particle);
  if( iextrafloat > 0) printf(" EXTRA FLOATs:");
  for(j=0;j<iextrafloat;j++) printf(" F%i %f",j+1,extrafloat[j]);
  if( iextralong > 0)  printf(" EXTRA LONGs:");
  for(j=0;j<iextralong;j++) printf(" L%i %d", j+1,extralong[j]);
  printf("\n");
  #endif

  return(OK);
}

short iaea_record_type::read_particle()
{
  float floatArray[NUM_EXTRA_FLOAT+7];
  IAEA_I32 longArray[NUM_EXTRA_LONG+7];
  //(MACG) -Wshadow int i,j,is,reclength;
  int i,j,l,is,reclength;
  char ctmp;

  // IAEA_I32 pos = ftell(p_file); // To check file position

  if( fread(&ctmp, sizeof(char),   1, p_file) != 1) // particle type is always read
  {
    fprintf(stderr, "\n ERROR: read_particle: Failed to read particle type\n");
    return (FAIL);;
  }

  particle = (short) ctmp;

  is = 1; // getting sign of Z director cosine w
  if(particle < 0) {is = -1; particle = -particle;}

  reclength = sizeof(char);      // particle type is always read
  unsigned int rec_to_read = 1;    // energy is always read

  if(ix > 0) rec_to_read++;
  if(iy > 0) rec_to_read++;
  if(iz > 0) rec_to_read++;
  if(iu > 0) rec_to_read++;
  if(iv > 0) rec_to_read++;
  if(iweight > 0) rec_to_read++;
  if(iextrafloat>0) rec_to_read += iextrafloat;

  if( fread(floatArray, sizeof(float), rec_to_read, p_file) != rec_to_read)
  {
    fprintf(stderr, "\n ERROR: read_particle: Failed to read FLOATs \n");
    return (FAIL);;
  }

  reclength += rec_to_read*sizeof(float);


  IsNewHistory = 0;
  if(floatArray[0]<0) IsNewHistory = 1; // like egsnrc
  energy = fabs(floatArray[0]);

  i = 0;
  if(ix > 0) x = floatArray[++i];
  if(iy > 0) y = floatArray[++i];
  if(iz > 0) z = floatArray[++i];
  if(iu > 0) u = floatArray[++i];
  if(iv > 0) v = floatArray[++i];
  if(iweight > 0) weight = floatArray[++i];
  for(j=0;j<iextrafloat;j++) extrafloat[j] = floatArray[++i];

  if(iw > 0)
  {
      w = 0.f;
      double aux = (u*u + v*v);
      if (aux<=1.0) w = (float) (is * sqrt((float)(1.0 - aux)));
      else
      {
            aux = sqrt((float)aux);
            u /= (float)aux;
            v /= (float)aux;
      }
  }

  if(iextralong > 0)
  {
     if( fread(longArray, sizeof(IAEA_I32), (size_t)iextralong, p_file) != (size_t)iextralong)
     {
       fprintf(stderr, "\n ERROR: read_particle: Failed to read LONGS\n");
       return (FAIL);
     }
     //(MACG) -Wshadow for(int l=0,j=0;j<iextralong;j++) extralong[j] = longArray[l++];
     for(l=0,j=0;j<iextralong;j++) extralong[j] = longArray[l++];
     reclength += (iextralong)*sizeof(IAEA_I32);
  }

  #ifdef DEBUG
  // charge defined
  int iaea_charge[MAX_NUM_PARTICLES]={0,-1,+1,0,+1};
  int charge = iaea_charge[particle - 1];

  printf("\n Read a particle with a record lenght %d (New History: %d)",
               reclength,IsNewHistory);
  printf("\n Q %d E %f X %f Y %f Z %f \n\t u %f v %f w %f W %f Part %d \n",
  charge, energy, x, y, z, u, v, w, weight, particle);
  if( iextrafloat > 0) printf(" EXTRA FLOATs:");
  for(j=0;j<iextrafloat;j++) printf(" F%i %f",j+1,extrafloat[j]);
  if( iextralong > 0)  printf(" EXTRA LONGs:");
  for(j=0;j<iextralong;j++)  printf(" L%i %d",j+1,extralong[j]);
  printf("\n");
  #endif
  return(reclength);
}
