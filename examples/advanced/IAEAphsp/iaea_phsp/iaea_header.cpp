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
//
// compile in Linux (copy all sources to a separate directory)
// cc *.cpp -lm -lstdc++ -o test_IAEAphsp
//
// Tested in RED HAT LINUX (gcc,cc,g95,icc) and WINDOWS XP (Microsoft C++ v4.2)
//
#if (defined WIN32) || (defined WIN64)
#include <iostream>  // so that namespace std becomes defined
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#if !(defined WIN32) && !(defined WIN64)
using namespace std;
#endif

#include "utilities.h"
#include "iaea_header.h"

int iaea_header_type::read_header ()
{
    char line[MAX_STR_LEN];

    if(fheader==NULL)
    {
      printf("\n ERROR: Unable to open header file \n");
        return(FAIL);
    }

  // ******************************************************************************
  // 1. PHSP format

      /*********************************************/
      if ( read_block(line,"FILE_TYPE") == FAIL )
      {
            printf("\nMandatory keyword FILE_TYPE is not defined in input\n");
            return FAIL;
      }
      else file_type = atoi(line);

      /*********************************************/
      if ( read_block(line,"CHECKSUM") == FAIL )
      {
            printf("\nMandatory keyword CHECKSUM is not defined in input\n");
            return FAIL;
      }
      else {
		  checksum = (IAEA_I64)atof(line);   //atol(line);  //atol() was not working properly in Windows
	  }

      /*********************************************/
      if ( read_block(line,"RECORD_LENGTH") == FAIL )
      {
            printf("\nMandatory keyword RECORD_LENGTH is not defined in input\n");
            return FAIL;
      }
      else record_length = atoi(line);

      /*********************************************/
      if ( read_block(line,"BYTE_ORDER") == FAIL )
      {
            printf("\nMandatory keyword BYTE_ORDER is not defined in input\n");
            return FAIL;
      }
      else byte_order = atoi(line);

      /*********************************************/
    if( get_blockname(line,"RECORD_CONTENTS") == FAIL)
    {
       printf("\nMandatory keyword RECORD_CONTENTS is not defined in input\n");
         return FAIL;
    }

    int i;
    for (i=0;i<9;i++) record_contents[i] = 0;

    for (i=0;i<9;i++)
    {
      if( get_string(fheader,line) == FAIL ) return FAIL;
      if( *line == SEGMENT_BEG_TOKEN ) break;
      record_contents[i] = atoi(line);
    };

    for(i=0;i<record_contents[7];i++)
    {
      if( get_string(fheader,line) == FAIL ) return FAIL;
      if( *line == SEGMENT_BEG_TOKEN ) break;
      extrafloat_contents[i] = atoi(line);
    }

    for(i=0;i<record_contents[8];i++)
    {
      if( get_string(fheader,line) == FAIL ) return FAIL;
      if( *line == SEGMENT_BEG_TOKEN ) break;
      extralong_contents[i] = atoi(line);
    }

    /*********************************************/
    if( get_blockname(line,"RECORD_CONSTANT") == FAIL)
    {
       printf("\nMandatory keyword RECORD_CONSTANT is not defined in input\n");
         return FAIL;
    }

    for (i=0;i<7;i++)
    {
        record_constant[i] = 32000.f;
        if(record_contents[i] > 0) continue;
        if( get_string(fheader,line) == FAIL ) return FAIL;
        if( *line == SEGMENT_BEG_TOKEN ) break;
        record_constant[i] = (float)atof(line);
    };

// ******************************************************************************
// 2. Mandatory description of the phsp

      if ( read_block(coordinate_system_description,"COORDINATE_SYSTEM_DESCRIPTION") == FAIL)
      {
            printf("\nMandatory keyword COORDINATE_SYSTEM_DESCRIPTION is not defined in input\n");
            return FAIL;
      }

  if(file_type == 1) // For event generators
  {
    /*********************************************/
      if ( read_block(line,"INPUT_FILE_FOR_EVENT_GENERATOR") == FAIL )
      {
            printf("\nMandatory keyword INPUT_FILE_FOR_EVENT_GENERATOR is not defined in input\n");
            return FAIL;
      }
  }

  if(file_type == 0) // for phsp files
  {
      /*********************************************/
      if ( read_block(line,"ORIG_HISTORIES") == FAIL )
      {
            printf("\nMandatory keyword ORIG_HISTORIES is not defined in input\n");
            return FAIL;
      }
      else {
          orig_histories = (IAEA_I64)atof(line); //atol(line);  //atol(0 was not working properly in Windows
          if( orig_histories == 0) printf(
           "\n The number of primary particles (ORIG_HISTORIES) is zero in the HEADER !\n");
      }

      /*********************************************/
      if ( read_block(line,"PARTICLES") == FAIL )
      {
            printf("\nMandatory keyword PARTICLES is not defined in input\n");
            return FAIL;
      }
      else nParticles = (IAEA_I64)atof(line);   //atol(line);  //atol(0 was not working properly in Windows

      /*********************************************/
      for(int itmp=0;itmp<MAX_NUM_PARTICLES;itmp++) particle_number[itmp]=0;
      IAEA_I64 npart;
	  //atol(line); commented  //atol(0 was not working properly in Windows
      if ( read_block(line,"PHOTONS") == OK )
            {npart = (IAEA_I64)atof(line); particle_number[0] = npart;}
      if ( read_block(line,"ELECTRONS") == OK )
            {npart = (IAEA_I64)atof(line); particle_number[1] = npart;}
      if ( read_block(line,"POSITRONS") == OK )
            {npart = (IAEA_I64)atof(line); particle_number[2] = npart;}
      if ( read_block(line,"NEUTRONS") == OK )
            {npart = (IAEA_I64)atof(line); particle_number[3] = npart;}
      if ( read_block(line,"PROTONS") == OK )
            {npart = (IAEA_I64)atof(line); particle_number[4] = npart;}
  }
// ******************************************************************************
// 3. Mandatory additional information
      if ( read_block(line,"IAEA_INDEX") == FAIL )
      {
            printf("\nMandatory keyword IAEA_INDEX is not defined in input\n");
            return FAIL;
      }
      else iaea_index = atoi(line);

      /*********************************************/
      if ( read_block(title,"TITLE") == FAIL )
      {
            printf("\nMandatory keyword TITLE is not defined in input\n");
            return FAIL;
      }

      /*********************************************/
      if ( read_block(machine_type,"MACHINE_TYPE") == FAIL)
      {
            printf("\nMandatory keyword MACHINE_TYPE is not defined in input\n");
            return FAIL;
      }

      /*********************************************/
      if ( read_block(MC_code_and_version,"MONTE_CARLO_CODE_VERSION") == FAIL )
      {
            printf("\nMandatory keyword MONTE_CARLO_CODE_VERSION is not defined in input\n");
            return FAIL;
      }

      /*********************************************/
      if ( read_block(line,"GLOBAL_PHOTON_ENERGY_CUTOFF") == FAIL )
      {
            printf("\nMandatory keyword GLOBAL_PHOTON_ENERGY_CUTOFF is not defined in input\n");
            return FAIL;
      }
      else global_photon_energy_cutoff = (float)atof(line);

      /*********************************************/
      if ( read_block(line,"GLOBAL_PARTICLE_ENERGY_CUTOFF") == FAIL )
      {
            printf("\nMandatory keyword GLOBAL_PARTICLE_ENERGY_CUTOFF is not defined in input\n");
            return FAIL;
      }
      else global_particle_energy_cutoff = (float)atof(line);

      /*********************************************/
      if ( read_block(transport_parameters,"TRANSPORT_PARAMETERS") == FAIL )
      {
            printf("\nMandatory keyword TRANSPORT_PARAMETERS is not defined in input\n");
            return FAIL;
      }

// ******************************************************************************
// 4. Optional description

      if ( read_block(beam_name,"BEAM_NAME") == FAIL )
            printf("ERROR reading BEAM_NAME\n");
      if ( read_block(field_size,"FIELD_SIZE") == FAIL )
            printf("ERROR reading FIELD_SIZE\n");
      if ( read_block(nominal_SSD,"NOMINAL_SSD") == FAIL )
            printf("ERROR reading NOMINAL_SSD\n");
      if ( read_block(variance_reduction_techniques,
                        "VARIANCE_REDUCTION_TECHNIQUES") == FAIL )
            printf("VARIANCE_REDUCTION_TECHNIQUES\n");
      if ( read_block(initial_source_description,"INITIAL_SOURCE_DESCRIPTION") == FAIL )
            printf("INITIAL_SOURCE_DESCRIPTION:\n");

      // Documentation sub-section
      /*********************************************/
      if ( read_block(MC_input_filename,"MC_INPUT_FILENAME") == FAIL )
            printf("MC_INPUT_FILENAME\n");
      if ( read_block(published_reference,"PUBLISHED_REFERENCE") == FAIL )
            printf("PUBLISHED_REFERENCE\n");
      if ( read_block(authors,"AUTHORS") == FAIL ) printf("AUTHORS\n");
      if ( read_block(institution,"INSTITUTION") == FAIL ) printf("INSTITUTION\n");
      if ( read_block(link_validation,"LINK_VALIDATION") == FAIL )
            printf("LINK_VALIDATION\n");
      if ( read_block(additional_notes,"ADDITIONAL_NOTES") == FAIL )
            printf("ADDITIONAL_NOTES\n");

// ******************************************************************************
// 5. Statistical information
      /*********************************************/
    if( get_blockname(line,"STATISTICAL_INFORMATION_PARTICLES") == OK)
    {
        for(i=0;i<MAX_NUM_PARTICLES;i++)
        {
              if( get_string(fheader,line) == FAIL ) return FAIL;
              if( *line == SEGMENT_BEG_TOKEN ) break;

              if(particle_number[i] == 0) continue;
              // -------------------------------------------------------
              // The fragment below replaces buggy sscanf() function
              int index0 =0, index1 =0, len = strlen(line), icnt =0;
              float fbuff[6];
              do
              {   if(advance(line, &index1, len) != OK) break; // Finding first number
                    fbuff[icnt++] = (float)atof(line+index1);
              if(advance(line, &index0, len) != OK) break; // Skipping spaces
              }   while (icnt < 6);
              // -------------------------------------------------------
              sumParticleWeight[i] = fbuff[0];
              minimumWeight[i]         = fbuff[1];
            maximumWeight[i]     = fbuff[2];
          averageKineticEnergy[i] = fbuff[3];
              minimumKineticEnergy[i] = fbuff[4];
              maximumKineticEnergy[i] = fbuff[5];
        }
    }

    if( get_blockname(line,"STATISTICAL_INFORMATION_GEOMETRY") == OK)
      {
        minimumX = minimumY = minimumZ = 32000.f;
        maximumX = maximumY = maximumZ = 0.f;

        for(i=0;i<3;i++)
        {
            if(record_contents[i] == 1)
            {
                  if( get_string(fheader,line) == FAIL ) return FAIL;
                  if( *line == SEGMENT_BEG_TOKEN ) break;

                // -------------------------------------------------------
                // The fragment below replaces buggy sscanf() function
                int index0 =0, index1 =0, len = strlen(line), icnt =0;
                float fbuff[2];
                do
                {   if(advance(line, &index1, len) != OK) break; // Finding first number
                      fbuff[icnt++] = (float)atof(line+index1);
                if(advance(line, &index0, len) != OK) break; // Skipping spaces
                } while (icnt < 2);
                // -------------------------------------------------------

                  if( i == 0) {minimumX = fbuff[0]; maximumX  = fbuff[1];}
                  if( i == 1) {minimumY = fbuff[0]; maximumY  = fbuff[1];}
                  if( i == 2) {minimumZ = fbuff[0]; maximumZ  = fbuff[1];}
            }
        }
      }

    return(OK);
}

int iaea_header_type::write_blockname(const char *blockname)
{
  return
  (fprintf(fheader,"%c%s%c\n",SEGMENT_BEG_TOKEN,blockname,SEGMENT_END_TOKEN));
}

int iaea_header_type::get_blockname(char *line, const char *blockname)
{
  char *begptr, *endptr;

  // printf("                       Reading block: %s ...\n",blockname);

  if(fheader==NULL)
  {
    printf("\n ERROR: Opening header file to Get Block \n"); return(FAIL);
  }
  rewind(fheader) ;

  while( get_string(fheader,line) == OK )
  {
    if( *line == SEGMENT_BEG_TOKEN )
      {
         begptr = line + 1 ;
             // printf("                       Read line %s ...\n",line);
         endptr = strchr(begptr,SEGMENT_END_TOKEN) ;
         if( endptr != NULL ) {
           *endptr = '\0' ;
               // printf("                       Comparing %s ...\n",begptr);
           if( strcmp(begptr,blockname) == 0 ) return(OK) ;
         }
      }
  };
  return(FAIL) ;
}

int iaea_header_type::get_block(char *lineread)
{
      int read = FAIL, count = 0;

      char line[MAX_STR_LEN];

      strcpy (lineread,""); // Deleting lineread contents
      while( get_string(fheader,line) == OK )
    {
        if( *line == SEGMENT_BEG_TOKEN ) break;
        strcat(lineread+count*MAX_NUMB_LINES,line); count++;
        read = OK;
    };
      return (read);
}

int iaea_header_type::read_block(char *lineread,const char *blockname)
{
    char line[MAX_STR_LEN];
    if( get_blockname(line,blockname) != OK) return FAIL;
      if( get_block(lineread) != OK) return FAIL;
      return OK;
}

int iaea_header_type::set_record_contents(iaea_record_type *p_iaea_record)
{
   int i;
   for(i=0;i<9;i++) record_contents[i] = 0;
   for(i=0;i<7;i++) record_constant[i] = 32000.f;

   if(p_iaea_record->ix > 0) record_contents[0] = 1;
   else record_constant[0] = p_iaea_record->x;

   if(p_iaea_record->iy > 0) record_contents[1] = 1;
   else record_constant[1] = p_iaea_record->y;

   if(p_iaea_record->iz > 0) record_contents[2] = 1;
   else record_constant[2] = p_iaea_record->z;

   if(p_iaea_record->iu > 0) record_contents[3] = 1;
   else record_constant[3] = p_iaea_record->u;

   if(p_iaea_record->iv > 0) record_contents[4] = 1;
   else record_constant[4] = p_iaea_record->v;

   if(p_iaea_record->iw > 0) record_contents[5] = 1;
   else record_constant[5] = p_iaea_record->w;

   if(p_iaea_record->iweight > 0) record_contents[6] = 1;
   else record_constant[6] = p_iaea_record->weight;

   if(p_iaea_record->iextrafloat>0) record_contents[7] = p_iaea_record->iextrafloat;
   if(p_iaea_record->iextralong>0) record_contents[8] = p_iaea_record->iextralong;

   record_length = 5; // To consider for particle type (1 byte) and energy (4 bytes)
   for(i=0;i<8;i++) record_length += record_contents[i]*sizeof(float);
   record_length -= 4; // 4 bytes substracted as w is not stored, just his sign
   record_length += record_contents[8]*sizeof(IAEA_I32);
   if(record_length > 0) return OK;
   else
   {
         printf("\nRECORD LENGTH IS ZERO, CHECK HEADER BLOCK RECORD_CONTENTS\n");
         return FAIL;
   }
}

int iaea_header_type::get_record_contents(iaea_record_type *p_iaea_record)
{
   int i;

   p_iaea_record->ix = 0;
   if(record_contents[0] == 1) p_iaea_record->ix = 1;
   else p_iaea_record->x = record_constant[0];

   p_iaea_record->iy = 0;
   if(record_contents[1] == 1) p_iaea_record->iy = 1;
   else p_iaea_record->y = record_constant[1];

   p_iaea_record->iz = 0;
   if(record_contents[2] == 1) p_iaea_record->iz = 1;
   else p_iaea_record->z = record_constant[2];

   p_iaea_record->iu = 0;
   if(record_contents[3] == 1) p_iaea_record->iu = 1;
   else p_iaea_record->u = record_constant[3];

   p_iaea_record->iv = 0;
   if(record_contents[4] == 1) p_iaea_record->iv = 1;
   else p_iaea_record->v = record_constant[4];

   p_iaea_record->iw = 0;
   if(record_contents[5] == 1) p_iaea_record->iw = 1;
   else p_iaea_record->w = record_constant[5];

   p_iaea_record->iweight = 0;
   if(record_contents[6] == 1) p_iaea_record->iweight = 1;
   else p_iaea_record->weight = record_constant[6];

   p_iaea_record->iextrafloat = 0;
   if(record_contents[7] > 0) p_iaea_record->iextrafloat = record_contents[7];

   p_iaea_record->iextralong = 0;
   if(record_contents[8] > 0) p_iaea_record->iextralong = record_contents[8];

   record_length = 5; // To consider for particle type (1 bytes) and energy (4 bytes)
   for(i=0;i<8;i++) record_length += record_contents[i]*sizeof(float);
   record_length -= 4; // 4 bytes substracted as w is not stored, just his sign
   record_length += record_contents[8]*sizeof(IAEA_I32);
   if(record_length > 0) return OK;
   else
   {
         printf("\nWRONG DEFINED HEADER BLOCK RECORD_CONTENTS\n");
         return FAIL;
   }
}

void iaea_header_type::initialize_counters()
{
  nParticles = read_indep_histories = 0;

  for(int i=0;i<MAX_NUM_PARTICLES;i++)
  {
        particle_number[i] = 0;
        sumParticleWeight[i] = 0.;
        maximumKineticEnergy[i] = 0.;
        minimumKineticEnergy[i] = 32000.;
        averageKineticEnergy[i] = 0.;
        minimumWeight[i] = 32000.;
        maximumWeight[i] = 0.;
  }
  minimumX = minimumY = minimumZ = 32000.f;
  maximumX = maximumY = maximumZ = -32000.f;

}

void iaea_header_type::update_counters(iaea_record_type *p_iaea_record)
{
  if (p_iaea_record->x > maximumX )  maximumX = p_iaea_record->x;
  if (p_iaea_record->x < minimumX )  minimumX = p_iaea_record->x;

  if (p_iaea_record->y > maximumY )  maximumY = p_iaea_record->y;
  if (p_iaea_record->y < minimumY )  minimumY = p_iaea_record->y;

  if (p_iaea_record->z > maximumZ )  maximumZ = p_iaea_record->z;
  if (p_iaea_record->z < minimumZ )  minimumZ = p_iaea_record->z;


  nParticles++;

  if ( p_iaea_record->IsNewHistory > 0 )
      read_indep_histories += p_iaea_record->IsNewHistory;

  int i = p_iaea_record->particle-1;
  if( i >= 0 && i < MAX_NUM_PARTICLES ) {
      particle_number[i]++;
      sumParticleWeight[i] +=  p_iaea_record->weight;
      averageKineticEnergy[i] += p_iaea_record->weight*
                                 fabs(p_iaea_record->energy);
      if (p_iaea_record->weight > maximumWeight[i] )
            maximumWeight[i] = p_iaea_record->weight;
      if (p_iaea_record->weight < minimumWeight[i] )
            minimumWeight[i] = p_iaea_record->weight;

      if (fabs(p_iaea_record->energy) > maximumKineticEnergy[i] )
         maximumKineticEnergy[i] = fabs(p_iaea_record->energy);
      if (fabs(p_iaea_record->energy) < minimumKineticEnergy[i] )
         minimumKineticEnergy[i] = fabs(p_iaea_record->energy);
  }

}

void iaea_header_type::print_statistics()
{
   printf("\n *************************************** \n");
   printf("           IAEA PHSP STATISTICS          \n");
   printf(" *************************************** \n");


   printf("\n Number of primary particles: %llu\n",orig_histories);
   printf(" Number of statistically independent histories read so far %llu\n",read_indep_histories);
   printf(" Total number of particles: %llu\n\n",nParticles);

   if(particle_number[0]>0)
   {
         printf(" PHOTONS: %llu,",particle_number[0]);
         if(sumParticleWeight[0]>0) printf(" <E> = %10.4G,",

             averageKineticEnergy[0]/sumParticleWeight[0]);

         printf("  Min.KE = %10.4G, Max.KE = %10.4G\n",

             minimumKineticEnergy[0],maximumKineticEnergy[0]);

         printf("  Total Weight = %15.6G \n",sumParticleWeight[0]);

         printf("  Min.Weight = %12.6G, Max.Weight = %12.6G\n\n",

             minimumWeight[0],maximumWeight[0]);

   }
   if(particle_number[1]>0)
   {
         printf(" ELECTRONS: %llu,",particle_number[1]);
         if(sumParticleWeight[1]>0) printf(" <E> = %10.4G,",

             averageKineticEnergy[1]/sumParticleWeight[1]);

         printf("  Min.KE = %10.4G, Max.KE = %10.4G\n",

             minimumKineticEnergy[1],maximumKineticEnergy[1]);

         printf("  Total Weight = %15.6G \n",sumParticleWeight[1]);

         printf("  Min.Weight = %12.6G, Max.Weight = %12.6G\n\n",

             minimumWeight[1],maximumWeight[1]);

   }
   if(particle_number[2]>0)
   {
         printf(" POSITRONS: %llu,",particle_number[2]);
         if(sumParticleWeight[2]>0) printf(" <E> = %10.4G,",

             averageKineticEnergy[2]/sumParticleWeight[2]);

         printf("  Min.KE = %10.4G, Max.KE = %10.4G\n",

             minimumKineticEnergy[2],maximumKineticEnergy[2]);

         printf("  Total Weight = %15.6G \n",sumParticleWeight[2]);

         printf("  Min.Weight = %12.6G, Max.Weight = %12.6G\n\n",

             minimumWeight[2],maximumWeight[2]);

   }
   if(particle_number[3]>0)
   {
         printf(" NEUTRONS: %llu,",particle_number[3]);
         if(sumParticleWeight[3]>0) printf(" <E> = %10.4G,",

             averageKineticEnergy[3]/sumParticleWeight[3]);

         printf("  Min.KE = %10.4G, Max.KE = %10.4G\n",

             minimumKineticEnergy[3],maximumKineticEnergy[3]);

         printf("  Total Weight = %15.6G \n",sumParticleWeight[3]);

         printf("  Min.Weight = %12.6G, Max.Weight = %12.6G\n\n",

             minimumWeight[3],maximumWeight[3]);

   }
   if(particle_number[4]>0)
   {
         printf(" PROTONS: %llu,",particle_number[4]);
         if(sumParticleWeight[4]>0) printf(" <E> = %10.4G,",

             averageKineticEnergy[4]/sumParticleWeight[4]);

         printf("  Min.KE = %10.4G, Max.KE = %10.4G\n",

             minimumKineticEnergy[4],maximumKineticEnergy[4]);

         printf("  Total Weight = %15.6G \n",sumParticleWeight[4]);

         printf("  Min.Weight = %12.6G, Max.Weight = %12.6G\n\n",

             minimumWeight[4],maximumWeight[4]);

   }
   printf(" GEOMETRY STATISTICS\n");
   if(record_contents[0] == 1)
         printf("  %7.3f < X coordinate < %7.3f\n",minimumX,maximumX);
   else printf("  X (constant) = %7.3f",record_constant[0]);
   if(record_contents[1] == 1)
         printf("  %7.3f < Y coordinate < %7.3f\n",minimumY,maximumY);
   else printf("  Y (constant) = %7.3f",record_constant[1]);
   if(record_contents[2] == 1)
         printf("  %7.3f < Z coordinate < %7.3f\n",minimumZ,maximumZ);
   else printf("  Z (constant) = %7.3f",record_constant[2]);
   printf("\n *************************************** \n\n");
}

int iaea_header_type::check_byte_order()
{
  /* Determine the byte order on this machine */
  float ftest=1.0f; /* assign a float to 1.0 */
  char *pf = (char *) &ftest;
  // printf("\n \t %x %x %x %x", pf[0],pf[1],pf[2],pf[3]);
  if(pf[0] == 0 && pf[3] != 0)
  {
    // printf("\nByte order: INTEL / ALPHA,LINUX -> LITLE_ENDIAN \n");
    return(LITTLE_ENDIAN);
  }else if(pf[0] != 0 && pf[3] == 0)
  {
    // printf("\nByte order: OTHER (SGI,SUN-SOLARIS) -> BIG_ENDIAN \n ");
    return(BIG_ENDIAN);
  }
  else
  {
    printf("\nERROR: indeterminate byte order");
    printf("\n \t %x %x %x %x", pf[0],pf[1],pf[2],pf[3]);
    return(UNKNOWN_ENDIAN);
  }
}

int iaea_header_type::write_header()
{
  if(fheader==NULL)
  {
      printf("\n ERROR: Opening header file to write \n"); return(FAIL);
  }

  rewind(fheader);

  if( write_blockname("IAEA_INDEX") == FAIL ) return(FAIL);

  fprintf(fheader,"%i   // Test header\n\n",iaea_index);

  write_blockname("TITLE");fprintf(fheader,"%s \n\n",title);

  write_blockname("FILE_TYPE");fprintf(fheader,"0\n\n"); // phasespace is assumed

  checksum = (IAEA_I64)record_length * nParticles;

  write_blockname("CHECKSUM");fprintf(fheader,"%llu \n\n",checksum);

  write_blockname("RECORD_CONTENTS");
  fprintf(fheader,"   %2i     // X is stored ?\n",record_contents[0]);
  fprintf(fheader,"   %2i     // Y is stored ?\n",record_contents[1]);
  fprintf(fheader,"   %2i     // Z is stored ?\n",record_contents[2]);

  fprintf(fheader,"   %2i     // U is stored ?\n",record_contents[3]);
  fprintf(fheader,"   %2i     // V is stored ?\n",record_contents[4]);
  fprintf(fheader,"   %2i     // W is stored ?\n",record_contents[5]);

  fprintf(fheader,"   %2i     // Weight is stored ?\n",record_contents[6]);
  fprintf(fheader,"   %2i     // Extra floats stored ?\n",record_contents[7]);
  fprintf(fheader,"   %2i     // Extra longs stored ?\n",record_contents[8]);

  int i;
  for(i=0;i<record_contents[7];i++) {
      if(extrafloat_contents[i] == 0) fprintf(fheader,
        "   %2i     // Generic float variable stored in the extrafloat array [%2i] \n",
        extrafloat_contents[i], i);
      if(extrafloat_contents[i] == 1) fprintf(fheader,
        "   %2i     // XLAST variable stored in the extrafloat array [%2i] \n",
        extrafloat_contents[i], i);
      if(extrafloat_contents[i] == 2) fprintf(fheader,
        "   %2i     // YLAST variable stored in the extrafloat array [%2i] \n",
        extrafloat_contents[i], i);
      if(extrafloat_contents[i] == 3) fprintf(fheader,
        "   %2i     // ZLAST variable stored in the extrafloat array [%2i] \n",
        extrafloat_contents[i], i);
  }

  for(i=0;i<record_contents[8];i++) {
      if(extralong_contents[i] == 0) fprintf(fheader,
        "   %2i     // Generic integer variable stored in the extralong array [%2i] \n",
        extralong_contents[i], i);

      if(extralong_contents[i] == 1) {
        fprintf(fheader,
        "   %2i     // Incremental history number stored in the extralong array [%2i] \n",
        extralong_contents[i], i);
        // orig_histories = read_indep_histories;
      }
      if(extralong_contents[i] == 2) fprintf(fheader,
        "   %2i     // LATCH EGS variable stored in the extralong array [%2i] \n",
        extralong_contents[i], i);

      if(extralong_contents[i] == 3) fprintf(fheader,
        "   %2i     // ILB5 PENELOPE variable stored in the extralong array [%2i] \n",
        extralong_contents[i], i);

      if(extralong_contents[i] == 4) fprintf(fheader,
        "   %2i     // ILB4 PENELOPE variable stored in the extralong array [%2i] \n",
        extralong_contents[i], i);

      if(extralong_contents[i] == 5) fprintf(fheader,
        "   %2i     // ILB3 PENELOPE variable stored in the extralong array [%2i] \n",
        extralong_contents[i], i);

      if(extralong_contents[i] == 6) fprintf(fheader,
        "   %2i     // ILB2 PENELOPE variable stored in the extralong array [%2i] \n",
        extralong_contents[i], i);

      if(extralong_contents[i] == 7) fprintf(fheader,
        "   %2i     // ILB1 PENELOPE variable stored in the extralong array [%2i] \n",
        extralong_contents[i], i);
  }

  fprintf(fheader,"\n");

  write_blockname("RECORD_CONSTANT");
  if(record_contents[0]==0)
    fprintf(fheader,"   %8.4f     // Constant X\n",record_constant[0]);
  if(record_contents[1]==0)
    fprintf(fheader,"   %8.4f     // Constant Y\n",record_constant[1]);
  if(record_contents[2]==0)
    fprintf(fheader,"   %8.4f     // Constant Z\n",record_constant[2]);
  if(record_contents[3]==0)
    fprintf(fheader,"   %8.5f     // Constant U\n",record_constant[3]);
  if(record_contents[4]==0)
    fprintf(fheader,"   %8.5f     // Constant V\n",record_constant[4]);
  if(record_contents[5]==0)
    fprintf(fheader,"   %8.5f     // Constant W\n",record_constant[5]);
  if(record_contents[6]==0)
    fprintf(fheader,"   %8.4f     // Constant Weight\n",record_constant[6]);

  fprintf(fheader,"\n");

  write_blockname("RECORD_LENGTH");fprintf(fheader,"%i\n\n",record_length);

  //(MACG) -Wshadow int byte_order = check_byte_order();
  byte_order = check_byte_order();
  write_blockname("BYTE_ORDER");fprintf(fheader,"%i\n\n",byte_order);

  write_blockname("ORIG_HISTORIES");
  if( orig_histories == 0) printf(

     "\n The number of primary particles (ORIG_HISTORIES) is zero in the HEADER !\n");


  fprintf(fheader,"%llu\n\n",orig_histories);

  write_blockname("PARTICLES");fprintf(fheader,"%llu\n\n",nParticles);

  if(particle_number[0]>0) {
        write_blockname("PHOTONS");fprintf(fheader,"%llu\n\n",particle_number[0]);}
  if(particle_number[1]>0) {
        write_blockname("ELECTRONS");fprintf(fheader,"%llu\n\n",particle_number[1]);}
  if(particle_number[2]>0) {
        write_blockname("POSITRONS");fprintf(fheader,"%llu\n\n",particle_number[2]);}
  if(particle_number[3]>0) {
        write_blockname("NEUTRONS");fprintf(fheader,"%llu\n\n",particle_number[3]);}
  if(particle_number[4]>0) {
        write_blockname("PROTONS");  fprintf(fheader,"%llu\n\n",particle_number[4]);}

  write_blockname("TRANSPORT_PARAMETERS"); fprintf(fheader,"\n");

  // 3. Mandatory additional information
  write_blockname("MACHINE_TYPE");  fprintf(fheader,"\n");

  write_blockname("MONTE_CARLO_CODE_VERSION"); fprintf(fheader,"\n");

  write_blockname("GLOBAL_PHOTON_ENERGY_CUTOFF");
  fprintf(fheader," %8.5f \n",global_photon_energy_cutoff);

  write_blockname("GLOBAL_PARTICLE_ENERGY_CUTOFF");
  fprintf(fheader," %8.5f \n",global_particle_energy_cutoff);

  write_blockname("COORDINATE_SYSTEM_DESCRIPTION"); fprintf(fheader,"\n");

  // 4. Optional information
  fprintf(fheader,"//  OPTIONAL INFORMATION\n\n");

  write_blockname("BEAM_NAME"); fprintf(fheader,"\n");

  write_blockname("FIELD_SIZE"); fprintf(fheader,"\n");

  write_blockname("NOMINAL_SSD"); fprintf(fheader,"\n");

  write_blockname("MC_INPUT_FILENAME"); fprintf(fheader,"\n");

  write_blockname("VARIANCE_REDUCTION_TECHNIQUES"); fprintf(fheader,"\n");

  write_blockname("INITIAL_SOURCE_DESCRIPTION"); fprintf(fheader,"\n");

  write_blockname("PUBLISHED_REFERENCE"); fprintf(fheader,"\n");

  write_blockname("AUTHORS"); fprintf(fheader,"\n");

  write_blockname("INSTITUTION"); fprintf(fheader,"\n");

  write_blockname("LINK_VALIDATION"); fprintf(fheader,"\n");

  write_blockname("ADDITIONAL_NOTES");
  fprintf(fheader,"%s\n","This is IAEA header as defined in the technical ");
  fprintf(fheader,"%s\n","report IAEA(NDS)-0484, Vienna, 2006");

  fprintf(fheader,"\n");
  // 5. Statistical information
  write_blockname("STATISTICAL_INFORMATION_PARTICLES");
  fprintf(fheader,

  "//        Weight        Wmin       Wmax       <E>         Emin         Emax    Particle\n");
  double eaver; char buffer[15];
  for(i=0;i<MAX_NUM_PARTICLES;i++)
  {
        if( particle_number[i] == 0 ) continue;

        switch (i)
        {
              case 0:
                    strcpy(buffer," PHOTONS");
                    break;
              case 1:
                    strcpy(buffer," ELECTRONS");
                    break;
              case 2:
                    strcpy(buffer," POSITRONS");
                    break;
              case 3:
                    strcpy(buffer," NEUTRONS");
                    break;
              case 4:
                    strcpy(buffer," PROTONS");
                    break;
        }

        eaver = 0.;
        if(sumParticleWeight[i]>0)
            eaver = averageKineticEnergy[i]/sumParticleWeight[i];
        fprintf(fheader,"  %15.6G %10.4G %10.4G %10.4G    %10.4G  %10.4G  %s\n",
        sumParticleWeight[i],minimumWeight[i],maximumWeight[i],
        eaver,minimumKineticEnergy[i],maximumKineticEnergy[i],buffer);
  }
  fprintf(fheader,"\n");

  write_blockname("STATISTICAL_INFORMATION_GEOMETRY");
  if(record_contents[0] == 1) fprintf(fheader," %G  %G\n",minimumX,maximumX);
  if(record_contents[1] == 1) fprintf(fheader," %G  %G\n",minimumY,maximumY);
  if(record_contents[2] == 1) fprintf(fheader," %G  %G\n\n",minimumZ,maximumZ);

  return(OK);

}

int iaea_header_type::print_header ()
{

    if(checksum == 0) printf("\n NEW PHASE SPACE FILE WILL BE CREATED\n");

    printf("\n\nIAEA_INDEX: %i\n",iaea_index);
    printf("TITLE: %s \n",title);

  // ******************************************************************************
  // 1. PHSP format

      /*********************************************/
    if(file_type == 0) printf("FILE TYPE: PHASESPACE \n");
    if(file_type == 1) printf("FILE TYPE: GENERATOR \n");

    if(checksum>0) printf("CHECKSUM: %llu\n",checksum);

    printf("RECORD LENGTH: %i\n",record_length);

      if(byte_order > 0) printf("BYTE ORDER: %i\n",byte_order);

    int i;
    printf("\nRECORD_CONTENTS:\n");
    for (i=0;i<7;i++)
    {
            if(record_contents[i] == 0)
            printf(" // Variable %1i is constant\n",i+1);
    }

    if(record_contents[7] > 0)
        printf(" // %1i extra FLOAT variable(s) defined\n",record_contents[7]);

    for(i=0;i<record_contents[7];i++)
    {
      if(extrafloat_contents[i] == 0) printf(
        " // Generic float variable stored in the extrafloat array [%2i] \n",i);
      if(extrafloat_contents[i] == 1) printf(
        " // XLAST variable stored in the extrafloat array [%2i] \n",i);
      if(extrafloat_contents[i] == 2) printf(
        " // YLAST variable stored in the extrafloat array [%2i] \n",i);
      if(extrafloat_contents[i] == 3) printf(
        " // ZLAST variable stored in the extrafloat array [%2i] \n",i);
    }

    if(record_contents[8] > 0)
        printf(" // %1i extra LONG variable(s) defined\n",record_contents[8]);

    for(i=0;i<record_contents[8];i++)
    {
      if(extralong_contents[i] == 0) printf(
        " // Generic integer variable stored in the extralong array [%2i] \n",i);
      if(extralong_contents[i] == 1) printf(
        " // Incremental history number stored in the extralong array [%2i] \n",i);
      if(extralong_contents[i] == 2) printf(
        " // LATCH EGS variable stored in the extralong array [%2i] \n",i);
      if(extralong_contents[i] == 3) printf(
        " // ILB5 PENELOPE variable stored in the extralong array [%2i] \n",i);
      if(extralong_contents[i] == 4) printf(
        " // ILB4 PENELOPE variable stored in the extralong array [%2i] \n",i);
      if(extralong_contents[i] == 5) printf(
        " // ILB3 PENELOPE variable stored in the extralong array [%2i] \n",i);
      if(extralong_contents[i] == 6) printf(
        " // ILB2 PENELOPE variable stored in the extralong array [%2i] \n",i);
      if(extralong_contents[i] == 7) printf(
        " // ILB1 PENELOPE variable stored in the extralong array [%2i] \n",i);
	}

    /*********************************************/
    printf("\nRECORD_CONSTANT:\n");
    for (i=0;i<7;i++)
    {
        if(record_contents[i] > 0) continue;
        printf(" %8.4f // Constant variable # %1i\n",record_constant[i],i+1);
    };
    printf("\n");
// ******************************************************************************
// 2. Mandatory description of the phsp

      if( strncmp(coordinate_system_description,"                ",15) > 0 )
          printf("\nCOORDINATE_SYSTEM_DESCRIPTION: \n%s\n",
            coordinate_system_description);

      if(file_type == 1)
      {
            // For event generators
          printf("INPUT FILE for event generator: %s \n",input_file_for_event_generator);
            return OK;
      }
      printf("\n");
      printf("NUMBER OF PRIMARY PARTICLES: %llu \n",orig_histories);

      printf("PARTICLES: %llu \n",nParticles);

      if(particle_number[0] > 0) printf("PHOTONS: %llu \n",particle_number[0]);
      if(particle_number[1] > 0) printf("ELECTRONS: %llu \n",particle_number[1]);
      if(particle_number[2] > 0) printf("POSITRONS: %llu \n",particle_number[2]);
      if(particle_number[3] > 0) printf("NEUTRONS: %llu \n",particle_number[3]);
      if(particle_number[4] > 0) printf("PROTONS: %llu \n",particle_number[4]);
      printf("\n");

// ******************************************************************************
// 3. Mandatory additional information
      /*********************************************/
      if( strncmp(machine_type,"                ",15) > 0 )
        printf("MACHINE_TYPE: %s\n",machine_type);

      if( strncmp(MC_code_and_version,"                ",15) > 0 )
        printf("MONTE_CARLO_CODE_VERSION: %s \n",MC_code_and_version);

      printf("GLOBAL_PHOTON_ENERGY_CUTOFF: %8.5f \n",global_photon_energy_cutoff);
      printf("GLOBAL_PARTICLE_ENERGY_CUTOFF: %8.5f \n",global_particle_energy_cutoff);
      printf("\n");

      if( strncmp(transport_parameters,"                ",15) > 0 )
        printf("\nTRANSPORT_PARAMETERS:\n%s\n",transport_parameters);

// ******************************************************************************
// 4. Optional description
      if( strncmp(beam_name,"                ",15) > 0 )
        printf("BEAM_NAME: %s\n",beam_name);
      if( strncmp(field_size,"                ",15) > 0 )
        printf("FIELD_SIZE: %s\n",field_size);
      if( strncmp(nominal_SSD,"                ",15) > 0 )
        printf("NOMINAL_SSD: %s\n",nominal_SSD);
      if( strncmp(variance_reduction_techniques,"                ",15) > 0 )
        printf("VARIANCE_REDUCTION_TECHNIQUES:\n%s\n",
        variance_reduction_techniques);
      if( strncmp(initial_source_description,"                ",15) > 0 )
        printf("INITIAL_SOURCE_DESCRIPTION: \n%s\n",
        initial_source_description);

      // Documentation sub-section
      /*********************************************/
      if( strncmp(MC_input_filename,"                ",15) > 0 )
        printf("MC_INPUT_FILENAME: %s\n",MC_input_filename);
      if( strncmp(published_reference,"                ",15) > 0 )
        printf("PUBLISHED_REFERENCE: \n%s\n",published_reference);
      if( strncmp(authors,"                ",15) > 0 )
        printf("AUTHORS: \n%s\n",authors);
      if( strncmp(institution,"                ",15) > 0 )
        printf("INSTITUTION: \n%s\n",institution);
      if( strncmp(link_validation,"                ",15) > 0 )
        printf("LINK_VALIDATION: \n%s\n",link_validation);
      if( strncmp(additional_notes,"                ",15) > 0 )
        printf("ADDITIONAL_NOTES: \n%s\n",additional_notes);

// ******************************************************************************
// 5. Optional statistical information

      print_statistics();

    return(OK);
}
