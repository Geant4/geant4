
* 
*             TEST PROGRAM FOR GEANT/SHOWER STUDIES
*                                                                               
*             Authors    R.Brun, M.Maire  *********
*                             
* This program generates showers within a cylinder made of an
* homogeneous material.The statistics of shower profiles and
* particle's flux are computed and plotted.
* It has been used to do the GEANT/EGS comparison.(CERN/DD/85/1)
*                                               
*     The present release has been modified :            
*                                          
* - to increase the flexibility of the volume definition from
*   data cards.(see description of the common /PVOLUM/ )
*                  
* - to permit to define all kinds of primary particles,i.e electrons 
*   as well pions ( ===> hadronic showers ).
*   see data card KINE.         
*
*     
* History                  
* -------                                            
*
*  27-08-01 : added histo 21 and 22 : edep distribution
*  08-02-99 : adapted for Geant4 comparison: electromagnetic/test/TestEm2
*  17-02-98 : adapted to Unix + interactive + graphic
*  21-01-90 : Simplified version (i.e. near gexam1)
*  06-04-89 : Add Hydrogen in the material library (imate=2)
*  10-02-89 : new histgrams for the resolution of the energy profile
*                                      
