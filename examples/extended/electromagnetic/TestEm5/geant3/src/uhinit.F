
      SUBROUTINE UHINIT
*
*
#include "histo.inc"
*
      CHARACTER*50 title(50)
*
      data title /
     1    'energy deposit in absorber',					 
     2    'energy of charged secondaries at creation',
     3    'energy of gammas at creation (log10(ekin/MeV))',
     C    'dummy','dummy','dummy','dummy','dummy','dummy',
     +	  '(transmit, charged) : kinetic energy at exit',
     1	  '(transmit, charged) : energy fluence dE(MeV)/dOmega',
     2	  '(transmit, charged) : space angle dN/dOmega',	 
     3	  '(transmit, charged) : projected angle at exit',		 
     4	  '(transmit, charged) : projected position at exit',
     5	  '(transmit, charged) : radius at exit',
     C    'dummy','dummy','dummy','dummy',          		 
     +	  '(transmit, neutral) : kinetic energy at exit',
     1	  '(transmit, neutral) : energy fluence dE(MeV)/dOmega',	 
     2	  '(transmit, neutral) : space angle dN/dOmega',
     3	  '(transmit, neutral) : projected angle at exit',
     C    'dummy','dummy','dummy','dummy','dummy','dummy', 
     +	  '(reflect , charged) : kinetic energy at exit',
     1	  '(reflect , charged) : energy fluence dE(MeV)/dOmega',
     2	  '(reflect , charged) : space angle dN/dOmega',
     3	  '(reflect , charged) : projected angle at exit',
     C    'dummy','dummy','dummy','dummy','dummy','dummy', 
     +	  '(reflect , neutral) : kinetic energy at exit',
     1	  '(reflect , neutral) : energy fluence dE(MeV)/dOmega',	 
     2	  '(reflect , neutral) : space angle dN/dOmega', 
     3	  '(reflect , neutral) : projected angle at exit',
     C    'dummy','dummy','dummy','dummy','dummy','dummy','dummy' /
* 
      if (histo(idhist)) call hdelet(idhist)
*      
      vmin = valmin
      vmax = valmax      
      if (idhist.eq.3) then
        vmin = log10(1000*valmin)
	vmax = log10(1000*valmax)
      endif		
      call hbook1(idhist,title(idhist),nbBins,vmin,vmax,0.)
*      
      histo   (idhist) = .true.
      binWidth(idhist) = (valmax-valmin)/nbBins
      if (valunit.le.0.) valunit = 1.                
      histUnit(idhist) = valunit        
*
      END
