
      SUBROUTINE UHINIT
*
*
#include "histo.inc"
*
      CHARACTER*50 title(16)
*
      data title /
     1    'energy deposit in absorber',					 
     2    'energy of charged secondaries at creation',			 
     3    'energy of gammas at creation (log10(ekin/MeV))',		 
     4	  '(transmit, charged) : kinetic energy at exit',		 
     5	  '(transmit, charged) : space angle at exit',			 
     6	  '(transmit, charged) : projected angle at exit',		 
     7	  '(transmit, charged) : projected position at exit',		 
     8	  '(transmit, neutral) : kinetic energy at exit',		 
     9	  '(transmit, neutral) : space angle at exit',
     +	  '(transmit, neutral) : projected angle at exit',		 
     1	  '(reflect , charged) : kinetic energy at exit',		 
     2	  '(reflect , charged) : space angle at exit',			 
     3	  '(reflect , charged) : projected angle at exit',		 
     4	  '(reflect , neutral) : kinetic energy at exit',		 
     5	  '(reflect , neutral) : space angle at exit',			 
     6	  '(reflect , neutral) : projected angle at exit'   /		 
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
