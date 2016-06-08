
      SUBROUTINE FFUSER(ikey)
*
*     The routine is called when a *key card is read
*
#include "calor.inc"
#include "histo.inc"
*
      CHARACTER*4  keyw
      CHARACTER*25 title      
*
      call UHTOC(ikey,4,keyw,4)
      if (keyw.eq.'HIST') then
        if (iexist(idhist).eq.1) call hdelet(idhist)
	title = 'Edep/Ebeam in absorber '//char(ichar('0')+(idhist-1))	
	call hbook1(idhist,title,nbinE,Emin,Emax,0.)
	iexist(idhist) = 1
      endif			   
*                                                                               
      END                                                                       
