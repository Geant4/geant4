
      SUBROUTINE FFUSER(ikey)
*
*     The routine is called when a *key card is read
*
      CHARACTER*4  keyw
*
      call UHTOC(ikey,4,keyw,4)
      if (keyw.eq.'HIST') call uhinit		   
*                                                                               
      END                                                                       
