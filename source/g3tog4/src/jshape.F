*     This code implementation is the intellectual property of
*     the RD44 GEANT4 collaboration.
*    
*     By copying, distributing or modifying the Program (or any work
*     based on the Program) you indicate your acceptance of this statement,
*     and all its terms.
*    
*     $Id: jshape.F,v 1.1 1999-01-07 16:06:54 gunter Exp $
*     GEANT4 tag $Name: not supported by cvs2svn $
*    
*     
*--   Author : Jouko Vuoskoski, CERN, Jouko.Vuoskoski@cern.ch
*     
      SUBROUTINE JSHAPE(RSHAPE,SHAPE)
      implicit none
      real rshape
      integer ishape
      CHARACTER*4      SHAPE

      ishape = rshape
      IF     (ISHAPE.EQ.1) THEN
         SHAPE='BOX '
      ELSEIF (ISHAPE.EQ.2) THEN
         SHAPE='TRD1'
      ELSEIF (ISHAPE.EQ.3) THEN
         SHAPE='TRD2'
      ELSEIF (ISHAPE.EQ.4) THEN
         SHAPE='TRAP'
      ELSEIF (ISHAPE.EQ.5) THEN
         SHAPE='TUBE'
      ELSEIF (ISHAPE.EQ.6) THEN
         SHAPE='TUBS'
      ELSEIF (ISHAPE.EQ.7) THEN
         SHAPE='CONE'
      ELSEIF (ISHAPE.EQ.8) THEN
         SHAPE='CONS'
      ELSEIF (ISHAPE.EQ.9) THEN
         SHAPE='SPHE'
      ELSEIF (ISHAPE.EQ.10) THEN
         SHAPE='PARA'
      ELSEIF (ISHAPE.EQ.11) THEN
         SHAPE='PGON'
      ELSEIF (ISHAPE.EQ.12) THEN
         SHAPE='PCON'
      ELSEIF (ISHAPE.EQ.13) THEN
         SHAPE='ELTU'
      ELSEIF (ISHAPE.EQ.14) THEN
         SHAPE='HYPE'
      ELSEIF (ISHAPE.EQ.28) THEN
         SHAPE='GTRA'
      ELSEIF (ISHAPE.EQ.29) THEN
         SHAPE='CTUB'
      ELSE
         write(6,*)'(E) jshape: shape number ',ishape,' is unknown'
         SHAPE=' '
      ENDIF
      END
