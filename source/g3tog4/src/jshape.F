*    
*     ********************************************************************
*     * DISCLAIMER                                                       *
*     *                                                                  *
*     * The following disclaimer summarizes all the specific disclaimers *
*     * of contributors to this software. The specific disclaimers,which *
*     * govern, are listed with their locations in:                      *
*     *   http://cern.ch/geant4/license                                  *
*     *                                                                  *
*     * Neither the authors of this software system, nor their employing *
*     * institutes,nor the agencies providing financial support for this *
*     * work  make  any representation or  warranty, express or implied, *
*     * regarding  this  software system or assume any liability for its *
*     * use.                                                             *
*     *                                                                  *
*     * This  code  implementation is the  intellectual property  of the *
*     * GEANT4 collaboration.                                            *
*     * By copying,  distributing  or modifying the Program (or any work *
*     * based  on  the Program)  you indicate  your  acceptance of  this *
*     * statement, and all its terms.                                    *
*     ********************************************************************
*    
*    
*     $Id: jshape.F,v 1.3 2001-07-11 09:59:03 gunter Exp $
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
