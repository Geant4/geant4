*    
*     ********************************************************************
*     * License and Disclaimer                                           *
*     *                                                                  *
*     * The  Geant4 software  is  copyright of the Copyright Holders  of *
*     * the Geant4 Collaboration.  It is provided  under  the terms  and *
*     * conditions of the Geant4 Software License,  included in the file *
*     * LICENSE and available at  http://cern.ch/geant4/license .  These *
*     * include a list of copyright holders.                             *
*     *                                                                  *
*     * Neither the authors of this software system, nor their employing *
*     * institutes,nor the agencies providing financial support for this *
*     * work  make  any representation or  warranty, express or implied, *
*     * regarding  this  software system or assume any liability for its *
*     * use.  Please see the license in the file  LICENSE  and URL above *
*     * for the full disclaimer and the limitation of liability.         *
*     *                                                                  *
*     * This  code  implementation is the result of  the  scientific and *
*     * technical work of the GEANT4 collaboration.                      *
*     * By using,  copying,  modifying or  distributing the software (or *
*     * any work based  on the software)  you  agree  to acknowledge its *
*     * use  in  resulting  scientific  publications,  and indicate your *
*     * acceptance of all terms of the Geant4 Software license.          *
*     ********************************************************************
*    
*    
*     $Id: jshape.F 67982 2013-03-13 10:36:03Z gcosmo $
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
