
      SUBROUTINE GUTRAK
*
*     User routine to control tracking of one track
*     Called by : GTREVE
*
#include "geant321/gckine.inc"
#include "geant321/gconst.inc"
#include "geant321/gctrak.inc"
*
#include "geomate.inc"
#include "runinfo.inc"
#include "histo.inc"
*
      logical transmit, reflect, charged, neutral 
*
* *** beginning of the track
*
      charged  = (charge.ne.0.)
      neutral  = .not.charged
*
      kflag = 0
      if (charge.eq.pkine(10))          kflag = 1
      if ((itra.eq.1).and.(istak.eq.0)) kflag = 2      
*
* *** perform tracking
*      
      CALL GTRACK
*
* *** end of the track
*
      transmit = (vect(1).ge. 0.49*xworld)
      reflect  = (vect(1).le.-0.49*xworld)
*
      if (transmit) itransmi = max(itransmi,kflag)
      if (reflect)  ireflect = max(ireflect,kflag)           
*
* *** histograms
*
      id = 0
*
*     energy spectrum at exit
*
      if (transmit.and.charged) id = 10
      if (transmit.and.neutral) id = 20
      if (reflect .and.charged) id = 30
      if (reflect .and.neutral) id = 40
      if (id.ne.0.and.histo(id)) then
        call hfill (id,gekin/histUnit(id),0.,1.)
      endif
*
*     space angle distribution at exit
*
      if (transmit.and.charged) id = 12
      if (transmit.and.neutral) id = 22
      if (reflect .and.charged) id = 32
      if (reflect .and.neutral) id = 42
      if (id.ne.0.and.histo(id)) then
        unit   = histUnit(id)
        theta  = acos(abs(vect(4)))
	dteta  = binWidth(id)*unit
	dOmega = twopi*sin(theta)*dteta
	weight = 0.
	if (dOmega.gt.(1./big)) weight = 1./dOmega	
        call hfill (id,theta/unit,0.,weight*unit*unit)
      endif
*
*     projected angle distribution at exit
*
      if (transmit.and.charged) id = 13
      if (transmit.and.neutral) id = 23
      if (reflect .and.charged) id = 33
      if (reflect .and.neutral) id = 43
      if (id.ne.0.and.histo(id)) then
        thety  = atan(vect(5)/abs(vect(4)))
        thetz  = atan(vect(6)/abs(vect(4)))
        call hfill (id,thety/histUnit(id),0.,1.)
        call hfill (id,thetz/histUnit(id),0.,1.)	
      endif
*
*     projected position at exit
*
      if (transmit.and.charged) id = 14
      if (id.ne.0.and.histo(id)) then
        call hfill (id,vect(2)/histUnit(id),0.,1.)
        call hfill (id,vect(3)/histUnit(id),0.,1.)	
      endif                        	                  
*
      END
