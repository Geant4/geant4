      subroutine plexk(ndots,ipart)
      vector vM(3,9)
      vector vdM(3,9)
      vector vT(3,9)
      vector vdT(3,9)
      vector ZcpZA(3,9)
      vector xlimits(3,9)
      vector vmass(9)
      vector xval
      vector yval
      dimension xM(3),xT(3)
      data idrn/99/
      data pi/3.141593/
C
      xmass = vmass(ipart)
      xmin = (xlimits(1,ipart))
      xmax = (xlimits(2,ipart))
      Vc = 1.44*ZcpZA(1,ipart)*(ZcpZA(2,ipart)-ZcpZA(1,ipart)-1.)/
     &                         (1.2*ZcpZA(3,ipart)**(1./3.))
C      print *,' %%% Vc =',Vc
      do i = 1, ndots
         x = xmin + RNDM(.5)*(xmax-xmin)
         E0 = x
         vp = sqrt(E0*(2.*xmass+E0))
         vk = (E0+vp)/2.
         E  = max(E0-Vc,0.)
         if(xmass.lt.900) then
            vk = vk + xmass
         endif
         y = 0.
         do k = 1, 3
            if(vM(k,ipart).gt.0.) then
C -- 1/2 of sigmas!! --
               xM(k) = vM(k,ipart) + 0.5*HRNDM1(idrn)*vdM(k,ipart)
               xT(k) = vT(k,ipart) + 0.5*HRNDM1(idrn)*vdT(k,ipart)
               if(ipart.lt.5) then
                  y=y+2.*xM(k)/sqrt(pi*xT(k)**3)*sqrt(E)*exp(-E/xT(k))
               elseif(ipart.eq.5) then
                  y=y+xM(k)*sqrt(E0*(E0+493.646*2.))*
     &                   exp(-E0/xT(k))
               elseif(ipart.eq.6) then
                  y=y+xM(k)/(xmax-xmin)*vp
               endif
            endif
         enddo
C#T         xval(i) = x
C#T         yval(i) = y
C#k
         xval(i) = vk
         yval(i) = y/vp
      enddo
      return
      end
