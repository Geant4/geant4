      subroutine plex(ndots,ipart)
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
      xmin = log10(xlimits(1,ipart))
      xmax = log10(xlimits(2,ipart))
      Vc = 1.44*ZcpZA(1,ipart)*(ZcpZA(2,ipart)-ZcpZA(1,ipart)-1.)/
     &                         (1.2*ZcpZA(3,ipart)**(1./3.))
C      print *,' %%% Vc =',Vc
      do i = 1, ndots
         x = xmin + RNDM(.5)*(xmax-xmin)
         E = max(10.**x-Vc,0.)
         vk = sqrt(E*(2.*xmass+E))
         y = 0.
         do k = 1, 3
            if(vM(k,ipart).gt.0.) then
C -- 1/2 of sigmas!! --
               xM(k) = vM(k,ipart) + 0.5*HRNDM1(idrn)*vdM(k,ipart)
               xT(k) = vT(k,ipart) + 0.5*HRNDM1(idrn)*vdT(k,ipart)
               if(ipart.lt.5) then
                  y=y+2.*xM(k)/sqrt(pi*xT(k)**3)*sqrt(E)*exp(-E/xT(k))
               elseif(ipart.eq.5) then
                  y=y+xM(k)*sqrt(10.**x*(10.**x+493.646*2.))*
     &                   exp(-10.**x/xT(k))
               elseif(ipart.eq.6) then
                  y=y+xM(k)/(10.**xmax-10.**xmin)
               endif
            endif
         enddo
         xval(i) = x
         yval(i) = y
      enddo
      return
      end
