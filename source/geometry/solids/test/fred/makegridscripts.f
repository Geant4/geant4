      PROGRAM MakeGridScripts
C     ------------------------------
C     Produces a set of GEANT4 scripts using
C     the "Fred" test program "shadow" feature
C
      IMPLICIT NONE
      Integer Icount, idx, idy, idxy, idz
      Real dx, dy, dz, theta, phi
C
      Real ddxy(0:6)/0.0, 0.1, -0.1, 0.2, -0.2, 0.4, 0.6/
C
      Real zddx(14)/  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  0.001, 
     +               -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -0.001  /
      Real zddy(14)/  0.0, -0.1, -0.2, -0.4, -0.6, -1.0, -100.0,
     +                0.0,  0.1,  0.2,  0.4,  0.6,  1.0,  100.0  /
      Real zddz(4)/ 0.0, 0.1, 0.2, 0.4 /
C
100   FORMAT( '/vis~/camera/viewpoint ', F13.9, ' ', F14.9 )
101   FORMAT( '/gridGun/origin ', 3F13.8 )
102   FORMAT( '/gridGun/direction ', 3F13.8 )
103   FORMAT( '/run/beamOn 1' )
104   FORMAT( '# view number ', I4 )
105   FORMAT( '/fred/pause' )
106   FORMAT( '/gridGun/grid1', 3F13.8 )
107   FORMAT( '/gridGun/grid2', 3F13.8 )
C
      Icount = 0
C
C     --- Scripts where the origin is a z-plane of -2
C
      open(unit=10,file='grids1.script',form='FORMATTED')
      write(10,'(''/fred/gun GRID'')')
      write(10,'(''/fred/draw SHADOW'')')
C
      write(10,106) 4.0, 0.0, 0.0
      write(10,107) 0.0, 4.0, 0.0
C
      do idy = 0, 6
        do idx = 0, 6
C
	  dx = ddxy(idx)
	  dy = ddxy(idy)
C
          theta = atan( sqrt(dx**2+dy**2) )*180/3.14159265
	  if (idx.eq.0.and.idy.eq.0) then
	    phi = 0.0
	  else
  	    phi = atan2(dy,dx)*180/3.14159265
	  endif
C
	  write(10,100) theta, phi
	  write(10,101) -2.0-2.0*dx, -2.0-2.0*dy, -2.0
	  write(10,102) dx, dy, 1.0
	  write(10,103)
	  write(10,104) icount
	  write(10,105)
	  icount = icount + 1
C
        enddo
      enddo
      close(10)
C
C     --- Scripts where the origin is a z-plane of +2
C
      open(unit=10,file='grids2.script',form='FORMATTED')
      write(10,'(''/fred/gun GRID'')')
      write(10,'(''/fred/draw SHADOW'')')
C
      write(10,106) 2.0, 0.0, 0.0
      write(10,107) 0.0, 2.0, 0.0
C
      do idy = 0, 6
        do idx = 0, 6
C
	  dx = ddxy(idx)
	  dy = ddxy(idy)
C
          theta = atan( sqrt(dx**2+dy**2) )*180/3.14159265
	  if (idx.eq.0.and.idy.eq.0) then
	    phi = 0.0
	  else
  	    phi = atan2(dy,dx)*180/3.14159265
	  endif
C
	  write(10,100) theta, phi
	  write(10,101) -2.0+2.0*dx, -2.0+2.0*dy, +2.0
	  write(10,102) -dx, -dy, -1.0
	  write(10,103)
	  write(10,104) icount
	  write(10,105)
	  icount = icount + 1
C
        enddo
      enddo
      close(10)
C
C     --- Okay: origin is a plane parallel to z-axis
C
      open(unit=10,file='grids3.script',form='FORMATTED')
      write(10,'(''/fred/gun GRID'')')
      write(10,'(''/fred/draw SHADOW'')')
C
      write(10,107) 0.0, 0.0, 4.0
C
      do idz = 1, 4
        dz = zddz(idz)
        do idxy = 1, 14 
C
          dx = zddx(idxy)
	  dy = zddy(idxy)
	  if (abs(dx).lt.0.01) then
C
C           --- Special case: origin is y-plane
C
            write(10,106) 4.0, 0.0, 0.0
C
            theta = atan( 1.0/dz )*180/3.14159265
	    if (dy.lt.0) then
	      phi = 270
	    else
  	      phi = 90
	    endif
            write(10,100) theta, phi
	    write(10,101) -2.0, sign(2.0,dy), -2.0+2.0*dz
	    write(10,102) 0.0, sign(1.0,-dy), -dz
	    write(10,103)
	    write(10,104) icount
	    write(10,105)
	    icount = icount + 1
	  else
C
            write(10,106) 0.0, 4.0, 0.0
C
            if (dz.eq.0) then
	      theta = 90.0
	    else
              theta = atan( sqrt(dx**2+dy**2)/dz )*180/3.14159265
            endif
	    phi = atan2(dy,dx)*180/3.14159265
C
            write(10,100) theta, phi
	    write(10,101) sign(2.0,-dx), -2.0-2.0*dy, -2.0-2.0*dz
	    write(10,102) dx, dy, dz
	    write(10,103)
	    write(10,104) icount
	    write(10,105)
	    icount = icount + 1
          endif
	enddo
      enddo
C
      close(10)
      STOP
      END
