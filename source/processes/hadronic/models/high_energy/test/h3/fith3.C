*** Begin of history file: 28/02/92 - 11.58.47
call fith3.ftn01
set ndvx 506
set ndvy 506
set ygti 0.2
set ymgu 1.
set ymgl 1.
set xmgr 0.5
set xmgl 1.5
set ywin 0.8
set xwin 1.2
set gsiz 0.3
set asiz 0.25
set vsiz 0.25
set ylab 0.6
set xlab 1.2
title_gl 'pA at 360 GeV/c' 
vect/create vx01(30) r 0. 1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12. _
  13. 14. 15. 16. 17. 18. 19. 20. 21. 22. 23. 24. 25. 26. 27. 28. 29. 30.
vect/create vxe01(30) r 
vect/create vy01(19) r 0.43 0.21 0.11 0.064 0.084 0.022 0.022 100. _
  0.022 0.022 100. 100. 100. 100. 100. 100. 100. _
  100.  0.022
vect/create vye01(19) r 0.05 0.07 0.05 0.035 0.04 0.02 0.02 0. 0.02 0.02 _
  0. 0. 0. 0. 0. 0. 0. 0. 0.02
vect/create vy02(30) r 0.24 0.12 0.087 0.077 0.11 0.063 0.047 0.047 _
  0.041 0.036 0.038 0.01 0.0052 0.01 0.02 0.0052 100. 0.01 0.01 0.0052 _
  100. 0.0052 0.0052 100. 0.0052 100. 100. 100. 100. 0.0052 
vect/create vye02(30) r 0. 0.03 0.02 0.02 0.03 0.02 0.015 0.015 0.015 _
  0.015 0.015 0.008 0.005 0.008 0.01 0.005 0. 0.008 0.008 0.005 0. _
  0.005 0.005 0. 0.005 0. 0. 0. 0. 0.005
vect/create vx01m(14) r 0. 2. 4. 6. 8. 10. 12. 14. 16. 18. 20. 22. _
  24. 26.
vect/create vy01m(8) r 0.29 0.16 0.063 0.025 0.013 0.0058 0.0024 0.001
vect/create vy01h(7) r 0.34 0.16 0.063 0.014 0.0052 0.0022 0.00075
vect/create vy02m(14) r 0.16 0.11 0.073 0.053 0.04 0.031 0.023 0.016 _
  0.011 0.0067 0.0039 0.0026 0.0016 0.00084
vect/create vy02h(13) r 0.2 0.11 0.073 0.053 0.04 0.031 0.021 0.013 _
  0.0085 0.0049 0.0028 0.0016 0.00078 
vect/create vx03(28) r 0.075 0.125 0.175 0.225 0.275 0.325 0.375 0.425 _
  0.475 0.525 0.575 0.625 0.675 0.725 0.775 0.825 0.875 0.925 0.975 _
  1.025 1.075 1.125 1.175 1.225 1.275 1.325 1.375 1.425
vect/create vxe03(28) r
vect/create vy03(28) r 1. 22. 55. 69. 64. 82. 57. 65. 46. 45. 52. 24. _
  20. 22. 17. 4. 8. 6. 6. 2. 2. 3. 1. 2. 1. 1. 1. 1.
vect/create vye03(28) r 1. 5. 8. 8. 8. 9. 8. 8. 7. 7. 7. 5. 5. 5. 4. _
  2. 3. 2.5 2.5 1.4 1.4 1.7 1. 1.4 1. 1. 1. 1.
vect/create vy03m(28) r 5.87 21.9 41.3 59.2 71.9 76.5 74.0 67.3 61.2 _
  56.1 49.0 42.3 34.2 26.0 20.4 16.3 12.8 10.2 8.2 6.9 5.6 4.6 3.3 _
  2.3 1.0
vect/create vx04(18) r -1.3 -0.8 -0.3 0.2 0.7 1.2 1.7 2.2 2.7 3.2 3.7 _
  4.2 4.7 5.2 5.7 6.2 6.7 7.2
vect/create vxe04(18) r
vect/create vy04(18) r 100000. 0.26 0.76 1.52 1.87 2.9 2.95 3.11 2.74 _
  3.5 2.44 2.68 1.8 1.45 0.76 0.58 0.21 0.18
vect/create vye04(18) r 0. 0.09 0.18 0.27 0.3 0.34 0.37 0.37 0.34 0.4 _
  0.3 0.37 0.27 0.24 0.18 0.18 0.09 0.09
vect/create vy05(18) r 0.15 0.52 1.28 3.15 4.19 4.99 4.8 4.61 3.94 3.65 _
  3.43 2.97 2.06 1.58 1.19 0.52 0.11 100000.
vect/create vye05(18) r 0.06 0.06 0.12 0.18 0.21 0.24 0.24 0.24 0.21 0.21 _
  0.18 0.18 0.15 0.12 0.09 0.06 0.06 0.
vect/create vx04m(11) r -1.2 -0.4 0.4 1.2 2.0 2.8 3.6 4.4 5.2 6.0 6.8
vect/create vy04m(11) r 0.091 0.3 1.19 2.35 2.92 3.05 2.77 2.16 1.43 0.76 _
  0.21
vect/create vy05m(11) r 0.3 1.07 2.74 4.26 4.69 4.57 3.93 2.77 1.58 0.7 _
  0.15
vect/create vy04l(11) r 0.03 0.24 1.43 2.50 3.14 3.14 2.77 2.16 1.43 _
  0.76 0.09
vect/create vy05l(11) r 0.05 0.43 2.74 4.26 4.93 4.66 3.93 2.77 1.58 0.58 _
  0.09
vect/create vx06(18) r -0.8 -0.3 0.2 0.7 1.2 1.7 2.2 2.7 3.2 3.7 4.2 4.7 _
  5.2 5.7 6.2 6.7 7.2 7.7
vect/create vxe06(18) r
vect/create vy06(18) r 10.2 12.6 7.1 3.7 3.1 2.2 2.11 1.74 2.23 1.60 _
  1.87 1.46 1.49 1.02 1.16 0.76 0.85 0.50
vect/create vye06(18) r 6.2 4.6 2. 0.7 0.6 0.2 0.2 0.3 0.2 0.4 0.3 0.3 _
  0.3 0.25 0.3 0.4 0.5 0.45
vect/create vy07(18) r 20.3 21.4 15.1 8.6 5.3 3.7 3.1 2.5 2.3 2.3 2.1 _
  1.65 1.65 1.16 1.01 0.4 0.000001 0.65
vect/create vye07(18) r 10. 4. 3. 1. 0.5 0.5 0.4 0.3 0.3 0.3 0.2 0.2 0.2 _
  0.2 0.2 0.15 0. 0.35
vect/create vx06m(12) r -1.2 -0.4 0.4 1.2 2. 2.8 3.6 4.4 5.2 6. 6.8 7.6
vect/create vy06m(12) r 19.9 8.9 3. 2. 1.8 1.65 1.6 1.38 1.12 0.98 0.98 _
  0.5
vect/create vy07m(12) r 78. 29.3 8. 3.5 2.9 2.4 2.3 1.71 1.34 0.79 0.47 _
  0.44
vect/create vy06l(12) r 3.2 2.2 1.9 2. 2.1 1.77 1.6 1.39 1.16 0.74 0.6 _
  0.64
vect/create vy07l(12) r 6.5 3.6 3.5 3.5 3.2 2.7 2.3 1.77 1.16 0.58 0.39 _
  0.37
vect/create vx10(8) r -4.2 -3.2 -2.75 -2.22 -1.75 -1.25 -0.75 -0.25
vect/create vxe10(8) r
vect/create vy10(8) r 16.9 23.1 12.4 20.5 16.4 10.9 14.3 4.8
vect/create vye10(8) r 7. 5.2 1.4 2.8 2.8 1.7 2.1 1.7
vect/create vx11(11) r -5.25 -4.75 -4.25 -3.75 -3.25 -2.75 -2.25 _
  -1.75 -1.25 -0.75 -0.25
vect/create vxe11(11) r
vect/create vy11(11) r 15.6 20.1 29.1 23.9 25.5 28.1 25.1 19.8 17. 14.9 7.1
vect/create vye11(11) r 7. 7.3 3.5 5.2 3.8 2.8 2.1 2.1 1.7 2.8 2.1
vect/create vy12(8) r 1.86 2.45 2.75 2.71 2.98 2.82 2.68 1.57
vect/create vye12(8) r 0.39 0.14 0.18 0.11 0.11 0.32 0.14 0.25
vect/create vy13(11) r 1.3 1.32 1.79 2.14 2.2 2.27 2.29 2.54 2.39 2.29 1.46
vect/create vye13(11) r 0.11 0.14 0.14 0.07 0.07 0.04 0.04 0.04 0.04 0.14 _
  0.32
vect/create vx10m(15) r -0.2 -0.6 -1. -1.4 -1.8 -2.2 -2.6 -3. -3.4 -3.8 _
  -4.2 -4.6 -5. -5.4 -5.8
vect/create vy10m(15) r 6.6 10. 10.3 14.1 15.9 16.9 17.6 17.2 15.9 13.8 _
  10.7 8.6 6.9 5.2 3.8
vect/create vy11m(15) r 7.3 11.4 15.3 18.7 22.5 25.3 27.1 27.1 25. 20.8 _
  16.3 11.4 8. 5.5 3.5
vect/create vy12m(15) r 2.28 2.75 2.92 2.93 2.90 2.82 2.75 2.68 2.57 2.46 _
  2.32 2.21 2.07 1.82 1.54
vect/create vy13m(15) r 2.21 2.53 2.61 2.57 2.50 2.43 2.39 2.36 2.32 2.25 _
  2.18 2.14 2.07 1.96 1.82   
gon
clr
*fortran/file 60 disk$scratch:[xv.fesefeldt]fith3_01.ps
*meta 60 -113 
hi/set/maximum 1 1.0
hi/set/maximum 2 1.0
hi/set/minimum 1 1.e-5
hi/set/minimum 2 1.e-5
set xsiz 13.
set ysiz 8.
opt logy
zone 2 1
hi/pl 1
atitle '    ' 'P(n?g!)'
text 20. 0.2 '(a) Al' 0.25 0.
text 18. 0.00005 '"0# data' 0.20 0.
text 18. 0.00002 '- GHEISHA' 0.20 0.
graph/hplot/errors vx01 vy01 vxe01 vye01 19 20 0.15
set dmod 0
graph 8 vx01m vy01m c
set dmod 2
graph 7 vx01m vy01h c
set dmod 0
hi/pl 2
atitle 'Number of slow protons n?g!' '   '
text 20. 0.2 '(b) Au' 0.25 0.
text 4. 0.00007 'curves"j#' 0.20 0.
text 4. 0.00004 '-- MCMHA' 0.20 0.
text 4.5 0.00004 '-' 0.20 0.
text 4. 0.00002 '-- Hegab and Hufner' 0.20 0.
graph/hplot/errors vx01 vy02 vxe01 vye02 30 20 0.15
set dmod 0
graph 14 vx01m vy02m c
set dmod 2
graph 13 vx01m vy02h c
exec hardcopy
*meta 0
*fort/close 60
*fortran/file 60 disk$scratch:[xv.fesefeldt]fith3_03.ps
*meta 60 -113
clr
set dmod 0
add 3 3 3 165. 0.
set xsiz 8.
set ysiz 8.
opt liny
zone 1 1
hi/pl 3
atitle 'Momentum p "M#GeV/c"N#' 'Number of events/(0.05 GeV/c)'
text 0.75 80. '"0# data' 0.20 0.
text 0.75 70. '- GHEISHA' 0.20 0.
text 0.75 60. 'curve"j# MCMHA' 0.20 0. 
graph/hplot/errors vx03 vy03 vxe03 vye03 28 20 0.17
graph 28 vx03 vy03m c
exec hardcopy
*meta 0
*fort/close 60
*fortran/file 60 disk$scratch:[xv.fesefeldt]fith3_04.ps 
*meta 60 -113
clr
add 5 5 5 10. 0.
vscale vy05 10. vy05
vscale vye05 10. vye05
vscale vy05m 10. vy05m
vscale vy05l 10. vy05l
hi/set/minimum 5 0.01
hi/set/maximum 5 100.
set xwin 2.
set xsiz 15.5
set ysiz 10.5
opt logy
zone 2 1
set dmod 0
hi/pl 5
set dmod 0
hi/pl 4 s
atitle '  ' '(1/N?ev!) dN/dy'
text 5. 40. 'Au(*10)' 0.25 0.
text 3. 1.  'Al' 0.25 0.
text 1. 0.2 '"0# data' 0.20 0.
text 1. 0.1 '- GHEISHA' 0.20 0.
text 1. 0.05 'curves"j#' 0.20 0.
text 1. 0.03 '-- MCMHA' 0.20 0.
text 1.2 0.03 '-' 0.20 0.
text 1. 0.02 '-- Lund' 0.20 0.
text -1. 65. '(a)' 0.25 0.
graph/hplot/errors vx04 vy04 vxe04 vye04 18 20 0.15
graph/hplot/errors vx04 vy05 vxe04 vye05 18 20 0.15
set dmod 0
graph 11 vx04m vy04m c
set dmod 2
graph 11 vx04m vy04l c
set dmod 0
graph 11 vx04m vy05m c
set dmod 2
graph 11 vx04m vy05l c
set dmod 0
add 7 7 7 10. 0.
add 9 9 9 10. 0.
hi/set/maximum 9 1000.
hi/set/minimum 9 0.1
hi/set/minimum 7 0.1
vscale vy07 10. vy07
vscale vye07 10. vye07
vscale vy07m 10. vy07m
vscale vy07l 10. vy07l
hi/pl 9
hi/pl 8 s
set dmod 0
hi/pl 7 s
hi/pl 6 s
atitle 'Rapidity y' 'Ratio R(y)'
text 2. 4. 'Al' 0.25 0.
text 2. 70. 'Au(*10)' 0.25 0.
text 3. 650. '"0# data' 0.20 0.
text 3. 400. '- GHEISHA' 0.20 0.
text -1. 0.4 'curves"j#' 0.20 0.
text 1. 0.4 '-- MCMHA' 0.20 0.
text 1.15 0.4 '-' 0.20 0.
text 1. 0.2 '-- Lund' 0.20 0.
text 1. 650. '(b)' 0.25 0. 
graph/hplot/errors vx06 vy06 vxe06 vye06 18 20 0.15
graph/hplot/errors vx06 vy07 vxe06 vye07 18 20 0.15
set dmod 0
graph 12 vx06m vy06m c
set dmod 2
graph 12 vx06m vy06l c
set dmod 0
graph 12 vx06m vy07m c
set dmod 2
graph 12 vx06m vy07l c
exec hardcopy
*meta 0
*fort/close 60
*fortran/file 60 disk$scratch:[xv.fesefeldt]fith3_10.ps 
*meta 60 -113
clr
set xwin 1.2
hi/set/maximum 10 40.
hi/set/maximum 11 40.
hi/set/maximum 12 4.
hi/set/maximum 13 4. 
add 10 10 10 1.1 0.
add 11 11 11 1.1 0.
add 12 12 12 1.1 0.
add 13 13 13 1.1 0.
set xsiz 10.5
set ysiz 10.5
opt liny
zone 2 2
set dmod 0
hi/zoom 10 ! 1 12
atitle '  ' '"L#n?s!"G#'
text -5. 35. '(a) Al' 0.25 0.
graph/hplot/errors vx10 vy10 vxe10 vye10 8 20 0.17
graph 15 vx10m vy10m c
hi/zoom 11 ! 1 12
text -5. 35. '(b) Au' 0.25 0.
graph/hplot/errors vx11 vy11 vxe11 vye11 11 20 0.17
graph 15 vx10m vy11m c
hi/zoom 12 ! 1 12
atitle '  ' '"L#y"G#'
text -5. 3.5 '(c) Al' 0.25 0.
text -4. 1. '- GHEISHA' 0.20 0.
text -4. 0.5 '"0# data' 0.20 0.
graph/hplot/errors vx10 vy12 vxe10 vye12 8 20 0.17
graph 15 vx10m vy12m c
hi/zoom 13 ! 1 12 
atitle '[D]y' '  '
text -5. 3.5 '(d) Au' 0.25 0.
text -5. 0.5 'curves"j# MCMHA' 0.20 0.
graph/hplot/errors vx11 vy13 vxe11 vye13 11 20 0.17
graph 15 vx10m vy13m c
exec hardcopy
*meta 0
*fort/close 60
*** End   of history file: 28/02/92 - 11.58.47
