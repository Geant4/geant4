set term postscript eps enhanced color "Helvetica" 46
set output "plot.eps"
set size 3.5,3.4

set style line 1 dt 1 pt 5 lc 0 lw 10.0 ps 3.0 #1.5
set style line 2 dt 2 pt 5 lc 0 lw 10.0 ps 3.0 #1.5
set style line 3 dt 3 pt 5 lc 0 lw 10.0 ps 3.0 #1.5
set style line 4 dt 4 pt 5 lc 0 lw 10.0 ps 3.0 #1.5
set style line 5 dt 4 pt 5 lc 3 lw 10.0 ps 3.0 #1.5

fileName1="Species.txt"

set logscale x
set format x "10^{%L}"
set xrange[1.01:0.999e6]

set xlabel "Time (ps)" font "Helvetica, 46" tc rgb "black" #white
set ylabel "G(Species/100 eV)" tc lt 0 font "Helvetica, 46" offset -0, -8

set multiplot
############### OH ##################
set size 1.4,1.7
set origin 0.0, 1.7
set rmargin 0
set lmargin 10
set bmargin 0
set tmargin 3

set format x ""
set xlabel ""
set yrange[2.51:6]
#^{\267}OH
plot fileName1 u 1:(stringcolumn(4) eq "OH^0" ? $2 : 1/0) w l ls 1 title "Default",\
     "data/OH.txt" index 0 u 1:2:3 with errorbars pt 4 ps 4 lc 0 notitle,\
     "data/OH.txt" index 1 u 1:2:3 with errorbars pt 5 ps 4 lc 0 notitle,\
     "data/OH.txt" index 2 u 1:2:3 with errorbars pt 6 ps 4 lc 0 notitle

############### e_aq #################
set size 1.4,1.7
set origin 1.4, 1.7
set lmargin 0
set rmargin 10
set bmargin 0
set tmargin 3

set format y ""
set ylabel ""

plot fileName1 u 1:(stringcolumn(4) eq "e_aq^-1" ? $2 : 1/0) w l ls 1 title "e-_{aq}",\
     "data/e_aq.txt" index 0 u 1:2:3 w errorbars pt 4 ps 4 lc 0 notitle,\
     "data/e_aq.txt" index 1 u 1:2:3 w errorbars pt 5 ps 4 lc 0 notitle,\
     "data/e_aq.txt" index 2 u 1:2:3 w errorbars pt 6 ps 4 lc 0 notitle,\
     "data/e_aq.txt" index 3 u 1:2:3 w errorbars pt 7 ps 4 lc 0 notitle
     

############### H3O #################

set size 1.4,1.7
set origin 2.415, 1.7
set lmargin 0
set rmargin 10
set bmargin 0
set tmargin 3

plot fileName1 u 1:(stringcolumn(4) eq "H3O^1" ? $2 : 1/0) w l ls 1 title "H_{3}O^{+}",\

############### H2O2 #################

set size 1.4,1.7
set origin 0.0, 0.0
set tmargin 0
set lmargin 10
set rmargin 0
unset bmargin

unset format y
set format x "10^{%L}"
set xlabel "Time (ps)" font "Helvetica, 46" tc rgb "black" #white

set yrange[0:0.99]
set ytics 0, 0.2, 1.

plot fileName1 u 1:(stringcolumn(4) eq "H2O2^0" ? $2 : 1/0) w l ls 1 title "H_{2}O_{2}",\

############### H2 #################

set size 1.4,1.7
set origin 1.4, 0.0
set tmargin 0
set lmargin 0
set rmargin 10
unset bmargin

set format y ""
set ylabel ""

plot fileName1 u 1:(stringcolumn(4) eq "H_2^0" ? $2 : 1/0) w l ls 1 title "H_{2}",\

############### H2 #################

set size 1.4,1.7
set origin 2.415, 0.0
set tmargin 0
set lmargin 0
set rmargin 10
unset bmargin

plot fileName1 u 1:(stringcolumn(4) eq "H^0" ? $2 : 1/0) w l ls 1 title "H^{\267}",\

set nomultiplot
