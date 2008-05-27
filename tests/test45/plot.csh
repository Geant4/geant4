#!/bin/csh -f
#trace on
#---------------------------------------------------------------
#
#  Author V.Ivanchenko 21 May 2005
#  test30
#
#----------------------------------------------------------------

cd $TARGET

rm -f exfor.root root.out *.gif *.png

echo "Start ROOT" > root.out

#if(-r $TEST45/$TARGET/exfor.root) then
#    ln -s $TEST45/$TARGET/exfor.root 
#endif

if(-r $TEST45/$TARGET/Plot.C) then
    echo 'make plots for ' $TARGET
    root  -b -q $TEST45/$TARGET/Plot.C >>& root.out
endif

#if(-r $TEST45/$TARGET/PlotLog.C) then
#    echo 'make plots for ' $TARGET
#    root  -b -q $TEST30/$TARGET/PlotLog.C >>& root.out
#endif

cd ../
