#!/bin/tcsh


echo " ######### This is a Gflash Test ! date: ######### " >test.log
date >> test.log
echo " ################################################# " >>test.log
../../bin/Linux-g++/systest test.mac  >>  test.log
less test.log | grep "####"  > test2.log
echo "ERRORS  :"
diff test2.log ref.log  
 rm test2.log
 rm test.log
exit
