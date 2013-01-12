echo 'Start run script'
cd ../../utils
gmake
echo 'utils was compiled'
cd ../standard
gmake
echo 'standard was compiled'
cd ../lowenergy
gmake
echo 'lowenergy was compiled'
cd test
gmake clean
gmake
echo 'executable is built'
