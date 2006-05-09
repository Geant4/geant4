echo "Start of run for all test of HARP"

setenv TARGET p_be_9gev
source run.csh
setenv TARGET p_al_13gev
source run.csh

echo "Test is done!" 
