#-------------------------------------------------------------------
# Last update: 20-Jan-2006
#
# This is a bash-shell script to execute a bunch of power jobs,
# each with a particular combination of  alpha , N  parameters.
# It uses the Python script:  replaceParameters.py .
#
#-------------------------------------------------------------------
#
echo "--- START jobs.sh --- "
date
#
. setup.sh
#
# --- 1000 ---
#
python replaceParameters.py 0.0 1000
gmake
echo " running: 0.0-1000 "
time power > output.log-0.0-1000 2>&1
echo " done"
#
python replaceParameters.py 0.1 1000
gmake
echo " running: 0.1-1000 "
time power > output.log-0.1-1000 2>&1
echo " done"
#
python replaceParameters.py 0.2 1000
gmake
echo " running: 0.2-1000 "
time power > output.log-0.2-1000 2>&1
echo " done"
#
python replaceParameters.py 0.3 1000
gmake
echo " running: 0.3-1000 "
time power > output.log-0.3-1000 2>&1
echo " done"
#
python replaceParameters.py 0.4 1000
gmake
echo " running: 0.4-1000 "
time power output.log-0.4-1000 2>&1
echo " done"
#
python replaceParameters.py 0.5 1000
gmake
echo " running: 0.5-1000 "
time power output.log-0.5-1000 2>&1
echo " done"
#
python replaceParameters.py 0.6 1000
gmake
echo " running: 0.6-1000 "
time power > output.log-0.6-1000 2>&1
echo " done"
#
python replaceParameters.py 0.7 1000
gmake
echo " running: 0.7-1000 "
time power > output.log-0.7-1000 2>&1
echo " done"
#
python replaceParameters.py 0.8 1000
gmake
echo " running: 0.8-1000 "
time power > output.log-0.8-1000 2>&1
echo " done"
#
python replaceParameters.py 0.9 1000
gmake
echo " running: 0.9-1000 "
time power > output.log-0.9-1000 2>&1
echo " done"
#
python replaceParameters.py 1.0 1000
gmake
echo " running: 1.0-1000 "
time power > output.log-1.0-1000 2>&1
echo " done"
#
#
# --- 2000 ---
#
python replaceParameters.py 0.0 2000
gmake
echo " running: 0.0-2000 "
time power > output.log-0.0-2000 2>&1
echo " done"
#
python replaceParameters.py 0.1 2000
gmake
echo " running: 0.1-2000 "
time power > output.log-0.1-2000 2>&1
echo " done"
#
python replaceParameters.py 0.2 2000
gmake
echo " running: 0.2-2000 "
time power > output.log-0.2-2000 2>&1
echo " done"
#
python replaceParameters.py 0.3 2000
gmake
echo " running: 0.3-2000 "
time power > output.log-0.3-2000 2>&1
echo " done"
#
python replaceParameters.py 0.4 2000
gmake
echo " running: 0.4-2000 "
time power > output.log-0.4-2000 2>&1
echo " done"
#
python replaceParameters.py 0.5 2000
gmake
echo " running: 0.5-2000 "
time power > output.log-0.5-2000 2>&1
echo " done"
#
python replaceParameters.py 0.6 2000
gmake
echo " running: 0.6-2000 "
time power > output.log-0.6-2000 2>&1
echo " done"
#
python replaceParameters.py 0.7 2000
gmake
echo " running: 0.7-2000 "
time power > output.log-0.7-2000 2>&1
echo " done"
#
python replaceParameters.py 0.8 2000
gmake
echo " running: 0.8-2000 "
time power > output.log-0.8-2000 2>&1
echo " done"
#
python replaceParameters.py 0.9 2000
gmake
echo " running: 0.9-2000 "
time power > output.log-0.9-2000 2>&1
echo " done"
#
python replaceParameters.py 1.0 2000
gmake
echo " running: 1.0-2000 "
time power > output.log-1.0-2000 2>&1
echo " done"
#
#
# --- 3000 ---
#
python replaceParameters.py 0.0 3000
gmake
echo " running: 0.0-3000 "
time power > output.log-0.0-3000 2>&1
echo " done"
#
python replaceParameters.py 0.1 3000
gmake
echo " running: 0.1-3000 "
time power > output.log-0.1-3000 2>&1
echo " done"
#
python replaceParameters.py 0.2 3000
gmake
echo " running: 0.2-3000 "
time power > output.log-0.2-3000 2>&1
echo " done"
#
python replaceParameters.py 0.3 3000
gmake
echo " running: 0.3-3000 "
time power > output.log-0.3-3000 2>&1
echo " done"
#
python replaceParameters.py 0.4 3000
gmake
echo " running: 0.4-3000 "
time power > output.log-0.4-3000 2>&1
echo " done"
#
python replaceParameters.py 0.5 3000
gmake
echo " running: 0.5-3000 "
time power > output.log-0.5-3000 2>&1
echo " done"
#
python replaceParameters.py 0.6 3000
gmake
echo " running: 0.6-3000 "
time power > output.log-0.6-3000 2>&1
echo " done"
#
python replaceParameters.py 0.7 3000
gmake
echo " running: 0.7-3000 "
time power > output.log-0.7-3000 2>&1
echo " done"
#
python replaceParameters.py 0.8 3000
gmake
echo " running: 0.8-3000 "
time power > output.log-0.8-3000 2>&1
echo " done"
#
python replaceParameters.py 0.9 3000
gmake
echo " running: 0.9-3000 "
time power > output.log-0.9-3000 2>&1
echo " done"
#
python replaceParameters.py 1.0 3000
gmake
echo " running: 1.0-3000 "
time power > output.log-1.0-3000 2>&1
echo " done"
#
#
echo "--- END jobs.sh --- "
date
#
