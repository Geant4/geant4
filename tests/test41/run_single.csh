#/bin/csh

setenv HISTFILENAME li_1_${1}.log
${2} ${3}li_1.in 

setenv HISTFILENAME li_2_${1}.log
${2} ${3}li_2.in 

setenv HISTFILENAME be_1_${1}.log
${2} ${3}be_1.in 

setenv HISTFILENAME be_2_${1}.log
${2} ${3}be_2.in 

setenv HISTFILENAME ch2_${1}.log
${2} ${3}ch2.in 

setenv HISTFILENAME c_${1}.log
${2} ${3}c.in 

setenv HISTFILENAME al_${1}.log
${2} ${3}al.in 

setenv HISTFILENAME fe_${1}.log
${2} ${3}fe.in 

setenv HISTFILENAME h2_1_${1}.log
${2} ${3}h2_1.in 

setenv HISTFILENAME h2_2_${1}.log
${2} ${3}h2_2.in 
