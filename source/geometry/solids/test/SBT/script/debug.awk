#!/bin/awk -f
$0 ~ /#/   { print }
$0 !~ /#/  { print "#" $0 ; print }
