##	FBtext.tcl
##	Window Manager
##       called by FBmain.proc
##
##	k.ohtubo(Tubocky)
##
##	1998.3.16
##	Tcl/Tk version 8.0
##      1998 July 5 GEANT4 Beta-01

wm title . "File Browser Text"

##	Set global variable
if [info exists env(MOMOPATH)] {
	set END [expr [string length $env(MOMOPATH)] - 1]
	if {[string index $env(MOMOPATH) $END] == "/"} {
		set HOME [strinf range $env(MOMOPATH) 0 [expr $END - 1]]
	} else {
		set HOME $env(MOMOPATH)
	}
} else {
	set HOME $env(HOME)/Momo/tcltk
}
set FILE $env(FB_Open_File_Name)
if {[lsearch [array names env] FONT] >= 0} {
	set FONT $env(FONT)
} else {
	set FONT ""
}
if {[lsearch [array names env] FORE_GROUND_COLOR] >= 0} {
	set F_COLOR $env(FORE_GROUND_COLOR)
} else {
	set F_COLOR #000000000
}
if {[lsearch [array names env] BACK_GROUND_COLOR] >= 0} {
	set B_COLOR $env(BACK_GROUND_COLOR)
} else {
	set B_COLOR #d90d90d90
}
set SEARCH Down
set REPLACE One
set errorCode NONE
set errorInfo ""


##	tclIndex
set auto_path [linsert $auto_path 0 $HOME/]
set auto_path [linsert $auto_path 0 $HOME/FBtext/]


##	Source
#source $HOME/FBtext/FBsearchstring.proc
#source $HOME/FBtext/FBclear.proc
#source $HOME/FBtext/FBsavelog.proc


##	Make window
frame .butt -highlightthickness 0
pack .butt -side top -fill x

menubutton .butt.func -text Function -bd 3 -menu .butt.func.menu \
	-relief raised -highlightthickness 0
pack .butt.func -side left -fill y

menu .butt.func.menu -tearoff 0
.butt.func.menu add command -label Clear -command Clear_Log
.butt.func.menu add command -label Save -command Before_Save
.butt.func.menu add command -label Search -command Search_String
.butt.func.menu add separator
.butt.func.menu add command -label "Exit File Browser Text" -command exit


frame .dir -highlightthickness 0
pack .dir -side top -fill x

label .dir.lbl -text Directory: -anchor w -highlightthickness 0
pack .dir.lbl -side left

set END [string length $env(HOME)]
set Path ~[string range [pwd] $END end]
label .dir.lbr -text $Path -anchor w -highlightthickness 0
pack .dir.lbr -side left -fill x -expand 1


frame .file -highlightthickness 0
pack .file -side top -fill x


label .file.lb -text "Load File:" -anchor w -highlightthickness 0
pack .file.lb -side left

entry .file.ent -highlightthickness 0
pack .file.ent -side left -fill x -expand 1


frame .fbt -highlightthickness 0
pack .fbt -side top -fill both -expand 1

text .fbt.text -width 80 -height 15 -bd 2 -yscrollcommand {.fbt.scroll set} \
	-highlightthickness 0
pack .fbt.text -side left -fill both -expand 1

scrollbar .fbt.scroll -command {.fbt.text yview} -highlightthickness 0
pack .fbt.scroll -side left -fill y


frame .info -highlightthickness 0
pack .info -side top -fill x


frame .info.lbl -highlightthickness 0
pack .info.lbl -side left

label .info.lbl.info -text Information -anchor w \
	-highlightthickness 0
pack .info.lbl.info -side top -fill x


frame .info.lbr -highlightthickness 0
pack .info.lbr -side left -fill x -expand 1

label .info.lbr.info -relief sunken -anchor w -highlightthickness 0
pack .info.lbr.info -side top -fill x -expand 1


if {$FONT != ""} {
	Change_Font $FONT .
}
if {$F_COLOR != "" && $B_COLOR != ""} {
	Change_Color $F_COLOR $B_COLOR .
}

Tab_off
Win_Size .
Del_Bind
Control


if {[file readable $FILE] && [file writable $FILE]} {
	.info.lbr.info configure -text "This file is readable & writable."
	set FILE_ID [open $FILE r]
	.fbt.text insert end [read $FILE_ID]
	close $FILE_ID
} elseif [file readable $FILE] {
	.info.lbr.info configure -text "This file is readable."
	set FILE_ID [open $FILE r]
	.fbt.text insert end [read $FILE_ID]

	.fbt.text configure -state disabled
	.butt.func.menu entryconfigure 0 -state disabled
	.butt.func.menu entryconfigure 1 -state disabled

	close $FILE_ID
} else {
	.info.lbr.info configure -text "This file is not readable."
}

.file.ent insert end $FILE
.file.ent configure -state disabled
