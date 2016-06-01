##	tMomo.tcl
##      This can invoke Tcl/Tk GAG, File borwser and Extensions.
##      This replaces Momo.tcl, an obsolete Momo
##
##	k.ohtubo(Tubocky)
##
##	1998.3.16
##	Tcl/Tk version 8.0
##      1998 July 5 GEANT4 Beta-01
##      1998 July 20 added "New Shell" and "continue" commands
##	Set global variable
#          MOMOPATH/tcltk/Momo/GAG hierarchy is required

if [info exists env(MOMOPATH)] {
	set END [expr [string length $env(MOMOPATH)] - 1]
	if {[string index $env(MOMOPATH) $END] == "/"} {
		set HOME [string range $env(MOMOPATH) 0 [expr $END - 1]]
	} else {
		set HOME $env(MOMOPATH)
	}
} else {
	set HOME $env(HOME)/Momo/tcltk
}
set WISH wish$tk_version
set env(G_PATH) $env(PWD)

set Env_Name [list FONT GUI_STYLE GUI_POSITION FORE_GROUND_COLOR BACK_GROUND_COLOR]

##     init files for Momo settings (font and color) and extension
set SAVE_FILE .Momorc
set ENV_FILE .Momoext

set errorCode NONE
set errorInfo ""

##	tclIndex
set auto_path [linsert $auto_path 0 $HOME/]

if [file exists $env(HOME)/$SAVE_FILE] {
	set file_ID [open $env(HOME)/$SAVE_FILE r]
	while {![eof $file_ID]} {
		set STRING [gets $file_ID]
		while {[string index $STRING 0] == " " || [string index $STRING 0 ] == "\t"} {
			set STRING [string range $STRING 1 end]
		}
		set LAST [string length $STRING]
		incr LAST -1
		while {[string index $STRING $LAST] == " " || [string index $STRING 0 ] == "\t"} {
			incr LAST -1
			set STRING [string range $STRING 0 $LAST]
		}
		set LAST [string first " " $STRING]
		set FIRST [expr $LAST + 1]
		incr LAST -1
		set env([string range $STRING 0 $LAST]) [string range $STRING $FIRST end]
	}
	unset env([string range $STRING 0 $LAST])
	close $file_ID
}

if {[lsearch [array names env] GUI_STYLE] < 0} {
	set env(GUI_STYLE) horizontal
}
if {[lsearch [array names env] GUI_POSITION] < 0} {
	set env(GUI_POSITION) left
}
if {[lsearch [array names env] FORE_GROUND_COLOR] < 0} {
	set env(FORE_GROUND_COLOR) #000000000
	set FC env(FORE_GROUND_COLOR)
	set FR 0
	set FG 0
	set FB 0
} else {
	set FC $env(FORE_GROUND_COLOR)
	if {[string index $FC 0] == "#"} {
		Reset_Color F
	} else {
		set FR 0
		set FG 0
		set FB 0
	}
}
if {[lsearch [array names env] BACK_GROUND_COLOR] < 0} {
	set env(BACK_GROUND_COLOR) #d90d90d90
	set BC env(BACK_GROUND_COLOR)
	set BR 0.848
	set BG 0.848
	set BB 0.848
} else {
	set BC $env(BACK_GROUND_COLOR)
	if {[string index $BC 0] == "#"} {
		Reset_Color B
	} else {
		set BR 0.848
		set BG 0.848
		set BB 0.848
	}
}


##	Make window
. configure -highlightthickness 0

menubutton .function -text Function -menu .function.menu -bd 3 \
	-relief raised -highlightthickness 0

menu .function.menu -tearoff 0
.function.menu add cascade -label Position \
	-menu .function.menu.position
.function.menu add cascade -label Style \
	-menu .function.menu.style
.function.menu add separator
.function.menu add command -label Exit -command {
	set file_ID [open $env(HOME)/$SAVE_FILE w]
	foreach VARIABLE $Env_Name {
		if [info exists env($VARIABLE)] {
			puts $file_ID "$VARIABLE $env($VARIABLE)"
		}
	}
	close $file_ID
	exit
}

menu .function.menu.position -tearoff 0
.function.menu.position add radiobutton -label Left \
	-variable env(GUI_POSITION) -value left -command Position_Pack
.function.menu.position add radiobutton -label Right \
	-variable env(GUI_POSITION) -value right -command Position_Pack

menu .function.menu.style -tearoff 0
.function.menu.style add radiobutton -label Horizontal \
	-variable env(GUI_STYLE) -value horizontal -command Style_Pack
.function.menu.style add radiobutton -label Vertical \
	-variable env(GUI_STYLE) -value vertical -command Style_Pack


menubutton .env -text Environment -menu .env.menu -bd 3 \
	-relief raised -highlightthickness 0

menu .env.menu -tearoff 0
.env.menu add command -label Color -command Define_Color
.env.menu add command -label Font -command Select_Font
.env.menu add separator
.env.menu add command -label "Env Variables" -command List_Env
#.env.menu entryconfigure Font -state disabled

button .gag -text "GAG" -highlightthickness 0 \
	-command {
		raise .
		exec $WISH $HOME/GAG/GAG.tcl &
	} -bd 3
#menubutton .gag -text "GAG(GEANT4 Adaptive GUI)" -bd 3 -menu .gag.menu \
#	-relief raised -highlightthickness 0

#menu .gag.menu -tearoff 0
#.gag.menu add command -label "Tcl/Tk(8.0) version" \
#	-command {
#		raise .
#		exec $WISH $HOME/GAG/GAG.tcl &
#	}
###### JAVA GAG removed
#.gag.menu add command -label "Java version" \
#	-command {
#		raise .
#		exec $JAVA $HOME/GAG/GAG &
#	}
##### GGE removed for Beta-01
#menubutton .gge -text "GGE(GEANT4 Geometry Editor)" -bd 3 \
#	-relief raised -menu .gge.menu -highlightthickness 0

#menu .gge.menu -tearoff 0
#.gge.menu add command -label "Tcl/Tk(8.0) version" \
#	-command {
#		raise .
#		exec $WISH $HOME/GGE/GGE.tcl &
#	}
#.gge.menu add command -label "Java version" \
#	-command {
#		raise .
#		exec $JAVA $HOME/GGE/GAG &
#	}


button .fb -text "File Browser" -bd 3 -highlightthickness 0 \
	-command {
		raise .
		FBmain
}

# added 1998 July 20

button .shel -text "New Shell" -bd 3 -highlightthickness 0 \
	-command {
		raise .
		exec xterm &
}



if [file exists $env(HOME)/$ENV_FILE] {
	menubutton .ext -text Plugins -bd 3 -menu .ext.menu \
		-relief raised -highlightthickness 0
	menu .ext.menu -tearoff 0

	set file_ID [open $env(HOME)/$ENV_FILE r]

	while {![eof $file_ID]} {
		set STRING [gets $file_ID]
		while {[string index $STRING 0] == " " || [string index $STRING 0 ] == "\t"} {
			set STRING [string range $STRING 1 end]
		}
		set LAST [string length $STRING]
		incr LAST -1
		while {[string index $STRING $LAST] == " " || [string index $STRING 0 ] == "\t"} {
			incr LAST -1
			set STRING [string range $STRING 0 $LAST]
		}
		if {$STRING == "Momoseparator"} {
			.ext.menu add separator
		} elseif {$STRING != "" && [string index $STRING 0] != "#"} {
			.ext.menu add command -label $STRING -command "exec $STRING &"
		}
	}
	close $file_ID

	bind .ext <ButtonRelease> {
		raise .
		raise .ext.menu
	}
}


##	Procedures
#	Change Style
proc Style_Pack {} {
	global env

	set WIDGET [winfo children .]
	pack forget .
	switch $env(GUI_STYLE) {
		horizontal {eval pack $WIDGET -side left -fill both
			.gag configure -text "GAG"
#			.gge configure -text "GGE(GEANT4 Geometry Editor)"
		}
		vertical {eval pack $WIDGET -side top -fill both
			.gag configure -text "GAG\n"
#			.gge configure -text "GGE\n(GEANT4 Geometry Editor)"
		}
	}
}

#	Change Position
proc Position_Pack {} {
	global env

	switch $env(GUI_POSITION) {
		left {wm geometry . +0+0}
		right {wm geometry . -0+0}
	}
}

wm deiconify .
wm geometry . +0+0
wm overrideredirect . 1
wm resizable . 0 0

#	Binding
bind .function	<ButtonRelease>	{
	raise .
	raise .function.menu
}
bind .env <ButtonRelease> {
	raise .
	raise .env.menu
}
#bind .gag <ButtonRelease> {
#	raise .
#	raise .gag.menu
#}
#bind .gge <ButtonRelease> {
#	raise .
#	raise .gge.menu
#}
bind .fb <ButtonRelease> {
	raise .
}

if {[lsearch [array names env] FONT] >= 0} {
	Change_Font $env(FONT) .
}

Change_Color $env(FORE_GROUND_COLOR) $env(BACK_GROUND_COLOR) .
Tab_off

Style_Pack
Position_Pack
Del_Bind
Control

#bind Menubutton <Leave> {
#	foreach w [winfo children .] {
#		grab release $w
#	}
#}
#bind Menubutton <Motion> {}
