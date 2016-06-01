##	GAG.tcl
##	GAG Window Manager
##
##	k.ohtubo(Tubocky)
##
##	1998.3.16

##      1998 July 5  Hajime
#         added "Exit GEANT4" and "Kill GEANT4" menus
#         yet to debug for invokation of non-GEANT4 executables
#                which hangs the pipe.
#         Non GEANT4 command which is suspended can be killed by the
#                new "Kill GEANT4" menu.
##      1998 July 20 Hajime
#         renamed Exec GEANT4 => Run GEANT4 
#         added Continue GEANT4
##	Tcl/Tk version 8.0
##      1998 July 5 GEANT4 Beta-01
##      1998 December 3 Beta-03 Added protection against G4 binaries with
#                               G4UIterminal session
#  
#       MOMOPATH/tcltk/Momo/GAG

wm title . GAG
wm protocol . WM_DELETE_WINDOW {
	.mbutt.function.command invoke last
}

##	Set global variable
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
if {[lsearch [array names env] FONT] >= 0} {
	set FONT $env(FONT)
} else {
	set FONT ""
}
set PROMPT "GAG> "
set Protocol "T1.0a"
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
if {![info exists env(G_PATH)]} {
	set env(G_PATH) $env(PWD)
}
set SEARCH Down
set REPLACE One
set HISTORY ""
set errorCode NONE
set errorInfo ""

##	Source
#source $HOME/GAG/GAGclearlog.proc
#source $HOME/GAG/GAGconnect.proc
#source $HOME/GAG/GAGdirectory.proc
#source $HOME/GAG/GAGparam.proc
#source $HOME/GAG/GAGsavelog.proc
#source $HOME/GAG/GAGsearchstring.proc
#source $HOME/Public.proc

##	tclIndex
set auto_path [linsert $auto_path 0 $HOME/]
set auto_path [linsert $auto_path 0 $HOME/GAG/]

##	Make window
. configure -highlightthickness 0

frame .mbutt -relief raised -bd 1 -highlightthickness 0
pack .mbutt -side top -fill x

menubutton .mbutt.function -text Function -relief raised -bd 3 \
	-menu .mbutt.function.command -highlightthickness 0
pack .mbutt.function -side left

menu .mbutt.function.command -tearoff 0
.mbutt.function.command add command -label "Run GEANT4" -command Directory
.mbutt.function.command add command -label "Continue" -command {
    puts $GEANT_ID continue
    .log.text insert end continue\n}
## GEANT4_ID is not yet checked 720

.mbutt.function.command add command -label "Exit GEANT4" -command {
	Message_Skeleton -icon question -button {Yes No} \
		-message "Do you exit GEANT4?"
	if {$FONT != ""} {
		Change_Font $FONT .msgskeleton
	}
	if {$F_COLOR != "" && $B_COLOR != ""} {
		Change_Color $F_COLOR $B_COLOR .msgskeleton
	}
	.msgskeleton.butt.0 configure -command {
		catch {puts $GEANT_ID exit}
		destroy .msgskeleton
##		exit
	}
	.msgskeleton.butt.1 configure -command {
		destroy .msgskeleton
	}
}
.mbutt.function.command add command -label "Kill GEANT4" -command {
	Message_Skeleton -icon question -button {Yes No} \
		-message "Do you really kill GEANT4?"
	if {$FONT != ""} {
		Change_Font $FONT .msgskeleton
	}
	if {$F_COLOR != "" && $B_COLOR != ""} {
		Change_Color $F_COLOR $B_COLOR .msgskeleton
	}
	.msgskeleton.butt.0 configure -command {
##test
                set PID [pid $GEANT_ID]
		catch {exec kill $PID}
		destroy .msgskeleton
##		exit
	}
	.msgskeleton.butt.1 configure -command {
		destroy .msgskeleton
	}
}

.mbutt.function.command add separator

.mbutt.function.command add command -label "Command History" -command Command_History
.mbutt.function.command add separator

.mbutt.function.command add command -label "Save Log" -command Before_Save
.mbutt.function.command add command -label "Search Log" \
	-command Search_String
.mbutt.function.command add command -label "Clear Log" -command Clear_Log

.mbutt.function.command add separator
.mbutt.function.command add command -label "Exit GAG" -command {
	Message_Skeleton -icon question -button {Yes No} \
		-message "Do you exit GAG?"
	if {$FONT != ""} {
		Change_Font $FONT .msgskeleton
	}
	if {$F_COLOR != "" && $B_COLOR != ""} {
		Change_Color $F_COLOR $B_COLOR .msgskeleton
	}
	.msgskeleton.butt.0 configure -command {
		catch {puts $GEANT_ID exit}
		catch {close $GEANT_ID}
		destroy .msgskeleton
		exit
	}
	.msgskeleton.butt.1 configure -command {
		destroy .msgskeleton
	}
}

frame .space0 -height 15 -relief flat -highlightthickness 0
pack .space0 -side top -fill x


frame .log -highlightthickness 0
pack .log -side top -fill both -expand 1

text .log.text -width 85 -height 15 -bd 2 -relief raised \
	-yscrollcommand {.log.scroll set} -highlightthickness 0
pack .log.text -side left -fill both -expand 1

scrollbar .log.scroll -command {.log.text yview} -highlightthickness 0
pack .log.scroll -side left -fill y


frame .space1 -height 5 -relief flat -highlightthickness 0
pack .space1 -side top -fill x


frame .comm -highlightthickness 0
pack .comm -side top -fill x

label .comm.label -text Command -highlightthickness 0
pack .comm.label -side left

entry .comm.ent -highlightthickness 0
pack .comm.ent -side left -fill x -expand 1


label .help -relief sunken -anchor w -highlightthickness 0
pack .help -side top -fill x

##	Dummy Button
button .dummy -command {
	set PARA [.comm.ent get]
	if {$PARA != ""} {
		GAG_Connect $PARA
	} else {
		.log.text insert end \n>
		.log.text see end
	}
}

##	Binding & Focus
#	Binding on command line
bind .comm.ent <Return> {.dummy invoke}
bind .comm.ent <Control-c> {}
bind .comm.ent <Enter>	{.help configure -text "direct command typein"}
bind .comm.ent <Leave>	{.help configure -text "setup phase"}
focus .comm.ent

#	Binding on Log board
bind .log <Enter>	{
	.help configure -text "editable log window"
}
bind .log <Leave>	{
	.help configure -text "setup phase"
	.log.text configure -relief raised
	focus .comm.ent
}

bind .log.text <Button> {
	.log.text configure -relief sunken
	focus .log.text
}

if {$FONT != ""} {
	Change_Font $FONT .
}
if {$F_COLOR != "" && $B_COLOR != ""} {
	Change_Color $F_COLOR $B_COLOR .
}
Win_Size .
Tab_off
Del_Bind
Control

#	Out put on Logboard
.log.text insert end "GEANT4 Adaptive GUI (GAG) : 1998 December \nGAG  protocol version $Protocol.\n####### GOOD LUCK TO YOUR  SIMULATION! #######\n> "

#bind Menubutton <Leave> {
#	foreach w [winfo children .mbutt] {
#		grab release $w
#	}
#}
#bind Menubutton <Motion> {}





