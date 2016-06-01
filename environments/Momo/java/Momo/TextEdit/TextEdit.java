/*
 * Copyright (c) 1997 John Jensen. All rights reserved.
 *
 * This software is FREE FOR COMMERCIAL AND NON-COMMERCIAL USE,
 * provided the following condition is met.
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose and without fee is hereby granted,
 * provided that any copy or derivative of this software or documentation
 * retaining the name "John Jensen" also retains this condition and the
 * following disclaimer.
 *
 * CopyrightVersion 1.0
 */

package TextEdit;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.awt.event.*;

import GAG.*;

public class TextEdit extends Frame implements GAGtextEditor, NoteListener, WindowListener
{
  // The NoteListen member fuction of NoteListener is overridded to support
  // the following actions.

  final private int NOACTION	= 0;
  final private int NEWFILE	= 1;
  final private int OPENFILE	= 2;
  final private int QUITAPP	= 3;

  // allocate standard components of TextEdit as part of initialization:

  public File file			= new File("GEANT4 LOG");

  private String openPrompt	= new String("Open text file");
  private String closePrompt	= new String("Save file as");

  // initialize some state flags

  private boolean neverNamed	= true;
  private boolean menuBusy	= false;

  // these three .. no four! .. four objects form the main window

  private Scrollbar	horiz;
  private Scrollbar	vert;
  private TextCanvas	textCanvas;
  private TextMenu	textMenu;

  // font

  private String		fontString = new String("Serif");
  private int			fontSize = 14;
  private int			tabSize = 4;

  // constructor puts together UI

  public TextEdit() {
    this.setTitle();

    vert = new Scrollbar(Scrollbar.VERTICAL);
    add("East", vert);
    horiz = new Scrollbar(Scrollbar.HORIZONTAL);
    add("South", horiz);

    textCanvas = new TextCanvas(this,horiz,vert);
    add("Center", textCanvas);
    textCanvas.setFont(fontString,fontSize);
    addKeyListener(textCanvas);

    textMenu = new TextMenu(this);
    setMenuBar(textMenu);

    addWindowListener(this);

    setSize(600, 600);
    setVisible(true);
    try { Thread.sleep(400); } catch(InterruptedException e) {}
  }
  public void writeString(String s){
    textCanvas.writeString(s);
  }
  public void writeLine(String s){
    textCanvas.writeLine(s);
  }

  // allow others to go straight to the canvas
  public TextCanvas getCanvas() {
    return textCanvas;
  }

  // add current filename to window title

  private void setTitle(){
    super.setTitle("TextEdit - " + file.getName());
  }
	
  // add the 1.1 WindowListener stuff

  public void windowDeiconified(WindowEvent event) {}
  public void windowIconified(WindowEvent event) {}
  public void windowActivated(WindowEvent event) {}
  public void windowDeactivated(WindowEvent event) {}
  public void windowOpened(WindowEvent event) {}
  public void windowClosed(WindowEvent event) {}

  public void windowClosing(WindowEvent event) {
    if (isItSafe(QUITAPP))
      System.exit(0);
  }

  // high level command processors (invoked my menu)

  public void setFont(String fs) {
    fontString = fs;
    textCanvas.setFont(fontString, fontSize);
  }

  public void setFontSize(int s) {
    fontSize = s;
    textCanvas.setFont(fontString, fontSize);
  }

  public void setTabSize(int s) {
    tabSize = s;
    textCanvas.setTab(tabSize);
  }

  public int getTabSize() {
    return tabSize;
  }

  public void cmdFileNew() {
    if (isItSafe(NEWFILE))
      newFile();
  }
	
  public void cmdFileOpen() {
    if (isItSafe(OPENFILE))
      openFile();
  }

  public void cmdFileSave() {
    saveFile(false);
  }

  public void cmdFileSaveAs() {
    saveFile(true);
  }

  public void cmdPrint() {
    textCanvas.print(file.getName());
  }

  public void cmdFileQuit() {
    if (isItSafe(QUITAPP))
      System.exit(0);
  }

  // if the file has been modified (dirty), fire off a notify dialog
  // and pass it the "action" you want returned to noteListen.  This
  // dataflow beats the "modal dialog problem" on Windows

  // I believe isItSafe is what Laurence Olivier said repeatedly to
  // Dustan Hoffman in Marathon Man.
  
  private boolean isItSafe(int action) {
    if (textCanvas.isDirty())
      {
	notify( "The current file has not been saved. ", action );
	return false;
      }
    else
      return true;
  }

  // implementation of the NoteListener interface - it does whatever
  // action the user agreed to by pressing "Continue" on a modal dialog

  public void noteListen(int action)
  {
    switch (action) 
      {
      case NOACTION:
	break;
      case NEWFILE:
	newFile();
	break;
      case OPENFILE:
	openFile();
	break;
      case QUITAPP:
	System.exit(0);
	break;
      }
  }

  // clear out the document

  private void newFile()
  {
    neverNamed = true;
    file = new File("Untitled");
    this.setTitle();
    textCanvas.clear();
    textCanvas.validate();
  }

  // call getFile to get a file path, and if you get one open it
  
  private void openFile()
  {
    if (getFile(FileDialog.LOAD))
      openTheFile();
  }

  public void openTheFile()
  {
    textCanvas.read(file);
    textCanvas.validate();
    neverNamed = false;
    this.setTitle();
  }

  // call getFile to get a file path, and if you get one save to it

  private void saveFile(boolean As)
  {
    boolean ok;

    if (neverNamed || As)
      {
	ok = getFile(FileDialog.SAVE);
	neverNamed = false;
      }
    else
      ok = true;
    
    if (ok)
      {
	textCanvas.write(file);
	neverNamed = false;
	this.setTitle();
      }
  }

  // use FileDialog to get a filepath

  private boolean getFile(int mode)
  {
    String prompt;
    String filename;
    String pathname;

    if (mode == FileDialog.LOAD)
      prompt = openPrompt;
    else
      prompt = closePrompt;
    
    FileDialog d = new FileDialog(this, prompt, mode);

    if (mode == FileDialog.LOAD)
      d.setFile("*");
    else
      d.setFile(file.getName());

    d.setDirectory(".");
    d.setVisible(true);

    filename = d.getFile(); 

    pathname = d.getDirectory() + filename; 

    d.dispose();

    if (filename != null)
      {
	file = new File(pathname);
	return true;
      }
    else
      return false;
  }

  // launch a note

  private void notify(String s, int action) {
    Note note = new Note(this,s,action);
    note.setVisible(true);
  }

  // main declares the TextEdit frame which does everything

  public static void main(String[] args) {
    TextEdit f = new TextEdit();
    if (args.length > 0) {
      f.file = new File(args[0]);
      if (f.file.isFile())
	f.openTheFile();
      else
	f.notify("File not found: "+args[0],0);
    }
  }
}
