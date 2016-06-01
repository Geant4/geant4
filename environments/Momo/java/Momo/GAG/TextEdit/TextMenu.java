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

import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class TextMenu extends MenuBar implements ActionListener {
  // menus
  private Menu fileMenu;
  private Menu editMenu;
  private Menu fontMenu;
  private Menu sizeMenu;
  private Menu helpMenu;
	
  // adaptors
  private EditMan	editMan;
  private SizeMan	sizeMan;
  private FontMan	fontMan;

  // allocate standard components of TextEdit as part of initialization:
  private String newString	= new String("New");
  private String openString	= new String("Open");
  private String saveString	= new String("Save");
  private String saveAsString = new String("Save As");
  private String printString	= new String("Print");
  private String quitString	= new String("Quit");
  private String aboutString	= new String("About");

  // remember who our parents are

  private TextEdit textEdit;
  private TextCanvas textCanvas;

  // constructor puts together UI

  public TextMenu(TextEdit parent) {
    super();

    textEdit = parent;
    textCanvas = textEdit.getCanvas();

    CheckboxMenuItem cbi;
						
    editMan = new EditMan(parent);
    sizeMan = new SizeMan(parent);
    fontMan = new FontMan(parent);

    fileMenu = new Menu("File");
    fileMenu.add(new MenuItem(newString));
    fileMenu.add(new MenuItem(openString));
    fileMenu.add(new MenuItem(saveString));
    fileMenu.add(new MenuItem(saveAsString));
    fileMenu.addSeparator();
    fileMenu.add(new MenuItem(printString));
    //fileMenu.addSeparator();
    //fileMenu.add(new MenuItem(quitString));    
    fileMenu.addActionListener(this);        
    add(fileMenu);

    editMenu = new Menu("Edit");
    editMan.addItems(editMenu);
    add(editMenu);

    fontMenu = new Menu("Font");
    String[] fontNames = parent.getToolkit().getFontList();
    for (int i = 0; i < fontNames.length; i++) {
      if (fontNames[i].equals("Courier")){
	cbi = add_checked(fontMenu, fontMan, fontNames[i]);
	fontMan.remember(cbi);
      }else{
	add_unchecked(fontMenu, fontMan, fontNames[i]);
      }
    }
    add(fontMenu);

    sizeMenu = new Menu("Size");
    add_unchecked(sizeMenu, sizeMan, new String("10"));
    add_unchecked(sizeMenu, sizeMan, new String("11"));
    cbi = add_checked(sizeMenu, sizeMan, new String("12"));
    sizeMan.remember(cbi);
    add_unchecked(sizeMenu, sizeMan, new String("14"));
    add_unchecked(sizeMenu, sizeMan, new String("16"));
    add_unchecked(sizeMenu, sizeMan, new String("18"));
    add_unchecked(sizeMenu, sizeMan, new String("36"));
    add(sizeMenu);

    helpMenu = new Menu("Help");
    helpMenu.add(new MenuItem(aboutString));
    helpMenu.addActionListener(this);        
    setHelpMenu(helpMenu);
  }
	
  private void add_unchecked( Menu menu, ItemListener adaptor, String string ) {
    CheckboxMenuItem cbi;
		
    cbi = new CheckboxMenuItem(string);
    cbi.addItemListener(adaptor);
    menu.add(cbi);
  }

  private CheckboxMenuItem add_checked( Menu menu, ItemListener adaptor, String string ) {
    CheckboxMenuItem cbi;
    
    cbi = new CheckboxMenuItem(string);
    cbi.addItemListener(adaptor);
    cbi.setState(true);
    menu.add(cbi);
    return cbi;
  }

  // action handler, part of std java event handling

  public void actionPerformed(ActionEvent evt) {
    int i;

    Object source = evt.getSource();
    String cmd = evt.getActionCommand();

    if (source == fileMenu)
      doFile(cmd);
    else
      if (source == helpMenu) {
	AboutDialog ab = new AboutDialog(textEdit);
	ab.show();
      }
  }

  private void doFile(String cmd){
    if (cmd.equals(newString))
      textEdit.cmdFileNew();
    else
      if (cmd.equals(openString))
	textEdit.cmdFileOpen();
      else
	if (cmd.equals(saveString))
	  textEdit.cmdFileSave();
	else
	  if (cmd.equals(saveAsString))
	    textEdit.cmdFileSaveAs();
	  else
	    if (cmd.equals(printString))
	      textEdit.cmdPrint();
	    else
	      if (cmd.equals(quitString))
		textEdit.cmdFileQuit();
  }
}

class AboutDialog extends Dialog implements WindowListener, ActionListener {
  Button ok;
  public AboutDialog(Frame parent) {
    super(parent,"About TextEdit",true);

    Panel center = new Panel(new GridLayout(0,1));
    center.add(new Label("GAG v0.1alpha", Label.CENTER));
    center.add(new Label("Copyright (c) 1997 Toshiaki Kodama"));
    center.add(new Label("Contact: toshiaki@naruto-u.ac.jp"));
    center.add(new Label("This GAG includs TextEditi. (c) 1997 John Jensen"));
    add("Center", center);

    Panel p = new Panel();
    ok = new Button("Ok");
    ok.addActionListener(this);
    p.add(ok);
    add("South",p);
    //setSize(240,240);
    pack();
    addWindowListener(this);
  }

  // dispose on 'ok'
  public void actionPerformed(ActionEvent event) {
    if (event.getSource() == ok)
      dispose();
  }

  // add the 1.1 WindowListener stuff

  public void windowDeiconified(WindowEvent event) {}
  public void windowIconified(WindowEvent event) {}
  public void windowActivated(WindowEvent event) {}
  public void windowDeactivated(WindowEvent event) {}
  public void windowOpened(WindowEvent event) {}
  public void windowClosed(WindowEvent event) {}
  public void windowClosing(WindowEvent event) {
    dispose();
  }
}

class SizeMan extends Object implements ItemListener {
  // remember stuff

  private TextEdit textEdit;
  private CheckboxMenuItem lastSize = null;
	
  // constructor puts together UI

  public SizeMan(TextEdit parent) {
    textEdit = parent;
  }
	
  public void remember( CheckboxMenuItem item )	{
    if (lastSize != null)
      lastSize.setState(false);
    lastSize = item;
  }
	
  public void itemStateChanged(ItemEvent evt) {
    remember((CheckboxMenuItem)evt.getSource());
    String cmd = (String)evt.getItem();
    Integer size = new Integer(cmd);
    textEdit.setFontSize(size.intValue());
  }
}

final class FontMan extends Object implements ItemListener {
  // remember stuff

  private TextEdit textEdit;
  private CheckboxMenuItem lastTab = null;
	
  // constructor puts together UI

  public FontMan(TextEdit parent) {
    textEdit = parent;
  }
	
  public void remember( CheckboxMenuItem item )	{
    if (lastTab != null)
      lastTab.setState(false);
    lastTab = item;
  }

  public void itemStateChanged(ItemEvent evt) {
    String cmd = (String)evt.getItem();
    remember((CheckboxMenuItem)evt.getSource());
    textEdit.setFont(cmd);
  }
}
