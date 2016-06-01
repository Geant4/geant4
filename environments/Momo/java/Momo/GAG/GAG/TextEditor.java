// jdk-1.1.3 (or later) + swing-1.0.2
// JFileChooser is to be used soon.
// 
package GAG;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.net.URL;
import java.util.*;

import com.sun.java.swing.text.*;
import com.sun.java.swing.*;

class TextEditor extends JFrame implements ActionListener {
  private JMenuItem newf, open, save, cut, copy, paste;
  private JTextArea editor;
  private JToolBar toolbar;
  protected FileDialog fileDialog;

  TextEditor() {
    super("GAG TextEditor");
    getContentPane().setLayout(new BorderLayout());
    editor = new JTextArea();
    editor.setFont(new Font("Serif", Font.BOLD, 12));
    JScrollPane scroller = new JScrollPane();
    //scroller.setVerticalScrollBarPolicy(JScrollPane.SCROLLBARS_ALWAYS);
    JViewport port = scroller.getViewport();
    port.add(editor);
    port.setBackingStoreEnabled(false);
    setJMenuBar( createMenubar() );
    getContentPane().add("Center", scroller);
    pack();
    setSize(500, 600);
    show();
  }
  JMenuBar createMenubar(){
    JMenu file, edit;
    JMenuBar mb = new JMenuBar();
    file = new JMenu("File");
    file.add(newf = new JMenuItem("New"));
    file.add(open = new JMenuItem("Open"));
    file.add(save = new JMenuItem("Save"));
    edit = new JMenu("Edit"); 
    edit.add(cut = new JMenuItem("Cut"));
    edit.add(copy = new JMenuItem("Copy"));
    edit.add(paste = new JMenuItem("Paste"));
    mb.add(file);
    mb.add(edit);
    newf.addActionListener(this);
    open.addActionListener(this);
    save.addActionListener(this);
    cut.addActionListener(this);
    copy.addActionListener(this);
    paste.addActionListener(this);
    return mb;
  }
  void writeString(String s){
    editor.append(s);
  }
  void writeLine(String s){
    editor.append(s+"\n");
  }

  public static void main(String[] args) {
    new TextEditor();
  }
  public void actionPerformed(ActionEvent e) {
    Object o = e.getSource();
    if (o == newf){newAction(); return;}
    if (o == open){openAction(); return;}
    if (o == save){saveAction(); return;}
    if (o == cut){editor.cut(); return;}
    if (o == copy){editor.copy(); return;}
    if (o == paste){editor.paste(); return;}
  }
  void newAction(){
    editor.setDocument(new PlainDocument());
    validate();
  }
  void openAction(){
    if (fileDialog == null) {
      fileDialog = new FileDialog(this);
    }
    fileDialog.setMode(FileDialog.LOAD);
    fileDialog.show();
    
    String file = fileDialog.getFile();
    if (file == null) {
      return;
    }
    newAction();
    String directory = fileDialog.getDirectory();
    File f = new File(directory, file);
    StringBuffer sb = new StringBuffer();
    if (f.exists()) {
      setTitle(file);
      try{
	Reader in = new FileReader(f);
	char[] buff = new char[4096];
	int nch;
	while ((nch = in.read(buff, 0, buff.length)) != -1) {
	  sb.append(new String(buff, 0, nch));
	}
      }catch(IOException e){ System.err.println(e.getMessage()); }
      editor.setText(sb.toString());
      validate();
    }
  }
  void saveAction(){
    String saveText;
    if (fileDialog == null) {
      fileDialog = new FileDialog(this);
    }
    fileDialog.setMode(FileDialog.LOAD);
    fileDialog.show();
    
    String file = fileDialog.getFile();
    if (file == null) {
      return;
    }
    String directory = fileDialog.getDirectory();
    File f = new File(directory, file);
    try{
      FileOutputStream fo = new FileOutputStream(f);
      PrintWriter outf = new PrintWriter(new DataOutputStream(fo));
      outf.println(editor.getText());
      outf.close();
    }catch(IOException e){ System.err.println(e.getMessage()); }
  }
}
