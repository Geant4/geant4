// GAG (Geant4 Adaptive GUI)
// Requires : JDK 1.1.3 (or later) + swing-1.0.2
// Version : alpha  1998/Mar/31
// GEANT4 GUI Group, Toshiaki Kodama
// H Yoshida 1998 June , for Beta release
//  Now, use swing 1.0.2  JFileChooser to select by file extensions
//  Now displayed labels change color to notify busy state.
// July 1, add @@Ask to accept users key-in with showInputDialog (imcomplete)
//  July 5 now use more readable fonts.
	
package GAG;

import TextEdit.TextEdit;
// Swing 1.0.2 ExampleFileFilter.class is in the same dir as gag.class
import ExampleFileFilter;

import com.sun.java.swing.preview.*;
import com.sun.java.swing.preview.filechooser.*;
// now use Swing-1.0.2 FileChooser, yet in the preview stage !!

import java.awt.*;
//import java.awt.Cursor.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import com.sun.java.swing.tree.*;

public class GAG extends JFrame implements TreeSelectionListener, ActionListener {
  static final String gagVer = "J1.0a";
  private String g4Ver, command;
  private GAGpipe geant4;
  private ErrStream errStream;
  private GAGtextEditor textEdit;
  private Hashtable paramsHash;
  private Vector disableCommands;
  private ParamPanel paramPanel;
  private JPanel geantPanel;
  private JSplitPane splitPane;
  private JFileChooser geant4Chooser, historyChooser;
  private JMenuItem execG4, exitG4, killG4, exitGAG, clearOneHist, clearAllHist, saveHist, htmlHelp, about;
  private JTextField directField;
  private JLabel promptLabel;
  private String fileName;
  private CardLayout cardLayout;
  private boolean endGAG;
//use p1' background color to show GEANT4's busy states. Not Yet!
  private JPanel p1;

// FileFilter in Swing-1.0.2 test
 private ExampleFileFilter exeFilter, outFilter, g4hFilter,
				g4histFilter, txtFilter, histFilter;

  public GAG(){
    super("GAG");
    getContentPane().setLayout(cardLayout = new CardLayout());
    geantPanel = new JPanel(new BorderLayout());

    JMenu geantMenu = new JMenu("GEANT4");

	geantMenu.setBackground(Color.orange);
	geantMenu.setFont(new Font("Serif", Font.BOLD, 12));

	execG4 = new JMenuItem("Exec GEANT4");
	exitG4 = new JMenuItem("Exit GEANT4");
	killG4 = new JMenuItem("Kill GEANT4");
	exitGAG = new JMenuItem("Exit GAG");
	execG4.setBackground(Color.green);
	exitG4.setBackground(Color.green);
	killG4.setBackground(Color.orange);
	exitGAG.setFont(new Font("Serif", Font.BOLD, 12));
	execG4.setFont(new Font("Serif", Font.BOLD, 12));
	exitG4.setFont(new Font("Serif", Font.BOLD, 12));
	killG4.setFont(new Font("Serif", Font.BOLD, 12));
	exitGAG.setFont(new Font("Serif", Font.BOLD, 12));

    geantMenu.add(execG4);
    geantMenu.add(exitG4);
    geantMenu.add(killG4);
    geantMenu.addSeparator();
    geantMenu.add(exitGAG);

    JMenu histMenu = new JMenu("History");

	histMenu.setBackground(Color.yellow);
	histMenu.setFont(new Font("Serif", Font.BOLD, 12));

	clearOneHist = new JMenuItem("Clear One");
	clearAllHist = new JMenuItem("Clear All");
	saveHist = new JMenuItem("Save history");
	clearOneHist.setBackground(Color.yellow);
	clearAllHist.setBackground(Color.yellow);
	saveHist.setBackground(Color.green);
	clearOneHist.setFont(new Font("Serif", Font.BOLD, 12));
	clearAllHist.setFont(new Font("Serif", Font.BOLD, 12));
	saveHist.setFont(new Font("Serif", Font.BOLD, 12));


    histMenu.add(clearOneHist);
    histMenu.add(clearAllHist);
    histMenu.add(saveHist);

    JMenu helpMenu = new JMenu("Help");
	helpMenu.setBackground(Color.white);
	helpMenu.setFont(new Font("Serif", Font.BOLD, 12));
    helpMenu.add(htmlHelp = new JMenuItem("Help"));
    helpMenu.add(about = new JMenuItem("About"));

    execG4.addActionListener(this);
    exitG4.addActionListener(this);
    killG4.addActionListener(this);
    exitGAG.addActionListener(this);
    clearOneHist.addActionListener(this);
    clearAllHist.addActionListener(this);
    saveHist.addActionListener(this);
    htmlHelp.addActionListener(this);
    about.addActionListener(this);;

    JMenuBar mb = new JMenuBar();
    mb.add(geantMenu);
    mb.add(histMenu);
    mb.add(helpMenu);    // mb.setHelpMenu(helpMenu);

//to the last position    executingMenu(false);
    setJMenuBar(mb);
    directField = new JTextField();
	directField.setFont(new Font("Serif", Font.BOLD, 12));
    directField.addActionListener(this);

    p1 = new JPanel(new BorderLayout());
    promptLabel = new JLabel("****** PROMPT ******", JLabel.CENTER);
    promptLabel.setForeground(Color.black);
	promptLabel.setFont(new Font("Serif", Font.BOLD, 12));
    p1.add("Center", directField);
    p1.add("West", promptLabel);

    geantPanel.add("South", p1);
	JLabel idleLabel =new JLabel("Welcome to GEANT4",JLabel.CENTER);
	idleLabel.setFont(new Font("Serif", Font.BOLD, 18));
    getContentPane().add("idle", idleLabel);
	JLabel treeLabel = new JLabel("Making command tree");
	treeLabel.setFont(new Font("Serif", Font.BOLD, 18));
    getContentPane().add("tree", treeLabel);
    getContentPane().add("geant",geantPanel);

    cardLayout.show(getContentPane(), "idle");
    setSize(600,400);
    setVisible(true);
    geant4Chooser = new JFileChooser(".");
// Swing 1.0.2 FileFilter  OK!! 1998 June 25
	geant4Chooser.setDialogTitle("GEANT4 Executable");

	exeFilter = new ExampleFileFilter("exe", "DOS executable");
	outFilter = new ExampleFileFilter("out", "Unix executable");

	geant4Chooser.addChoosableFileFilter(exeFilter);
	geant4Chooser.addChoosableFileFilter(outFilter);
// 

    historyChooser = new JFileChooser(".");

// Swing 1.0.2 FileFilter  OK!! 1998 June 25
	historyChooser.setDialogTitle("Command History");
	g4hFilter = new ExampleFileFilter("g4h", "Geant4 History");
	g4histFilter = new ExampleFileFilter("g4hist", "Geant4 History");
	txtFilter = new ExampleFileFilter("txt", "Geant4 History");
	histFilter = new ExampleFileFilter("hist", "Geant4 History");

	historyChooser.addChoosableFileFilter(g4hFilter);
	historyChooser.addChoosableFileFilter(g4histFilter);
	historyChooser.addChoosableFileFilter(txtFilter);
	historyChooser.addChoosableFileFilter(histFilter);

// end of Swing test

    executingMenu(false);

  }
  private void executingMenu(boolean b){
    execG4.setEnabled(!b);
//originally     exitG4.setEnabled(b);
	    exitG4.setEnabled(!b);
    killG4.setEnabled(b);
//yoshida
//	exitGAG.setEnabled(b);
    clearOneHist.setEnabled(b);
    clearAllHist.setEnabled(b);
    saveHist.setEnabled(b);

//
//this hide allthe time	p1.setVisible(!b);
	if (!b) {
		p1.setBackground(Color.green); repaint();
	}else{
		p1.setBackground(Color.red); repaint();
	}
  }

  public void mainLoop(String args[]){
    endGAG = false;
    if (args.length > 0) {
      fileName = args[0];
      geantStart(args);
    }
    String[] fn = new String[1];
    while(!endGAG){
      waiting();
      if (endGAG) return;
      fn[0] = fileName;
      geantStart(fn);
    }
  }
  private synchronized void waiting(){
    try{
      wait();
    }catch(InterruptedException e){
    }
  }
  private synchronized void resume(){
    notify();
  }
  private void geantStart(String[] args){
    String line, prompt = "G4";
    geant4 = new GAGpipe(args);
    boolean paramChanged = false;
    if (geant4.isError()){
      JOptionPane.showMessageDialog(this, geant4.getErrorMsg(),
				    "Warning", JOptionPane.WARNING_MESSAGE);
      geant4 = null;
      return;
    }
    executingMenu(true);
    if (textEdit == null) textEdit = new TextEdit();
    errStream = new ErrStream(geant4);
    errStream.start();
    geant4.writeLine("@@GAGmodeJAVA");
    while ((line = geant4.readStdLine()) != null){
      if (line.startsWith("@@")){
	if (line.equals("@@JTreeBegin")){
	  if (makeCommandTree()) break;
	  continue;
	}
	if (line.equals("@@JParamBegin")){
	  if (receiveParams()) break;
	  paramChanged = true;
	  continue;
	}
	if (line.equals("@@JDirGuideBegin")){
	  if (receiveCascade()) break;
	  continue;
	}
	if (line.equals("@@DisableListBegin")){
	  if (receiveDisables()) break;
	  paramChanged = true;
	  continue;
	}
	if (line.equals("@@Ready")){
//yoshida effective even if paramPanel is null. But PROMPT case is enough.
	p1.setVisible(true);

	p1.setBackground(Color.green); repaint();
	  if (paramPanel != null) {
	    paramPanel.toReady();
	    if (paramChanged && !command.equals("")) setParamPanel(command);
	  }
	  paramChanged = false;
	  textEdit.writeString(prompt + ">");
	  exitG4.setEnabled(true);
	  continue;
	}
	if (line.startsWith("@@PROMPT")){
	  prompt = line.substring(line.indexOf("\"")+1, line.lastIndexOf("\""));
	  promptLabel.setText(prompt);
//yoshida effective even if paramPanel is null
	p1.setVisible(true);
	p1.setBackground(Color.green); repaint();
	setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
	show();

	  continue;
	}
	if (line.startsWith("@@State")){
	  continue;
	}
// direct interaction with users
	if (line.startsWith("@@Ask")){
	String answer = JOptionPane.showInputDialog(this,line);
	geant4.writeLine(answer); 
	continue;
	}
//direct questionnaire

	if (line.startsWith("@@CurrentValue")){
	  paramPanel.loadCurrentValues(line);
	  textEdit.writeLine(line.substring(line.indexOf(" ")+1));
	  continue;
	}
	if (line.startsWith("@@Version")){
	  StringTokenizer st = new StringTokenizer( line );
	  st.nextToken();
	  g4Ver = st.nextToken();
	  continue;
	}
	if (line.startsWith("@@ErrResult")){
	  JOptionPane.showMessageDialog(this, 
		line.substring(line.indexOf("\"")+1, line.lastIndexOf("\"")),
		"Warning", JOptionPane.WARNING_MESSAGE);
	  continue;
	}
      }
      textEdit.writeLine(line);
    }
    geantTerminate();
  }
  private void geantTerminate(){
    executingMenu(false);
    if (geant4 != null){
      geant4.processKill();
      geant4 = null;
    }
    paramsHash = null;
    if (errStream != null){
      errStream.close();
      errStream = null;
    }
    cardLayout.show(getContentPane(), "idle");
    paramPanel = null;
  }
  private boolean makeCommandTree() {
    String pLine, label;
    command = "";
    DefaultMutableTreeNode topTree, nowTree, childTree;
    boolean found;
    cardLayout.show(getContentPane(), "tree");
    if (splitPane != null) geantPanel.remove(splitPane);
    int lastIndex = fileName.lastIndexOf("/");
    topTree = new DefaultMutableTreeNode( (lastIndex<0) ? fileName :  fileName.substring(lastIndex+1) );
    while( !(pLine = geant4.readStdLine()).equals("@@JTreeEnd") ){
      if (pLine == null){ return true;}
      nowTree = topTree; childTree = null; found = true;
      StringTokenizer stCo = new StringTokenizer(pLine,"/");
      while(stCo.hasMoreTokens()){
	label = stCo.nextToken();
	if (found){
	  found = false;
	  for (int i=0; i<nowTree.getChildCount(); i++){
	    childTree = (DefaultMutableTreeNode)nowTree.getChildAt(i);
	    if (label.equals( (String)childTree.getUserObject() )){
	      found = true;
	      break;
	    }
	  }
	}
	if (found == false) nowTree.add(childTree = new DefaultMutableTreeNode(label));
	nowTree = childTree;
      }
    }
    JTree tree = new JTree(topTree);
    tree.addTreeSelectionListener(this);
    paramPanel = new ParamPanel(this);
    splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, new JScrollPane(tree), paramPanel );

// test 1998 June 25 to keep left space
	splitPane.setDividerLocation(200);

    geantPanel.add("Center", splitPane);
    cardLayout.show(getContentPane(), "geant");
    paramsHash = new Hashtable();
    disableCommands = new Vector();
    return false;
  }
  private boolean receiveParams() {
    String lines, name, guide[], range;
    int guideEnts, paramEnts;
    name = geant4.readStdLine();
    lines = geant4.readStdLine();
    try{
      guideEnts = Integer.parseInt(lines);
    }catch(NumberFormatException e){
      return true;
    }
    guide = new String[guideEnts];
    for (int i=0; i<guideEnts; i++){
      guide[i] = geant4.readStdLine();
    }
    range = geant4.readStdLine();
    if (geant4.isError()) return true;
    lines = geant4.readStdLine();
    try{
      paramEnts = Integer.parseInt(lines);
    }catch(NumberFormatException e){
      return true;
    }
    GAGcommandItem ci = new GAGcommandItem(name, guide, range, paramEnts);
    for (int i=0; i<ci.getParamEntries(); i++){
      for (int j=0; j<ci.MAX; j++){
	lines = geant4.readStdLine();
	if (lines == null) return true;
	ci.params[i][j]=lines;
      }
    }
    if ( !(geant4.readStdLine()).equals("@@JParamEnd") ) return true;
    paramsHash.put(name, ci);
    return false;
  }
  private boolean receiveCascade() {
    String name, guide;
    name = geant4.readStdLine();
    guide = geant4.readStdLine();
    // System.out.println("CASCADE: "+name+"  "+guide);
    if ( !(geant4.readStdLine()).equals("@@JDirGuideEnd") ) return true;
    paramsHash.put(name, guide);
    return false;
  }
  private boolean receiveDisables(){
    String line;
    disableCommands.removeAllElements();
    while(!(line = geant4.readStdLine()).equals("@@DisableListEnd")){
      if (line == null) return true;
      disableCommands.addElement(line);
    }
    return false;
  }
  void reqCurrent(String cmd){
    sendCommand("?"+cmd);
    return;
  }
  void sendCommand(String com){
    textEdit.writeLine(com);
    geant4.writeLine(com);
   exitG4.setEnabled(false);
//yoshida cursor
	setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
	p1.setBackground(Color.red); repaint();
	p1.setVisible(true);
//	executingMenu(true);
// too much?
  }
  boolean setParamPanel(String key){
    GAGcommandItem ci = (GAGcommandItem)paramsHash.get(key);
    if (ci == null) return false;
    command = key;
    paramPanel.setParamBox(ci, isEnableCommand(key));
    return true;
  }
  public void valueChanged(TreeSelectionEvent te){
    if ( te.isAddedPath() == false ) return;
    String dirHelp, tempPath;
    StringBuffer com = new StringBuffer();
    Object path[] = te.getPath().getPath();
    if (path.length > 1) {
      for (int i=1; i<path.length; i++) com.append(("/"+path[i]));
      command = com.toString();
    }else{
      command = "";
    }

    DefaultMutableTreeNode nowNode = (DefaultMutableTreeNode)path[path.length-1];
    int childCount = nowNode.getChildCount();
    if ( childCount > 0 ){
      DefaultMutableTreeNode tempNode;
      if (path.length > 1) {
	dirHelp = (String)paramsHash.get(command+"/");
      }else{
	dirHelp = new String(path[0].toString());
      }
      JLabel[][] helpLabel = new JLabel[2][childCount];
      for (int i=0; childCount > i; i++){
	tempNode = (DefaultMutableTreeNode)nowNode.getChildAt(i);
	tempPath = command + "/" + tempNode;
	if (tempNode.getChildCount() == 0){
	  helpLabel[0][i] = new JLabel(tempNode+" ");
	  helpLabel[0][i].setForeground(isEnableCommand(tempPath) ? Color.red : Color.lightGray);
	  helpLabel[1][i] = new JLabel(((GAGcommandItem)(paramsHash.get(tempPath))).getTitle());
	}else{
	  helpLabel[0][i] = new JLabel(tempNode+"/ ");
	  helpLabel[0][i].setForeground(isEnableCommand(tempPath+"/") ? Color.blue : Color.lightGray);
	  helpLabel[1][i] = new JLabel((String)paramsHash.get(tempPath +"/"));
	}
      }
      paramPanel.setTreeHelp(dirHelp, helpLabel);
    }else{
      GAGcommandItem ci = (GAGcommandItem)paramsHash.get(command);
      if (ci == null) return;
      paramPanel.setParamBox(ci, isEnableCommand(command));
    }
    return;
  }
  boolean isEnableCommand(String command){
    return !disableCommands.contains(command);
  }
  public void actionPerformed(ActionEvent ae){


    Object o = ae.getSource();
    if ( o == directField){
      sendCommand(directField.getText());
      paramPanel.addHistory(directField.getText());
      directField.setText("");
      return;
    }
    if (o == exitGAG){
      if (geant4 != null) {
	if(JOptionPane.showConfirmDialog(this, "Kill GEANT4 automatically.", "Exit GAG", JOptionPane.YES_NO_OPTION) == 1) {return;}
	geantTerminate();
	endGAG = true;
      }else{
	resume();
	endGAG = true;
      }
      return;
    }
    if (o == exitG4){
      sendCommand("exit");
      return;
    }
    if (o == execG4){
      if (geant4Chooser.showOpenDialog(this) == -1 )return;
// Swing1.0.2
      File f = geant4Chooser.getSelectedFile();
      if (f.isFile()) {
	fileName = f.getPath();
	resume();
      }
      return;
    }
    if (o == killG4){
      geantTerminate();
      return;
    }
    if (o == htmlHelp){
      return;
    }
    if (o == about){
      new AboutDialog(this, gagVer, g4Ver);
      return;
    }
    if (o == clearAllHist){
      paramPanel.clearAllHistory();
      return;
    }
    if (o == clearOneHist){
      paramPanel.clearOneHistory();
      return;
    }
    if (o == saveHist){
      if (historyChooser.showSaveDialog(this) == -1 ) return;
// Swing1.0.2
     try{
	FileOutputStream fo = new FileOutputStream(historyChooser.getSelectedFile());
	PrintWriter outf = new PrintWriter(new DataOutputStream(fo));
	String[] items = paramPanel.getHistoryItems();
	for (int i=0; i<items.length; i++) outf.println(items[i]);
	outf.close();
      }catch(IOException e){ System.out.println(e.getMessage()); }
      return;
    }
  }
}

class ErrStream extends Thread implements ActionListener {
  private GAGpipe geant4;
  private Frame errFrame;
  private TextArea errArea;
  private Button close;
  ErrStream(GAGpipe geant4){
    this.geant4 = geant4;
    errFrame = new Frame("Std Error");
    errFrame.setLayout(new BorderLayout());
    errFrame.add("Center", errArea = new TextArea());
    errFrame.add("South", close = new Button("Close"));
    close.addActionListener(this);
    errFrame.setSize(400, 200);
  }
  public void run(){
    String line;
    while((line=geant4.readErrLine())!=null){
      errArea.append(line+"\n");
      if (!errFrame.isVisible()) errFrame.setVisible(true);
    }
  }
  public void actionPerformed(ActionEvent ae){
    errFrame.setVisible(false);
    errArea.setText("");
  }
  void close(){
    errFrame.dispose();
    stop();
  }
}

class AboutDialog extends Dialog implements ActionListener {
  Button ok;
  public AboutDialog(Frame parent, String ver1, String ver2){
    super(parent, "About GAG", true);
    setForeground(Color.black);
    Panel center = new Panel(new GridLayout(0,1));
    center.add(new Label("GAG Version "+ ver1));
    center.add(new Label("G4UIGAG Version "+ ver2));
    center.add(new Label("Toshiaki Kodama, GEANT4 GUI GROUP"));
    add("Center", center);

    Panel p = new Panel();
    ok = new Button("Ok");
    ok.addActionListener(this);
    p.add(ok);
    add("South",p);
    pack();
    setVisible(true);
  }
  public void actionPerformed(ActionEvent event) {
    if (event.getSource() == ok)
      dispose();
  }
}


