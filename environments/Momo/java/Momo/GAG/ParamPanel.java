// GAG (Geant4 Adaptive GUI)
// Requires : JDK 1.1.3 (or later) + swing-1.0.2
// Version : alpha  1998/Mar/31
// GEANT4 GUI Group, Toshiaki Kodama
// modified by  H Yoshida, 1998 June 25
// use of FileFilter of Swing 1.0.2, to select after file extensions
// improvement of notification method of GEANT4 states
// by using setVisible and setBackground of relevant buttons
// specially "ok" button.

package GAG;

import com.sun.java.swing.preview.*;
//FileFilter of Swing1.0.2
import ExampleFileFilter;

import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.preview.filechooser.*;
// for Swing-1.0.2 FileChooser

class ParamPanel extends JPanel implements ActionListener, ItemListener {
  private GAG gag;
  private GridBagLayout gbl; GridBagConstraints gbc;
  private JPanel paramArea;
  private JButton def, current, clear, ok, fileBrowser;
  private JComboBox history;
  private int posY;
  private GAGcommandItem ci;
  private Component paraBox[];
  private boolean isExecuting;
  private JFileChooser fileChooser;
// FileFilter in Swing-1.0.2 test
 private ExampleFileFilter g4mFilter, macroFilter, g4macroFilter,
 txtFilter;



  ParamPanel(GAG gag){
    isExecuting = false;
    this.gag = gag;
    setLayout( new BorderLayout() );

    gbc = new GridBagConstraints();
    gbl = new GridBagLayout();

    paramArea = new JPanel();
    paramArea.setLayout(gbl);
    paramArea.setBackground(new Color(255, 250, 230));
    add("Center", new JScrollPane(paramArea));

    JPanel buttons = new JPanel();
    buttons.setLayout(new GridLayout( 1, 0, 2, 2));
//  Yoshida color
	def  = new JButton("Default");
	current = new JButton("Current");
	clear = new JButton("Clear");
	ok = new JButton("OK");
	ok.setBackground(Color.green);
//no good??	ok.setFont(new Font("OK",1,12));
	def.setBackground(Color.blue);
	current.setBackground(Color.orange);
	clear.setBackground(SystemColor.white);
//	ok.setBackground(SystemColor.control);
// color

    buttons.add( def );
    buttons.add( current );
    buttons.add( clear );
    buttons.add( ok );
    add("South", buttons);

    history = new JComboBox();
    history.addItem("# GEANT4 history macro.");
    history.addItemListener(this);
    add("North", history);
//yoshida setinvisible and setdisable
	def.setVisible(false);
    def.setEnabled(false);
    def.setMargin(new Insets(1,3,1,3));
    def.addActionListener(this);

//yoshida
	current.setVisible(false);
    current.setEnabled(false);
    current.setMargin(new Insets(1,3,1,3));
    current.addActionListener(this);
//yoshida
	clear.setVisible(false);
    clear.setEnabled(false);
    clear.setMargin(new Insets(1,3,1,3));
    clear.addActionListener(this);

//yoshida
	ok.setVisible(false);
    ok.setEnabled(false);
    ok.setMargin(new Insets(1,3,1,3));
    ok.addActionListener(this);

    fileBrowser = new JButton("Open a Macro File");
    fileBrowser.addActionListener(this);
    fileChooser = new JFileChooser(".");
// Swing 1.0.2 FileFilter OK!! 1998 June 25
	fileChooser.setDialogTitle("Macro File");
        g4mFilter = new ExampleFileFilter("g4m", "Geant4 macro");
        g4macroFilter = new ExampleFileFilter("g4macro", "Geant4 macro");
        txtFilter = new ExampleFileFilter("txt", "Geant4 macro");
        macroFilter = new ExampleFileFilter("macro", "Geant4 macro");

        fileChooser.addChoosableFileFilter(g4mFilter);
        fileChooser.addChoosableFileFilter(g4macroFilter);
        fileChooser.addChoosableFileFilter(txtFilter);
        fileChooser.addChoosableFileFilter(macroFilter);

// end of Swing-1.0.2 FileChooser


  }
  void setParamBox(GAGcommandItem ci, boolean enable){
    paramArea.removeAll();
    this.ci = ci;
    posY = 0;
    JLabel l;
    if (!enable){
      l = new JLabel("This command is not available.");
      l.setForeground(Color.red);
      addComponent(l);
    }
    addComponent(new JLabel(ci.getCommandName()));
    String guide[] = ci.getCommandGuide(); 
    for (int i=0; i<guide.length; i++){
      addComponent(new JLabel(guide[i]));
    }
    if (!ci.getCommandRange().equals("")){
      l = new JLabel(ci.getCommandRange());
      l.setForeground(Color.green);
      addComponent(l);
    }
    if (ci.getParamEntries() > 0){
      def.setEnabled(true);
//yoshida
	def.setVisible(true);
      current.setEnabled(!isExecuting && enable);
//yoshida
      current.setVisible(!isExecuting && enable);
      clear.setEnabled(true);
//yoshida
	clear.setVisible(true);
      paraBox = new JComponent[ci.getParamEntries()];
      String candidateList;
      for (int i=0; i<ci.getParamEntries(); i++){
	candidateList = ci.getParamCandidate(i);
	if (candidateList.equals("")){
	  if (ci.getParamType(i).equals("b")){
	    paraBox[i] = new BooleanCombo(ci.getParamDefault(i));
	  }else{
	    paraBox[i] = new JTextField(20);
	    ((JTextField)paraBox[i]).setText(ci.getParamDefault(i));
	  }  
	}else{
	  paraBox[i] = new JComboBox();
	  StringTokenizer st = new StringTokenizer(candidateList);
	  while(st.hasMoreTokens()){
	    ((JComboBox)paraBox[i]).addItem(st.nextToken());
	  }
	  if (!ci.getParamDefault(i).equals("")){
	    ((JComboBox)paraBox[i]).setSelectedItem(ci.getParamDefault(i));
	  }
	}
	addComponent(new JLabel(ci.getParamName(i)), paraBox[i], new JLabel("("+ci.getParamType(i) +") "+ci.getParamGuide(i) ));
	if (ci.getParamName(i).equals("fileName")) addComponent(new JLabel(), fileBrowser);
      }
    }else{
//yoshida setinvisible and then setdisable
	def.setVisible(false);
      def.setEnabled(false);
	current.setVisible(false);
      current.setEnabled(false);
	clear.setVisible(false);
      clear.setEnabled(false);
    }
    ok.setEnabled(!isExecuting && enable);
//yoshida OK button order of enable and visible is important
    ok.setVisible(!isExecuting && enable);

    validate();
    repaint();
    return;
  }

  void setTreeHelp(String dirHelp, JLabel helps[][]){
//yoshida setinvisible then setdisable
	ok.setVisible(false);
    ok.setEnabled(false);
	def.setVisible(false);
    def.setEnabled(false);
	current.setVisible(false);
    current.setEnabled(false);
	clear.setVisible(false);
    clear.setEnabled(false);
//
    paramArea.removeAll();
    posY = 0; ci = null;
    JLabel title = new JLabel(dirHelp);
    title.setForeground(Color.green);
    addComponent(title);
    for (int i=0; i<helps[0].length; i++){
      addComponent(helps[0][i], helps[1][i]);
    }
    validate();
    repaint();
    return;
  }
  private String makeCommandLine(){
    StringBuffer param = new StringBuffer(ci.getCommandName());
    if (ci.getParamEntries() > 0){
      String val;
      for (int i=0; i<ci.getParamEntries(); i++){
	if (paraBox[i] instanceof JTextField){
	  val = (((JTextField)paraBox[i]).getText()).trim();
	  if (val.equals("")){
	    if (!ci.isOmittable(i)){
	      JOptionPane.showMessageDialog(this, ci.getParamName(i)+"\nThis parameter cannot omitted.", "Error", JOptionPane.WARNING_MESSAGE);
	      return null;
	    }else{
	      if (!ci.getParamDefault(i).equals(""))
		val = ci.getParamDefault(i);
	      else{
		JOptionPane.showMessageDialog(this, ci.getParamName(i)+"\nThis omittable parameter has not default value.", "Error", JOptionPane.WARNING_MESSAGE);
		return null;
	      }
	    }
	  }else if ( val.indexOf(" ") >= 0){
	    JOptionPane.showMessageDialog(this, ci.getParamName(i)+"\nThis parameter cannot include SPACE character", "Error", JOptionPane.WARNING_MESSAGE);
	    return null;
	  }
	}else{
	  val = (String)((JComboBox)paraBox[i]).getSelectedItem(); 
	}
	param.append(" "+val);
      }
//yoshida
	current.setVisible(false);
      current.setEnabled(false);
    }
//yoshida
	ok.setVisible(false);
    ok.setEnabled(false);
    history.addItem(param.toString());
    return param.toString();
  }
  private void addComponent(Component compo){
    gbc.gridx = 0; gbc.gridy = posY++;
    gbc.gridwidth = 3; gbc.gridheight = 1;
    gbc.weightx = 1.0; gbc.weighty = 0.0;
    gbc.anchor = GridBagConstraints.WEST;
    gbl.setConstraints(compo, gbc);
    paramArea.add(compo);
  }
  private void addComponent(Component compo1, Component compo2){
    gbc.gridx = 0; gbc.gridy = posY;
    gbc.gridwidth = 1; gbc.gridheight = 1;
    gbc.weightx = 0.0; gbc.weighty = 0.0;
    gbc.anchor = GridBagConstraints.EAST;
    gbl.setConstraints(compo1, gbc);
    paramArea.add(compo1);
    gbc.gridx = 1;
    gbc.anchor = GridBagConstraints.WEST;
    gbl.setConstraints(compo2, gbc);
    paramArea.add(compo2);
    posY++;
  }
  private void addComponent(Component title, Component compo1, Component compo2){
    gbc.gridx = 0; gbc.gridy = posY;
    gbc.gridwidth = 1; gbc.gridheight = 1;
    gbc.weightx = 0.0; gbc.weighty = 0.0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.EAST;
    gbl.setConstraints(title, gbc);
    paramArea.add(title);
    gbc.gridx = 1;
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbl.setConstraints(compo1, gbc);
    paramArea.add(compo1);
    gbc.gridx = 2;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.WEST;
    gbl.setConstraints(compo2, gbc);
    paramArea.add(compo2);
    posY++;
  }
  private void loadValue(int i, String str){
    if (paraBox[i] instanceof JTextField){
      ((JTextField)paraBox[i]).setText(str);
    }
    if (paraBox[i] instanceof JComboBox){
      ((JComboBox)paraBox[i]).setSelectedItem(str);
    }
  }
  void addHistory(String line){
    if (line.startsWith("/")) history.addItem(line);
  }
  void clearAllHistory(){
    history.removeAllItems();
  }
  void clearOneHistory(){
    int ix = history.getSelectedIndex();
    if (ix >= 0) history.removeItemAt(ix);
  }
  String[] getHistoryItems(){
    String[] item = new String[history.getItemCount()];
    for (int i=0; i>history.getItemCount(); i++){
      item[i] = (String)history.getItemAt(i);
    }
    return item;
  }
  void loadCurrentValues(String line){
    StringTokenizer st = new StringTokenizer( line );
    st.nextToken();
    int i=0;
    while(st.hasMoreTokens()){
      loadValue(i++, st.nextToken());
    }
  }
  public void itemStateChanged(ItemEvent ie){
    if (isExecuting) return;
    if (ie.getStateChange() != ItemEvent.SELECTED) return;
    String str = (String)ie.getItem();
    if ( !str.startsWith("/") ) return;
    StringTokenizer st = new StringTokenizer( str );
    if (!gag.setParamPanel(st.nextToken())) return;
    int i=0;
    while(st.hasMoreTokens()){
      loadValue(i++, st.nextToken());
      if (i >= ci.getParamEntries()) return;
    }
  }
  public void actionPerformed(ActionEvent ev){
    Object source = ev.getSource();
    if ( source == ok ){
      String p = makeCommandLine();
      if (p == null) return;
      isExecuting = true;
      gag.sendCommand(p);
//yoshida  instead of sandwatch, ok button becomes invisible
	ok.setVisible(false);
	ok.setEnabled(false);
// Be careful!  mouse may freeze, if setinvisible but not setdisabled!! 
      return;
    }else if (source == def){
      for (int i=0; ci.getParamEntries()>i; i++){
	loadValue(i, ci.getParamDefault(i));
      }
      return;
    }else if (source == current){
      gag.reqCurrent(ci.getCommandName());
      return;
    }else if (source == clear){
      for (int i=0; ci.getParamEntries()>i; i++){
	if (paraBox[i] instanceof JTextField){
	  ((JTextField)paraBox[i]).setText("");
	}
      }
      return;
    }else if (source == fileBrowser){
//swing102 filechooser showDialog(this, null) or the next
      if (fileChooser.showOpenDialog(this) == -1 ) return;
      File f = fileChooser.getSelectedFile();
      if (f.isFile()) {
	if (paraBox[0] instanceof JTextField){
	  ((JTextField)paraBox[0]).setText(f.getPath());
	}
      }
      return;
    }
  }
  void toReady(){
    isExecuting = false;
    if (ci == null) return;
    ok.setEnabled(true);
//yoshida setenable then setvisible
	ok.setVisible(true);
    if (ci.getParamEntries() == 0) return;
    current.setEnabled(true);
//yoshida
	current.setVisible(true);
  }
}

class BooleanCombo extends JComboBox {
  BooleanCombo(String str){
    String def = str.toUpperCase();
    if (def.equals("0") || def.equals("1")){
      addItem("0");
      addItem("1");
    }else if (def.equals("Y") || def.equals("N")){
      addItem("Y");
      addItem("N");
    }else if (def.equals("YES") || def.equals("NO")){
      addItem("YES");
      addItem("NO");
    }else if (def.equals("T") || def.equals("F")){
      addItem("T");
      addItem("F");
    } else {
      addItem("FALSE");
      addItem("TRUE");
    }
    if (!def.equals("")){
      setSelectedItem(def);
    }
  }
}
