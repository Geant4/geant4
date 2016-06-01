/*
 *
 */

import com.sun.java.swing.*;
import com.sun.java.swing.plaf.basic.*;
import com.sun.java.swing.border.*;
import com.sun.java.swing.event.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;

class VisAtbPanel extends JPanel implements ActionListener, ListSelectionListener {
  private JButton addVis, delVis;
  JList visList;
  VisAtbItem visItem;
  DefaultListModel listModel;
  private VisAtbDialog visDialog;
  VisAtbPanel(VolumesFrame parent){
    setLayout(new BorderLayout());
    setBorder(new TitledBorder(LineBorder.createBlackLineBorder(),
                                                      " VisAttrib "));
    add("Center",visList = new JList(listModel = new DefaultListModel()));
    JPanel buttons = new JPanel(new BorderLayout());
    buttons.add("Center", addVis = new JButton(" Create "));
    buttons.add("East", delVis = new JButton("Del"));
    delVis.setEnabled(false);
    add("North", buttons);
    visDialog = new VisAtbDialog(parent);
    addVis.addActionListener(this);
    delVis.addActionListener(this);
    visList.addListSelectionListener(this);
    visList.addMouseListener(new editVis());
    visList.setCellRenderer(new VisCellRenderer());
  }
  VisAtbItem getSelectedVis(){
    return (VisAtbItem)visList.getSelectedValue();
  }
  public void actionPerformed(ActionEvent ae){
    Object o = ae.getSource();
    if (o == addVis){

      addVis.setForeground(Color.red);
      VisAtbItem visItem = visDialog.getNewAtb();
      if (visItem == null) return;
      listModel.addElement(visItem);
      addVis.setForeground(Color.black);
      repaint();
    }else if (o == delVis){
      delVis.setForeground(Color.red);
      int index = visList.getSelectedIndex();
      if (index < 0) return;
      listModel.removeElementAt(index);
    }
  }
 

  public void visCloseAct(){addVis.setForeground(Color.black);}  

  public void visClear(){
    visList.setModel(listModel = new DefaultListModel());
    repaint();
  }
 
  public void valueChanged(ListSelectionEvent le){
    if (le.getValueIsAdjusting()) return;
    System.out.println(le);
    delVis.setEnabled(!(visList.getSelectedIndex() < 0));
  }
  class editVis extends MouseAdapter {
    public void mouseClicked(MouseEvent e){
      if (e.getClickCount() == 2){
	int index = visList.locationToIndex(e.getPoint());
	if (index < 0) return;
	VisAtbItem newItem = visDialog.getChangedAtb((VisAtbItem)listModel.getElementAt(index));
	if (newItem == null) return;
	listModel.setElementAt(newItem, index);
      }
    }
  }
  String getCPP(){
    StringBuffer cpp = new StringBuffer("\n// Visualization attributes\n\n");
    if (listModel.getSize()==0){
      return cpp.toString();
    }else{
     for(int i=0; i<listModel.getSize(); i++){
      VisAtbItem vObject = (VisAtbItem)listModel.elementAt(i);
      cpp.append("\nG4VisAttributes * ");
      cpp.append(vObject.name+"= new G4VisAttributes( G4Colour(");
      cpp.append(vObject.color.getRed()+"/255. ,"+vObject.color.getGreen()+"/255. ,"+vObject.color.getBlue()+"/255. ));\n");
      for (int j=0; j<vObject.atbsIsValid.length; j++){
        if (vObject.atbsIsValid[j]){
          cpp.append(vObject.name+ "->Set"+vObject.atbMethodName[j]+"("+vObject.atbsValue[j]+");\n");
        }
      }
     }
     return cpp.toString();
    }
  }
}

class VisCellRenderer extends BasicListCellRenderer {
    //class VisCellRenderer extends JLabel implements ListCellRenderer {
  public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus){
    setIcon(new VisIcon(((VisAtbItem)value).color));
    return super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
    //    return this;
  }
  class VisIcon implements Icon {
    Color color;
    VisIcon(Color color){
      this.color = color;
    }
    public void paintIcon(Component c, Graphics g, int x, int y){
      g.setColor(color);
      g.fillOval(0,0,10,10);
    }
    public int getIconWidth(){
      return 8;
    }
    public int getIconHeight(){
      return 8;
    }
  }
}





