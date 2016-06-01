/* Changes of Physical Panel(Repeated vols Button) 9,13,1998
 *
 */

import com.sun.java.swing.*;
import com.sun.java.swing.border.*;
import com.sun.java.swing.event.*;
//import FileChooserDemo.src.ExampleFileFilter;
import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

class PhysicalPanel extends JPanel implements ActionListener{

  private JButton single, replicas, trans, rotate, transSlice, 
                  rotateSlice, parame;
  private JLabel phyLabel,repeatedLabel,singleLabel,sliceLabel,parameLabel;
  PhysicalItem phyItem;

  VolumesFrame parent;
  PhysicalPanel(VolumesFrame parent){
    this.parent = parent;
    GridBagLayout gbl = new GridBagLayout();
    GridBagConstraints gbc = new GridBagConstraints();
    setLayout(gbl);
    setBorder(new TitledBorder(LineBorder.createBlackLineBorder()," Physical Volume "));
    
      gbc.gridx = 0; gbc.gridy = 0;
      gbc.weightx = 1; gbc.weighty = 1;
      gbc.fill = GridBagConstraints.BOTH;
      gbc.gridwidth = GridBagConstraints.REMAINDER;
      single   = new JButton(" Single Placement");
      gbl.setConstraints(single, gbc);
      add(single);

      JPanel repPanel = new JPanel();
      repPanel.setLayout(new BorderLayout());
      repeatedLabel = new JLabel("  Repeted Placement ",SwingConstants.LEFT);
      repeatedLabel.setForeground(Color.blue);
      repPanel.add("North",repeatedLabel);      
      
      JPanel axialortrans = new JPanel();
       axialortrans.setLayout(new GridLayout(1,2));

        trans = new JButton("Translation");
         axialortrans.add("Center",trans);      

        rotate = new JButton("AxialSymmetric");
         axialortrans.add("East",rotate);
      
      repPanel.add("Center",axialortrans);      

      gbc.gridx = 0; gbc.gridy = 1;
      gbc.weightx = 1; gbc.weighty = 1;
      gbc.fill = GridBagConstraints.BOTH;
      gbc.gridwidth = GridBagConstraints.REMAINDER;
      gbl.setConstraints(repPanel, gbc);
      add(repPanel);
      
      JPanel slicePanel = new JPanel();
      slicePanel.setLayout(new BorderLayout());
      sliceLabel = new JLabel("  Slicing ",SwingConstants.LEFT);
      sliceLabel.setForeground(Color.blue);
      slicePanel.add("North",sliceLabel);

      JPanel axiatransSlice = new JPanel();
       axiatransSlice.setLayout(new GridLayout(1,2));

        transSlice = new JButton("Translation");
         axiatransSlice.add("Center",transSlice);

        rotateSlice = new JButton("AxialSymmetric");
         axiatransSlice.add("East",rotateSlice);

      slicePanel.add("Center",axiatransSlice);



      gbc.gridx = 0; gbc.gridy = 2;
      gbc.weightx = 1; gbc.weighty = 1;
      gbc.fill = GridBagConstraints.BOTH;
      gbc.gridwidth = GridBagConstraints.REMAINDER;
      gbl.setConstraints(slicePanel, gbc);
      add(slicePanel);

      gbc.gridx = 0; gbc.gridy = 3;
      gbc.weightx = 1; gbc.weighty = 1;
      gbc.fill = GridBagConstraints.BOTH;
      gbc.gridwidth = GridBagConstraints.REMAINDER;
      parame = new JButton("Parametrisation");
      gbl.setConstraints(parame, gbc);
      add(parame);

    single.addActionListener(this);
    single.setForeground(Color.red);
    repeatedLabel.setForeground(Color.blue);
    trans.addActionListener(this);
    rotate.addActionListener(this);
    trans.setForeground(Color.red);
    rotate.setForeground(Color.red);
    transSlice.addActionListener(this);
    transSlice.setForeground(Color.red);
    rotateSlice.addActionListener(this);
    rotateSlice.setForeground(Color.red);
    parame.addActionListener(this);
    parame.setForeground(Color.red);
//    parame.setEnabled(false);
    setSize(280, 400);
  }

  public void actionPerformed(ActionEvent ae){
    Object o = ae.getSource();
    if (o == single){ 
      parent.singleOpen();
      repaint();
    }else if (o == trans){
      parent.repeatOpen();
    }else if (o == rotate){
      parent.repeatRotateOpen();
    }else if (o == transSlice){
      parent.replicaOpen();
    }else if (o == rotateSlice){
      parent.slicingOpen();
    }else if (o == parame){
      parent.parameOpen();
    }
  }


//C++ source code
  String getCPP(){
    StringBuffer cpp = new StringBuffer("\n// Physical Volumes ----  Single Positioned Placement,   Repeated Placement,   Slicing  --------------------------- \n\n");
    return cpp.toString();
  }

}















