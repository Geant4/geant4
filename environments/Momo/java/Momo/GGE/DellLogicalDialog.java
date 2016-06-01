
//example
//8.5

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class DellLogicalDialog extends JDialog implements ActionListener{
   private JButton b1, b2;
   private JLabel label1;
   private VolumesFrame volFrame;
   public DellLogicalDialog(VolumesFrame volFrame){
      super(volFrame,"Delete Logical Volume",false);
      this.volFrame = volFrame;
      JPanel createLogPanel = new JPanel();
       createLogPanel.setLayout(new BorderLayout());
       b1 = new JButton("OK");
       b2 = new JButton("CANCEL");
       b1.addActionListener(this);
       b2.addActionListener(this);

        label1 = new JLabel("   Click the pName to delete");
        label1.setFont(new Font("Serif",Font.BOLD, 13));

       createLogPanel.add(label1, BorderLayout.NORTH);
       createLogPanel.add(b1, BorderLayout.CENTER);
       createLogPanel.add(b2, BorderLayout.EAST);
      getContentPane().add(createLogPanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 80);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        volFrame.logPanel.dell();
//        this.setVisible(false); 
      }else if(e.getSource() == b2){
        this.setVisible(false);       
        volFrame.logPanel.dellLogCloseAct();
      }
   }
 
}






