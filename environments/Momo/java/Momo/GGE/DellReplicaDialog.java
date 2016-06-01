
//example
//8.5

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class DellReplicaDialog extends JDialog implements ActionListener{
   private JButton b1, b2;
   private JLabel label1;
   private ReplicasFrame repFrame;
   public DellReplicaDialog(ReplicasFrame repFrame){
      super(repFrame,"Delete Slicing Translation",false);
      this.repFrame = repFrame;
      JPanel createRepPanel = new JPanel();
       createRepPanel.setLayout(new BorderLayout());
       b1 = new JButton("OK");
       b2 = new JButton("CANCEL");
       b1.addActionListener(this);
       b2.addActionListener(this);

        label1 = new JLabel("   Click the pName to delete");
        label1.setFont(new Font("Serif",Font.BOLD, 13));

       createRepPanel.add(label1, BorderLayout.NORTH);
       createRepPanel.add(b1, BorderLayout.CENTER);
       createRepPanel.add(b2, BorderLayout.EAST);
      getContentPane().add(createRepPanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 80);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        repFrame.dellRep();
//        this.setVisible(false); 
      }else if(e.getSource() == b2){
        this.setVisible(false);       
        repFrame.dellRepCloseAct();
      }
   }
 
}






