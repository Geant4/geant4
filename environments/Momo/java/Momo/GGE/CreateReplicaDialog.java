
//example
//8.5

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class CreateReplicaDialog extends JDialog implements ActionListener{
   private JButton b1, b2, b3;
   private JLabel label1, label2;
   private ReplicasFrame repFrame;
   private InsertReplicaDialog insertReplicaDialog;
   public CreateReplicaDialog(ReplicasFrame repFrame){
      super(repFrame,"Create Slicing Translation",false);
      this.repFrame = repFrame;
      JPanel createRepPanel = new JPanel();
       createRepPanel.setLayout(new BorderLayout());
       b1 = new JButton("APPEND");
       b2 = new JButton("INSERT");
       b3 = new JButton("END");
       b1.addActionListener(this);
       b2.addActionListener(this);
       b3.addActionListener(this);
       JPanel labelPanel = new JPanel();
        labelPanel.setLayout(new BorderLayout());
        label1 = new JLabel("  Choose Append or Insert");
        label1.setFont(new Font("Serif",Font.BOLD, 12));
//        label2 = new JLabel(" ");
//        label2.setFont(new Font("Serif",Font.BOLD, 12));
        labelPanel.add(label1, BorderLayout.NORTH);
//        labelPanel.add(label2, BorderLayout.CENTER);
       createRepPanel.add(labelPanel, BorderLayout.NORTH);
       createRepPanel.add(b3, BorderLayout.EAST);
       createRepPanel.add(b2, BorderLayout.CENTER);
       createRepPanel.add(b1, BorderLayout.WEST);
      getContentPane().add(createRepPanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 100);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        repFrame.appendRep(); 
      }else if(e.getSource() == b2){
        insertReplicaDialog = new InsertReplicaDialog(repFrame);
        insertReplicaDialog.setVisible(true);
      }else if(e.getSource() == b3){
        this.setVisible(false);      
        repFrame.createRepCloseAct();
      } 
   }    
 
}





