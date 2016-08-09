/*
 * Main.java
 *
 * Created on 17 de agosto de 2006, 12:29
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package bluecaterpillar;

import gui.*;

import javax.swing.JFrame;
import javax.swing.UIManager;

/**
 *
 * @author julian
 */
public class Main {
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            UIManager.setLookAndFeel(
                UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) { 
            System.out.println(e.getMessage());
        }        
        
        // TODO code application logic here
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {

                BlueCaterpillarUI ui = new BlueCaterpillarUI();
                ui.pack();
                //ui.setResizable(false);
                ui.setVisible(true);
            }
        });        
    }
}
