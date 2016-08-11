/*
 * Results.java
 *
 * Created on 26 de agosto de 2006, 18:22
 */

package gui;

import bluecaterpillar.*;

import javax.swing.*;
import javax.swing.text.html.HTMLEditorKit;
import java.awt.datatransfer.*;
import java.awt.BorderLayout;
import javax.swing.text.*;
import java.io.StringReader;

/**
 *
 * @author  julian
 */
public class Results extends javax.swing.JDialog {
    
    PopulationManager pm = null;
    String html = null;
    
    
    /** Creates new form Results */
    public Results(JFrame parent, PopulationManager pm) {
        super(parent);
       
        this.pm = pm;
       
        initComponents();
        
        //output.setContentType( "text/html" );
        //output.setEditable( false );
    }

    public void setMaximum(int max){
        progressBar.setMaximum(max);
        progressBar.setValue(0);
        progressBar.setStringPainted(true);   
    }
    
    public javax.swing.JProgressBar getBar() { return progressBar; }
    
    private JPanel newPanel(Object [][] data) {
        
        String [] names = new String [1 + pm.populations_X.length * 2];
        
        names [0] = new String(" ");
        for (int i = 1; i <= pm.populations_X.length; i++){
            names[i] = pm.populations_X[i-1];
            names[i + pm.populations_X.length] = pm.populations_X[i-1];
        }          
                
        JPanel panel = new JPanel();
        JTable table = new JTable(data, names);

        panel.add(table, BorderLayout.CENTER);

        return panel;
    }
    
    public void addResults()
    {
        pauseButton.setEnabled(false);
        progressBar.setEnabled(false);
        stateText.setText("Done");
        this.setTitle("Blue Caterpillar");
        cancelButton.setText("Exit");
        
        String [] names = new String [1 + pm.populations_X.length * 2];
        
        names [0] = new String(" ");
        for (int i = 1; i <= pm.populations_X.length; i++){
            names[i] = pm.populations_X[i-1];
            names[i + pm.populations_X.length] = pm.populations_X[i-1];
        }            

        if (pm.goldstein)
        {
            tabs.addTab("Goldstein", new JTable(pm.goldsteinTable, names));
        }
        if (pm.slambda)
            tabs.addTab("S-Lambda", new JTable(pm.slambdaTable, names));            
        if (pm.likelihood)
            tabs.addTab("Likelihood", new JTable(pm.likelihoodTable, names));
        if (pm.conformity)
            tabs.addTab("Conformity", new JTable(pm.conformityTable, names)); 
        if (pm.distances){
            tabs.addTab("Nei's standard distance", new JTable(pm.neisStandarTable, names)); 
            tabs.addTab("Nei's minimum distance", new JTable(pm.neisMinimumTable, names)); 
            tabs.addTab("Nei's Da distance", new JTable(pm.neisDaTable, names)); 
            tabs.addTab("Chord distance", new JTable(pm.chordTable, names)); 
        }
        
        copyButton.setEnabled(true);
    }
    /*public void addResults()
    {
        if (html == null){
            html = new String("<b>Bootstrap: " + pm.times);
            html = new String(html + " Rarity: " + pm.rarity + "</b><br>");
        }
        else
            html = new String(html + "<br>jojo");
        
        output.setText("<html><head></head>" + html + "<body></body></html>");
    }
    
    private String newHtmlTable()
    {
        String table = new String("");
        return new String();
    }*/
    
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc=" Generated Code ">//GEN-BEGIN:initComponents
    private void initComponents() {
        progressBar = new javax.swing.JProgressBar();
        stateText = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        tabs = new javax.swing.JTabbedPane();
        cancelButton = new javax.swing.JButton();
        pauseButton = new javax.swing.JButton();
        copyButton = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DO_NOTHING_ON_CLOSE);
        setTitle("Proccesing...");

        stateText.setText("Please wait");

        jLabel2.setText("Output");

        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cancelButtonActionPerformed(evt);
            }
        });

        pauseButton.setText("  Pause  ");
        pauseButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                pauseButtonActionPerformed(evt);
            }
        });

        copyButton.setText("Copy");
        copyButton.setEnabled(false);
        copyButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                copyButtonActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(stateText)
                    .add(progressBar, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 496, Short.MAX_VALUE)
                    .add(layout.createSequentialGroup()
                        .add(jLabel2)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, 292, Short.MAX_VALUE)
                        .add(pauseButton)
                        .add(16, 16, 16)
                        .add(cancelButton)
                        .add(12, 12, 12))
                    .add(tabs, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 496, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.TRAILING, copyButton))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(stateText)
                .add(7, 7, 7)
                .add(progressBar, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(layout.createSequentialGroup()
                        .add(35, 35, 35)
                        .add(jLabel2))
                    .add(layout.createSequentialGroup()
                        .add(15, 15, 15)
                        .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(pauseButton)
                            .add(cancelButton))))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(tabs, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 189, Short.MAX_VALUE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(copyButton, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 23, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void copyButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_copyButtonActionPerformed

        JTable table = (JTable) tabs.getSelectedComponent();
        
        int rows = table.getRowCount();
        int columns = table.getColumnCount();
        String data = new String();
        
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                String value = (String) table.getValueAt(i,j);
                
                if (value == null) value = new String(" ");
                
                data = new String(data + value  + "\t");
            }
            data = new String(data + "\n");
        }
        
        //System.out.println(data);
        
        Clipboard clipboard = getToolkit().getSystemClipboard();
        StringSelection data_to_copy = new StringSelection(data);
        clipboard.setContents(data_to_copy, data_to_copy);
    }//GEN-LAST:event_copyButtonActionPerformed

    private void pauseButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_pauseButtonActionPerformed
// TODO add your handling code here:
        if (pauseButton.getText().equals("  Pause  "))
        {
            pm.pauseThread();
            pauseButton.setText("Continue ");
        }else{
            pm.continueThread();
            pauseButton.setText("  Pause  ");
        }            
    }//GEN-LAST:event_pauseButtonActionPerformed

    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cancelButtonActionPerformed
// TODO add your handling code here:
        pm.stopThread();
        this.setVisible(false);        
    }//GEN-LAST:event_cancelButtonActionPerformed
    
   
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton cancelButton;
    private javax.swing.JButton copyButton;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JButton pauseButton;
    private javax.swing.JProgressBar progressBar;
    private javax.swing.JLabel stateText;
    private javax.swing.JTabbedPane tabs;
    // End of variables declaration//GEN-END:variables
    
}
