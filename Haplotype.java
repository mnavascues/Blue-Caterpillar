/*
 * Haplotype.java
 *
 * Created on 14 de agosto de 2006, 16:40
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package bluecaterpillar;

/**
 *
 * @author julian
 */
public class Haplotype {
    
    double [] haplotype_freq_in_pop = null;
    double haplotype_freq_in_ref = 0;
    private String haplo = null;
    int nPopulations = 0;
    boolean simulated = false;
    
    /** Creates a new instance of Haplotype */
    public Haplotype(String haplo, int nPopulations) {
        this.haplo = haplo;
        haplotype_freq_in_pop = new double [nPopulations + 1];
        this.nPopulations = nPopulations;
    }

    public double getFreqTotalRef() { return haplotype_freq_in_ref; }
    
    public double getFreqInPop(int population) { return haplotype_freq_in_pop[population]; }
    
    public void setFreqRef(double freq) { this.haplotype_freq_in_ref = freq; }
    
    public void increment(int population){
        haplotype_freq_in_pop[population]++;
    }

    public void setZero(int population){
        haplotype_freq_in_pop[population] = 0;
    }
    
    public void normalize (int pop_size){
        haplotype_freq_in_pop[nPopulations] = haplotype_freq_in_pop[nPopulations] / pop_size;
    }
    
    public void normalize (int [] population_sizes, boolean [] is_reference, int haplo_ref_counter){
        int ref_size = 0;
        
        for (int i = 0; i < nPopulations; i++){
            haplotype_freq_in_pop[i] = haplotype_freq_in_pop[i] / population_sizes[i];

            if (haplo_ref_counter > 0 && is_reference[i])
                ref_size += population_sizes[i];

        }
        
        if (haplo_ref_counter > 0){
            System.out.println("Dividiendo: " + haplo_ref_counter + "/" + ref_size);
            haplotype_freq_in_ref = haplo_ref_counter * 1.0 / ref_size;
        }            
    }
    
    public String toString()
    {
        String pop_freqs = new String();
        for (int i=0; i < nPopulations + 1; i++)
            pop_freqs = new String(pop_freqs + " " + haplotype_freq_in_pop[i]);
            
        return new String(getHaplo() + " " + pop_freqs + " Ref: " + haplotype_freq_in_ref);
    }
    
    public void setSimulated(boolean value){
        simulated = value;
    }
    
    public boolean getSimulated(){
        return simulated;
    }

    public String getHaplo() {
        return haplo;
    }


}
