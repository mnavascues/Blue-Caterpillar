/*
 * PopulationManager.java
 *
 * Created on 2 de agosto de 2006, 19:23
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package bluecaterpillar;


import gui.Results;
import java.io.*;
import java.util.*;



/**
 * A class that load input data froma a file and makes the statics
 * @author julian
 */
public class PopulationManager {
    
    /**
     * Contains the population, name and haplotype of each individual expressed in integers
     */
    
    private Thread wt = null;
    
    private int [][] data = null;
    
    private int [] population_sizes = null;
    
    private boolean [] is_reference = null;
    
    private double [][][] freqs = null;
    
    private double [][] meanSizes = null;
    
    int rows = 0, columns = 0, max_allele = 0, nlocus = 0;    
    
    /**
     * Counts the different populations to assign a unique value to each (id)
     */
    private int population_counter = 0;

    /**
     * Counts the different individuals to assign a unique value to each (id)
     */
    private int individual_counter = 0;

    /**
     * Hash to get a "population name" from a id
     */
    private Hashtable p2n;

    /**
     * Hash to get a id from a "population name"
     */
    private Hashtable n2p;

    /**
     * Hash to get a "individual name" from a id
     */
    private Hashtable i2n;

    /**
     * Hash to get a id from a "individual name"
     */
    private Hashtable n2i;
    
    Hashtable haplotype_freqs;    

    
    private double absentAlleleFreq;
    private double rate;
    private boolean useConstant;
    
    public int times;
    public double rarity;
    public boolean distances, goldstein, slambda, likelihood, conformity;
    public String [] populations_X, populations_Y;    
    Results dialogResults;

    public Object [][] goldsteinTable = null;
    public Object [][] slambdaTable = null;
    public Object [][] likelihoodTable = null;
    public Object [][] conformityTable = null;
    public Object [][] neisStandarTable = null;
    public Object [][] neisMinimumTable = null;
    public Object [][] neisDaTable = null;
    public Object [][] chordTable = null;
    
    /**
     * Creates a new instance of PopulationManager
     * Loads all data from a file
     */

    public PopulationManager(String filename) throws Exception {
        
        p2n = new Hashtable();
        n2p = new Hashtable();
        i2n = new Hashtable();
        n2i = new Hashtable();
        haplotype_freqs = new Hashtable(); 
        
        int population_counter_prev = 0;
        boolean counted = false;
        String str;

        try {
            File f = new File(filename);
            BufferedReader br = new BufferedReader (new FileReader(f));
            Hashtable h = new Hashtable();
            
            while ((str = br.readLine()) != null){
                StringTokenizer tok = new StringTokenizer(str);
                String population_name = new String("");
                
                if (str.length() > 1)
                    population_name = new String(tok.nextToken());

                
                if (h.get(population_name) == null){
                    population_counter_prev++;
                    h.put(population_name,population_name);
                }
                
                if (!counted){
                    columns = Utils.getNumberTokens(str);
                    counted = true;
                }
                
                if (columns == Utils.getNumberTokens(str))
                    rows++;                
            }
            
        } catch (Exception e) {
            throw (new Exception("Error calculating data matrix size: " + e.getMessage()));
        }        
        
        // nuevo !!
        data = new int [rows * 2][columns];            
        
        population_sizes = new int [population_counter_prev + 1];
        is_reference = new boolean [population_counter_prev + 1]; 
        
        for(int k = 0; k < is_reference.length; k++)
            is_reference[k] = true;
        
        int i = 0, j = 0;
        
        try {
            File f = new File(filename);
            BufferedReader br = new BufferedReader (new FileReader(f));

            while ((str = br.readLine()) != null){

                int contadero = Utils.getNumberTokens(str);
                
                if (contadero != columns){
                    System.out.println("Ignoring line: " + str);
                }else {
               
                    j = 0;
                    StringTokenizer tokens = new StringTokenizer(str);
                    String population_name = tokens.nextToken();
                    String individual_name = tokens.nextToken();

                    // if population exists takes it, if not, creates new one
                    if (p2n.get(population_name) == null){
                        p2n.put(population_name, population_counter);
                        n2p.put(population_counter, population_name);       
                        population_sizes[population_counter]++; // population size inc.
                        population_counter++;                    
                    }else
                        population_sizes[((Integer) p2n.get(population_name))]++;

                    i2n.put(individual_name, individual_counter);
                    n2i.put(individual_counter, individual_name);                    

                    data[i][j] = ((Integer) p2n.get(population_name));
                    data[i][++j] = individual_counter;

                    while (tokens.hasMoreTokens()){ 
                        String locus = tokens.nextToken();

                        if (max_allele < Integer.valueOf(locus))
                            max_allele = Integer.valueOf(locus);

                        data[i][++j] = Integer.valueOf(locus);
                    }
                    individual_counter++;

                    System.out.println("individual_counter es " + individual_counter);
                    i++;
                }
            }



            // number of locus
            nlocus = columns - 2;
            // nuevo !!
            freqs = new double [population_counter + 1][nlocus][max_allele];
            meanSizes = new double [population_counter + 1][nlocus];           
            
        } catch (Exception e) {
            throw (new Exception("Wrong format? Line: "+ (i+1) + " Column: "+ j + " " + e.getMessage()));
        }
        
        for (i = 0; i < population_counter + 1; i++){
            System.out.println("Tamaño: " + population_sizes[i]);
        }
        
        this.calculateAlleleFreqs(0, rows);
        this.meanSize(0, rows);        
    }
    
    public String [] getPopulationNames(){
        String [] names = new String[population_counter];
       
        for (int  i= 0; i < population_counter; i++)
            names[i] = (String) n2p.get(i);
        
        return names;
    }
    
    public void setAsReference(String name){
        int id = (Integer) p2n.get(name);
        
        is_reference[id] = true;
        
        System.out.println("Poblacion " + id + " ahora es de referencia");
    }
    
    public void unsetAsReference(String name){
        int id = (Integer) p2n.get(name);
        
        is_reference[id] = false;
        
        System.out.println("Poblacion " + id + " YA NO ES de referencia");
    }    
    
    public void calculateHaplotypeFreqs2()
    {
        /*System.out.println(" ");
        for (int i = rows; i < rows + population_sizes[population_counter]; i++){
            for (int j =0; j < columns; j++)
                System.out.print(data[i][j] + " ");
            System.out.println(" ");
        }*/
        
        for (Enumeration e = haplotype_freqs.elements() ; e.hasMoreElements() ;) {
            Haplotype h = (Haplotype) e.nextElement();
            h.setZero(population_counter);
        }        
        
        Hashtable temp = new Hashtable ();

        for (int individual = rows; individual < population_sizes[population_counter] + rows; individual++){
    
            // constructs string haplotype
            String haplotype = new String("");
            for (int locus = 2; locus < columns; locus++)  // iteration over each locus
                haplotype = new String(String.valueOf(data[individual][locus]) + "." + haplotype);
    
            if (temp.get(haplotype) == null) {

                Haplotype h = ((Haplotype) haplotype_freqs.get(haplotype));

                for (int i = individual; i < population_sizes[population_counter] + rows; i++){
                    boolean equal = true;

                    for (int locus = 2; locus < columns; locus++){
                        if (data[individual][locus] != data[i][locus]){  
                            equal = false;
                            break;
                        }
                    }
                    if (equal)
                        h.increment(data[i][0]);
                }

                h.normalize(population_sizes[population_counter]);
                haplotype_freqs.put(haplotype, h);                
                temp.put(haplotype, haplotype);
                
                //System.out.println(h.toString());                
            }
        }    
        //System.out.println("***********************************************");
    }
    
    public void calculateHaplotypeFreqs()
    {
        for (int individual = 0; individual < rows; individual++){
            
            // constructs string haplotype
            String haplotype = new String("");
            for (int locus = 2; locus < columns; locus++)  // iteration over each locus
                haplotype = new String(String.valueOf(data[individual][locus]) + "." + haplotype);
    
            if (haplotype_freqs.get(haplotype) == null){
                
                Haplotype h = new Haplotype(haplotype, population_counter);

                //int haplo_pop_counter = 0;
                int haplo_ref_counter = 0;                
                
                System.out.println("Encontrado el haplo: " + haplotype);
                
                for (int i = individual; i < rows; i ++){
                    boolean equal = true;
                    
                    for (int locus = 2; locus < columns; locus++){
                        if (data[individual][locus] != data[i][locus]){  
                            equal = false;
                            break;
                        }
                    }
                    if (equal){
                        h.increment(data[i][0]);
                        if (is_reference[ data[i][0] ]){
                            System.out.println("Icrementando. Pobacion: " + data[i][0] + "individuo " + i);
                            for(int k=0;k<is_reference.length;k++)
                                System.out.println(is_reference[k] + " ");
                            haplo_ref_counter++;
                        }
                    }
                }
                
                h.normalize(population_sizes, is_reference, haplo_ref_counter);
                
                //System.out.println(h.toString());
                
                haplotype_freqs.put(haplotype, h);
            
            }    

        }
    }

    public int getRows() { return rows; }
    
    public void calculateAlleleFreqs(int begin, int end)
    {
        /*System.out.println(" ");
        for (int i = begin; i < end; i++){
            for (int j =0; j < columns; j++)
                System.out.print(data[i][j] + " ");
            System.out.println(" ");
        }*/
   
        boolean [][][] checked = new boolean [population_counter + 1][nlocus][max_allele];
        
        for (int locus = 2; locus < columns; locus++){ // iteration over each locus
            for (int individual = begin; individual < end; individual++){
                int population = data[individual][0];
                int locus_value = data[individual][locus];
                int locus_counter = 1;
                
                // if we havent calculated freq yet, calcule it for this value of locus
                if (checked[population][locus-2][locus_value-1] == false){
                    for(int i = individual+1; i < end; i++){
                        // if they are in the same population and have the same locus
                        if (data[i][0] == population && data[i][locus] == locus_value) { 
                            locus_counter++;
                        }
                    }
        
                    freqs[population][locus-2][locus_value-1] = locus_counter * 1.0 / population_sizes[population];
                    checked[population][locus-2][locus_value-1] = true;
                }
                
            }
        }

           /* 
        //for (int i = 0; i < population_counter + 1; i++){
            System.out.println("\n" + "POBLACION "+ population_counter + "------------------------\n");
            for (int j = 0; j < nlocus; j++){
                for (int k=0; k < max_allele; k++){
                    System.out.print("Locus " +j+ " Valor " + (k+1) + " " + freqs[population_counter][j][k] + " *** ");
                }
                System.out.println(" ");
            }
        //} */
    }
    
    public void calculateAlleleFreqs()
    {
        /*System.out.println(" ");
        for (int i = rows; i < rows + population_sizes[population_counter]; i++){
            for (int j =0; j < columns; j++)
                System.out.print(data[i][j] + " ");
            System.out.println(" ");
        }*/
   
        for (int locus_value = 1; locus_value <= max_allele; locus_value++){
            for (int locus = 2; locus < columns; locus++) {
                int locus_counter = 0;
                for (int individual = rows; individual < rows + population_sizes[population_counter]; individual++){
                    if (data[individual][locus] == locus_value)
                        locus_counter++;
                }
                    
                freqs[population_counter][locus-2][locus_value-1] = locus_counter * 1.0 / population_sizes[population_counter];
                   
            }
        }
           /* 
        //for (int i = 0; i < population_counter + 1; i++){
            System.out.println("\n" + "POBLACION "+ population_counter + "------------------------\n");
            for (int j = 0; j < nlocus; j++){
                for (int k=0; k < max_allele; k++){
                    System.out.print("Locus " +j+ " Valor " + (k+1) + " " + freqs[population_counter][j][k] + " *** ");
                }
                System.out.println(" ");
            }
*/
    }    
    
    
    
    public void meanSize(int begin, int end) {
  /*      
        System.out.println(" ");
        for (int i = rows; i < rows + population_sizes[population_counter]; i++){
            for (int j =0; j < columns; j++)
                System.out.print(data[i][j] + " ");
            System.out.println(" ");
        }
    */    
        for (int locus = 2; locus < columns; locus++) {        

            boolean [] checked = new boolean [population_counter + 1];
            
            for (int i = begin; i < end; i++){
                
                int counter = 0;
                int pop = data[i][0];

                if (! checked[pop]){
                    for (int j = i; j < end; j++){
                        if (data[j][0] == pop){
                            counter += data[j][locus];
                        }
                    }

                    checked[pop] = true;
                    meanSizes[pop][locus-2] = counter * 1.0 / population_sizes[pop];
                }
            }
        }
      
        /*
       for (int i = 0; i < population_counter + 1; i++){
            for (int j = 0; j < nlocus; j++){
                System.out.println("poblacion " + i + " locus " + j + " media " + meanSizes[i][j]);
            }
       }
*/
    }
    
    public double [] distances (int population_X, int population_Y) {
        
        // Allele based ***********************************************
        
        double Jxy = 0, Jx = 0, Jy = 0, SqrtXY = 0;
        double Jxy_total = 0, Jx_total = 0, Jy_total = 0, SqrtXY_total = 0, chord = 0;
        
        for (int locus = 0; locus < nlocus; locus++) {        
            Jxy = Jx = Jy = SqrtXY = 0;
            for (int allele = 0; allele < max_allele; allele++){
                Jxy += freqs[population_X][locus][allele] * freqs[population_Y][locus][allele];
                Jx += freqs[population_X][locus][allele] * freqs[population_X][locus][allele];
                Jy += freqs[population_Y][locus][allele] * freqs[population_Y][locus][allele];
                SqrtXY += Math.sqrt(freqs[population_X][locus][allele] * freqs[population_Y][locus][allele]);
            }
            
            Jxy_total += Jxy;
            Jx_total += Jx;
            Jy_total += Jy;
            SqrtXY_total += SqrtXY;
            
            chord += Math.sqrt(2 * (1.0 - SqrtXY));
                        
            //System.out.println("Locus " + locus + "  Jxy= " + Jxy + " Jx= " + Jx + " Jy= " + Jy + " SqrtXY= " + SqrtXY);
        }
        
        Jxy_total = Jxy_total / nlocus; 
        Jx_total = Jx_total / nlocus;
        Jy_total = Jy_total / nlocus;
        SqrtXY_total = SqrtXY_total / nlocus;
        chord = chord / nlocus;
        
        double Ds = -1.0 * Math.log(Jxy_total / Math.sqrt(Jx_total * Jy_total));
        double Dm = ((Jx_total + Jy_total) / 2.0)  - Jxy_total;
        double Da = 1.0 - SqrtXY_total;
        double Dc = 2.0 / 3.1416 * chord; //Math.PI  * chord;

        //System.out.println("Ds= " + Ds + " Dm= " + Dm + " Da= " + Da + " Dc= " + Dc );        

        double [] return_value = new double [4];
        return_value[0] = Ds; return_value[1] = Dm; return_value[2] = Da; return_value[3] = Dc;
        
        return return_value;
    }
    
    
    public double likelihood(int population_X, int population_Y, double rate){
        
        int end = rows;
        
        if (population_Y == population_counter)
            end = rows + population_sizes[population_counter];
        
        double likelihood = 1;
        for (int individual = 0; individual < end; individual++){
         
            if (data[individual][0] == population_Y){
                // constructs string haplotype
                String haplotype = new String("");
                for (int locus = 2; locus < columns; locus++)  // iteration over each locus
                    haplotype = new String(String.valueOf(data[individual][locus]) + "." + haplotype);                        
                
                Haplotype h = ((Haplotype) haplotype_freqs.get(haplotype));
                
                double freq_h_X = h.getFreqInPop(population_X);
                double freq_h_TOT = h.getFreqTotalRef();
                String haplotype_sup = new String("");
                String haplotype_inf = new String("");                
                String temp_sup, temp_inf;
                Haplotype h_sup, h_inf;
                double sum_x = 0, sum_tot = 0;
                for (int locus_m = 2; locus_m < columns; locus_m++){
                    for (int locus = 2; locus < columns; locus++){
                        
                        if (locus == locus_m){
                            temp_sup = String.valueOf(data[individual][locus] + 1);
                            temp_inf = String.valueOf(data[individual][locus] - 1);
                        }else{
                            temp_sup = String.valueOf(data[individual][locus]);
                            temp_inf = temp_sup;
                        }
                        
                        haplotype_sup = new String(temp_sup + "." + haplotype_sup);                                        
                        haplotype_inf = new String(temp_inf + "." + haplotype_inf);                                        
                    }
                    
                    h_sup = ((Haplotype) haplotype_freqs.get(haplotype_sup));
                    h_inf = ((Haplotype) haplotype_freqs.get(haplotype_inf));
                    
                    if (h_sup != null){ 
                        sum_x += h_sup.getFreqInPop(population_X);
                        sum_tot += h_sup.getFreqTotalRef();        
                    }
                    if (h_inf != null){ 
                        sum_x +=  h_inf.getFreqInPop(population_X);
                        sum_tot += h_inf.getFreqTotalRef();        
                    }
                }
                
                
                
                double numerator = (1-rate)*freq_h_X + rate*(sum_x);
                double denominator = (1-rate)*freq_h_TOT + rate*(sum_tot);
                
                
                System.out.println("Haplo: " + haplotype + " Frecuencias: " + h.getFreqInPop(population_X) + " " + h.getFreqTotalRef());
                
                if (numerator == 0){
                    if (useConstant)
                        numerator = absentAlleleFreq;
                    else
                        numerator = 1 / population_sizes[population_X];
                }
                
                if (denominator == 0){
                    if (useConstant)
                        denominator = absentAlleleFreq;
                    else
                        denominator = 1 / getRefsSize();
                }                
                
                if (numerator == 0 || denominator  == 0) 
                    return 0; // likelihood = 0 !!
        
                
                /* MODIFICACION AQUI ¿? Y NO OLVIDES EL COPY */
                 
                
                likelihood *= numerator / denominator;
            }
        }
        
        //System.out.println("Likelihood: "+ likelihood);
        
        return likelihood;
    }

    public int getRefsSize(){
    
        int size = 0;
        for (int i = 0; i < is_reference.length; i++){
            if (is_reference[i]){
                size += population_sizes[i];
            }
        }
        
        return size;
    }
    
    public double conformity(int population_X, int population_Y, double rarity)
    {
        int end = rows;
        
        if (population_Y == population_counter)
            end = rows + population_sizes[population_counter];        
        
        Hashtable tempura = new Hashtable();
        double conformity = 1;
        double value;
        for (int individual = 0; individual < end; individual++){
            
            if (data[individual][0] == population_Y){
                String haplotype = new String("");         
                for (int locus = 2; locus < columns; locus++)  // iteration over each locus
                    haplotype = new String(String.valueOf(data[individual][locus]) + "." + haplotype);                        
            
                if (tempura.get(haplotype) == null){
                    Haplotype h = ((Haplotype) haplotype_freqs.get(haplotype));
                    double freq = h.getFreqInPop(population_X);
                    
                    if (freq < rarity)
                        value = 1 - Math.pow((1 - freq), population_sizes[population_Y]);
                    else
                        value = 1;
                    
                    conformity *= value;
                    tempura.put(haplotype, true);
                }
            }                
        }   
        
        //System.out.println(" Conformity: " + conformity);        
        
        return conformity;
    }
    
    
    public double SLambda(int population_X, int population_Y){
        int end = rows;
        
        if (population_Y == population_counter)
            end = rows + population_sizes[population_counter];             
        
        
        double Slambda = 0;
        Hashtable temp = new Hashtable();
        for (int individual = 0; individual < end; individual++){
            
            if (data[individual][0] == population_Y || data[individual][0] == population_X){
                String haplotype = new String("");         
                for (int locus = 2; locus < columns; locus++)  // iteration over each locus
                    haplotype = new String(String.valueOf(data[individual][locus]) + "." + haplotype);                        
            
                if (temp.get(haplotype) == null){
                    Haplotype h = ((Haplotype) haplotype_freqs.get(haplotype));
                    Slambda += Math.pow(h.getFreqInPop(population_X) - h.getFreqInPop(population_Y), 2);
                    temp.put(haplotype, true);
                }
            }                
        }        
        //System.out.println(" Slambda: " + Slambda);
        
        return Slambda;
    }
    
    public double goldstein(int population_X, int population_Y)
    {
        double Goldstein = 0;

        for (int locus = 0; locus < nlocus; locus++)
            Goldstein += Math.pow(meanSizes[population_X][locus] - meanSizes[population_Y][locus], 2.0);

        
        Goldstein /= nlocus;
        
        
        //System.out.println(" Goldstein " + Goldstein);
        
        return Goldstein;
    }
    
    
    public void resampling(int population_X, int population_Y, double rate){

        // elimina haplotipos simulados 
        for (Enumeration e = haplotype_freqs.elements() ; e.hasMoreElements() ;) {
            Haplotype h = (Haplotype) e.nextElement();
            if (h.getSimulated()){
                haplotype_freqs.remove(h.getHaplo());
                h = null;
            }
        } 
        
        
        int [] individuals = new int [population_sizes[population_X]];
        for (int i = 0, j = 0; i < rows; i++){
            if (data[i][0] == population_X){
                individuals[j] = data[i][1];
                j++;
            }
        }

        int min = 0, max = population_sizes[population_X] - 1;
        
        for (int i = 0; i < population_sizes[population_Y]; i++){
            double aleat = Math.random() * ((max + 1) - min);
            aleat = Math.floor(aleat);
            aleat += min;                    
            
            int id = individuals[((int) aleat)];
            
            data[rows + i][0] = population_counter;
            data[rows + i][1] = rows + i;
            
            String haplotype = new String("");
                
            for (int locus = 2; locus < columns; locus++){
                aleat = Math.random();
                
                if (aleat >= rate){
                    data[rows + i][locus] = data[id][locus];
                }else{ 
                    
                    if (aleat < rate && aleat >= rate/2)
                        data[rows + i][locus] = data[id][locus] + 1;
                    else
                        data[rows + i][locus] = data[id][locus] - 1;
                }
                
                haplotype = String.valueOf(data[rows + i][locus]) + "." + haplotype;
            }
            
            Haplotype h = ((Haplotype) haplotype_freqs.get(haplotype));
            
            if (h == null){
                h = new Haplotype(haplotype, population_counter);
                h.setSimulated(true);
                haplotype_freqs.put(haplotype, h); 
            }
        /* MODIFICACION AQUI ¿? Y NO OLVIDES EL COPY */
        
        }
        
        
        population_sizes[population_counter] = population_sizes[population_Y];
        
        /*System.out.println(" ");
        for (int i = 0; i < rows + population_sizes[population_Y]; i++){
            for (int j =0; j < columns; j++)
                System.out.print(data[i][j] + " ");
            System.out.println(" ");
        } */       
    }
    
    
    
    class WorkingThread extends Thread {

        int progress = 0;

        
        private Object [][] makeTable()
        {
                int resultRows = populations_Y.length + 2;
                int resultCols =  1 + (populations_X.length * 2);        
            
                Object [][] table = new Object [resultRows][resultCols];
                
                for (int i = 0; i < populations_Y.length; i++)
                    table[i + 2][0] = populations_Y[i];

                for (int i = 1; i <= populations_X.length; i++){
                    table[0][i] = populations_X[i-1];
                    table[0][i + populations_X.length] = populations_X[i-1];
                    table[1][i] = new String(" Value");
                    table[1][i + populations_X.length] = new String (" %");
                }           
                
                return table;
        }

        public void run ()
        {
            calculateHaplotypeFreqs();        
            
            goldsteinTable = makeTable();
            slambdaTable = makeTable();
            likelihoodTable = makeTable();
            conformityTable = makeTable();
            neisStandarTable = makeTable();
            neisMinimumTable = makeTable();
            neisDaTable = makeTable();
            chordTable = makeTable();
            
            
            for (int i = 0; i < populations_Y.length; i++){
                for (int j = 0; j < populations_X.length; j++) {

                    int poblacion_Y = (Integer) p2n.get(populations_Y[i]);
                    int poblacion_X = (Integer) p2n.get(populations_X[j]);

                    System.out.println("************************************");
                    System.out.println("X: " + populations_X[j] + " Y: " + populations_Y[i] );                    
                    
                    double [] dist = null;
                    double g = 0, s = 0, l = 0, c = 0;

                    double [] dist2 = null;
                    double g2 = 0, s2 = 0, l2 = 0, c2 = 0;

                    double counter_Ds = 0, counter_Dm = 0, counter_Da = 0, counter_Dc = 0;
                    double counter_g = 0, counter_s = 0, counter_l = 0, counter_c = 0;

                    if (distances)  dist = distances(poblacion_X, poblacion_Y);                
                    if (goldstein)  g = goldstein(poblacion_X, poblacion_Y);                
                    if (slambda)    s = SLambda(poblacion_X,  poblacion_Y);
                    if (likelihood) l = likelihood(poblacion_X, poblacion_Y, rate);
                    if (conformity) c = conformity(poblacion_X, poblacion_Y, rarity);
                    

                    //System.out.println(l + " " + g + " " + s + " " + c);
                    

                    int iterations = times;
                    //iterations = 0;
                    while (iterations > 0) {
                        //System.out.println(iterations + " iterations remaining...");
                        g2 = s2 = l2 = c2 = 0;

                        resampling(poblacion_X, poblacion_Y, rate);

                        if (distances) {
                            calculateAlleleFreqs();

                            dist2 = distances(poblacion_X, population_counter);

                            if (dist2[0] >= dist[0]) counter_Ds++;
                            if (dist2[1] >= dist[1]) counter_Dm++;
                            if (dist2[2] >= dist[2]) counter_Da++;
                            if (dist2[3] >= dist[3]) counter_Dc++;
                        }

                        if (slambda || likelihood || conformity){
                            calculateHaplotypeFreqs2();
                        }

                        if (goldstein) {
                            meanSize(rows, rows + population_sizes[population_counter]); 
                            g2 = goldstein(poblacion_X, population_counter);
                            if (g2 >= g) counter_g++;
                        }                    

                        if (slambda)  {
                            s2 = SLambda(poblacion_X, population_counter);
                            if (s2 >= s) counter_s++;
                        }

                        if (likelihood){
                            l2 = likelihood(poblacion_X, population_counter, rate);
                            if (l2 <= l) counter_l++;
                        }                


                        if (conformity) {
                            c2 = conformity(poblacion_X, population_counter, rarity);
                            if (c2 <= c) counter_c++;
                        }
                        iterations--;

                        dialogResults.getBar().setValue(++progress);
                    }            

                    counter_Ds /= times;
                    counter_Dm /= times;
                    counter_Da /= times;
                    counter_Dc /= times;
                    counter_g /= times;
                    counter_s /= times;
                    counter_l /= times;
                    counter_c /= times;        

                    System.out.println("Ds " + counter_Ds);
                    System.out.println("Dm " + counter_Dm);        
                    System.out.println("Da " + counter_Da);
                    System.out.println("Dc " + counter_Dc);                
                    System.out.println("Goldstein " + counter_g);                        
                    System.out.println("S-Lambda " + counter_s);                                
                    System.out.println("Likelihood " + counter_l);                                
                    System.out.println("Conformity " + counter_c);

                    
                    
                    goldsteinTable[i+2][j+1] = String.format( "%g", g);
                    goldsteinTable[i+2][j+1+populations_X.length] = String.format( "%g", counter_g);
                    
                    slambdaTable[i+2][j+1] = String.format( "%g", s);;
                    slambdaTable[i+2][j+1+populations_X.length] = String.format( "%g", counter_s);
                    
                    likelihoodTable[i+2][j+1] = String.format( "%g", l);
                    likelihoodTable[i+2][j+1+populations_X.length] = String.format( "%g", counter_l);

                    conformityTable[i+2][j+1] = String.format( "%g", c);
                    conformityTable[i+2][j+1+populations_X.length] = String.format( "%g", counter_c);                    

                    if (distances) {
                        neisStandarTable[i+2][j+1] = String.format( "%g", dist[0]);
                        neisStandarTable[i+2][j+1+populations_X.length] = String.format( "%g", counter_Ds);                    

                        neisMinimumTable[i+2][j+1] = String.format( "%g", dist[1]);
                        neisMinimumTable[i+2][j+1+populations_X.length] = String.format( "%g", counter_Dm);                       

                        neisDaTable[i+2][j+1] = String.format( "%g", dist[2]);
                        neisDaTable[i+2][j+1+populations_X.length] = String.format( "%g", counter_Da);                      

                        chordTable[i+2][j+1] = String.format( "%g", dist[3]);
                        chordTable[i+2][j+1+populations_X.length] = String.format( "%g", counter_Dc);  
                    }
                  
                    
                    
                    /*System.out.println("************************************************************************");
                    System.out.println("Pareja: " + populations_Y[i] + " - " + populations_X[j]);
                    System.out.println("Posicion: " + (i+2) + " - " + (j+1) + " Valor: " + g);
                    System.out.println("Posicion: " + (i+2) + " - " + (j+1+populations_X.length) + " Valor: " + counter_g);
                    */

                    
                    
                    
                    
                    
                    
                    
                }
            }                
            dialogResults.addResults();
            wt = null;
        }
    }
    
    
    public void configure (int times, String []  populations_X, String [] populations_Y, double rarity,
            boolean distances, boolean goldstein, boolean slambda, boolean likelihood, boolean conformity,
            boolean useConstant, double absentAlleleFreq, double rate, Results dialogResults)
    {
        this.times = times;
        this.populations_X = populations_X;
        this.populations_Y = populations_Y;
        this.rarity = rarity;
        this.distances = distances;
        this.goldstein = goldstein;
        this.slambda = slambda;
        this.likelihood = likelihood;
        this.conformity = conformity;      
        this.dialogResults = dialogResults;
        this.useConstant = useConstant;
        this.absentAlleleFreq = absentAlleleFreq;
        this.rate = rate;
        
        if (wt == null){
            wt = new WorkingThread();
            wt.start();
        }
    }
    
    public void stopThread()
    {
        if (wt != null){
            wt.stop();
            wt = null;
        }            
    }        
    
    public void pauseThread()
    {
        if (wt != null)
            wt.suspend();
    }      
    
    public void continueThread()
    {
        if (wt != null)
            wt.resume();
    }      

}
    
    


