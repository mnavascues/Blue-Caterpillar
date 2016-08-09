/*
 * Utils.java
 *
 * Created on 16 de agosto de 2006, 13:43
 *
 * To change this template, choose Tools | Options and locate the template under
 * the Source Creation and Management node. Right-click the template and choose
 * Open. You can then make changes to the template in the Source Editor.
 */

package bluecaterpillar;

import java.io.*;
import java.util.*;
/**
 *
 * @author julian
 */
public class Utils {
    
    /** Creates a new instance of Utils */
    public Utils() {}
    
    static int getNumberTokens(String str){
        StringTokenizer toffee = new StringTokenizer(str);
        int contadero = 0;
        while (toffee.hasMoreTokens())
        {
            toffee.nextToken();
            contadero++;
        }
        return contadero;
    }
}