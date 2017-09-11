/*
# GAPE - Genetic Algorithm for Parameter Estimation of SIFTER tool
#
# Created by Eng. (C) Julian Camilo Castañeda Alonso and Msc. Tania Andrea Rodriguez Quiñones on August 2017.
# Copyright (c) 2017. Eng. (C) Julian Camilo Castañeda Alonso and Msc. Tania Andrea Rodriguez Quiñones. Universidad Antonio Narino. All rights reserved.
#
# GAPE is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License v3.0 found in the LICENSE file in the root directory of this project.
*/

package gape.output; 

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Eng. (C) Julian Camilo Castañeda Alonso - Msc. Tania Andrea Rodriguez Quiñones
 *
 */

public class PrintFiles {
	
	/**
	 * Default constructor
	 */
	public PrintFiles() {

	}
	
	/**
	 * This method print in file the Transition Matrix
	 * @param output
	 * @param transitionMatrix
	 * @param rowNames
	 */
	public void printTransitionMatrix(String output, double[][] transitionMatrix, int[] rowNames) {
		
		File o = new File(output);
		
		FileWriter fw;
		
		if (!o.exists()) 
        {
            try 
            {
                o.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(PrintFiles.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + o);
                System.exit(0);
            }
        }
		
		BufferedWriter bw = null;
        
        try 
    	{
    		fw = new FileWriter(o.getAbsoluteFile(), true);
    		bw = new BufferedWriter(fw);
    		
    		bw.write("# scale parameter: " + 20.0 + "\n");
    		
    		for (int i = 0; i < transitionMatrix.length; i++) {
				
    			bw.write("\t" + String.valueOf(rowNames[i]));
    			
			}
    		
    		bw.write("\n");
    		
    		for (int i = 0; i < transitionMatrix.length; i++) {
    			
    			bw.write(String.valueOf(rowNames[i]));
    			
    			for (int j = 0; j < transitionMatrix[i].length; j++) {
					
					bw.write("\t" + transitionMatrix[i][j]);
					
				}
    			
    			bw.write("\n");
			}
    		
			bw.flush();
		} 
    	catch (IOException e) 
    	{
			e.printStackTrace();
		}
	}
	
	/**
	 * This method print in file the Alpha Matrix
	 * @param output
	 * @param alpha
	 */
	public void printAlpha(String output, double[] alpha) {
		
		File o = new File(output);
		
		FileWriter fw;
		
		if (!o.exists()) 
        {
            try 
            {
                o.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(PrintFiles.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + o);
                System.exit(0);
            }
        }
		
		BufferedWriter bw = null;
		
        try 
    	{
    		fw = new FileWriter(o.getAbsoluteFile(), true);
    		bw = new BufferedWriter(fw);
    		
    		for (int i = 0; i < alpha.length; i++) {
				bw.write(alpha[i] + "\n");
			}
    		
			bw.flush();
		} 
    	catch (IOException e) 
    	{
			e.printStackTrace();
		}

	}
	
	/**
	 * This method print in file the Sigma Matrix
	 * @param output
	 * @param sigma
	 */
	public void printSigma(String output, double[] sigma) {
		
		File o = new File(output);
		
		FileWriter fw;
		
		if (!o.exists()) 
        {
            try 
            {
                o.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(PrintFiles.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + o);
                System.exit(0);
            }
        }
		
		BufferedWriter bw = null;
		
        try 
    	{
    		fw = new FileWriter(o.getAbsoluteFile(), true);
    		bw = new BufferedWriter(fw);
    		
    		bw.write("species\t" + sigma[0] + "\n");
    		bw.write("duplication\t" + sigma[1] + "\n");
    		
			bw.flush();
		} 
    	catch (IOException e) 
    	{
			e.printStackTrace();
		}
        
	}
	
	public void printBestIndividual(String output, String bestIndividual) {
		
		File o = new File(output);
		
		FileWriter fw;
		
		if (!o.exists()) 
        {
            try 
            {
                o.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(PrintFiles.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + o);
                System.exit(0);
            }
        }
		
		BufferedWriter bw = null;
		
        try 
    	{
    		fw = new FileWriter(o.getAbsoluteFile(), true);
    		bw = new BufferedWriter(fw);
    		
    		bw.write(bestIndividual);    		
			bw.flush();
		} 
    	catch (IOException e) 
    	{
			e.printStackTrace();
		}
        
	}
	
	public void printTimeExecution(String output, String timeExecution) {
		
		File o = new File(output);
		
		FileWriter fw;
		
		if (!o.exists()) 
        {
            try 
            {
                o.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(PrintFiles.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + o);
                System.exit(0);
            }
        }
		
		BufferedWriter bw = null;
		
        try 
    	{
    		fw = new FileWriter(o.getAbsoluteFile(), true);
    		bw = new BufferedWriter(fw);
    		
    		bw.write("\nExecution Time: " + timeExecution);    		
			bw.flush();
		} 
    	catch (IOException e) 
    	{
			e.printStackTrace();
		}
        
	}
	
}
