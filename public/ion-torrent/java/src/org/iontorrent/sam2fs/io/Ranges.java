/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2fs.io;

import java.io.*;
import java.util.*;
import net.sf.picard.reference.*;
import net.sf.samtools.*;

/**
 * Stores a list of genomic intervals (ranges).
 *
 * Note: ranges should be specified in the following format:
 * <name>:<start>:<end>
 * where <name> is the sequence name, <start> is the one-based
 * start position, and <end> is the one-based end position (inclusive).
 *
 * @author nils.homer@lifetech.com
 */
public class Ranges {

    /**
     * The list of genomic ranges, in genomic sorted order.
     */
    private LinkedList<Range> ranges = null;
    
    /**
     * 
     */
    public Ranges() {
        
    }
    /**
     * Creates a list of ranges for all sequences in the reference dictionary.
     * @param referenceDictionary the reference dictionary.
     */
    public Ranges(SAMSequenceDictionary referenceDictionary)
    {
        int i;
        this.ranges = new LinkedList<Range>();
        for(i=0;i<referenceDictionary.size();i++) {
            this.ranges.add(new Range(i, 1, referenceDictionary.getSequence(i).getSequenceLength()));
        }
    }

    /**
     * Creates a list of ranges given the intervals in the given file.  The intervals are extended
     * on either end given the offset.
     * @param file the input file.
     * @param referenceDictionary the reference dictionary.
     * @param offset the offset.
     */
    public Ranges(File file, SAMSequenceDictionary referenceDictionary, int offset)
    {
        BufferedReader br = null;
        String line = null;
        int i, lineNumber = 1;
        Map<String, Integer> hm = new HashMap<String, Integer>();

        try {
            // open
            br = new BufferedReader(new FileReader(file));

            // init
            this.ranges = new LinkedList<Range>();
            for(i=0;i<referenceDictionary.size();i++) {
                hm.put(referenceDictionary.getSequence(i).getSequenceName(), new Integer(i));
            }

            // read the file
            while(null != (line = br.readLine())) {
                this.addRange(line, lineNumber, hm, referenceDictionary, offset);
                lineNumber++;
            }

            // close
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Creates a list of one range given the input range.  The interval is extended
     * on either end given the offset.
     * @param range the input range.
     * @param referenceDictionary the reference dictionary.
     * @param offset the offset.
     */
    public Ranges(String range, SAMSequenceDictionary referenceDictionary, int offset)
        throws Exception
    {
        Map<String, Integer> hm = new HashMap<String, Integer>();
        int i;
        int startPosition, endPosition;

        // init
        this.ranges = new LinkedList<Range>();
        for(i=0;i<referenceDictionary.size();i++) {
            hm.put(referenceDictionary.getSequence(i).getSequenceName(), new Integer(i));
        }
                
        this.addRange(range, -1, hm, referenceDictionary, offset);

    }

    /**
     * Creates a list of one range given the input range.  
     * @param range the input range.
     * @param referenceDictionary the reference dictionary.
     */
    public Ranges(String range, SAMSequenceDictionary referenceDictionary)
        throws Exception
    {
        this(range, referenceDictionary, 0);
    }

    /**
     * Parses the given input string and adds it to the list of ranges.
     * @param line the input range string.
     * @param lineNumber the line if the range string is from a file, negative otherwise.
     * @param hm the map of sequence names in the dictionary to give the genomic ordering.
     * @param referenceDictionary the reference dictionary.
     * @param offset the offset.
     */
    private void addRange(String line, int lineNumber, Map<String, Integer> hm, SAMSequenceDictionary referenceDictionary, int offset)
        throws Exception
    {
        int referenceIndex=-1;
        int startPosition=-1, endPosition=-1;
        int i=0;
        int colon, dash;

        String exceptionEndMessage = null;
        if(lineNumber < 0) {
            exceptionEndMessage = new String(".");
        }
        else {
            exceptionEndMessage = new String(" on line " + lineNumber + ".");
        }
        
        // get delimiters
        colon = line.indexOf(':');
        dash = line.indexOf('-', colon+1);

        // Only chromosome name
        if(colon < 0) { // referenceIndex only
            referenceIndex = (int)hm.get(line);
            startPosition = 1;
            endPosition = referenceDictionary.getSequence(referenceIndex).getSequenceLength();
        }
        else {
            if(colon <= 0 || (0 < dash && (dash - colon) <= 1) || (0 < dash && (line.length() - dash) <= 1)) {
                throw new Exception("RANGE was improperly specified" + exceptionEndMessage);
            }   
                
            String chrName = line.substring(0, colon);
            if(null == hm.get(chrName)) {
                throw new Exception("Could not find reference name ["+chrName+"] in RANGES" + exceptionEndMessage);
            }
            referenceIndex = (int)hm.get(chrName);

            if(dash <= 0) {
                startPosition = Integer.parseInt(line.substring(colon+1));
                endPosition = referenceDictionary.getSequence(referenceIndex).getSequenceLength();
            }
            else {
                startPosition = Integer.parseInt(line.substring(colon+1, dash));
                endPosition = Integer.parseInt(line.substring(dash+1));
            }
        }

        if(startPosition <= 0 || referenceDictionary.getSequence(referenceIndex).getSequenceLength() < startPosition) {
            throw new Exception("startPosition was out of bounds in RANGE" + exceptionEndMessage);
        }
        else if(endPosition <= 0 || referenceDictionary.getSequence(referenceIndex).getSequenceLength() < endPosition) {
            throw new Exception("endPosition was out of bounds in RANGE" + exceptionEndMessage);
        }
        else if(endPosition < startPosition) {
            throw new Exception("endPosition < startPosition in RANGE" + exceptionEndMessage);
        }

        startPosition -= offset;
        if(startPosition <= 0) {
            startPosition = 1;
        }
        endPosition += offset;
        if(referenceDictionary.getSequence(referenceIndex).getSequenceLength() < endPosition) {
            endPosition = referenceDictionary.getSequence(referenceIndex).getSequenceLength();
        }

        if(null != this.ranges.peek()) {
            Range last = this.ranges.getLast();
            if(referenceIndex < last.referenceIndex) {
                throw new Exception("Ranges must be in sorted order in RANGES (line numbers " + (lineNumber-1) + "-" + lineNumber);
            }
            else if(referenceIndex == last.referenceIndex) { // sam reference
                if(startPosition < last.startPosition) {
                    throw new Exception("Ranges must be in sorted order in RANGES (line numbers " + (lineNumber-1) + "-" + lineNumber);
                }
                else if(last.endPosition + 1 < startPosition) {
                    // just add
                    this.ranges.add(new Range(referenceIndex, startPosition, endPosition));
                }
                else if(last.endPosition < endPosition) {
                    // merge over-lapping
                    last.endPosition = endPosition;
                }
            }
            else {
                // just add
                this.ranges.add(new Range(referenceIndex, startPosition, endPosition));
            }
        }
        else {
            // just add
            this.ranges.add(new Range(referenceIndex, startPosition, endPosition));
        }
    }

    /**
     * @return an iterator over the ranges.
     */
    public Iterator<Range> iterator() 
    {
        return this.ranges.iterator();
    }

    /**
     * @return the ith (zero-based) range.
     */
    public Range get(int i)
    {
        return this.ranges.get(i);
    }

    /**
     * @return the number of ranges.
     */
    public int size()
    {
        return this.ranges.size();
    }
}
