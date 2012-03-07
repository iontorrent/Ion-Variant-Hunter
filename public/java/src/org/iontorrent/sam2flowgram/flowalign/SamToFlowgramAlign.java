/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/*
 * LICENSE to be determined
 */
package org.iontorrent.sam2flowgram.flowalign;

import org.iontorrent.sam2flowgram.util.*;
import org.iontorrent.sam2flowgram.io.*;

import net.sf.samtools.*;
import net.sf.picard.cmdline.*;
import net.sf.picard.io.IoUtil;

import java.io.*;
import java.util.*;
import java.lang.Math;

/**
 * TODO
 *
 * @author nils.homer@lifetech.com
 */
public class SamToFlowgramAlign extends CommandLineProgram {
    
    public static final String program_version = SamToFlowgramAlign.class.getPackage().getImplementationVersion();
    
    public static final String PROGRAM_VERSION = "0.0.4";
    @Usage (programVersion=PROGRAM_VERSION)
        public final String USAGE = getStandardUsagePreamble() + "Sam to Flowgram Alignment Utility.";
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM or BAM file.")
        public List<File> INPUT = new ArrayList<File>();
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference FASTA file.")
        public File REFERENCE=null;
    @Option(doc="The phase penalty for the flowgram alignment.", optional=true)
        public int PHASE_PENALTY=1;
    @Option(doc="Whether to suppress job-progress info on System.err.", optional=true)
        public boolean QUIET_STDERR=false;
    @Option(doc="The maximum size of the SAM records queue.", optional=true)
        public int MAX_QUEUE_SIZE = 65536;
    @Option(doc="The number of threads for parallel processing.", optional=true)
        public int NUM_THREADS = 1;
    @Option(doc="A range to examine.", optional=true)
        protected String RANGE=null;
    @Option(doc="Extra reference bases to consider during alignment.", optional=true)
        public int OFFSET = 10;
    @Option(doc="Validate alignments.", optional=true)
        public boolean VALIDATE_ALIGNMENTS = false;
    @Option(doc="prints program version to stdout", optional=true)
        private boolean printVersion = false;
    
    // Variables to track the run time of the program.
    private long startTime;
    private long endTime;

    // Variables to track the progress of the program.
    private int SamToFlowSpace_OUTPUT_CTR = 0; // the number of records processed.
    private int maxOutputStringLength = 0;
    
    // An added offset for flushing the graph/re-alignment.
    private final static int SamToFlowSpace_OFFSET_PROCESS = 100; 

    // Variables to handle the input reference sequence(s).
    private ReferenceSequence referenceSequence = null;

    // records ready to be added to the phred table
    private ThreadPoolLinkedList recordQueue = null; 

    // handles the input/output of SAM/BAM records.
    private SAMRecordIO io = null;

    // for tracking with the RANGES option.
    // for inputting within RANGES.
    private Ranges inputRanges = null; // the input ranges to consider.
    private Iterator<Range> inputRangesIterator = null; // an iterator over the input ranges to consider.
    // for outputting within RANGES.
    private Ranges outputRanges = null; // the output ranges to consider.
    private Iterator<Range> outputRangesIterator = null; // an iterator over the output ranges to consider.
    private Range outputRange = null; // the current output range.
    /**
     * @param args the command line arguments
     */
    public static void main(final String[] args) {
        
        
        SamToFlowgramAlign s = new SamToFlowgramAlign();
        
        new SamToFlowgramAlign().instanceMain(args);
    }

    /**
     * The default command line validation.
     */
    protected String[] customCommandLineValidation()
    {
     if (printVersion)
            System.out.println("program version: " + program_version);
        return super.customCommandLineValidation();
    }

    /**
     * The main work-horse function.
     */
    protected int doWork() 
    {
        FlowAlignRecord rec = null;

        try { 
            this.startTime = System.nanoTime();

            // Check some parameters
            IoUtil.assertFileIsReadable(REFERENCE);

            if(!QUIET) {
                QUIET=true; // Over-ride
                //throw new Exception("Please use option 'QUIET' when outputting to stdout.");
            }

            // Initialize the reference
            this.referenceSequence = new ReferenceSequence(REFERENCE);

            // Get ranges
            if (RANGE == null) {
                this.inputRanges = new Ranges(this.referenceSequence.getDictionary());
                this.outputRanges = new Ranges(this.referenceSequence.getDictionary());
                this.io = new SAMRecordIO(INPUT, null, null, PROGRAM_VERSION, false, this.referenceSequence.getDictionary());
            }
            else {
                 this.inputRanges = new Ranges(RANGE, this.referenceSequence.getDictionary(), 0);
                 this.outputRanges = new Ranges(RANGE, this.referenceSequence.getDictionary(), 0);
                 this.io = new SAMRecordIO(INPUT, null, null, PROGRAM_VERSION, true, this.referenceSequence.getDictionary());
            }
            
            this.inputRangesIterator = this.inputRanges.iterator();
            this.outputRangesIterator = this.outputRanges.iterator();
            if(!this.inputRangesIterator.hasNext()) {
                return 0;
            }
            this.outputRange = this.outputRangesIterator.next();
    
            // Go through each input range
            while(this.inputRangesIterator.hasNext()) {
                Range inputRange = this.inputRangesIterator.next();

                int prevReferenceIndex=-1;
                int prevAlignmentStart=-1;

                this.recordQueue = new ThreadPoolLinkedList();
                this.io.query(this.referenceSequence.getDictionary().getSequence(inputRange.referenceIndex).getSequenceName(),
                        inputRange.startPosition, 
                        inputRange.endPosition);

                // Get first record
                if(null == rec) {
                    rec = this.getNextAlignRecord();
                }

                // skip this range
                if(null != rec &&
                        rec.referenceIndex != inputRange.referenceIndex) 
                {
                    continue;
                }

                // Get the reference bases
                // TODO: could retrieve just the specified bases
                this.referenceSequence.moveTo(inputRange.referenceIndex);

                while(null != rec) {
                    if(rec.record.getReadUnmappedFlag()) { 
                        // ignore
                    }
                    else {
                        int curReferenceIndex = rec.referenceIndex;
                        int curAlignmentStart = rec.positionStart;

                        // Make sure that it is sorted
                        if(rec.referenceIndex < prevReferenceIndex 
                                || (rec.referenceIndex == prevReferenceIndex && rec.positionStart < prevAlignmentStart)) 
                        {
                            throw new Exception("SAM/BAM file is not co-ordinate sorted.");
                        }

                        if(curReferenceIndex != inputRange.referenceIndex) {
                            break;
                        }

                        prevReferenceIndex = rec.referenceIndex;
                        prevAlignmentStart = rec.positionStart;

                        // add to the table
                        if(this.MAX_QUEUE_SIZE <= this.recordQueue.size()) {
                            this.realign();
                        }

                        // Add the current record to the graph list
                        this.recordQueue.add(rec);
                    }

                    // get new record
                    rec = this.getNextAlignRecord();
                }
                // add to the table
                this.realign();
                // nullify
                this.recordQueue = null;
            }

            // Close input/output files
            this.io.closeAll();

            this.endTime = System.nanoTime();

            // to end it all
            if(!QUIET_STDERR) {
                System.err.println("");
                System.err.println("Processing complete");
                // Memory
                double totalMemory = (double)Runtime.getRuntime().totalMemory();
                double totalMemoryLog2 = Math.log(totalMemory) / Math.log(2.0);
                if(totalMemoryLog2 < 10) {
                    System.err.println("Total memory usage: " + (int)totalMemory + "B");
                } 
                else if(totalMemoryLog2 < 20) {
                    System.err.println("Total memory usage: " + (Math.round(100 * totalMemory / Math.pow(2, 10)) / 100) + "KB");
                }
                else {
                    System.err.println("Total memory usage: " + (Math.round(100 * totalMemory / Math.pow(2, 20)) / 100) + "MB");
                }
                // Run time
                long seconds = (this.endTime - this.startTime) / 1000000000;
                long hours = seconds / 3600; seconds -= hours * 3600; 
                long minutes = seconds / 60; seconds -= minutes* 60; 
                System.err.println("Total execution time: " + hours + "h : " + minutes + "m : " + seconds + "s");
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("Please report bugs to nils.homer@lifetech.com");
            System.exit(1);
        }

        // this is annoying
        QUIET = true;

        return 0;
    }

    /**
     * @return the next record to process.
     */
    private FlowAlignRecord getNextAlignRecord()
        throws Exception
    {
        if(this.io.hasNextAlignRecord()) {
            return this.io.getNextAlignRecord();
        }
        else {
            return null;
        }
    }

    /*
     * Outputs a progress message to stderr.
     * @param rec the last SAM record processed.
     */
    private void outputProgress(SAMRecord rec)
    {
        if(QUIET_STDERR) {
            return;
        }
        else {
            // TODO: enforce column width ?
            int i;
            String outputString = new String("Records processsed: " + this.SamToFlowSpace_OUTPUT_CTR + " (last " + rec.getReferenceName() + ":" + rec.getAlignmentStart() + "-" + rec.getAlignmentEnd() + ")");
            int outputStringLength = outputString.length();
            if(this.maxOutputStringLength < outputStringLength) {
                this.maxOutputStringLength = outputStringLength;
            }
            System.err.print("\r" + outputString);
            for(i=outputStringLength;i<this.maxOutputStringLength;i++) {
                System.err.print(" ");
            }
        }
    }

    /**
     * Processes the records in the scaffolding queue via multiple threads.  The 
     * processed records are added to a queue for signal aggregation.
     */
    private void realign()
        throws Exception
    {
        // Process alignments
        if(0 < this.recordQueue.size()) { 

            int i, size;
            LinkedList<Thread> threads = null;
            LinkedList<LinkedList<FlowAlignRecord>> inputThreadLists = null;
            LinkedList<LinkedList<FlowAlignRecord>> outputThreadLists = null;
            SAMRecord lastRecord = null;

            // Get the records for the threads 
            lastRecord = recordQueue.getLast().record;
            size = recordQueue.size();
            inputThreadLists = recordQueue.getThreadLists(this.NUM_THREADS, -1);
            outputThreadLists = recordQueue.getThreadLists(this.NUM_THREADS, -1);

            // Create threads
            threads = new LinkedList<Thread>();
            for(i=0;i<this.NUM_THREADS;i++) {
                threads.add(new RealignThread(i, 
                            inputThreadLists.get(i).listIterator(),
                            outputThreadLists.get(i).listIterator()));
            }

            // Start
            for(i=0;i<this.NUM_THREADS;i++) {
                threads.get(i).start();
            }

            // Join
            for(i=0;i<this.NUM_THREADS;i++) {
                threads.get(i).join();
            }

            // Output
            List<ListIterator<FlowAlignRecord>> iters = new LinkedList<ListIterator<FlowAlignRecord>>();
            for(i=size=0;i<this.NUM_THREADS;i++) {
                size +=  outputThreadLists.get(i).size();
                iters.add(outputThreadLists.get(i).listIterator());
            }
            for(i=0;0<size;i++) {
                if(this.NUM_THREADS <= i) {
                    i=0;
                }
                if(iters.get(i).hasNext()) {
                    FlowAlignRecord rec = iters.get(i).next();
                    lastRecord = rec.record;

                    // output
                    //System.out.println("" + rec.record.getReadName().toString() + "\t" + rec.alignment.getScore());
                      System.out.print( rec.getSAMString() );
                      System.out.println(rec.alignment.getAlignmentString());
                      System.out.println(rec.positionStart+","+rec.positionEnd+","+rec.alignment.getScore());
                    //rec.alignment.print(System.out, Integer.MAX_VALUE);

                    this.SamToFlowSpace_OUTPUT_CTR++;
                    size--;
                }
            }

            // Output progress
            this.outputProgress(lastRecord);
        }
    }

    /**
     * Handles multi-threading the graph scaffolding and signal aggregation.
     */
    private class RealignThread extends Thread {
        /**
         * The thread identifier.
         */
        private int threadID;

        /**
         * The records to process.
         */
        private ListIterator<FlowAlignRecord> arIn;
          
        /**
         * The records to output.
         */
        private ListIterator<FlowAlignRecord> arOut;
        
        /**
         * @param threadID the thread identifier.
         * @param arIn the records to process.
         */
        public RealignThread(int threadID,
                ListIterator<FlowAlignRecord> arIn,
                ListIterator<FlowAlignRecord> arOut)
        {
            this.threadID = threadID;
            this.arIn = arIn;
            this.arOut = arOut;
        }

        /**
         * Starts the thread and process the records.
         */
        public void run()
        {
            FlowAlignRecord rec = null;
            int i;

            try {
                while(this.arIn.hasNext()) {
                    // Get record
                    rec = this.arIn.next();

                    // set the alignment
                    rec.setAlignment(referenceSequence, PHASE_PENALTY, OFFSET, VALIDATE_ALIGNMENTS); 
                    if(rec.strand) { // reverse compliment
                        rec.alignment.reverseCompliment();
                    }

                    // add to the output
                    this.arOut.add(rec);
                }
            } catch (Exception e) {
                if(null != rec) {
                    System.err.println("On read: " + rec.record.getReadName().toString());
                }
                e.printStackTrace();
                System.err.println("Please report bugs to nils.homer@lifetech.com");
                System.exit(1);
            }
        }
    }
}
