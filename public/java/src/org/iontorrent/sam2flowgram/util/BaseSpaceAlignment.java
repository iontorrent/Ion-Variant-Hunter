/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2flowgram.util;

import java.util.*;
import java.io.PrintStream;
import net.sf.samtools.*;

// TODO: banded version for speed
/**
 * A simple global Smith Waterman alignment class.
 *
 * @author nils.homer@lifetech.com
 */
public class BaseSpaceAlignment 
{
    /**
     * The score for a base match.
     */
    private static int SCORE_MATCH = 1;

    /**
     * The penalty for a base mismatch.
     */
    private static int PEN_MISMATCH = 3;

    /**
     * The penalty for a gap open.
     */
    private static int PEN_GAP_OPEN = 5;

    /**
     * The penalty for a gap extension.
     */
    private static int PEN_GAP_EXT = 2;

    /**
     * The minimum score (negative infinity).
     */
    public static int MINOR_INF = -1000000; // how to set this?

    /**
     * Indicates that the previous location in the path is from a match.
     */
    public final static int FROM_M = 0;

    /**
     * Indicates that the previous location in the path is from a insertion.
     */
    public final static int FROM_I = 1;

    /**
     * Indicates that the previous location in the path is from a deletion.
     */
    public final static int FROM_D = 2;

    /**
     * Indicates that the previous location in the path is from the start.
     */
    public final static int FROM_S = 3;


    /**
     * The zero-based index in the input tseq where the alignment starts.
     */
    public int tseqStart;

    /**
     * The zero-based index in the input tseq where the alignment ends.
     */
    public int tseqEnd;

    /**
     * A class to store the dp in the dynamic programming matrix.
     */
    private static class BaseSpaceAlignmentCell {
        /**
         * Stores the score for extending with a match.
         */
        public int matchScore;

        /**
         * Stores the score for extending with a insertion.
         */
        public int insScore;

        /**
         * Stores the score for extending with a deletion.
         */
        public int delScore;

        /**
         * Stores the previous cell in the path to a match.
         */
        public int matchFrom;

        /**
         * Stores the previous cell in the path to a insertion.
         */
        public int insFrom;

        /**
         * Stores the previous cell in the path to a deletion.
         */
        public int delFrom;

        /**
         * Creates a new cell.
         */
        public BaseSpaceAlignmentCell() 
        {
            this.matchScore = MINOR_INF; 
            this.matchFrom = FROM_S; 
            this.insScore = MINOR_INF; 
            this.insFrom = FROM_S; 
            this.delScore = MINOR_INF; 
            this.delFrom = FROM_S; 
        }

        /**
         * Debugging print function.
         * @param out the output stream.
         */
        public void print(PrintStream out)
        {
            out.println("[" + this.matchScore + "," + this.matchFrom
                    + ":" + this.insScore + "," + this.insFrom
                    + ":" + this.delScore + "," + this.delFrom
                    + "]");
        }

    }

    /**
     * Default constructor.
     */
    public BaseSpaceAlignment()
    {
    }
    
    /**
     * Performs a local Smith Waterman on the two DNA sequences.
     * @param qseq the query DNA sequence.
     * @param tseq the target DNA sequence.
     * @return a list of cigar elements representing the alignment.
     */
    public List<CigarElement> align(char qseq[], char tseq[])
        throws Exception

    {
        return align(qseq, tseq, true, true);
    }

    /**
     * Performs a global Smith Waterman on the two DNA sequences.
     * @param qseq the query DNA sequence.
     * @param tseq the target DNA sequence.
     * @return a list of cigar elements representing the alignment.
     */
    public List<CigarElement> align(char qseq[], 
            char tseq[],
            boolean startLocal,
            boolean endLocal)
        throws Exception
    {
        int i, j;
        BaseSpaceAlignmentCell dp[][] = null;
        int bestScore = MINOR_INF, bestFrom = FROM_S;
        int bestI = -1, bestJ = -1;

        // DEBUG
        /*
        System.err.println("startLocal=" + startLocal + " endLocal=" + endLocal);
        System.err.println("Aligning\n[" + new String(qseq) + "]\n[" + new String(tseq) + "]");
        */

        dp = new BaseSpaceAlignmentCell[qseq.length+1][tseq.length+1];
        for(i=0;i<qseq.length+1;i++) {
            for(j=0;j<tseq.length+1;j++) {
                dp[i][j] = new BaseSpaceAlignmentCell();
            }
        }

        // start cell
        dp[0][0].matchScore = 0; 
        // rows
        for(i=1;i<qseq.length+1;i++) {
            // allow for beginning insertion
            if(1 == i) {
                dp[i][0].insScore = dp[i-1][0].matchScore - PEN_GAP_OPEN - PEN_GAP_EXT;
                dp[i][0].insFrom = FROM_M; 
                //dp[i][0].insLength = dp[i-1][0].matchLength + 1;
            }
            else {
                dp[i][0].insScore = dp[i-1][0].insScore - PEN_GAP_EXT;
                dp[i][0].insFrom = FROM_I; 
                //dp[i][0].insLength = dp[i-1][0].matchLength + 1;
            }
        }
        // cols
        if(startLocal) {
            for(j=1;j<tseq.length+1;j++) {
                dp[0][j].matchScore = 0; // allow to begin anywhere in the reference 
                dp[0][j].matchFrom = FROM_S;
            }
        }

        // fill in the rest of the matrix
        for(i=1;i<qseq.length+1;i++) {
            for(j=1;j<tseq.length+1;j++) {
                // ins
                if(dp[i-1][j].matchScore - PEN_GAP_OPEN < dp[i-1][j].insScore) {
                    dp[i][j].insScore = dp[i-1][j].insScore - PEN_GAP_EXT;
                    dp[i][j].insFrom = FROM_I;
                }
                else {
                    dp[i][j].insScore = dp[i-1][j].matchScore - PEN_GAP_OPEN - PEN_GAP_EXT;
                    dp[i][j].insFrom = FROM_M;
                }

                // del
                if(dp[i][j-1].matchScore - PEN_GAP_OPEN < dp[i][j-1].delScore) {
                    dp[i][j].delScore = dp[i][j-1].delScore - PEN_GAP_EXT;
                    dp[i][j].delFrom = FROM_D;
                }
                else {
                    dp[i][j].delScore = dp[i][j-1].matchScore - PEN_GAP_OPEN - PEN_GAP_EXT;
                    dp[i][j].delFrom = FROM_M;
                }

                // match/mismatch
                int score = (qseq[i-1] == tseq[j-1]) ? SCORE_MATCH : (-1 * PEN_MISMATCH);
                if(dp[i-1][j-1].matchScore < dp[i-1][j-1].insScore) {
                    if(dp[i-1][j-1].delScore <= dp[i-1][j-1].insScore) {
                        dp[i][j].matchScore = dp[i-1][j-1].insScore + score;
                        dp[i][j].matchFrom = FROM_I;
                    }
                    else {
                        dp[i][j].matchScore = dp[i-1][j-1].delScore + score;
                        dp[i][j].matchFrom = FROM_D;
                    }
                }
                else {
                    if(dp[i-1][j-1].delScore <= dp[i-1][j-1].matchScore) {
                        dp[i][j].matchScore = dp[i-1][j-1].matchScore + score;
                        dp[i][j].matchFrom = FROM_M;
                    }
                    else {
                        dp[i][j].matchScore = dp[i-1][j-1].delScore + score;
                        dp[i][j].matchFrom = FROM_D;
                    }
                }

                /*
                System.err.print("i=" + i + " j=" + j);
                dp[i][j].print(System.err);
                */
            }
        }

        // Get the best path
        // Note: we want to align the full read, but not necessarily the full reference
        if(endLocal) {
            for(j=1;j<tseq.length+1;j++) {
                if(bestScore < dp[qseq.length][j].matchScore) {
                    bestScore = dp[qseq.length][j].matchScore;
                    bestI = qseq.length;
                    bestJ = j;
                    bestFrom = FROM_M;
                }
                if(bestScore < dp[qseq.length][j].insScore) {
                    bestScore = dp[qseq.length][j].insScore;
                    bestI = qseq.length;
                    bestJ = j;
                    bestFrom = FROM_I;
                }
                if(bestScore < dp[qseq.length][j].delScore) {
                    bestScore = dp[qseq.length][j].delScore;
                    bestI = qseq.length;
                    bestJ = j;
                    bestFrom = FROM_D;
                }
            }
        }
        else {
            // Note: we want to align the full read against the full reference
            if(bestScore < dp[qseq.length][tseq.length].matchScore) {
                bestScore = dp[qseq.length][tseq.length].matchScore;
                bestI = qseq.length;
                bestJ = tseq.length;
                bestFrom = FROM_M;
            }
            if(bestScore < dp[qseq.length][tseq.length].insScore) {
                bestScore = dp[qseq.length][tseq.length].insScore;
                bestI = qseq.length;
                bestJ = tseq.length;
                bestFrom = FROM_I;
            }
            if(bestScore < dp[qseq.length][tseq.length].delScore) {
                bestScore = dp[qseq.length][tseq.length].delScore;
                bestI = qseq.length;
                bestJ = tseq.length;
                bestFrom = FROM_D;
            }
        }

        // DEBUG
        /*
        System.err.println("bestScore=" + bestScore 
                + " bestI=" + bestI
                + " bestJ=" + bestJ
                + " bestFrom=" + bestFrom);
                */

        // Recover the path
        List<CigarElement> elements;
        CigarOperator prevOperator = null, curOperator = null;
        int prevOperatorLength = 0;
        int curI, curJ, curFrom, nextFrom;

        // initialize start and end offsets
        this.tseqStart = 0;
        this.tseqEnd = bestJ-1;

        elements = new LinkedList<CigarElement>();
        curI = bestI;
        curJ = bestJ;
        curFrom = bestFrom;
        while(FROM_S != curFrom && 0 < curI) { // read bases left

            // get operator
            switch(curFrom) {
                case FROM_M:
                    curOperator = CigarOperator.M;
                    break;
                case FROM_I:
                    curOperator = CigarOperator.I;
                    break;
                case FROM_D:
                    curOperator = CigarOperator.D;
                    break;
                default:
                    throw new Exception("Unknown from");
            }

            // add
            if(null != prevOperator && curOperator == prevOperator) {
                prevOperatorLength += 1;
                elements.set(elements.size()-1, new CigarElement(prevOperatorLength, prevOperator));
            }
            else {
                elements.add(new CigarElement(1, curOperator));
                prevOperator = curOperator;
                prevOperatorLength = 1;
            }

            this.tseqStart = curJ-1; // since curJ is one-based

            // update from
            // update row/col
            switch(curFrom) {
                case FROM_M:
                    curFrom = dp[curI][curJ].matchFrom;
                    curI--;
                    curJ--;
                    break;
                case FROM_I:
                    curFrom = dp[curI][curJ].insFrom;
                    curI--;
                    break;
                case FROM_D:
                    curFrom = dp[curI][curJ].delFrom;
                    curJ--;
                    break;
                case FROM_S:
                    curI--;
                    curJ--;
                    break;
                default:
                    throw new Exception("Unknown from");
            }

            // DEBUG
            //System.err.println("Next curFrom=" + curFrom + "\tcurI=" + curI + "\tcurJ=" + curJ);
        }

        if(0 == elements.size()) {
            throw new Exception("Zero elements added");
        }

        // reverse, since we added in post order
        Collections.reverse(elements);

        return elements;
    }
}
