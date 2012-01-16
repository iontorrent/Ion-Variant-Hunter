/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2fs.flowspace; 

import java.io.PrintStream;
import java.lang.Math;
import java.lang.Integer;
import org.iontorrent.sam2fs.util.SamToFlowSpaceUtil;

/**
 * Represents an alignment in flow space.
 *
 * @author nils.homer@lifetech.com
 */
public class FlowSpaceAlignment 
{
    /**
     * The alignment was extended from a match. 
     */
    public static final int FROM_M = 0;

    /**
     * The alignment was extended from an insertion.
     */
    public static final int FROM_I = 1;

    /**
     * The alignment was extended from an deletion.
     */
    public static final int FROM_D = 2;

    /**
     * The alignment was extended from a match with an empty flow.
     */
    public static final int FROM_ME = 3;

    /**
     * The alignment was extended from an insertion with an empty flow.
     */
    public static final int FROM_IE = 4;

    /**
     * The alignment was extended from a phased match (skipping a phase).
     */
    public static final int FROM_MP = 5;

    /**
     * The alignment was extended from a phased insertion (skipping a phase).
     */
    public static final int FROM_IP = 6;

    /**
     * The alignment was extended from an insertion.
     */
    public static final int FROM_S = 7;

    /**
     * The lower bound on the alignment score, or negative infinity.
     */
    public static final int MINOR_INF = -1000000; // how to set this?

    /**
     * A flow deletion in the alignment string.
     */
    public static final char ALN_DEL = '-';

    /**
     * A flow insertion in the alignment string.
     */
    public static final char ALN_INS = '+';

    /**
     * A flow match in the alignment string.
     */
    public static final char ALN_MATCH = '|';

    /**
     * A flow mismatch in the alignment string.
     */
    public static final char ALN_MISMATCH = ' ';

    /**
     * The alignment length.
     */
    public int length = 0; // the alignment length

    /**
     * The memory allocated for the alignment.
     */
    private int mem = 0; // the mememory allocated for the alignment

    /**
     * The best alignment score found so far.
     */
    private int score;
    
    // Depend on the alignment order 
    /**
     * The flow order for the alignment, including deleted reference bases.
     */
    public byte flowOrder[] = null;

    /**
     * The query or read sequence in the alignment, including gaps.
     */
    public int qseq[] = null; // read

    /**
     * The target or reference sequence in the alignment, including gaps.
     */
    public int tseq[] = null; // reference

    /**
     * The alignment string.
     */
    public char aln[] = null; // alignment string

    /**
     * The index of the first non-empty query flow.
     */
    public int nonEmptyFlowFirst; // the index of the first non-empty read flow

    /**
     * The index of the last non-empty query flow.
     */
    public int nonEmptyFlowLast; // the index of the last non-empty read flow

    /**
     * The zero-based index in the input tseq where the alignment starts.
     */
    public int tseqStart;
    
    /**
     * The zero-based index in the input tseq where the alignment ends.
     */
    public int tseqEnd;

    /**
     * The tseq length, for this.reverseCompliment().
     */
    public int tseqLength;

    /**
     * The penalty for phasing a flow.
     */
    public static final int FLOW_SPACE_PHASE_PENALTY = 1;

    // qseq - query - read
    // tseq - target - reference
    /**
     * Represents an alignment in flow space.
     *
     * Notes: we want to align the flow flowQseq to a subsequence of tseq
     *
     * @param flowQseq the query's flow sequence.
     * @param tseq the target base in integer format.
     * @param qseqFlowOrder the flow order of the query flow sequence.
     */
    public FlowSpaceAlignment(FlowSeq flowQseq, byte tseq[], 
            FlowOrder qseqFlowOrder)
        throws Exception
    {
        this(flowQseq, tseq, qseqFlowOrder, false, false, 
                FLOW_SPACE_PHASE_PENALTY);
    }
    
    /**
     * Represents an alignment in flow space.
     *
     * Notes: we want to align the flow flowQseq to a subsequence of tseq
     *
     * @param flowQseq the query's flow sequence.
     * @param tseq the target base in integer format.
     * @param qseqFlowOrder the flow order of the query flow sequence.
     * @param startLocal false if the we must begin the alignment at the start of the target, true otherwise 
     * @param endLocal false if the we must end the alignment at the end of the target, true otherwise 
     * @param phasePenalty the penalty for phasing in the alignment.
     */
    public FlowSpaceAlignment(FlowSeq flowQseq, byte tseq[], FlowOrder qseqFlowOrder,
            boolean startLocal, boolean endLocal, int phasePenalty)
        throws Exception
    {
        this.init(flowQseq, tseq, qseqFlowOrder, startLocal, endLocal, phasePenalty);
    }
    
    // qseq - query - read
    // tseq - target - reference
    /**
     * Represents an alignment in flow space.
     *
     * Notes: we want to align the flow flowQseq to a subsequence of tseq
     *
     * @param flowQseq the query's flow sequence.
     * @param tseq the target base in integer format.
     * @param qseqFlowOrder the flow order of the query flow sequence.
     * @param startLocal false if the we must begin the alignment at the start of the target, true otherwise 
     * @param endLocal false if the we must end the alignment at the end of the target, true otherwise 
     * @param phasePenalty the penalty for phasing in the alignment.
     */
    public void init(FlowSeq flowQseq, byte tseq[], FlowOrder qseqFlowOrder,
            boolean startLocal, boolean endLocal, int phasePenalty)
        throws Exception
    {
        int i, j, k, l;
        int vScoreP, vScoreE, vFromP, vFromE;
        int cType, iFrom;
        int bestI, bestJ, bestCType;
        int gapSumsI[] = null;
        FlowSeq flowTseq = null;
        FlowSpaceAlignmentCell dp[][] = null;
        FlowOrder tseqFlowOrder = null;

        // create tseq flow order
        tseqFlowOrder = new FlowOrder(tseq, true);

        // for this.reverseCompliment().
        tseqLength = tseq.length;

        // HERE
        /*
        for(i=0;i<tseqFlowOrder.length;i++) {
            System.err.println("i=" + i + " tseqFlowOrder.flowOrder[i]=" + tseqFlowOrder.flowOrder[i]);
        }
        */

        // convert bases to flow space
        flowTseq = new FlowSeq(tseq, tseqFlowOrder.flowOrder);

        // HERE
        /*
        System.err.println("startLocal=" + startLocal + " endLocal=" + endLocal);
        flowQseq.print(System.err);
        qseqFlowOrder.print(System.err);
        flowTseq.print(System.err);
        */
        
        // init
        this.mem = flowQseq.length + flowTseq.length;
        this.flowOrder = new byte[mem];
        this.qseq = new int[mem];
        this.tseq = new int[mem];
        this.aln= new char[mem];
        this.length = 0;

        //init gap sums & dp matrix 
        gapSumsI = new int[flowQseq.length];
        dp = new FlowSpaceAlignmentCell[1+flowQseq.length][1+flowTseq.length];

        for(i=0;i<=flowQseq.length;i++) {
            if (i < flowQseq.length) {
                k = i % qseqFlowOrder.length;
                //j = (i < qseqFlowOrder.jumpRev[k]) ? 0 : (i - qseqFlowOrder.jumpRev[k]);
                //j = (i < qseqFlowOrder.jumpRev[k]) ? 0 : (i - qseqFlowOrder.jumpRev[k] + 1);
                j = (i < qseqFlowOrder.jumpRev[k]) ? 0 : (i - qseqFlowOrder.jumpRev[k] + 1);
                gapSumsI[i] = phasePenalty;
                while(j <= i) {
                    gapSumsI[i] += flowQseq.flow[j];
                    j++;
                }
                /*
                System.err.println("i=" + i
                        + " flowQseq.flow[i]=" + flowQseq.flow[i]
                        + " qseqFlowOrder.flowOrder[i]=" + qseqFlowOrder.flowOrder[i]
                        + " qseqFlowOrder.jumpRev[i]=" + qseqFlowOrder.jumpRev[i]
                        + " gapSumsI[i]=" + gapSumsI[i]); 
                        */
            }
            // dp matrix init
            for(j=0;j<=flowTseq.length;j++) {
                dp[i][j] = new FlowSpaceAlignmentCell();
                dp[i][j].matchScore = dp[i][j].insScore = dp[i][j].delScore = MINOR_INF;
                dp[i][j].matchFrom = dp[i][j].insFrom = dp[i][j].delFrom = FROM_S;
            }
            
            if (i > 0) {
                k = (i-1) % qseqFlowOrder.length;
                iFrom = ((i < qseqFlowOrder.jumpRev[k]) ? 0 : (i - qseqFlowOrder.jumpRev[k])); 

                // vertical
                // only allow phasing from an insertion
                if(0 == iFrom) {
                    dp[i][0].insScore = 0 - gapSumsI[i-1];
                    dp[i][0].insFrom = FROM_IP;
                }
                else {
                    dp[i][0].insScore = dp[iFrom][0].insScore - gapSumsI[i-1];
                    dp[i][0].insFrom = FROM_IP;
                }
            }
        }

        // init start cells
        dp[0][0].matchScore = 0; 
        // align
        for(i=1;i<=flowQseq.length;i++) { // query
            k = (i-1) % qseqFlowOrder.length;
            iFrom = ((i < qseqFlowOrder.jumpRev[k]) ? 0 : (i - qseqFlowOrder.jumpRev[k])); 
            for(j=1;j<=flowTseq.length;j++) { // target
                // horizontal
                if(dp[i][j-1].delScore < dp[i][j-1].matchScore) {
                    if(dp[i][j-1].insScore <= dp[i][j-1].matchScore) {
                        dp[i][j].delScore = dp[i][j-1].matchScore - flowTseq.flow[j-1];
                        dp[i][j].delFrom = FROM_M;
                    }
                    else {
                        dp[i][j].delScore = dp[i][j-1].insScore - flowTseq.flow[j-1];
                        dp[i][j].delFrom = FROM_I;
                    }
                }
                else {
                    if(dp[i][j-1].insScore <= dp[i][j-1].delScore) {
                        dp[i][j].delScore = dp[i][j-1].delScore - flowTseq.flow[j-1];
                        dp[i][j].delFrom = FROM_D;
                    }
                    else {
                        dp[i][j].delScore = dp[i][j-1].insScore - flowTseq.flow[j-1];
                        dp[i][j].delFrom = FROM_I;
                    }
                }

                // vertical
                // four moves:
                // 1. phased from match
                // 2. phased from ins
                // 3. empty from match
                // 4. empth from ins
                // Note: use the NEXT reference base for flow order matching
                if(j == flowTseq.length // no next reference base
                        || (1 == i) // always start with leading phasing
                        || (qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length] == tseqFlowOrder.flowOrder[j % tseqFlowOrder.length])) 
                {
                    vScoreE = MINOR_INF; 
                    vFromE = FROM_ME;
                }
                else {
                    if(dp[i-1][j].insScore <= dp[i-1][j].matchScore) {
                        vScoreE = dp[i-1][j].matchScore - flowQseq.flow[i-1];
                        vFromE = FROM_ME;
                    }
                    else {
                        vScoreE = dp[i-1][j].insScore - flowQseq.flow[i-1];
                        vFromE = FROM_IE;
                    }
                    // Start anywhere in tseq
                    if(i == 1 && vScoreE + flowQseq.flow[i-1] < 0) {
                        vScoreE = 0 - flowQseq.flow[i-1];
                        vFromE = FROM_S;
                    }
                }
                // phased from ...
                if(dp[iFrom][j].insScore <= dp[iFrom][j].matchScore) {
                    vScoreP = dp[iFrom][j].matchScore - gapSumsI[i-1];
                    vFromP = FROM_MP;
                }
                else {
                    vScoreP = dp[iFrom][j].insScore - gapSumsI[i-1];
                    vFromP = FROM_IP;
                }
                // compare empty vs. phased
                if(vScoreP <= vScoreE) { // Note: always choose empty over phased
                    dp[i][j].insScore = vScoreE;
                    dp[i][j].insFrom = vFromE;
                }
                else {
                    dp[i][j].insScore = vScoreP;
                    dp[i][j].insFrom = vFromP;
                }

                // diagonal
                if(qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length] != tseqFlowOrder.flowOrder[(j-1) % tseqFlowOrder.length]) {
                    // out of phase, do not want
                    dp[i][j].matchScore = MINOR_INF;
                    dp[i][j].matchFrom = FROM_S;
                }
                else {
                    int s = ((flowQseq.flow[i-1] < flowTseq.flow[j-1]) ? (flowTseq.flow[j-1]-flowQseq.flow[i-1]) : (flowQseq.flow[i-1]-flowTseq.flow[j-1]));
                    // NB: do not penalize it full on the first or last flow
                    if(i == 1 || i == flowQseq.length) {
                        s = SamToFlowSpaceUtil.getFlowSignalFromBaseCall(SamToFlowSpaceUtil.getBaseCallFromFlowSignal(flowQseq.flow[i-1]));
                        s = (flowQseq.flow[i-1] < s) ? (s - flowQseq.flow[i-1]) : (flowQseq.flow[i-1] - s);
                    }
                    if(dp[i-1][j-1].insScore <= dp[i-1][j-1].matchScore) {
                        if(dp[i-1][j-1].delScore <= dp[i-1][j-1].matchScore) {
                            dp[i][j].matchScore = dp[i-1][j-1].matchScore - s;
                            dp[i][j].matchFrom = FROM_M;
                        }
                        else {
                            dp[i][j].matchScore = dp[i-1][j-1].delScore - s;
                            dp[i][j].matchFrom = FROM_D;
                        }
                    }
                    else {
                        if(dp[i-1][j-1].delScore <= dp[i-1][j-1].insScore) {
                            dp[i][j].matchScore = dp[i-1][j-1].insScore - s;
                            dp[i][j].matchFrom = FROM_I;
                        }
                        else {
                            dp[i][j].matchScore = dp[i-1][j-1].delScore - s;
                            dp[i][j].matchFrom = FROM_D;
                        }
                    }

                    // Start anywhere in tseq
                    if(startLocal && 1 == i && dp[i][j].matchScore + s < 0) {
                        dp[i][j].matchScore = 0 - s;
                        dp[i][j].matchFrom = FROM_S;
                    }
                }
                
                // HERE
                /*
                System.err.print("i=" + i + " j=" + j 
                        + " qseq=" + qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length]
                        + " tseq=" + tseqFlowOrder.flowOrder[(j-1) % tseqFlowOrder.length]
                        + " ");
                dp[i][j].print(System.err);
                */
            }
        }

        // Get best scoring cell
        this.score = MINOR_INF-1; 
        bestCType = FROM_S;
        bestI = -1;
        bestJ = -1;

        // TODO: want to map the query into a sub-sequence of the target
        // We can end anywhere in the target, but we haven't done the beginning.
        // We also need to return where the start end in the target to update start/end position(s).
        if(endLocal) {
            for(j=1;j<=flowTseq.length;j++) { // target
                /*
                System.err.println("j=" + j
                        + " " + tseqFlowOrder.flowOrder[(j-1) % tseqFlowOrder.length]
                        + " " + flowTseq.flow[j-1]
                        + " " + dp[flowQseq.length][j].delScore
                        + " " + dp[flowQseq.length][j].insScore
                        + " " + dp[flowQseq.length][j].matchScore
                        + " " + this.score
                        );
                */
                //if(this.score <= dp[flowQseq.length][j].delScore) {
                if(0 < flowTseq.flow[j-1] && this.score <= dp[flowQseq.length][j].delScore) {
                    bestI = flowQseq.length;
                    bestJ = j;
                    this.score = dp[flowQseq.length][j].delScore;
                    bestCType = FROM_D;
                }
                if(this.score <= dp[flowQseq.length][j].insScore) {
                    bestI = flowQseq.length;
                    bestJ = j;
                    this.score = dp[flowQseq.length][j].insScore;
                    bestCType = FROM_I;
                }
                //if(0 < flowTseq.flow[j-1] && this.score <= dp[flowQseq.length][j].matchScore) {
                if(this.score <= dp[flowQseq.length][j].matchScore) {
                    bestI = flowQseq.length;
                    bestJ = j;
                    this.score = dp[flowQseq.length][j].matchScore;
                    bestCType = FROM_M;
                }
            }
        }
        else {
            if(this.score <= dp[flowQseq.length][flowTseq.length].delScore) {
                bestI = flowQseq.length;
                bestJ = flowTseq.length;
                this.score = dp[flowQseq.length][flowTseq.length].delScore;
                bestCType = FROM_D;
            }
            if(this.score <= dp[flowQseq.length][flowTseq.length].insScore) {
                bestI = flowQseq.length;
                bestJ = flowTseq.length;
                this.score = dp[flowQseq.length][flowTseq.length].insScore;
                bestCType = FROM_I;
            }
            if(this.score <= dp[flowQseq.length][flowTseq.length].matchScore) {
                bestI = flowQseq.length;
                bestJ = flowTseq.length;
                this.score = dp[flowQseq.length][flowTseq.length].matchScore;
                bestCType = FROM_M;
            }
        }

        /*
        System.err.println("flowQseq.length=" + flowQseq.length
                + " flowTseq.length=" + flowTseq.length);
        System.err.println("bestI=" + bestI
                + " bestJ=" + bestJ
                + " score=" + score
                + " bestCType=" + bestCType);
        */
        
        // Calculate tseqEnd
        this.tseqEnd = 0;
        for(j=0;j<bestJ;j++) {
            this.tseqEnd += SamToFlowSpaceUtil.getBaseCallFromFlowSignal(flowTseq.flow[j]);
            //System.err.println(SamToFlowSpaceUtil.DNA[tseqFlowOrder.flowOrder[j]] + " " + flowTseq.flow[j] + " " + this.tseqEnd);
        }
        this.tseqEnd--; 

        i = bestI;
        j = bestJ;
        cType = bestCType;


        // trace path back
        while(0 < i) { // qseq flows left
            int nextCType = -1;

            // HERE
            //System.err.println("i=" + i + " j=" + j + " cType=" + cType);
        
            if(FROM_M == cType) {
                nextCType = dp[i][j].matchFrom;
                this.add(flowQseq.flow[i-1], flowTseq.flow[j-1], qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length]);
                i--;
                j--;
            }
            else if(FROM_I == cType) {
                nextCType = dp[i][j].insFrom;
                if(dp[i][j].insFrom == FROM_ME || dp[i][j].insFrom == FROM_IE) {
                    this.add(flowQseq.flow[i-1], 0, qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length]);
                    i--;
                }
                else if(dp[i][j].insFrom == FROM_MP || dp[i][j].insFrom == FROM_IP) {
                    k = (i-1) % qseqFlowOrder.length;
                    iFrom = ((i < qseqFlowOrder.jumpRev[k]) ? 0 : (i - qseqFlowOrder.jumpRev[k])); 
                    while(iFrom < i) {
                        k = (i-1) % qseqFlowOrder.length;
                        this.add(flowQseq.flow[i-1], -1, qseqFlowOrder.flowOrder[k]);
                        i--;
                    }
                }
                else if(dp[i][j].insFrom == FROM_S) {
                    while(0 < i) {
                        // always a start insertion
                        this.add(flowQseq.flow[i-1], -1, qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length]);
                        //this.add(flowQseq.flow[i-1], 0, qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length]);
                        i--;
                    }
                }
                else {
                    System.err.println("dp[i][j].insFrom=" + dp[i][j].insFrom);
                    throw new Exception("bug encountered");
                }
            }
            else if(FROM_D == cType) {
                nextCType = dp[i][j].delFrom;
                //System.err.println("FOUND!");
                this.add(-1, flowTseq.flow[j-1], tseqFlowOrder.flowOrder[(j-1) % tseqFlowOrder.length]);
                j--;
            }
            else {
                System.err.println("cType=" + cType);
                System.err.println("i=" + i + " j=" + j);
                qseqFlowOrder.print(System.err);
                flowQseq.print(System.err);
                tseqFlowOrder.print(System.err);
                flowTseq.print(System.err);
                throw new Exception("bug encountered");
            }

            // HERE
            //System.err.println("nextCType=" + nextCType);
            switch(nextCType) {
                case FROM_M:
                case FROM_I:
                case FROM_D:
                case FROM_S:
                    cType = nextCType;
                    break;
                case FROM_ME:
                case FROM_MP:
                    cType = FROM_M;
                    break;
                case FROM_IE:
                case FROM_IP:
                    cType = FROM_I;
                    break;
                default:
                    throw new Exception("bug encountered");
            }
        }
        // Calculate tseqStart
        this.tseqStart = 0;
        for(i=0;i<j;i++) {
            this.tseqStart += SamToFlowSpaceUtil.getBaseCallFromFlowSignal(flowTseq.flow[i]);
        }
        
        // reverse the arrays tseq, qseq, aln, flowOrder
        this.reverse();
        
        // TODO: are these needed?
        this.nonEmptyFlowFirst = 0;
        for(i=0;i<this.length;i++) {
            if(0 < SamToFlowSpaceUtil.getBaseCallFromFlowSignal(this.qseq[i])) {
                this.nonEmptyFlowFirst = i;
                break;
            }
        }
        this.nonEmptyFlowLast = 0;
        for(i=this.length-1;0<=i;i--) {
            if(0 < SamToFlowSpaceUtil.getBaseCallFromFlowSignal(this.qseq[i])) {
                this.nonEmptyFlowLast = i;
                break;
            }
        }

        // HERE
        //this.print(System.err);
    }

    /** 
     * Reverse the alignment.
     *
     * Note: this does not reverse non-empty first/last index.
     */
    private void reverse()
    {
        int i;
        for(i=0;i<this.length/2;i++) {
            int b;
            char c;
            byte by;
            
            b = this.qseq[i];
            this.qseq[i] = this.qseq[this.length-i-1];
            this.qseq[this.length-i-1] = b;
            
            c = this.aln[i];
            this.aln[i] = this.aln[this.length-i-1];
            this.aln[this.length-i-1] = c;
           
            b = this.tseq[i];
            this.tseq[i] = this.tseq[this.length-i-1];
            this.tseq[this.length-i-1] = b;
            
            by = this.flowOrder[i];
            this.flowOrder[i] = this.flowOrder[this.length-i-1];
            this.flowOrder[this.length-i-1] = by;
        }
    }
    
    /** 
     * Reverse compliments this alignment, to faciliate its representation
     * on the forward genomic strand for mapped reverse strand queries.
     */
    public void reverseCompliment()
    {
        int i;
        for(i=0;i<this.length/2;i++) {
            int b;
            char c;
            byte by;
            
            b = this.qseq[i];
            this.qseq[i] = this.qseq[this.length-i-1];
            this.qseq[this.length-i-1] = b;
            
            c = this.aln[i];
            this.aln[i] = this.aln[this.length-i-1];
            this.aln[this.length-i-1] = c;
           
            b = this.tseq[i];
            this.tseq[i] = this.tseq[this.length-i-1];
            this.tseq[this.length-i-1] = b;
            
            by = SamToFlowSpaceUtil.NTINT2COMP[(int)this.flowOrder[i]];
            this.flowOrder[i] = SamToFlowSpaceUtil.NTINT2COMP[(int)this.flowOrder[this.length-i-1]];
            this.flowOrder[this.length-i-1] = by;
        }
        if(1 == (this.length % 2)) {
            this.flowOrder[i] = SamToFlowSpaceUtil.NTINT2COMP[(int)this.flowOrder[i]];
        }

        // reverse tseq bounds
        i = this.tseqStart;
        this.tseqStart = this.tseqLength - this.tseqEnd - 1;
        this.tseqEnd = this.tseqLength - i - 1;

        // convert leading empty flows to insertions
        for(i=0;i<this.length;i++) {
            if(0 == this.tseq[i] && (ALN_MATCH == this.aln[i] || ALN_MISMATCH == this.aln[i])) {
                this.tseq[i] = 0;
                this.aln[i] = ALN_INS;
            }
            else {
                break;
            }
        }

        // Reverse non-empty first/last
        //nonEmptyFlowFirst = this.length - nonEmptyFlowFirst - 1;
        //nonEmptyFlowLast = this.length - nonEmptyFlowLast - 1;
    }

    /**
     * Cells in the dynamic programming matrix.
     */
    private class FlowSpaceAlignmentCell 
    {
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
        public FlowSpaceAlignmentCell() 
        {
            // do nothing
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
     * Adds the given flow to the alignment.
     * @param qseqN the number of query bases.
     * @param tseqN the number of target bases.
     * @param base the flow base.
     */
    private void add(int qseqN, int tseqN, byte base) 
        throws Exception
    {
        // not enough memory
        if(this.mem <= this.length) {
            throw new Exception("this.mem <= this.length ["+this.mem+"<="+this.length+"]");
        }
        // add in the alignment
        this.flowOrder[this.length] = base;
        this.qseq[this.length] = qseqN;
        this.tseq[this.length] = tseqN;
        if(-1 == qseqN) {
            this.aln[this.length] = ALN_DEL; 
            this.qseq[this.length] = 0;
        }
        else if(-1 == tseqN) {
            this.aln[this.length] = ALN_INS;
            this.tseq[this.length] = 0;
        }
        else if(SamToFlowSpaceUtil.getBaseCallFromFlowSignal(qseqN) == SamToFlowSpaceUtil.getBaseCallFromFlowSignal(tseqN)) {
            this.aln[this.length] = ALN_MATCH;
        }
        else {
            this.aln[this.length] = ALN_MISMATCH;
        }
        this.length++;
    }

    /**
     * Moves flow insertions and deletions.
     */
    public boolean rightAdjustIndels()
    {
        int i, j, k;
        boolean adjusted = false;
        int lower;

        lower = 0;
        // skip over start indels
        while(lower < this.length) {
            if(ALN_MATCH == this.aln[lower] || ALN_MISMATCH == this.aln[lower]) {
                break;
            }
            lower++;
        }

        i = this.length - 2; // since if we end with a deletion, we have no bases to shift
        while(lower < i) {
            // find the end of the indel
            if(ALN_INS == this.aln[i] || ALN_DEL == this.aln[i]) {
                // get the left-most indel base
                j = i; // start with theright-most indel base
                while(lower < j-1 && this.aln[i] == this.aln[j-1]) { 
                    j--;
                }
                // while we have a base to shift
                k = i + 1; // right non-indel base
                // shift over while valid shifts exist
                while(lower < j 
                        && (ALN_MATCH == this.aln[k] || ALN_MISMATCH == this.aln[k])
                        && this.tseq[k] == 0
                        && this.flowOrder[k] == this.flowOrder[j]) 
                {
                    int tmpInt;
                    byte tmpByte;
                    char tmpChar;

                    // swap
                    //System.err.println("SWAPING k=" + k + " j=" + j + " " + this.qseq[k] + " with " + this.qseq[j]);
                    //tmpInt = this.qseq[k]; this.qseq[k] = this.qseq[j]; this.qseq[j] = tmpInt;
                    tmpInt = this.tseq[k]; this.tseq[k] = this.tseq[j]; this.tseq[j] = tmpInt;
                    //tmpByte = this.flowOrder[k]; this.flowOrder[k] = this.flowOrder[j]; this.flowOrder[j] = tmpByte;
                    tmpChar = this.aln[k]; this.aln[k] = this.aln[j]; this.aln[j] = tmpChar;

                    j++;
                    k++;
                    adjusted = true;
                }
                // update past the indel
                i = j - 1;
            }
            else {
                i--;
            }
        }

        return adjusted;
    }
    
    /**
     * Debugging print function.
     * @param out the output stream.
     */
    public void print(PrintStream stream) 
    {
        this.print(stream, 150, true);
    }
    
    public void print(PrintStream stream, int colWidth, boolean align) 
    {
        stream.println(this.getAlignmentString(colWidth, align));
    }

    private int charWidth(int val)
    {
        if(val <= 0) {
            return 1;
        }
        else {
            return (1 + (int)Math.log10(val));
        }
    }
    
    public String getAlignmentString(int colWidth, boolean align) 
    {
        int i, j, numChars = 0;
        StringBuilder tseq = new StringBuilder();
        StringBuilder qseq = new StringBuilder();
        StringBuilder aln = new StringBuilder();
        StringBuilder flowOrder = new StringBuilder();
        StringBuilder string = new StringBuilder();

        // get the widths
        for(i=0;i<this.length;i++) {
            int w, maxWidth = 1;
            w = this.charWidth(this.qseq[i]); 
            if(maxWidth < w) {
                maxWidth = w;
            }
            w = this.charWidth(this.tseq[i]); 
            if(maxWidth < w) {
                maxWidth = w;
            }
            if(colWidth < numChars + maxWidth + 1) {
                string.append(qseq.toString() + "\n" + aln.toString() + "\n" + tseq.toString() + "\n" + flowOrder.toString() + "\n"); 
                // reset
                tseq = new StringBuilder();
                qseq = new StringBuilder();
                aln = new StringBuilder();
                flowOrder = new StringBuilder();
                numChars = 0;
            }
            else if(0 < i) {
                qseq.append(",");
                aln.append(",");
                tseq.append(",");
                flowOrder.append(",");
                numChars++;
            }
            // qseq
            if(align) {
                for(j=this.charWidth(this.qseq[i]);j<maxWidth;j++) {
                    qseq.append(" ");
                }
            }
            if(ALN_DEL == this.aln[i]) {
                qseq.append(ALN_DEL);
            }
            else {
                qseq.append(this.qseq[i]);
            }
            // aln
            if(align) {
                for(j=1;j<maxWidth;j++) {
                    aln.append(" ");
                }
            }
            aln.append(this.aln[i]);
            // tseq
            if(align) {
                for(j=this.charWidth(this.tseq[i]);j<maxWidth;j++) {
                    tseq.append(" ");
                }
            }
            if(ALN_INS == this.aln[i]) {
                tseq.append(ALN_INS);
            }
            else {
                tseq.append(this.tseq[i]);
            }
            // flow order
            if(align) {
                for(j=1;j<maxWidth;j++) {
                    flowOrder.append(" ");
                }
            }
            flowOrder.append(SamToFlowSpaceUtil.DNA[this.flowOrder[i % this.flowOrder.length]]);
            numChars += maxWidth;
        }
        string.append(qseq.toString() + "\n" + aln.toString() + "\n" + tseq.toString() + "\n" + flowOrder.toString()); 
        return string.toString();
    }

    public String getAlignmentString() 
    {
        return this.getAlignmentString(Integer.MAX_VALUE, false);
    }
    
    public int getScore() {
        return this.score;
    }
}
