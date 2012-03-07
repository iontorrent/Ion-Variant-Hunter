/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
package org.iontorrent.sam2flowgram.io;

/**
 * Stores a genomic interval.
 *
 * @author nils.homer@lifetech.com
 */
public class Range {

    /**
     * The zero-based reference index.
     */
	public int referenceIndex = -1;

    /**
     * The one-based start position of this interval.
     */
	public int startPosition = -1;

    /**
     * The one-based end position of this interval.
     */
	public int endPosition = -1;

    /**
     * Creates a new genomic interval.
     * @param referenceIndex the zero-based reference index.
     * @param startPosition the one-based start position of this interval.
     * @param endPosition the one-based end position of this interval.
     */
	public Range(int referenceIndex, int startPosition, int endPosition)
	{
		this.referenceIndex = referenceIndex;
		this.startPosition = startPosition;
		this.endPosition = endPosition;
	}
}

