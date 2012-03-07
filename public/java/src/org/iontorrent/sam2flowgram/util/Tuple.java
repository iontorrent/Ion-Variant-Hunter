/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/*
 * LICENSE to be determined
 */
package org.iontorrent.sam2flowgram.util;

/**
 *
 * @author nils.homer@lifetech.com
 */
public class Tuple<T> {
    public T one;
    public T two;

    /**
     * Creates an empty tuple.
     */
    public Tuple() {
    }

    /**
     * Creates a tuple.
     * @param one the first element of the tuple
     * @param two the second element of the tuple
     */
    public Tuple(T one, T two) {
        this.one = one;
        this.two = two;
    }
}
