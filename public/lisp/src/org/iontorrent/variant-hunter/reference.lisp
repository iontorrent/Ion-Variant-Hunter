;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

(in-package :cl-user)

;; requires parse-util.lisp

#|
.fai file is in the format of:
refName	lengSeq  	filePos      	??	??	
chr1	247249719	6       	50	51
chr2	242951149	252194726	50	51
chr3	199501827	500004904	50	51
|#

(defclass fai-entry ()
  ((contig-name :accessor contig-name :initarg :contig-name)
   (length-contig :accessor length-contig :initarg :length-contig)
   (file-pos :accessor file-pos :initarg :file-pos)
   (length-line :accessor length-line :initarg :length-line)
   (length-line-w-cr :accessor length-line-w-cr :initarg :length-line-w-cr)))

(defclass reference-reader ()
  ((fai-entries :accessor fai-entries :initarg :fai-entries)
   (fai-entry-hash :accessor fai-entry-hash :initarg :fai-entry-hash)
   (ref-file :accessor ref-file :initarg :ref-file)
   (fai-file :accessor fai-file :initarg :fai-file)))

(defun make-reference-reader (ref-file)
  (let ((reference-index
	 (make-instance 'reference-reader :ref-file ref-file))
	(idx (format nil "~a.fai" ref-file)))
    (with-slots (fai-file fai-entries fai-entry-hash)
	reference-index
      (setq fai-file idx)
      (setq fai-entries nil)
      (setq fai-entry-hash (make-hash-table :test #'equal))

      (with-delimited-file (fai-line idx)
	(let ((fai-entry (make-instance 'fai-entry
					:contig-name (car fai-line)
					:length-contig (parse-integer (second fai-line))
					:file-pos (parse-integer (third fai-line))
					:length-line (parse-integer (fourth fai-line))
					:length-line-w-cr (parse-integer (fifth fai-line)))))
	  (when (gethash (contig-name fai-entry) fai-entry-hash)
	    (error "Repeated contig name in ~a" idx))
	  (setf (gethash (contig-name fai-entry) fai-entry-hash)
		fai-entry)
	  (push fai-entry fai-entries)))
      (setq fai-entries (reverse fai-entries)))
    reference-index))

(defgeneric get-file-position (ref-reader contig startPos &optional endPos start-padding end-padding))

(defmethod get-file-position ((ref-reader reference-reader) contig startPos &optional endPos
			      (start-padding 0) (end-padding 0))
  (let ((fai-entry (gethash contig (fai-entry-hash ref-reader)))
	num-whitespace-chars)
    (with-slots (length-contig file-pos length-line length-line-w-cr)
	fai-entry
      (when (> (or endPos -1) (length-contig fai-entry))
	(error "ERROR: End position ~a is off '~a' which has a length of ~a.~%" endPos contig (length-contig fai-entry)))
      (when (< startPos 1)
	(error "ERROR: Start position ~a is off '~a' which has a length of ~a.~%" startPos contig (length-contig fai-entry)))
      
      ;; add padding if possible
      (setq startPos (- startPos start-padding))
      (when (< startPos 1)
	(setq startPos 1))
      (setq endPos (+ endPos end-padding))
      (when (> endPos (length-contig fai-entry))
	(setq endPos (length-contig fai-entry)))
	
      (let ((cr-length (- length-line-w-cr length-line)))
	(setq num-whitespace-chars (* (floor (1- startPos) length-line) cr-length)))
      (values
       (+ file-pos startPos num-whitespace-chars -1)
       startPos)
      )))

(defun extract-contig-start-end (position-search)
  (setq position-search (remove #\, position-search))
  (let ((contig-pos (parse-string position-search #\:))
	start-end
	contig
	start-pos
	end-pos)
    (setq contig (car contig-pos))
    (setq start-end (parse-string (second contig-pos) #\-))
    (setq start-pos (parse-integer (car start-end)))
    (setq end-pos (parse-integer (second start-end)))
    (values contig start-pos end-pos)))

(defgeneric extract-reference (ref-reader contig startPos endPos
					  &optional start-padding end-padding))
(defmethod extract-reference ((ref-reader reference-reader) contig startPos endPos
			      &optional (start-padding 0) (end-padding 0))
  (multiple-value-bind (file-pos start-pos-w-pad)
      (get-file-position ref-reader contig startPos endPos
			 start-padding end-padding)
    (with-open-file (refStrm (ref-file ref-reader) :direction :input)
      (file-position refStrm file-pos)
      (let (line
	    (extract-ref "")
	    (continue? t))
	  (while (and continue?
		      (not (eql :eof
				(setq line (read-line refStrm nil :eof)))))
	    (setq extract-ref
		  (concatenate 'string extract-ref line))
	    (when (> (length extract-ref)
		     (+ (- endPos startPos) start-padding end-padding))
	      (setq continue? nil)))
	  (values
	   (subseq extract-ref 0 (min (length extract-ref)
				      (- (+ endPos end-padding) (- startPos start-padding)-1)))
	   start-pos-w-pad
	   )))))

(defgeneric extract-reference-from-position-search (ref-reader position-search))
(defmethod extract-reference-from-position-search ((ref-reader reference-reader) position-search)
  (multiple-value-bind (contig start-pos end-pos)
      (extract-contig-start-end position-search)
    (extract-reference ref-reader contig start-pos end-pos)))

(defun rc-na (na)
  (cond ((eql na #\A) #\T)
	((eql na #\T) #\A)
	((eql na #\C) #\G)
	((eql na #\G) #\C)

	((eql na #\a) #\t)
	((eql na #\t) #\a)
	((eql na #\c) #\g)
	((eql na #\g) #\c)

	((eql na #\.) #\.)
	((eql na #\?) #\?)
	((eql na #\N) #\N)
	(t #\*)))

(defun reverse-complement (seq)
  (let (rc-list)
    (dotimes (pos (length seq))
      (push (rc-na (char seq pos)) rc-list))
    (format nil "~{~a~^~}" rc-list)))
