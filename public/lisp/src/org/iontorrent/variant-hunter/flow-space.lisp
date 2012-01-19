;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

(in-package :cl-user)

;; requires parse-util.lisp
;; requires sam-parse.lisp

(defvar *flow-order* "TACGTACGTCTGAGCATCGATCGATGTACAGC")

(defun make-base-calls-from-na-n-counts (nas counts)
  (let (seq)
    (mapcar #'(lambda (na count)
		(dotimes (i count)
		  (push na seq)))
	    nas counts)
    (format nil "~{~a~^~}" (reverse seq))))

(defun make-base-calls-from-flow-seq (flow-seq &optional (flow-order *flow-order*))
  (let ((flow-number 0)
	nas
	counts-of-each-na
	residuals-of-each)
  (dolist (flow-intensity (cdr flow-seq))
    (let ((na (char flow-order (mod flow-number (length flow-order)))))
      (push na nas)
      (multiple-value-bind (count residual)
	  (round-half-up flow-intensity 100)
	(push count counts-of-each-na)
	(push residual residuals-of-each)))
    (incf flow-number))
  (setq nas (reverse nas))
  (setq counts-of-each-na (reverse counts-of-each-na))
  (setq residuals-of-each (reverse residuals-of-each))
  (let ((base-calls (make-base-calls-from-na-n-counts nas counts-of-each-na)))
    (values base-calls
	    nas
	    counts-of-each-na
	    residuals-of-each))))

(defun check-and-incf-base (cur-base flow-calls)
  (if (eql (car (car flow-calls)) cur-base)
      (incf (cdr (car flow-calls)))
      (push (cons cur-base 1) flow-calls))
  flow-calls)

(defun convert-base-seq-into-flow-space (base-seq &optional (start-flow-number 0)
					 (flow-order *flow-order*))
  (setq base-seq (string-upcase base-seq))
  (let ((continue? t)
	(current-flow start-flow-number)
	(current-ref-num 0)
	flow-calls)
    (while continue?
      (let ((cur-flow-base
	     (char flow-order (mod current-flow (length flow-order))))
	    (current-ref-base
	     (char base-seq current-ref-num)))
	(if (eql cur-flow-base current-ref-base)
	    (progn  ;; match so consume a ref base
	      (if (eql (car (car flow-calls)) cur-flow-base)
		  (incf (cdr (car flow-calls)))
		  (push (cons cur-flow-base 1) flow-calls))
	      (incf current-ref-num))
	    (if (find current-ref-base "ATGC")
	      (progn  ;; no match, so consume a flow base
		(unless (eql cur-flow-base (car (car flow-calls)))
		  (push (cons cur-flow-base 0) flow-calls))
		(incf current-flow))
	      (incf current-ref-num))) ;; skip ref base
	(setq continue? (< current-ref-num (length base-seq)))))
    (reverse flow-calls)))
;; Example format of result is
;; ((#\T . 0) (#\A . 1) (#\C . 0) (#\G . 0) (#\T . 1) (#\A . 0) (#\C . 1) (#\G . 2))

(defun convert-base-seq-to-hp-count-seq (base-seq)
  (let (hp-count-seq)
    (dotimes (x (length base-seq))
      (let ((base
	     (char base-seq x)))
	(if hp-count-seq
	    (let ((last-hp-count (pop hp-count-seq)))
	      (if (eql (car last-hp-count) base)
		  (progn ;; same base as last
		    (incf (cdr last-hp-count))
		    (push last-hp-count hp-count-seq))
		  (progn
		    (push last-hp-count hp-count-seq)
		    (push (cons base 1) hp-count-seq))))
	    (push (cons base 1) hp-count-seq)
	    )))
    (reverse hp-count-seq)))

(defun get-flow-counts (flow-number flow-seq)
  ;; gets the nth flow counts, where n is 0-based
  (cdr (nth flow-number flow-seq)))


(defun combine-flow-signals-w-flow-base (flow-signals &optional (flow-order *flow-order*)
					 (offset 0))
  (let ((counter (- offset 1)))
    (mapcar #'(lambda (flow-signal)
		(incf counter)
		(let ((base (char flow-order (mod counter (length flow-order)))))
		  (cons base flow-signal)))
	    (cdr flow-signals))))

;; Examples of ref segments		    
;; Deletion
;; read aaaaatcgtcaat----tactggcaaaccaaatccagcaacacatca
;; ref  aaaaatcgtcaataaaatactggcaaaccaaatccagcaacacatca
;; Need to add remove flows from ref

;; Insertion
;; read aaaaatcgtcaataaaagacgtactggcaaaccaaatccagcaacacatca
;; ref  aaaaatcgtcaataaaa----tactggcaaaccaaatccagcaacacatca
;; Need to add flows into ref

(defun list-diff (list1 list2)
  ;; diff will be as long as the shortest list
  (mapcar #'(lambda (item1 item2)
	      (- item1 item2))
	  list1 list2))
;;  Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

;; Nil's flow space alignment from FlowSpaceAlignment.java
;; https://iontorrent.jira.com/wiki/display/TM/An+algorithm+for+converting+alignments+to+flow+space
;; java @ trap/src/trap/flowspace/FlowSpaceAlignment.java

#|  Using keywords for these, instead of these fake enums
(defvar *from-m* 0)
(defvar *from-i* 1)
(defvar *from-d* 2)
(defvar *from-me* 3)
(defvar *from-ie* 4)
(defvar *from-mp* 5)
(defvar *from-ip* 6)
(defvar *from-ip* 6)
(defvar *from-s* 7)
|#

(defvar *align-del*      #\- )
(defvar *align-ins*      #\+ )
(defvar *align-match*    #\| )
(defvar *align-mismatch* #\. )

(defvar *minor-inf* -1000000)

(defclass flow-seq ()
  ;;Represents a base sequence in flow space
  ((GAP :accessor GAP :initarg :GAP)
   (flow-seq :accessor flow-seq :initarg :flow-seq)))

(defun make-flow-seq (base-seq &optional
		      (flow-order *flow-order*)
		      (start-flow-index 0)
		      (gap #\- ))
  (make-instance
   'flow-seq
   :gap gap
   :flow-seq (convert-base-seq-into-flow-space base-seq start-flow-index flow-order)))

(defgeneric get-flow-length (flow-seq))
(defmethod get-flow-length((flow-seq-obj flow-seq))
  (length (flow-seq flow-seq-obj)))

(defclass flow-order ()
  ((flow-order :accessor flow-order :initarg :flow-order)
   (key-seq :accessor key-seq :initarg :key-seq)
   (jump-fwd :accessor jump-fwd :initarg :jump-fwd)
   (jump-rev :accessor jump-rev :initarg :jump-rev)))

(defgeneric rotate (flow-order start-flow-index))
(defmethod rotate ((flow-order-obj flow-order) start-flow-index)
  (with-slots (flow-order)
      flow-order-obj
    (unless (eql 0 start-flow-index)
      (let ((tmp-flow-order (copy-seq flow-order)))
	(dotimes (i (length flow-order))
	  (let ((j (mod (+ i start-flow-index) (length flow-order))))
	    (setf (char flow-order i) (char tmp-flow-order j)))))))
  flow-order-obj)

(defgeneric create-jump-tables (flow-order))
(defmethod create-jump-tables ((flow-order-obj flow-order))
  (with-slots (flow-order jump-fwd jump-rev)
      flow-order-obj
    (setq jump-fwd (make-array (length flow-order)))
    (setq jump-rev (make-array (length flow-order)))
    (dotimes (i (length flow-order))
      (let ((k 1)
	    (j (mod (1+ i) (length flow-order))))
	(while (and (not (eql (char flow-order i)
			      (char flow-order j)))
		    (not (eql i j)))
	  (setq j (mod (1+ j) (length flow-order)))
	  (incf k))
	(setf (aref jump-fwd i) k)
	(setf (aref jump-rev j) k)))))

(defun remove-homopolymers (seq)
  (let (last-na)
    (with-output-to-string (strm)
      (dotimes (x (length seq))
	(let* ((na (char seq x))
	       (na-up (char (string-upcase na) 0)))
	  (unless (eql na-up last-na)
	    (format strm "~a" na))
	  (setq last-na na-up))))))
	   
;; with key sequence, flow order is circular rotated, i.e.
;; TACGTACGTCTGAGCATCGATCGATGTACAGC =(key=TCAG)=>
;;        GTCTGAGCATCGATCGATGTACAGCTACGTAC"

(defun make-flow-order (flow-order-seq &key key-seq remove-hp)
  (when remove-hp
    (setq flow-order-seq (remove-homopolymers flow-order-seq)))
  (let ((j 0)
	(flow-order-obj
	 (make-instance 'flow-order
			:flow-order (copy-seq flow-order-seq)
			:key-seq (copy-seq key-seq))))
    (when (and key-seq (not (eql 0 (length key-seq))))
      (dotimes (i (length key-seq))
	(while (not (eql (char key-seq i)
			 (char flow-order-seq j)))
	  (setq j (mod (1+ j) (length flow-order-seq)))))
      (rotate flow-order-obj j))
    (create-jump-tables flow-order-obj)
    flow-order-obj))

(defgeneric make-gap-sums (flow-order flow-query-seq))
(defmethod make-gap-sums ((flow-order-obj flow-order) flow-query-seq)
  (with-slots (flow-order jump-rev)
      flow-order-obj
    (let ((gap-sums-i (make-list (length flow-query-seq)))
	  (i 0)
	  k)
      (mapcar #'(lambda (item)
		  item
		  (setq k (mod i (length flow-order)))
		  (let ((total 0)
			(j (if (< i (aref jump-rev k))
			       0
			       (- i  (aref jump-rev k)))))
		    (while (< j i)
		      (incf total (cdr (nth j flow-query-seq)))
		      (incf j))
		    (incf i)
		    total))
	      gap-sums-i))))

(defclass flowspace-alignment-cell ()
  ((match-score :accessor match-score :initarg :match-score)
   (ins-score :accessor ins-score :initarg :ins-score)
   (del-score :accessor del-score :initarg :del-score)
   (match-from :accessor match-from :initarg :match-from)
   (ins-from :accessor ins-from :initarg :ins-from)
   (del-from :accessor del-from :initarg :del-from)))

(defclass flow-space-aligner ()
  ((bam :accessor bam :initarg :bam)
   (dp-matrix :accessor dp-matrix :initarg :dp-matrix)
   (query-flow-seq :accessor query-flow-seq :initarg :query-flow-seq)
   (target-flow-seq :accessor target-flow-seq :initarg :target-flow-seq)
   ;;qseqFlowOrder
   (query-flow-order :accessor query-flow-order :initarg :query-flow-order)
   (target-flow-order :accessor target-flow-order :initarg :target-flow-order)
   (gap-sums-i :accessor gap-sums-i :initarg :gap-sums-i)
   (dimension1 :accessor dimension1 :initarg :dimension1)
   (dimension2 :accessor dimension2 :initarg :dimension2)

   ;; Score (sum of the residauls off from target)
   ;; Negative values, closer to 0 is closer to target
   ;; Only in the java version
   ;; (align-score :accessor align-score :initform nil)
   ;; This is stored in the bam record

   ;; These are lists produced during alignment (so no initarg)
   ;; Not really flow order, but rather the flow bases for the entire alignment
   (flow-order :accessor flow-order :initform nil)
   ;; Note, the above flow-order can have extra bases to represent deleted bases
   ;;       present in the reference, but missing from the read

   ;; Query counts is a list of intensities / 100.
   (query-counts :accessor query-counts :initform nil)
   ;; Target is the number of reference bases
   (target-counts :accessor target-counts :initform nil)
   (symbols :accessor symbols :initform nil) ;; gap/match/mismatch symbols

   ;; Query-counts and symbols both can be corrected.  If this happens, orig-query-counts
   ;; will be filled in with the original values.
   (orig-query-counts :accessor orig-query-counts :initform nil)
   (orig-symbols :accessor orig-symbols :initform nil)

   ;; This is calculated from the flow-order tag and the actual flow order
   (flow-numbers :accessor flow-numbers :initform nil)

   ;; Positions can actually shift during flow space alignment, so store resulting target positions
   ;; From the java version
   (target-start-pos :accessor target-start-pos :initform nil)
   (target-end-pos :accessor target-end-pos :initform nil)
   ;; Note, original java start pos is 0-based and includes that position and
   ;;       the original java end pos is 0-base but excludes that position.
   ;; Here, both the start and end pos are 1-base, inclusive.
   ))

;; Helper function for printing
(defun break-residuals-into-tenths-hundreths (query-counts target-counts)
  (let (rounded-units
	;;the rest represent the residuals
	signs
	tenths
	hundreths
	remainders)
    (mapc #'(lambda (query-count target-count)
	      (multiple-value-bind (factor1 remainder1)
		  (round-half-up (- target-count query-count) 1)
		(multiple-value-bind (factor2 remainder2)
		    (floor (abs remainder1) 1/10)
		  (multiple-value-bind (factor3 remainder3)
		      (floor remainder2 1/100)
		    (push factor1 rounded-units)
		    (push (cond ((eql remainder1 0) #\Space)
				((< remainder1 0)   :+ )
				(t                  :-))
			  signs)
		    (push factor2 tenths)
		    (push factor3 hundreths)
		    (push remainder3 remainders)
		    ))))
	  (reverse query-counts)
	  (reverse target-counts))
    (values rounded-units signs tenths hundreths remainders)))

(defgeneric print-alignment (aligner &optional residuals? strm))
(defmethod print-alignment ((aligner flow-space-aligner) &optional residuals?
			    (strm t))
  (with-slots (flow-order query-counts target-counts symbols)
      aligner
    (format strm "F: ~{~a~^ ~}~%" flow-order)
    (format strm "Q: ~{~a~^ ~}~%" (mapcar #'round query-counts))
    (format strm "   ~{~a~^ ~}~%" symbols)
    (format strm "T: ~{~a~^ ~}~%" target-counts)
    (with-slots (orig-query-counts orig-symbols)
	aligner
      (when orig-query-counts
	(format strm "Original alignment:~%")
	(format strm "Q: ~{~a~^ ~}~%" (mapcar #'round orig-query-counts))
	(format strm "   ~{~a~^ ~}~%" orig-symbols)
	(format strm "T: ~{~a~^ ~}~%" target-counts)))
    (when residuals?
      (multiple-value-bind (rounded-units signs tenths hundreths)
	  (break-residuals-into-tenths-hundreths query-counts target-counts)
	(format strm "   ~{~a~^ ~}~%" (mapcar #'(lambda (item)
					 (cond ((eql item 0) #\Space)
					       ((< item 0) :- )
					       ((> item 0) :+ )))
					      rounded-units))
	(format strm "   ~{~a~^ ~}~%" (mapcar #'abs rounded-units))
	(format strm "R: ~{~a~^ ~}~%" signs)
	(format strm "R: ~{~a~^ ~}~%" tenths)
	(format strm "R: ~{~a~^ ~}~%" hundreths)
	)
      )
    (when residuals?  ;; probably should be another option
      (format strm "ALIGN FLOW INTENSITIES: ")
      (let ((is-first? t))
	(mapc #'(lambda (query-count symbol)
		  (unless (eql symbol *align-del*)
		    (format strm (if is-first? "~a" ",~a")
			    (* 100 query-count)))
		  (setq is-first? nil))
	      query-counts symbols))
      (format strm "~%")
      )))

;; Working with large intensity values which are pretty inaccurate.  The
;; routine below takes in a maximum flow intensity value, checks the reference,
;; and makes a choice whether to correct it to reference.

;; Helper function to define max. intensities based on which base
(defun make-base-max-intensity-hash (max-intensity-per-base-string)
  "Parses string in format of A:500,T:600,G:650,C:700 and makes max-intensity-for-base-hash"
  (let ((bases-intensities (parse-string max-intensity-per-base-string #\,))
	(max-intensity-for-base-hash (make-hash-table)))
    (dolist (base-intensity bases-intensities)
      (let* ((parts (parse-string base-intensity #\:))
	     (base (char (car parts) 0))
	     (max-intensity (parse-integer (second parts))))
	(setf
	 (gethash base max-intensity-for-base-hash)
	 max-intensity)))
    max-intensity-for-base-hash))

;; Timing, took 10.3 sec to do a million non-altering iterations, & 14.3, altering
;; on a first gen i7 2.33GHz.  Most alignments will be non-altering.
(defgeneric alter-alignment-for-large-intensities (aligner max-intensity-for-base-hash))
(defmethod alter-alignment-for-large-intensities ((aligner flow-space-aligner) max-intensity-for-base-hash)
  (unless max-intensity-for-base-hash
    (return-from alter-alignment-for-large-intensities aligner))
  (with-slots (flow-order query-counts target-counts symbols orig-query-counts orig-symbols)
      aligner
    ;; Note, flow-order is really the bases in the alignment, see notes in defclass.
    (when orig-query-counts
      ;; Ensure always working with the original alignment values
      (setq query-counts orig-query-counts)
      (setq symbols orig-symbols)
      (setq orig-query-counts nil)
      (setq orig-symbols nil))
    (let (has-correction?
	  new-query-counts
	  new-symbols)
      (mapc #'(lambda (flow-base query-count target-count symbol)
		(let ((max-intensity (gethash flow-base max-intensity-for-base-hash)))
		;;(when (eql *align-mismatch* symbol)
		;;  (format t "Max: ~a, Flow: ~a ~a ~a ~a~%"
		;;	  max-intensity
		;;	  flow-base query-count target-count symbol))

		  ;; There's uncertanty in this count, and there's a mismatch
		  ;; between target and query, then take reference.
		(when (and (eql *align-mismatch* symbol)
			   (> (* 100 query-count) max-intensity)
			   (> (* 100 (+ target-count 2)) max-intensity))		   
		  (setq has-correction? t)

		  ;; set query count to reference, but with a very high deviation
		  (setq query-count (- target-count 49/100))
		  (setq symbol #\|)
		  )
		(push query-count new-query-counts)
		(push symbol new-symbols)))
	    flow-order query-counts target-counts symbols)
      (when has-correction?
	(setq orig-query-counts query-counts)
	(setq orig-symbols symbols)
	(setq query-counts (reverse new-query-counts))
	(setq symbols (reverse new-symbols)))))
  aligner)

;; Check to see if returned target sequence and positions for SamToFlowSpace.jar
;; are consistent to the actual reference.
;;
;; This check is pretty expensive, so probably should do only if
;; there are issues found in the flow space alignments.
(defgeneric print-n-check-target-seq-n-positions (aligner &optional strm warn-strm))
(defmethod print-n-check-target-seq-n-positions ((aligner flow-space-aligner) &optional (strm nil) (warn-strm t))
  (with-slots (flow-order target-counts target-start-pos target-end-pos bam)
      aligner
    (let (chrom
	  strand2
	  (number-align-target-bases 0)
	  align-target-seq
	  target-seq-ref
	  bam-read-name
	  (align-strm (make-string-output-stream :element-type 'character))
	  found-issue? did-print?
	  hp-bases-adjustment
	  )
      (when bam
	(with-slots (strand ref-name)
	    bam
	  (setq chrom ref-name)
	  (setq strand2 strand))
	(with-slots (ref-seq-pos ref-seq reverse-complement-ref read-name)
	    bam
	  (setq bam-read-name read-name)
	  ;;(format align-strm "~a~%" read-name)
	  (format align-strm "~a(~a):~a-~a FWD: ~a~%" chrom strand2 ref-seq-pos (+ ref-seq-pos (length ref-seq) -1) ref-seq)
	  (format align-strm "~a(~a):~a-~a REV: ~a~%" chrom strand2 ref-seq-pos (+ ref-seq-pos (length ref-seq) -1) reverse-complement-ref)  ))
      (setq align-target-seq
	    (with-output-to-string (t-strm)
	      (mapc #'(lambda (flow-base target-counts)
			(incf number-align-target-bases target-counts)
			(dotimes (x target-counts)
			  (format t-strm "~a" flow-base)))
		    flow-order target-counts)))
      (format align-strm "~a(~a):~a-~a JFS: ~a~%" chrom strand2 target-start-pos target-end-pos align-target-seq)
      (when bam
	(with-slots (ref-seq-pos ref-seq)
	    bam
	  (setq target-seq-ref (safe-subseq ref-seq (- target-start-pos ref-seq-pos)
					    (+ (- target-start-pos ref-seq-pos) (- target-end-pos target-start-pos -1))))
	  (when (eql strand2 :- )
	    (setq target-seq-ref (reverse-complement target-seq-ref)))
	  (if (eql strand2 :+ )
	      (format align-strm "~a(~a):~a-~a FWD: ~a~%" chrom strand2 target-start-pos target-end-pos target-seq-ref)
	      (format align-strm "~a(~a):~a-~a REV: ~a~%" chrom strand2 target-start-pos target-end-pos target-seq-ref))))
      (close align-strm)

      ;; Sometimes there's not enough reference stored in the bam-record orbject, so can't do 
      ;; a full comparison, but can still check the other bases.  So, to do this:
      ;; Replace ? in ref with bases in Flow space alignment
      (setq target-seq-ref (replace-end-?s-w-subst target-seq-ref align-target-seq))

      ;; Now some checks
      (unless (equal target-seq-ref align-target-seq) ;;(eql number-align-target-bases (1+ (- target-end-pos target-start-pos)))
	(setq found-issue? t)
	(format warn-strm "For ~a, target range ~a-~a seq w/ length of ~a is _not_the_same_ as alignment target seq w/ length of ~a.~%"
		(if bam
		    (get-basic-attrib-string bam)
		    bam-read-name)
		target-start-pos target-end-pos
		(1+ (- target-end-pos target-start-pos))
		number-align-target-bases)
	(setq hp-bases-adjustment (- number-align-target-bases
				     (1+ (- target-end-pos target-start-pos))
				     ))
	(when t ;;(eql 0 (random 10)) ;; limit output to 1 in 25.
	  (setq did-print? t)
	  (format warn-strm "Details:~%~a~%" (get-output-stream-string align-strm))
	  (format warn-strm "Alignment:~%")
	  (print-alignment aligner t warn-strm))
	)
      (when strm
	(format strm "~a" (get-output-stream-string align-strm)))
      (values
       found-issue? hp-bases-adjustment did-print?))))
;;(find bam-read-name '("XOU43:999:717" "XOU43:648:504" "XOU43:966:1209" "XOU43:93:201" "XOU43:966:804") :test #'string=)

#|
;; if got position correct, no need to run this function
(defmethod alter-alignment-target-positions ((aligner flow-space-aligner) hp-bases-adjustment)
  (with-slots (bam target-start-pos target-end-pos target-counts)
      aligner
    (when (eql :- (strand bam))
      (decf target-start-pos hp-bases-adjustment))
    ))
|#

;; The flowing is the LISP implementation of SamToFlowSpace.jar.
;;
;; This makes an aligner object for a particular read
;; which will do the alignment in flow space.
(defun make-flow-space-aligner (query-flow-seq target-base-seq query-flow-order-seq key-seq
				&optional bam)
  (let (fs-align-cell-array
	query-flow-order
	gap-sums
	query-length
	target-length ;; (length target-flow-seq))
	target-flow-order
	target-flow-seq )

    ;; Setup the query properties
    (setq query-flow-order (make-flow-order query-flow-order-seq :key-seq key-seq))
    (setq query-length (length query-flow-seq))

    ;; Setup the target properties
    (setq target-flow-order (make-flow-order target-base-seq
					     :remove-hp t))
    (setq target-flow-seq (flow-seq (make-flow-seq target-base-seq (flow-order target-flow-order))))
    (setq target-length (length target-flow-seq))

    ;; Setup the gap sums
    (setq gap-sums (make-gap-sums query-flow-order query-flow-seq))

    ;; Setup the dp matrix
    (setq fs-align-cell-array (make-array (list (1+ query-length) (1+ target-length))))
    (dotimes (i (1+ query-length))
      (dotimes (j (1+ target-length))
	(let ((fs-align-cell (make-instance 'flowspace-alignment-cell
					    :match-score *minor-inf*
					    :ins-score *minor-inf*
					    :del-score *minor-inf*
					    :match-from :from-s
					    :ins-from :from-s
					    :del-from :from-s)))
	  (setf (aref fs-align-cell-array i j) fs-align-cell))))

    ;; Make object
    (make-instance 'flow-space-aligner
		   :bam bam
		   :dp-matrix fs-align-cell-array
		   :query-flow-seq query-flow-seq
		   :target-flow-seq target-flow-seq
		   :query-flow-order query-flow-order
		   :target-flow-order target-flow-order
		   :gap-sums-i gap-sums
		   :dimension1 (car (array-dimensions fs-align-cell-array))
		   :dimension2 (second (array-dimensions fs-align-cell-array)))))

(defgeneric get-i-from (aligner i))
(defmethod get-i-from ((fs-align flow-space-aligner) i)
  (with-slots (query-flow-order)
      fs-align
    (let* ((k (mod (- i 1) (length (flow-order query-flow-order))))
	   (i-from (if (< i (aref (jump-rev query-flow-order) k))
		       0
		       (- i (aref (jump-rev query-flow-order) k)))))
      i-from)))

;;Commented as // init start cells
(defgeneric init-gap-penalties (aligner phase-penalty))
(defmethod init-gap-penalties ((aligner flow-space-aligner) phase-penalty)
  (with-slots (dp-matrix query-flow-order gap-sums-i dimension1 dimension2)
      aligner
    (setf (match-score (aref dp-matrix 0 0)) 0)
    (loop for i from 1 to (1- dimension1) do
	 (let* ((i-from (get-i-from aligner i)))
	   ;; vertical
	   ;; only allow phasing from an insertion
	   (setf (ins-from (aref dp-matrix i 0)) :from-ip)
	   (if (eql 0 i-from)
	       (setf (ins-score (aref dp-matrix i 0))
		     (- 0 (nth (1- i) gap-sums-i) phase-penalty))
	       (setf (ins-score (aref dp-matrix i 0))
		     (- (ins-score (aref dp-matrix i-from 0)) (nth (1- i) gap-sums-i) phase-penalty)))
	   ))))

;; Will need to update the entire dp matrix
;;         for(i=1;i<=flowQseq.length;i++) { // query
;;            k = (i-1) % qseqFlowOrder.length;
;;            iFrom = ((i < qseqFlowOrder.jumpRev[k]) ? 0 : (i - qseqFlowOrder.jumpRev[k])); 
;;            for(j=1;j<=flowTseq.length;j++) { // target

  

;; Different components of the alignment

(defgeneric calc-horizontal-scores (aligner i j))
(defmethod calc-horizontal-scores ((aligner flow-space-aligner) i j)
  (with-slots (dp-matrix target-flow-seq)
      aligner
    (with-slots (del-score ins-score match-score)
	(aref dp-matrix i (- j 1))
      (let ((a-target-flow (cdr (nth (- j 1) target-flow-seq )))) ;; flowTseq.flow[j-1]
	(if (< del-score match-score)
	    (if (<= ins-score match-score)
		(setf (del-score (aref dp-matrix i j))
		      (- match-score a-target-flow)
		      (del-from (aref dp-matrix i j))
		      :from-m)
		(setf (del-score (aref dp-matrix i j))
		      (- ins-score a-target-flow)
		      (del-from (aref dp-matrix i j))
		      :from-i))
	    (if (<= ins-score del-score)
		(setf (del-score (aref dp-matrix i j))
		      (- del-score a-target-flow)
		      (del-from (aref dp-matrix i j))
		      :from-d)
		(setf (del-score (aref dp-matrix i j))
		      (- ins-score a-target-flow)
		      (del-from (aref dp-matrix i j))
		      :from-i)))))))
#|
// vertical
// four moves:
// 1. phased from match
// 2. phased from ins
// 3. empty from match
// 4. empth from ins
// Note: use the NEXT reference base for flow order matching 
|#
(defgeneric calc-vertical-scores (aligner i j phase-penalty &optional i-from))
(defmethod calc-vertical-scores ((aligner flow-space-aligner) i j phase-penalty &optional i-from)
  (setq i-from (or i-from (get-i-from aligner i)))
  (with-slots (dp-matrix query-flow-seq target-flow-seq query-flow-order target-flow-order
			 gap-sums-i)
      aligner
    (let ((v-score-e *minor-inf*)
	  (v-from-e :from-me)
	  (a-query-flow (cdr (nth (- i 1) query-flow-seq)))) ;;flowQseq.flow[i-1]
      (unless (or (eql j (length target-flow-seq)) ;; no next reference base
		  (eql i 1)  ;;always start with leading phasing
		  ;;(qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length] == tseqFlowOrder.flowOrder[j % tseqFlowOrder.length])) 
		  (eql
		   (char (flow-order query-flow-order)
			 (mod (- i 1) (length (flow-order query-flow-order))))
		   (char (flow-order target-flow-order)
			 (mod j (length (flow-order target-flow-order)))))
		  )
	(with-slots (ins-score match-score)
	    (aref dp-matrix (1- i) j)
	  (if (<= ins-score match-score)
	      (setq v-score-e (- match-score a-query-flow) 
		    v-from-e :from-me)
	      (setq v-score-e (- ins-score a-query-flow)
		    v-from-e :from-ie))
	  ;; Start anywhere in tseq
	  (when (and (eql 1 i)
		     (< (+ v-score-e a-query-flow)
			0))
	    (setq v-score-e (- 0  a-query-flow))
	    (setq v-from-e :from-s))))
      ;;phased from . . .
      (let (v-score-p v-from-p)
	(with-slots (ins-score match-score)
	    (aref dp-matrix i-from j)
	  (if (<= ins-score match-score)
	      (setq v-score-p (- match-score (nth (1- i) gap-sums-i) phase-penalty)
		    v-from-p :from-mp)
	      (setq v-score-p (- ins-score (nth (1- i) gap-sums-i) phase-penalty)
		    v-from-p :from-ip)))

	;; compare empty vs. phased
	(with-slots (ins-score ins-from)
	    (aref dp-matrix i j)
	  (if (<= v-score-p v-score-e) ;; Note: always choose empty over phased
	      (setq ins-score v-score-e
		    ins-from v-from-e)
	      (setq ins-score v-score-p
		    ins-from v-from-p)))))))

;; diagonal
(defgeneric calc-diagonal-scores (aligner i j phase-penalty &optional start-local i-from))
(defmethod calc-diagonal-scores ((aligner flow-space-aligner) i j phase-penalty &optional start-local i-from)
  (setq i-from (or i-from (get-i-from aligner i)))
  (with-slots (dp-matrix query-flow-seq target-flow-seq
			 query-flow-order target-flow-order gap-sums-i)
      aligner
    (unless (eql (char (flow-order query-flow-order)
		       (mod (- i 1)
			    (length (flow-order query-flow-order))))
		 (char (flow-order target-flow-order)
		       (mod (- j 1)
			    (length (flow-order target-flow-order)))))
      ;; out of phase, so skip
      (with-slots (match-score match-from)
	  (aref dp-matrix i j)
	(setq match-score *minor-inf*)
	(setq match-from :from-nowhere))
      (return-from calc-diagonal-scores))
    (let ((score-delta
	   (if (< (cdr (nth (- i 1) query-flow-seq ))
		  (cdr (nth (- j 1) target-flow-seq)))
	       ;; (flowTseq.flow[j-1]-flowQseq.flow[i-1]) 
	       (- (cdr (nth (- j 1) target-flow-seq))
		  (cdr (nth (- i 1) query-flow-seq )))
	       ;; (flowQseq.flow[i-1]-flowTseq.flow[j-1]))
	       (- (cdr (nth (- i 1) query-flow-seq))
		  (cdr (nth (- j 1) target-flow-seq)))
	       ))
	  new-match-score new-match-from)
      ;; Get score from cell diagonal to i/j, and calc new score
      (with-slots (ins-score del-score match-score)
	  (aref dp-matrix (- i 1) (- j 1))
	(if (<= ins-score match-score)
	    (if (<= del-score match-score)
		(setq new-match-score (- match-score score-delta)
		      new-match-from :from-m)
		(setq new-match-score (- del-score score-delta)
		      new-match-from :from-d))
	    (if (<= del-score ins-score)
		(setq new-match-score (- ins-score score-delta)
		      new-match-from :from-i)
		(setq new-match-score (- del-score score-delta)
		      new-match-from :from-d))))
      ;; Update score
      (with-slots (match-score match-from)
	  (aref dp-matrix i j)
	(setq match-score new-match-score)
	(setq match-from new-match-from)

	;;// Start anywhere in tseq
	;;startLocal && 1 == i && dp[i][j].matchScore + s < 0)
	(when (and start-local ;;start-local
		   (eql i 1)
		   (< (+ match-score score-delta)
		      0))
	  (setq match-score (- score-delta))
	  (setq match-from :from-s))))))

(defgeneric calc-score-matrix (aligner phase-penalty &optional start-local?))
(defmethod calc-score-matrix ((aligner flow-space-aligner) phase-penalty 
			      &optional start-local?)
  (with-slots (dp-matrix query-flow-seq target-flow-seq
			 query-flow-order target-flow-order gap-sums-i)
      aligner
    ;; fill in the matrix
    (dotimes (i2 (length query-flow-seq))
      (let* ((i (1+ i2))
	     (i-from (get-i-from aligner i)))
	(dotimes (j2 (length target-flow-seq))
	  (let ((j (1+ j2)))
	    (calc-horizontal-scores aligner i j)
	    (calc-vertical-scores aligner i j phase-penalty i-from)
	    (calc-diagonal-scores aligner i j phase-penalty start-local? i-from)
	    ))))
    ))

(defun make-number-seq (min max &optional (incr 1))
  (let (seq
	(cur-val min))
    (while (<= cur-val max)
      (push cur-val seq)
      (incf cur-val incr))
    (reverse seq)))

;; Takes a look at the matrix's terminus to determine the best score.
;; Must be at the end of the query sequence (i) but can be either
;; end of target or anywhere in the target, if end-local? is t.
(defgeneric calc-best-score (aligner &optional end-local?))
(defmethod calc-best-score ((aligner flow-space-aligner) &optional end-local?)
  (with-slots (dp-matrix query-flow-seq target-flow-seq)
      aligner
    (let ((i-s (list (length query-flow-seq)))
	  (j-s (list (length target-flow-seq)))
	  (best (list (- *minor-inf* 1) :from-s -1 -1)))  ; best = (score, type, i, j)
      (when end-local?
	(setq j-s (make-number-seq 1 (length target-flow-seq))))
      (mapc #'(lambda (i)
		(mapc #'(lambda (j)
			  (with-slots (del-score ins-score match-score)
			      (aref dp-matrix i j)
			    (mapc #'(lambda (score type)
				      (when (<= (car best) score)
					(setf (car best)    score)
					(setf (second best) type )
					(setf (third best)  i    )
					(setf (fourth best) j    )))
				  (list del-score ins-score match-score)
				  (list :from-d :from-i :from-m))))
		      j-s))
	    i-s)
      best)))


(defgeneric add-to-alignment (aligner query-count target-count flow-base))
(defmethod add-to-alignment ((aligner flow-space-aligner) query-count target-count flow-base)
  (let (symbol)
    (with-slots (flow-order query-counts target-counts symbols)
	aligner
      (cond ((eql -1 query-count)  (setq query-count 0
					 symbol *align-del*))
	    ((eql -1 target-count) (setq target-count 0
					 symbol *align-ins*))
	    ((eql (round-half-up query-count)
		  target-count)   (setq symbol *align-match*))
	    (t
	     (setq symbol *align-mismatch*)))
      (push flow-base flow-order)
      (push query-count query-counts)
      (push target-count target-counts)
      (push symbol symbols)
      )
    ))

(defgeneric get-query-flow-order-na (aligner position))
(defmethod get-query-flow-order-na ((aligner flow-space-aligner) position)
  (with-slots (query-flow-order)
      aligner
    (char (flow-order query-flow-order )
	  (mod position
	       (length (flow-order query-flow-order))))))

(defgeneric get-target-flow-order-na (aligner position))
(defmethod get-target-flow-order-na ((aligner flow-space-aligner) position)
  (with-slots (target-flow-order)
      aligner
    (char (flow-order target-flow-order )
	  (mod position
	       (length (flow-order target-flow-order))))))

(defgeneric trace-path-back (aligner best))
(defmethod trace-path-back ((aligner flow-space-aligner) best)
  ;; best = (score, type, i, j
  (pop best) ;; score
  (let ((c-type (pop best))
	(i (pop best))
	(j (pop best))
	next-c-type)
  (with-slots (flow-order query-counts target-counts symbols)
      aligner
    (setq flow-order nil)
    (setq query-counts nil) 
    (setq target-counts nil)
    (setq symbols nil))
  (while (< 0 i)
    ;;(format t "c-type=~a, i=~a, j=~a, next-c-type=~a~%" c-type i j next-c-type)
    (setq next-c-type nil)
    (with-slots (match-from ins-from del-from)
	(aref (dp-matrix aligner) i j)
      (with-slots (query-flow-seq target-flow-seq query-flow-order)
	  aligner
	(cond
	  ((eql c-type :from-m)
	   (setq next-c-type match-from)
	   ;;this.add(flowQseq.flow[i-1], flowTseq.flow[j-1], qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length]);
	   (add-to-alignment aligner
			     (get-flow-counts (- i 1) query-flow-seq)
			     (get-flow-counts (- j 1) target-flow-seq)
			     (get-query-flow-order-na aligner (- i 1)))
	   (decf i)
	   (decf j))

	  ((eql c-type :from-i)
	   (setq next-c-type ins-from)
	   (cond ((or (eql ins-from :from-me)
		      (eql ins-from :from-ie))
		  ;; this.add(flowQseq.flow[i-1], 0, qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length]);
		  (add-to-alignment aligner
				    (get-flow-counts (- i 1) query-flow-seq)
				    0
				    (get-query-flow-order-na aligner (- i 1)))
				    
		  (decf i))
		 ((or (eql ins-from :from-mp)
		      (eql ins-from :from-ip))
		  (let ((i-from (get-i-from aligner i)))
		    (while (< i-from i)
		      (let ((k (mod (- i 1) (length (flow-order (query-flow-order aligner))))))
			k
			;; this.add(flowQseq.flow[i-1], -1, qseqFlowOrder.flowOrder[k]);
			(add-to-alignment aligner
					  (get-flow-counts (- i 1) query-flow-seq)
					  -1
					  (get-query-flow-order-na aligner k))
			(decf i)))))
		 ((eql ins-from :from-s)
		  (while (< 0 i)
		    ;; always a start insertion
		    ;; this.add(flowQseq.flow[i-1], -1, qseqFlowOrder.flowOrder[(i-1) % qseqFlowOrder.length]);
		    (add-to-alignment aligner
				      (get-flow-counts (- i 1) query-flow-seq)
				      -1
				      (get-query-flow-order-na aligner (- i 1)))
		    (decf i)))
		 (t
		  (error "ERROR: dp[i][j].insFrom=~a~%Value not allowed.~%" ins-from))))
	  ((eql c-type :from-d)
	   (setq next-c-type del-from)
	   ;;this.add(-1, flowTseq.flow[j-1], tseqFlowOrder.flowOrder[(j-1) % tseqFlowOrder.length]);
	   (add-to-alignment aligner
			     -1
			     (get-flow-counts (- j 1) target-flow-seq)
			     (get-target-flow-order-na aligner (- j 1)))
	   (decf j))
	  (t
	   (error "ERROR: c-type=~a, next-c-type=~a, i=~a, j=~a~%" c-type next-c-type i j))))
      (cond ((find next-c-type '(:from-m :from-i :from-d :from-s))
	     (setq c-type next-c-type))
	    ((find next-c-type '(:from-me :from-mp))
	     (setq c-type :from-m))
	    ((find next-c-type '(:from-ie :from-ip))
	     (setq c-type :from-i))
	    (t
	     (error "ERROR with next-c-type = ~a~%" next-c-type)))
      ))))

