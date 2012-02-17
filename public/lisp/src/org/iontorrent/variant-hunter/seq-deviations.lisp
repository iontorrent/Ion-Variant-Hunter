;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

(in-package :cl-user)

;; requires parse-util.lisp
;; requires reference.lisp
;; requires flow-space.lisp
;; requires sam-parse.lisp
;; requires stats.lisp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; On a read by read basis, find sequence deviations   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; find deviations between the query (read) and target (ref) sequences inside the alignment
(defclass seq-deviation ()
  ((ref-name :accessor ref-name :initarg :ref-name :initform nil)
   (bam :accessor bam :initarg :bam)
   ;; (align :accessor align :initarg :align)  ;; FOR DEBUG ONLY!!!!

   (absolute-target-pos :accessor absolute-target-pos :initarg :absolute-target-pos
			:initform nil)
   (rightmost-target-pos :accessor rightmost-target-pos :initarg :rightmost-target-pos)
   (relative-target-pos :accessor relative-target-pos :initarg :relative-target-pos)
   (query-seq :accessor query-seq :initarg :query-seq)
   (target-seq :accessor target-seq :initarg :target-seq)
   (dev-type :accessor dev-type :initarg :dev-type)
   (residual-infos :accessor residual-infos :initarg :residual-infos)
   (full-read-residuals :accessor full-read-residuals)
   (merged-residual-infos :accessor merged-residual-infos :initarg :merged-residual-infos)

   ;; the read's min and max flow number, related to length of key and aligned read
   (min-flow :accessor min-flow :initarg :min-flow)
   (max-flow :accessor max-flow :initarg :max-flow)

   ;; merged-from is a list of seq-deviations this was merged from
   (merged-from :accessor merged-from :initarg :merged-from :initform nil)

   ;; There may be some issue with the deviation, and if so, may not want to use it
   (issue :accessor issue :initarg :issue :initform nil)

   ;; original
   (orig-target-pos :accessor orig-target-pos)
   (orig-query-seq :accessor orig-query-seq)
   (orig-target-seq :accessor orig-target-seq)
   ))

(defun get-indel-type-from-querry-target-seq (current-query current-target)
  (let ((q-lng (length current-query))
	(t-lng (length current-target)))
    (cond ((eql q-lng t-lng) :snp)
	  ((and (eql q-lng 0)
		(< 0 t-lng)) :deletion)
	  ((and (eql t-lng 0)
		(> q-lng 0)) :insertion)
	  ((< q-lng t-lng)   :subst-deletion)
	  ((> q-lng t-lng)   :subst-insertion))))

(defun make-homopolymer-strech (nuc size previous-nucs)
  (with-output-to-string (strm)
    (format strm "~a" (or previous-nucs ""))
    (dotimes (x size)
      x
      (format strm "~a" nuc))))

(defgeneric find-seq-deviations (aligner read &optional aligner-method))
(defmethod find-seq-deviations ((aligner flow-space-aligner) (read null) &optional aligner-method)
  aligner-method ;;ignored
  (let (variants
	current-pos
	current-query
	current-target
	current-residuals
	(relative-target-pos 0)
	last-residual
	flow-range
	full-read-residuals)
    (with-slots (flow-order query-counts target-counts symbols flow-numbers)
	aligner
      (setq flow-range (cons
			(nth (position-if-not #'(lambda (item) (eql item #\-))
					      symbols :from-end nil)
			     flow-numbers)
			(nth (position-if-not #'(lambda (item) (eql item #\-))
					      symbols :from-end t)
			     flow-numbers)))
      (mapc #'(lambda (nuc query-count target-count symbol flow-number)
		(multiple-value-bind (round-query-count residual)
		    (round-half-up query-count)
		  ;; Have extra info associated w/ residual value
		  (setq residual
			(list nuc flow-number
			      symbol
			      (* 100 residual)
			      (* 100 query-count)
			      target-count
			      (* 100 (- query-count target-count))
			      (- round-query-count target-count)
			      (if (not (eql round-query-count target-count))
				  :deviation
				  )))
		  ;; (push residual full-read-residuals)

		  ;; check if there's a sequence deviation
		  (unless (eql round-query-count target-count)
		    ;;(push last-residual current-residuals) ;; get residual before variant
		    ;; Put in last residual before the deviation
		    (unless current-residuals
		      (push last-residual current-residuals))
		    (push residual current-residuals)
		    (setq current-pos (or current-pos relative-target-pos))
		    (setq current-query (make-homopolymer-strech nuc round-query-count current-query))
		    (setq current-target (make-homopolymer-strech nuc target-count current-target)))
		  (when (and (eql round-query-count target-count)
			     current-residuals)
		    (push residual current-residuals))

		  ;; Non-empty flow terminates the sequence deviation at this flow
		  (when (and current-pos
			     (not (eql round-query-count 0))
			     (eql round-query-count target-count))
		    ;; check to see if query and target are the same
		    ;; TODO, shouldn't need this, but workaround for
		    ;; a bug in the LISP realigner
		    (unless (string=  current-target current-query)
		      
		      (push (make-instance 'seq-deviation
					   :bam (bam aligner)
					   ;; :align aligner  ;; FOR DEBUG ONLY!!!!
					   :relative-target-pos current-pos
					   :query-seq current-query
					   :target-seq current-target
					   :dev-type (get-indel-type-from-querry-target-seq
						      current-query current-target)
					   :residual-infos (reverse current-residuals)
					   :min-flow (car flow-range)
					   :max-flow (cdr flow-range))
			    variants))
		    (setq current-pos nil)
		    (setq current-query nil)
		    (setq current-target nil)
		    (setq current-residuals nil))
		  (unless current-pos
		    (setq last-residual residual))
		  )
		(unless (eql 0 target-count)
		  (incf relative-target-pos target-count))
		)
	    flow-order query-counts target-counts symbols flow-numbers))
    (when current-pos
      (push (make-instance 'seq-deviation
			   :bam (bam aligner)
			   ;; :align aligner  ;; FOR DEBUG ONLY!!!!!!
			   :relative-target-pos current-pos
			   :query-seq current-query
			   :target-seq current-target
			   :dev-type (get-indel-type-from-querry-target-seq
				      current-query current-target)
			   :residual-infos (reverse current-residuals)
			   :min-flow (car flow-range)
			   :max-flow (cdr flow-range))
	    variants))
    (setq full-read-residuals (reverse full-read-residuals))
    (dolist (seq-dev variants)
      (setf (full-read-residuals seq-dev) full-read-residuals))
    (reverse variants)))
;; The above relative target positions is taken from the begining of the
;; read.  If the target has no sequence, it is the position immediately after
;; the query insertion.

;; The following function adds absolute position, taking strand into account.
;; For 0 length targets, it also uses the more standard representation of using 
;; the position of the base immediately _before_ the query insertion.
(defgeneric add-target-absolute-position (deviation target-start-pos
						    &optional ref-name-in strand))
(defmethod add-target-absolute-position ((deviation seq-deviation)
					 (target-start-pos number)
					 &optional ref-name-in strand)
  (with-slots  (;;with these and above arguments
		relative-target-pos target-seq
		;; will set these attributes
		ref-name absolute-target-pos
		orig-target-pos ;; preserve what position is before trimming
		)
      deviation
    (setq ref-name ref-name-in)
    (setq absolute-target-pos
	  (if (eql strand :-)
	      (- target-start-pos relative-target-pos (length target-seq) -1)
	      (+ target-start-pos relative-target-pos)))
    (when (and ;;(eql strand :+ )  ;; with neg. strand, will naturally be before insertion
	       (eql (length target-seq) 0))  ;; make it immediately before insertion
      (decf absolute-target-pos))
    (setq orig-target-pos absolute-target-pos)
    )
  deviation)

;; takes in a list of deviations
(defmethod add-target-absolute-position ((deviations cons)
					 (target-start-pos number)
					 &optional ref-name-in strand)
  (mapc #'(lambda (deviation)
	    (add-target-absolute-position deviation
					  target-start-pos
					  ref-name-in
					  strand))
	deviations)
  deviations)

;; Makes query/target sequence all relative to + strand
;; by reverse complementing the - strand hits.
;; Note, residual list is not changed.
(defgeneric make-relative-to-positive-strand (deviation))
(defmethod make-relative-to-positive-strand ((deviation seq-deviation))
  (with-slots (query-seq target-seq bam
			 orig-query-seq orig-target-seq)
      deviation
    (when (eql :- (strand bam))
      (setq query-seq (reverse-complement query-seq))
      (setq target-seq (reverse-complement target-seq))
      )
    (setq orig-query-seq (copy-seq query-seq))
    (setq orig-target-seq (copy-seq target-seq))
    )
  deviation)

;; Note, all deviations in this list need to be on the same strand
(defmethod make-relative-to-positive-strand ((deviations cons))
  (when (< 1 (length (remove-duplicates (mapcar #'(lambda (d) (strand (bam d))) deviations))))
    (error "ERROR: Multiple strands in make-relative-to-positive-strand within the read."))
  (let ((new-deviations
	 (mapcar #'make-relative-to-positive-strand deviations)))
    (if (eql :- (strand (bam (car deviations))))
	(reverse new-deviations)
	new-deviations
	)
    ))

;; Now algorithms that try to make a uniform leftmost call from any alignments
;; that occur.
#|
Insertions:
          12345--45678
Insertion TTAGATTGATTG  
Ref1      TTA----GATTG  GATT ins after TAA
Ref2      TTAGA----TTG  TTGA ins after AGA
ref3      TT--A--GATTG  AG ins after TT, TT insertion after A
|#

(defun merge-residual-infos (residual-infos1 residual-infos2)
  (setq residual-infos1 (remove nil residual-infos1))
  (setq residual-infos2 (remove nil residual-infos2))

  (let ((first-flow1 
	 (second ;; TODO, second representes flow number
	  (car residual-infos1)))
	(first-flow2
	 (second ;; TODO, second representes flow number
	  (car residual-infos2))))
    (if (< first-flow1 first-flow2)
	(append residual-infos1 residual-infos2)
	(append residual-infos2 residual-infos1))))

;; merging of deviations that are on the SAME read.
(defgeneric check-if-can-merge (dev1 dev2))       ;; Checks (close enough?) if can merge
(defgeneric merge-deviations-general (dev1 dev2)) ;; Gets merged sequences

(defgeneric merge-deviations (dev1 dev2))         ;; Makes new merged seq dev if possible
(defmethod  merge-deviations ((dev1 seq-deviation) (dev2 seq-deviation))

  (unless (check-if-can-merge dev1 dev2)
    (return-from  merge-deviations nil))  
  (labels ((make-merged-obj (new-pos new-query-seq new-target-seq dev-type issue)
	     (make-instance 'seq-deviation
			    :ref-name (ref-name dev1)
			    :bam (bam dev1)
			    ;; :align (align dev1) ;; FOR DEBUG ONLY!!!!!!
			    :dev-type dev-type
			    :absolute-target-pos new-pos
			    :query-seq new-query-seq
			    :target-seq new-target-seq
			    :residual-infos nil
			    :merged-residual-infos (merge-residual-infos
						    (or
						     (residual-infos dev1)
						     (merged-residual-infos dev1))
						    (or
						     (residual-infos dev2)
						     (merged-residual-infos dev2)))
			    :min-flow (min-flow dev1)
			    :max-flow (max-flow dev1)
			    :merged-from (append (or (merged-from dev1) (list dev1))
						 (or (merged-from dev2) (list dev2)))
			    :issue issue
			    )))
    #|
    (when (and (eql (dev-type dev1) :deletion)
	       (eql (dev-type dev2) :deletion))
      (return-from  merge-deviations
	(merge-deletions dev1 dev2 #'make-merged-obj)))
    (when (and (find (dev-type dev1) '(:deletion :subst-deletion))
	       (find (dev-type dev2) '(:deletion :subst-deletion)))
      (return-from  merge-deviations
	(merge-subst-deletions dev1 dev2 #'make-merged-obj)))
    |#
    (when (or (issue dev1) (issue dev2))  ;; one dev already had issue
      (return-from merge-deviations nil)) ;; so don't try to merge
    (multiple-value-bind (merged-pos merged-q-seq merged-t-seq issue)
	(merge-deviations-general dev1 dev2)
      (when merged-pos
	(make-merged-obj merged-pos merged-q-seq merged-t-seq :merged issue)))))

#|
For the homozygous indel, CFTR.6.210s 90  rs67140043  GATT        -/-, rs67140043, 
which is slightly incorrect, but there are several correct representations:
          888889999999
          567890123456
Ref:85-96 TTAGATTGATTG
Leftmost  TTA----GATTG  GATT del @ pos 88-91
Rightmost TTAGA----TTG  TTGA del @ pos 90-94
Alt:      TT--A--GATTG  AG and TT del @ pos 87-88 and 90-92, respectively
Alt:      TTAG----ATTG  ATTG del @ pos 89-92
|#

(defvar *allowed-merged-deletion-distance* 4)
(defmethod check-if-can-merge ((dev1 seq-deviation) (dev2 seq-deviation))
  (unless (eql (bam dev1) (bam dev2))
    (error "ERROR, cannot merge seq. deviations from two different bams.  BAM1: ~a, BAM2: ~a"
           (bam dev1) (bam dev2)))
  ;; Should come in order already, but just checking
  (unless (< (absolute-target-pos dev1)
             (absolute-target-pos dev2))
    (format *error-output* "Note, for ~a, seq. deviations did not come in order.~%" (bam dev1))
    (let ((temp-dev2 dev2))
      (setq dev2 dev1)
      (setq dev1 temp-dev2)))

  (when (or (eql (dev-type dev1) :snp)
	    (eql (dev-type dev2) :snp))
    ;; No merging if either one is a SNP
    (return-from check-if-can-merge nil))

  (let* ((end-pos-dev1 (+ (absolute-target-pos dev1)
			  (length (target-seq dev1))
			  -1))
	 ;; distance between the end of dev1 and the start of dev2
	 (space-size (- (absolute-target-pos dev2) end-pos-dev1)))
      ;; Find if two deletions are close enough for merging.
      (when (> space-size  ;; space size too big
               *allowed-merged-deletion-distance*)
	(return-from check-if-can-merge nil))
      t))

(defvar *max-non-empty-flows-for-merge* 3)
(defmethod merge-deviations-general ((dev1 seq-deviation) (dev2 seq-deviation))
  (let ((q-seq1 (query-seq dev1))
	(t-seq1 (target-seq dev1))
	(pos1 (absolute-target-pos dev1))
	(q-seq2 (query-seq dev2))
	(t-seq2 (target-seq dev2))
	(pos2 (absolute-target-pos dev2))
	(i-seq "")
	i-start  ;; Start and end postion of the
	i-end    ;; space between devs.  Fill with ref.
	issue
	)
    #|
    (format t "PB: ~a ~a~%" pos1 pos2)
    |#

    (when (eql (length t-seq1) 0)
      (incf pos1))
    (when (eql (length t-seq2) 0)
      (incf pos2))
    (setq i-start (+ pos1 (length t-seq1)))
    (setq i-end (- pos2 1))
    
    (setq i-seq (extract-reference-from-bam-obj (bam dev1)
						i-start i-end))
    #| DEBUG statements
    (format t "PA:  ~a ~a~%" pos1 pos2)
    (format t "IP: ~a-~a~%" i-start i-end)
    (format t "Q: ~a:~a:~a~%" q-seq1 i-seq q-seq2)
    (format t "T: ~a:~a:~a~%" t-seq1 i-seq t-seq2)
    |#

    ;; Check to see if goes off reference
    ;; Sometimes extraction of ref is inperfect, but also
    ;; happens in certain imposible merging situations.
    (when (eql i-seq :no-reference)
      (setq issue :cannot-merge)
      (setf (issue dev1) issue)
      (setf (issue dev2) issue)
      (return-from merge-deviations-general (values nil nil nil issue)))

    (when (> (length (convert-base-seq-to-hp-count-seq i-seq))
	     *max-non-empty-flows-for-merge*)
      ;; too far appart for merging
      (return-from merge-deviations-general (values nil nil nil issue)))

    (let ((merged-q-seq (format nil "~a~a~a" q-seq1 i-seq q-seq2))
	  (merged-t-seq (format nil "~a~a~a" t-seq1 i-seq t-seq2))
	  (merged-pos pos1))
      (when (eql 0 (length merged-t-seq))
	(decf merged-pos))  ;; empty target, pos is before insertion

      #| DEBUG statements
      (format t "Merged result:~%")
      (format t "Q: ~a~%" merged-q-seq)
      (format t "T: ~a~%" merged-t-seq)
      (format t "Target Position: ~a~%" merged-pos)
      |#
      (values merged-pos merged-q-seq merged-t-seq issue)
      )))

(defun merge-deviation-list (deviations)
  (let ((orig-deviations deviations)
	deviations-w-merged)
    (push (pop orig-deviations) deviations-w-merged)
    (while orig-deviations
      (let ((dev1 (car deviations-w-merged))
	    (dev2 (pop orig-deviations))
	    merged)
	(setq merged (merge-deviations dev1 dev2))
	(when merged
	  ;;(setq *merged-dev-debug* merged)
	  (pop deviations-w-merged) ;; will be replaced with newly merged version
	  (push merged deviations-w-merged))
	(unless merged
	  (unless (issue dev2) ;;dev1 and dev2 had issues, so don't do dev2 again
	    (push dev2 deviations-w-merged)))
	))
    (reverse deviations-w-merged)
    ))

;; Old leftmost function
(defun make-leftmost-old (deviation)
  ;;(setq *test-dev1* deviation)
  (with-slots (dev-type query-seq target-seq bam absolute-target-pos)
      deviation
    (let (del-size ins-size
	  ref-pos
	  (ref-seq (ref-seq bam))
	  (ref-seq-pos (ref-seq-pos bam)))
      (when (eql :deletion dev-type)
	(while (progn
		 (setq del-size (length target-seq))
		 (setq ref-pos (- absolute-target-pos ref-seq-pos 1))
		 (and
		  ;;(< ref-pos (length ref-seq))
		  (>= ref-pos 0)
		  (eql (char target-seq (- del-size 1))
		       (char ref-seq ref-pos))))
	  ;;(format t "~a ~a ~a ~a~%" del-size ref-pos (char target-seq (- del-size 1))
	  ;;	  (char ref-seq ref-pos))
	  (decf absolute-target-pos)
	  (setq target-seq (format nil "~a~a"
				   (char target-seq (1- del-size))
				   (subseq target-seq 0 (1- del-size))))))
      (when (eql :insertion dev-type)
	(while (progn
		 (setq ins-size (length query-seq))
		 (setq ref-pos (- absolute-target-pos ref-seq-pos )) ;; no -1 for insertion
		 ;;(format t "~a ~a ~a ~a~%" ins-size ref-pos (char query-seq (- ins-size 1))
		 ;; (char ref-seq ref-pos)
	  	 ;; )
		 (and (>= ref-pos 0)
		      (< ref-pos (length ref-seq))
		      (eql (char query-seq (- ins-size 1)) ;; last char
			   (char ref-seq ref-pos))))
	  
	  (decf absolute-target-pos)
	  (setq query-seq (format nil "~a~a"
				  (char query-seq (1- ins-size))
				  (subseq query-seq 0 (1- ins-size))))))
    ))
  deviation
  )

;; Trimming
(defun seq-trim-right (query-seq target-seq &optional (target-pos-shift 0))
  (if (and (> (length query-seq) 0)
	   (> (length target-seq) 0)
	   (eql (char query-seq  (- (length query-seq) 1))
		(char target-seq (- (length target-seq) 1))))
      (seq-trim-right (subseq query-seq  0 (- (length query-seq) 1))
		      (subseq target-seq 0 (- (length target-seq) 1)))
      (values query-seq
	      target-seq
	      (if (eql (length target-seq) 0)
		  (- target-pos-shift 1)
		  target-pos-shift))))

(defun seq-trim-left (query-seq target-seq &optional (target-pos-shift 0))
  (if (and (> (length query-seq) 0)  ; check to see if there's sequence left
	   (> (length target-seq) 0)
	   (eql (char query-seq  0)  ; check to see if base is the same
		(char target-seq 0)))
      (seq-trim-left (subseq query-seq 1)  ; now try next base
		     (subseq target-seq 1)
		     (1+ target-pos-shift))
      (values query-seq  ; terminate
	      target-seq
	      (if (eql (length target-seq) 0)
		  (- target-pos-shift 1)
		  target-pos-shift))))

(defun seq-trim (query-seq target-seq &optional (direction-first :right))
  "Trims one direction first, then the other, for a complete trimming"
  (let (first-func
	second-func
	(result-list (list query-seq target-seq 0)))
    (if (eql direction-first :right)
	(setq first-func #'seq-trim-right
	      second-func #'seq-trim-left)
	(setq first-func #'seq-trim-left
	      second-func #'seq-trim-right))
    (when (and (> (length query-seq) 0)     ;; Check if both query and
	       (> (length target-seq) 0)) ;; target have seq left 
      (setq result-list
	    (multiple-value-list (apply first-func result-list))))
    (when (and (> (length (car result-list)) 0)     ;; Check if both query and
	       (> (length (second result-list)) 0)) ;; target have seq left
      (setq result-list
	  (multiple-value-list (apply second-func result-list))))
    (values-list result-list)))

(defgeneric trim-seq-deviation (deviation))
(defmethod trim-seq-deviation  ((deviation seq-deviation))
  (with-slots (query-seq target-seq absolute-target-pos relative-target-pos)
      deviation
    (multiple-value-bind (new-q-seq new-t-seq target-pos-shift)
	(seq-trim (query-seq deviation) (target-seq deviation))
      (setq query-seq new-q-seq)
      (setq target-seq new-t-seq)
      (incf absolute-target-pos target-pos-shift)
      (when (slot-boundp deviation 'relative-target-pos)
	(incf (relative-target-pos deviation) target-pos-shift))))
  deviation)

;; Making lefmost/rightmost
#|
Making leftmost algorithm. . .
-/CGA  --(expand left)--> A/ACGA
A/ACGA --(trim right) --> -/ACG

-/ACG  --(expand left) --> G/GACG
G/GACG --(trim right) --> -/GAC

-/GAC  --(expand left) --> C/CGAC
C/CGAC --(trim right) --> -/CGA

;; Termination condition
-/CGA  --(expand left) --> Q/QCGA
Q/QCGA --(trim right) --> Q/QCGA

;; Cleanup
Q/QCGA --(trim left) --> -/CGA

** Terminate **
|#

(defgeneric make-leftmost (deviation))
(defgeneric set-rightmost-position (deviation))

(defun seq-shift (query-seq target-seq target-bases
		  &optional (direction :left) (number-shifts 0))
  ;; Target bases are a list of bases immediately before (left) ar after (right)
  ;; the current target-seq
  (let (expand-trim
	cleanup-trim
	expand-query-seq
	expand-target-seq)
    (labels ((combine-seq (seq1 seq2)
	       (format nil "~a~a" (or seq1 "") (or seq2 ""))))
      (when (eql direction :left)
	;; Expand left
	(setq expand-query-seq
	      (combine-seq (car target-bases) query-seq))
	(setq expand-target-seq
	      (combine-seq (car target-bases) target-seq))
	(setq expand-trim #'seq-trim-right)
	(setq cleanup-trim #'seq-trim-left)
	)
      (when (eql direction :right)
	;; Expand right
	(setq expand-query-seq
	      (combine-seq query-seq (car target-bases)))
	(setq expand-target-seq
	      (combine-seq target-seq (car target-bases)))
	(setq expand-trim #'seq-trim-left)
	(setq cleanup-trim #'seq-trim-right)
	)
      )

  ;; Trim
  (multiple-value-bind (trim-q-seq trim-t-seq)
      (funcall expand-trim expand-query-seq expand-target-seq)
    (if (or (not (cdr target-bases)) ;; no more target bases to consider
	    ;; or if no further trimming possible
	    (and (eql (length trim-q-seq) (length expand-query-seq))
		 (eql (length trim-t-seq) (length expand-target-seq))))
	;; termination
	(if (or (eql 0 (length query-seq))
		(eql 0 (length target-seq)))
	    (values expand-query-seq expand-target-seq
		    number-shifts t )  ;; made-non-blank? is true
	    (values query-seq target-seq number-shifts))
	;; recurse
	(seq-shift trim-q-seq trim-t-seq (cdr target-bases) direction
		   (1+ number-shifts))))))

(defun get-target-bases-before (ref-seq ref-seq-pos
				target-seq target-seq-pos)
  (let (target-bases)
    (dotimes (x (- target-seq-pos ref-seq-pos
		   (if (eql 0 (length target-seq))
		       -1 0)))
      (when (< x (length ref-seq))
	(push (char ref-seq x) target-bases)))
    target-bases))

(defun get-target-bases-after (ref-seq ref-seq-pos
				target-seq target-seq-pos &optional debug-strm)
  (let (target-bases
	ref-rel-pos)
    (dotimes (x (- (+ (length ref-seq) ref-seq-pos)
		   target-seq-pos
		   (length target-seq)
		   (if (eql 0 (length target-seq))
		       1 0) ;; don't include base before insertion
		   ))
      (setq ref-rel-pos (- (length ref-seq) x 1))
      (when (< ref-rel-pos 0)
	;; check to see if bases go off the availible sequence
	(format debug-strm "Note, requested unavailible reference. ")
	(format debug-strm "[length(ref-seq), ref-seq-pos, target-seq, target-seq-pos] = ")
	(format debug-strm "[~a, ~a, ~a, ~a]~%" (length ref-seq) ref-seq-pos target-seq target-seq-pos)
	(return-from get-target-bases-after nil))
      (push (char ref-seq ref-rel-pos) target-bases))
    target-bases))

(defmethod make-leftmost ((deviation seq-deviation))
  (with-slots (query-seq target-seq absolute-target-pos bam)
      deviation
    (let (num-shifts
	  made-non-blank?)
      (multiple-value-setq (query-seq target-seq num-shifts made-non-blank?)
	(seq-shift query-seq target-seq
		   (get-target-bases-before (ref-seq bam) (ref-seq-pos bam)
					    target-seq absolute-target-pos)
		   :left))
      (when (and made-non-blank?
		 (eql 1 (length query-seq)))
	(decf absolute-target-pos))
      (decf absolute-target-pos num-shifts)))
  deviation
  )

;;(defvar *debug-deviation*)

;; Set rightmost position in deviation (but doesn't change the actual representation)
(defmethod set-rightmost-position ((deviation seq-deviation))
  (with-slots (query-seq target-seq absolute-target-pos bam
			 rightmost-target-pos)
      deviation
    (setq rightmost-target-pos nil)
    (let ((target-bases-after
	   (get-target-bases-after (ref-seq bam) (ref-seq-pos bam)
				   target-seq absolute-target-pos)))
      ;; if goes off the availble reference, leave rightmost position at nil
      (unless target-bases-after
	;;(setq *debug-deviation* deviation)
	(return-from set-rightmost-position deviation))

      (multiple-value-bind (right-query-seq right-target-seq num-shifts made-non-blank?)
	  (seq-shift query-seq target-seq target-bases-after
		     :right)
	;; ignored for now
	right-query-seq right-target-seq

	;; set rightmost positions
	(setf rightmost-target-pos absolute-target-pos)
	(incf rightmost-target-pos (length target-seq)) ;; put it on other side of indel
	(when (eql 0 (length target-seq))
	  (incf rightmost-target-pos))  ;; put it on the right side of insertion
	(unless made-non-blank?
	  (decf rightmost-target-pos))
	(incf rightmost-target-pos num-shifts)
	)))
  deviation
  )

;; From a list of deviations that should have the same rightmost position
(defun get-rightmost-position-from-deviations (deviations)
  (let ((uniq-rightmost-positions
	 (remove nil
		 (remove-duplicates
		  (mapcar #'rightmost-target-pos deviations)))))
    (if (> (length uniq-rightmost-positions) 1)
	(format *error-output* "Note, multiple rightmost positions found from deviations coming from [with a limit of 20] ~{~a~^,~}.~%"
		(mapcar #'get-basic-attrib-string
			(let ((bams
			       (mapcar #'bam deviations)))
			  ;; limit list size to 20
			  (subseq bams 0
				  (min 20 (length deviations))))))
	;; found unique one, so return that
	(return-from get-rightmost-position-from-deviations
	  (car uniq-rightmost-positions)))

    ;; Return the most common read
    (car (car
	  (sort (mapcar #'(lambda (item)
			    (cons
			     item
			     (count item deviations :key #'rightmost-target-pos)))
			uniq-rightmost-positions) #'> :key #'cdr)))))

;; Printing
(defgeneric print-seq-deviation (deviation &optional
					   stream
					   merged-list?
					   print-flow-intensities? ))

(defun print-seq-deviations (seq-deviations &optional
			     stream
			     merged-list?
			     (header? t))
  (when header?
    (when (absolute-target-pos (car seq-deviations))
      (format stream "~a~a~a~a~a~a" 
	      (if merged-list? "Merged?" "")
	      (if merged-list? #\Tab "")
	      "ref" #\Tab "abs-tar-pos" #\Tab))
    (format stream "~a~a~a~a~a~a~a~a~a~%"
	    "rel-tar-pos" #\Tab
	    "strand"      #\Tab
	    "query"       #\Tab
	    "target"      #\Tab
	    "nuc:flow:align:resid:sig:target_bases:query-target:base_diff:has_deviation?"
	    ))
  (mapcar #'(lambda (deviation)
	      (print-seq-deviation deviation stream merged-list?)
	      )
	  seq-deviations)
  t)


(defmethod print-seq-deviation ((deviation seq-deviation) &optional
				stream
				merged-list?
				(print-flow-intensities? t))
  (let ((out-format (format nil "~~a~a~~a~a~~a~a~~a~~{~~{~~a~~^:~~}~~^~a~~}~%"
			    #\Tab  #\Tab  #\Tab  #\Tab)))
    (labels ((dash-for-blank (seq)
	       (if (eql (length seq) 0)
		   "-"
		   seq)))
      (with-slots (bam query-seq
		       target-seq
		       relative-target-pos
		       ref-name
		       absolute-target-pos
		       residual-infos merged-from)
	  deviation
	(when absolute-target-pos
	  (format stream "~a~a~a~a~a~a~a~a~a~a"
		  (if merged-list? (if residual-infos "NotMerged" "IsMerged") "")
		  (if merged-list? #\Tab "")
		  (read-name bam) #\Tab
		  ref-name #\Tab
		  absolute-target-pos #\Tab
		  (strand bam) #\Tab
		  ))
	(format stream out-format
		(if residual-infos
		    relative-target-pos "")
		(dash-for-blank query-seq)
		(dash-for-blank target-seq )
		(if residual-infos "" "MergedFrom:")
		(or residual-infos
		    (mapcar #'(lambda (child)
				(list
				 (dash-for-blank (query-seq child))
				 (dash-for-blank (target-seq child))))
			    merged-from))
		)
	;; recursively print out the children
	(when merged-from
	  (mapcar #'(lambda (child)
		      (print-seq-deviation child stream merged-list? nil))
		  merged-from))
	(when print-flow-intensities?
	  (format stream "FLOW INTENSITIES: ~{~a~^,~}~%" (flow-intensities bam)))
	)
      )))


(defmethod find-seq-deviations ((aligner flow-space-aligner) (bam bam-record) &optional (aligner-method :lisp))
  (let ((read-centric-deviations
	 (find-seq-deviations aligner nil)))
    (unless (eql (bam aligner) bam)
      (format *error-output* "Warning, bam record given to find-seq-deivations doesn't match bam from alignment.  Using bam found within the aligner object.~%"))
    (with-slots (ref-seq-pos ref-seq ref-name strand)
	(bam aligner)
      (if read-centric-deviations
	  (make-relative-to-positive-strand 
	   (add-target-absolute-position read-centric-deviations
					 (case aligner-method
					   (:lisp (if (eql strand :+)
						      ref-seq-pos
						      (+ ref-seq-pos (length ref-seq) -1)))
					   (:java (with-slots (target-start-pos target-end-pos)
						      aligner
						    (if (eql strand :+)
							;;ref-pos and ref-end-pos could be slightly offset during fs alignment
							;;so just use values reported by it
							target-start-pos
							target-end-pos)))
					   (t (error "Error in find seq deviations")))
					 ref-name
					 strand))
	  read-centric-deviations
	  )
      )))

;; Util functions for deviations
(defgeneric seq-deviation< (dev1 dev2))
(defmethod seq-deviation< ((dev1 seq-deviation)
			   (dev2 seq-deviation))
  (cond ((string= (ref-name dev1) (ref-name dev2))
	 (< (absolute-target-pos dev1) (absolute-target-pos dev2)))
	(t
	 (string< (ref-name dev1) (ref-name dev2)))))

(defgeneric sort-deviations (deviations))
(defmethod sort-deviations ((deviations cons))
  (when deviations
    (sort deviations #'seq-deviation<)))


(defgeneric remove-deviations-with-issues (deviations))
(defmethod remove-deviations-with-issues ((deviations cons))
  (remove-if
   #'(lambda (deviation)
       (when (issue deviation)
	 (format *error-output* "WARNING, sequence deviation has issue ~a which was found in ~a~%"
		 (issue deviation) (get-basic-attrib-string (bam deviation)))
	 t))
   deviations)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Variant candidate
;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Now take list of deviations and combine them together to make calls

(defclass variant-candidate ()
  ((ref-name :accessor ref-name :initarg :ref-name)
   (leftmost-pos :accessor leftmost-pos :initarg :leftmost-pos)
   (max-rightmost :accessor max-rightmost :initform nil)
   ;; sequence deviations that were pileuped by position
   (seq-deviations :accessor seq-deviations :initarg :seq-deviations)

   ;; Calculated score for variant candidate
   (score :accessor score :initarg :score)

   (overall-strand-prob :accessor overall-strand-prob :initform nil) ;; binomial probability
   (num-spanning-reads :accessor num-spanning-reads :initform nil)
   (num-spanning-ref-reads :accessor num-spanning-ref-reads :initform nil)
   (num-plus-spanning-reads :accessor num-plus-spanning-reads :initform nil)
   (num-variant-reads :accessor num-variant-reads :initform nil)
   (variant-freq :accessor variant-freq :initform nil)
   (zygosity-call :accessor zygosity-call :initform nil)

   (overall-filtered? :accessor overall-filtered? :initform nil)

   ;;;;;;;;;
   ;; Grouped by actual sequence deviation
   ;;;;;;;;;;;
   ;; hash
   (seq-dev-hash :accessor seq-dev-hash :initarg :seq-dev-hash)
   ;; list of cons(score q-t-seq), sorted by score
   (scores-n-q-t-seqs :accessor scores-n-q-t-seqs :initarg :scores-n-q-t-seqs)
   ;; list of rightmost positions in order of scores-n-q-t-seqs (score)
   (rightmost-positions :accessor rightmost-positions :initform nil)
   ;; list of variant frequences
   (variant-freqs-by-score :accessor variant-freqs-by-score :initform nil)
   ;; strand counts and probabilities
   (strand-counts :accessor strand-counts :initform nil)
   (strand-probs :accessor strand-probs :initform nil)
   (strand-bias-scores :accessor strand-bias-scores :initform nil)
   ;; is this particular variant filtered?
   ;;(variants-filtered? :accessor variants-filtered? :initform nil)

   ))

(defun make-variant-candidate-list-from-seq-deviations (sorted-seq-deviation-list)
  (let (seq-deviation-pileup
	last-ref
	last-pos)
    (dolist (seq-dev sorted-seq-deviation-list)

      ;; get last reference and position from top of pileup
      (when (car seq-deviation-pileup)
	(with-slots (ref-name absolute-target-pos)
	    (car (last (seq-deviations (car seq-deviation-pileup))))
	  (setq last-ref ref-name)
	  (setq last-pos absolute-target-pos)))

      ;; Now check current seq-dev to see if it can be included
      ;; into the current pileup.
      (with-slots (ref-name bam absolute-target-pos)
	  seq-dev
	(if (and (equal ref-name last-ref)
		 (eql absolute-target-pos last-pos))
	    ;; Yes, add a deviation to the current pileup
	    (push seq-dev (seq-deviations (car seq-deviation-pileup)))
	    ;; No, add to a new pileup
	    (push (make-instance 'variant-candidate
				 :ref-name ref-name
				 :leftmost-pos absolute-target-pos
				 :seq-deviations (list seq-dev))
		  seq-deviation-pileup)))
      )
    (reverse seq-deviation-pileup)))

;;;;;;;;;;;;;;;;;;;
;; Variant score/probability calculation
;;;;;;;;;;;;;;;;

;;; Determine a score/likelihood of a single sequence deviation


;; Get elements of the deviation required to find the score
(defgeneric get-signal-deviations (seq-dev))
(defmethod get-signal-deviations ((seq-dev seq-deviation))
  (let (dev-flows
	non-dev-flows
	;;deletion is represented by inserted ref bases into the floworder
	(num-deleted-flows 0)
	;; inserted bases, are skipped flows (aka impossible flows) used to 
	;; align the ref target to the read query sequence.
	(num-inserted-flows 0)
	(list-pos 0)
	)
    (with-slots (residual-infos merged-residual-infos)
	seq-dev
      (dolist (residual-info (or residual-infos merged-residual-infos))

	;; TODO, should figure out why NIL occurs at the begining
	(unless residual-info 
	  (when (eql 0 list-pos)
	    ;; expected
	    (format nil "Warning, a NIL residual info found.~%"))
	  (unless (eql 0 list-pos)
	    (error "A NIL residual info found other than at the begining.~%")))
	(incf list-pos)

	(when residual-info
	  ;; TODO, really should make residual-info an object
	  (let ((symbol (third residual-info))
		(query-target-diff (seventh residual-info))
		(is-dev-from-ref (car (last residual-info))))
	    (when (eql *align-ins* symbol)
	      (incf num-inserted-flows))
	    (if (eql *align-del* symbol) ;; Added base is artifical, so wouldn't want that
		(incf num-deleted-flows) ;;   as part of the distance calculation.
		(if is-dev-from-ref  ;; Separate diff from ref from agreement with ref.
		    (push query-target-diff dev-flows)
		    (push query-target-diff non-dev-flows)))))))
    (values
     num-deleted-flows
     num-inserted-flows
     (euclidean-distance dev-flows)
     (euclidean-distance non-dev-flows)
     (reverse dev-flows)
     (reverse non-dev-flows)
     )))

(defvar *score-coef-delete*)
(defvar *score-coef-insert*)
(defvar *score-coef-intensity*)

(setq *score-coef-delete* 10)
(setq *score-coef-insert* 10)
(setq *score-coef-intensity* 1/10)

;; Calculate the score
(defgeneric calc-deviation-score (num-deleted-flows &optional num-inserted-flows
						    distance-dev-flows distance-nondev-flows))
;; from values
(defmethod calc-deviation-score ((num-deleted-flows number)
				 &optional num-inserted-flows
				 distance-dev-flows distance-nondev-flows)
  (+
   1
   ;; Inserted and deleted items in flow space are strong evidence
   (* *score-coef-delete* num-deleted-flows)
   ;; inserted flows are more numerous, but on average ~ 4 flows to get to a base
   (* *score-coef-insert* (/ num-inserted-flows 4))
   ;; Intensity differences
   (/ (* *score-coef-intensity* (+ distance-dev-flows distance-nondev-flows)) 100)
   0d0))

(defgeneric get-flow-number-range (seq-dev))
(defmethod get-flow-number-range ((seq-dev seq-deviation))
  (let ((ri (or (residual-infos seq-dev)
		(merged-residual-infos seq-dev)))
	low-flow
	high-flow)
    (dolist (i-ri ri)
      (when (second i-ri)
	(unless low-flow
	  (setq low-flow (second i-ri)))
	(setq high-flow (second i-ri))))
    (when low-flow
      (cons low-flow high-flow))))
#|
    (cons
     (or (second (car ri))  ;;TODO bug to fix, somethings this is NIL
	 (second (second ri)))
     (second (find-if #'numberp ri :key #'second :from-end t))) ;;(car (last ri))
    ))
|#


(defgeneric calc-ends-of-read-factor (seq-dev))
(defmethod calc-ends-of-read-factor  ((seq-dev seq-deviation))
  (let ((bp-factor 0.1)
	(max-num-bases-for-factor 4)
	(range (get-flow-number-range seq-dev))
	num-bases-from-end
	)
    (unless range
      (setf (issue seq-dev) :no-flow-intensity-positions)
      (return-from calc-ends-of-read-factor nil))
    (with-slots (min-flow max-flow)
	seq-dev
      (setq num-bases-from-end
	    (min (- (car range) min-flow)
		 (- max-flow (cdr range))))
      (when (< num-bases-from-end 0)
	(error "ERROR: calc-ends-of-read-factor calc to be less than 0 in ~a.~%"
	       seq-dev))
      (if (< num-bases-from-end max-num-bases-for-factor)
	  (* bp-factor num-bases-from-end)
	  1))))

;; from a single deviation object
(defmethod calc-deviation-score ((seq-dev seq-deviation)
				  &optional num-inserted-flows
				 distance-dev-flows distance-nondev-flows)
  ;; ignored variables
  num-inserted-flows distance-dev-flows distance-nondev-flows

  (multiple-value-bind (num-deleted-flows num-inserted-flows
					  distance-dev-flows distance-nondev-flows)
      (get-signal-deviations seq-dev)
    (values
     (* (or (calc-ends-of-read-factor seq-dev) 0)
	(calc-deviation-score num-deleted-flows num-inserted-flows
			      distance-dev-flows distance-nondev-flows))
     num-deleted-flows num-inserted-flows
     distance-dev-flows distance-nondev-flows)))

(defun calc-unified-variant-score (scores)
  (* (average scores) (log (length scores)))
  )

;; from a list of devation objects
(defmethod calc-deviation-score ((seq-deviations cons)
				 &optional num-inserted-flows
				 distance-dev-flows distance-nondev-flows)
  ;; ignored variables
  num-inserted-flows distance-dev-flows distance-nondev-flows

  (let (scores
	num-deleted-flows-lst
	num-inserted-flows-lst
	distance-dev-flows-lst
	distance-nondev-flows-lst)
    (dolist (seq-dev seq-deviations)
      (multiple-value-bind (score num-deleted-flows num-inserted-flows
				  distance-dev-flows distance-nondev-flows)
	  (calc-deviation-score seq-dev)
	(push score scores)
	(push num-deleted-flows num-deleted-flows-lst)
	(push num-inserted-flows num-inserted-flows-lst)
	(push distance-dev-flows distance-dev-flows-lst)
	(push distance-nondev-flows distance-nondev-flows-lst)))
    (values
     (calc-unified-variant-score scores)
     scores
     num-deleted-flows-lst
     num-inserted-flows-lst
     distance-dev-flows-lst
     distance-nondev-flows-lst)))

(defmethod calc-deviation-score ((var-cand variant-candidate)
				  &optional num-inserted-flows
				 distance-dev-flows distance-nondev-flows)

  ;; ignored variables
  num-inserted-flows distance-dev-flows distance-nondev-flows

  (with-slots (seq-deviations ref-name leftmost-pos score)
      var-cand
    ;; debug
    ;;(format t "~a ~a~%" ref-name leftmost-pos)

    (let ((score-values
	   (multiple-value-list (calc-deviation-score seq-deviations))))
      (setq score (car score-values))
      (values-list 
       (reverse (cons leftmost-pos
		      (cons ref-name
			    (reverse
			     score-values))) )))))

;; The next several functions work on a list of sequence deviations, that
;; have uses '(query-seq target-seq) as the key, and has a list of
;; scores and '(query-seq target-seq).

;; Group by sequence, the sequence deviations within a variant candidate
(defgeneric group-sequence-deviations (var-cand))
(defmethod group-sequence-deviations ((var-cand variant-candidate))
  (with-slots (seq-deviations seq-dev-hash scores-n-q-t-seqs)
      var-cand
    (setq seq-dev-hash (make-hash-table :test #'equal))
    (setq scores-n-q-t-seqs nil)
    (dolist (seq-dev seq-deviations)
      (let ((q-t-seq (list (query-seq seq-dev)
			   (target-seq seq-dev))))
	(unless (gethash q-t-seq seq-dev-hash)
	  (push q-t-seq scores-n-q-t-seqs))
	(push seq-dev (gethash q-t-seq seq-dev-hash))
	;;(push (calc-deviation-score seq-dev) uniq-seq-scores)
	))
    (setq scores-n-q-t-seqs
	  (mapcar #'(lambda (q-t-seq)
		      (cons (calc-deviation-score
			     (gethash q-t-seq seq-dev-hash)) q-t-seq))
		  scores-n-q-t-seqs))
    (setq scores-n-q-t-seqs (sort scores-n-q-t-seqs #'> :key #'car))
    t))

;; Set rightmost for the deviations in a variant candidate
;; and then find a single rightmost position for each variant found.
;; Max is max rightmost positon of all unfiltered variants.
(defmethod set-rightmost-position ((var-cand variant-candidate))
  (with-slots (scores-n-q-t-seqs seq-dev-hash
				 seq-deviations
				 rightmost-positions max-rightmost)
      var-cand
    (mapcar #'set-rightmost-position seq-deviations)
    (setq max-rightmost nil)
    (setq rightmost-positions
	  (mapcar #'(lambda (score-n-q-t-seq)
		      ;; get sequence deviations
		      (let ((rightmost
			     (get-rightmost-position-from-deviations
			      (gethash (list (second score-n-q-t-seq)
					     (third score-n-q-t-seq)) seq-dev-hash))))
			(when rightmost
			  (unless max-rightmost
			    (setq max-rightmost rightmost))
			  (when ;;(and (not variant-filtered?)
			      (> rightmost max-rightmost)
			    (setq max-rightmost rightmost))
			  rightmost)))
		  scores-n-q-t-seqs))))

(defun determine-zygosity-call (freq min-het-freq min-hom-freq)
  (cond ((>= freq min-hom-freq)
	 :homozygous)
	((>= freq min-het-freq)
	 :heterozygous)
	(t
	 :no-call)))

(defun zygosity-to-genotype-call (zygosity-call)
  (case zygosity-call
    (:homozygous "1/1")
    (:heterozygous "0/1")
    (t ".")))

;; Set total read counts
(defgeneric calc-spanning-reads-n-freqs (var-cand &optional total-spanning-reads total-plus-spanning-reads spanning-read-names))
(defmethod calc-spanning-reads-n-freqs ((var-cand variant-candidate) &optional total-spanning-reads total-plus-spanning-reads spanning-read-names)
  ;; Number of reads that span both leftmost and rightmost positions
  ;; Optional parameters are used to set the spanning read counts from python call.
  (with-slots (scores-n-q-t-seqs seq-dev-hash
				 num-spanning-ref-reads zygosity-call
				 num-variant-reads
				 num-spanning-reads num-plus-spanning-reads variant-freq variant-freqs-by-score)
      var-cand

    (setq num-variant-reads nil)
    (when (and total-spanning-reads total-plus-spanning-reads)
      (setq num-spanning-reads total-spanning-reads)
      (setq num-plus-spanning-reads total-plus-spanning-reads))
    (let ((total-variant-reads 0)
	  read-name-sets)
      (setq variant-freqs-by-score
	    (mapcar #'(lambda (score-n-q-t-seq)
			(let* ((seq-devs
				(gethash (list (second score-n-q-t-seq)
					       (third score-n-q-t-seq)) seq-dev-hash))
			       (num-reads-for-seq (length seq-devs)))
			  (push (mapcar #'read-name
					(mapcar #'bam seq-devs))
				read-name-sets)
			  (incf total-variant-reads num-reads-for-seq)
			  (push num-reads-for-seq num-variant-reads)
			  num-reads-for-seq))
		    scores-n-q-t-seqs))
      (setq num-variant-reads (reverse num-variant-reads))

      ;; Have info to calc. num-spanning-ref-reads
      ;; Other case is just recalc. after filtering
      (when (and total-spanning-reads total-plus-spanning-reads)
	;; Note, spanning-read-names can be the null set, so can't use that
	;; for seeing which case you are in.
        (setq num-spanning-ref-reads
              (set-difference-size spanning-read-names (apply #'append read-name-sets)
                                   :test #'equal)))

      ;; ok, now have num-spanning-ref-reads set
       (let ((tot-num-reads (+ num-spanning-ref-reads total-variant-reads))
	     (min-het-freq 2/10)
	     (min-hom-freq 7/10)
	     (max-freq 0)
	     )
	 (setq variant-freqs-by-score
	       (mapcar #'(lambda (num-reads-for-seq)
			   (if (eql 0 tot-num-reads)
			       :infinity
			       (let ((freq
				      (/ num-reads-for-seq tot-num-reads)))
				 (when (> freq max-freq)
				   (setq max-freq freq))
				 freq )))
		       variant-freqs-by-score))
	 (setq zygosity-call (determine-zygosity-call max-freq min-het-freq min-hom-freq))
	 (setq variant-freq
	       (if (eql 0  tot-num-reads)
		   :infinity
		   (/ total-variant-reads tot-num-reads)))
	 ))))

;;(defvar *debug-variant-candidate* nil)

(defun calc-strand-bias (num-plus-spanning-reads num-spanning-reads
			 num-plus-variant num-variant)
  ;; very simple strand bias score calculation
  ;; 0 is no bias
  ;; 1 is bias on one, but not biased on other
  ;; 2 is each one is completely on different strands
  (when (and (> num-spanning-reads 10)
	     (> num-variant 10))
    (let ((overall-ratio
	   (/ num-plus-spanning-reads num-spanning-reads))
	  (variant-ratio
	   (/ num-plus-variant num-variant))
	  strand-bias-score)
      (setq strand-bias-score
	     (* 2 (abs (- variant-ratio overall-ratio)))))))

;; find strand counts and probability
(defgeneric determine-strand-counts-n-bias (var-cand))
(defmethod determine-strand-counts-n-bias ((var-cand variant-candidate))
  (let (binomial-prob
	(overall-pos 0)
	(total-reads 0)
	(prob-plus-strand 0.5d0))
    (with-slots (scores-n-q-t-seqs seq-dev-hash
                 strand-counts strand-probs
		 strand-bias-scores
		 num-plus-spanning-reads num-spanning-reads
		 overall-strand-prob)
	var-cand
      (values
       (setq
	strand-counts
	(mapcar #'(lambda (q-t-seq)
		    (let* ((strands
			    (mapcar #'strand
				    (mapcar #'bam (gethash q-t-seq seq-dev-hash))))
			   (pos-hits (count :+ strands))
			   (neg-hits (count :- strands)))
		      (incf overall-pos pos-hits)
		      (incf total-reads pos-hits)
		      (incf total-reads neg-hits)
		      (push (cumalative-binomial-probability
			     pos-hits
				  (+ pos-hits neg-hits)
				  prob-plus-strand ;; prob. of any read hitting pos strand
				  ) binomial-prob)
		      (list pos-hits
			    neg-hits)))
		(mapcar #'cdr scores-n-q-t-seqs)))
       (setq strand-probs (reverse binomial-prob))
       (setq strand-bias-scores
	     (mapcar #'(lambda (strand-count)
			 (when (and num-plus-spanning-reads num-spanning-reads)
			   (calc-strand-bias num-plus-spanning-reads num-spanning-reads (car strand-count) (+ (car strand-count) (second strand-count)))))
		     strand-counts))
       (setq overall-strand-prob
	     (cumalative-binomial-probability overall-pos total-reads prob-plus-strand))))))

;; Filter variants
(defun make-filtered-variant-count-statement (variants stream)
  (format stream "Total number of variants: ~a~%Total that are filtered: ~a~%" (length variants) (count t variants :key #'overall-filtered?)))

(defgeneric filter-variants (var-cand filter-hash))
(defmethod filter-variants ((var-cand variant-candidate) filter-hash)
  (with-slots (seq-dev-hash scores-n-q-t-seqs variant-freqs-by-score
			    strand-bias-scores
			    strand-counts strand-probs
			    rightmost-positions
			    overall-filtered?)
      var-cand
    (let (new-scores new-var-freqs new-strand-counts new-strand-probs
		     new-strand-bias-scores
		     new-rightmost-positions)
      ;;(setq variants-filtered?
	    (mapcar
	     #'(lambda (score-n-q-t-seq var-freq 
			strand-bias strand-count strand-prob rightmost-pos)
		 (let ((num-reads (length (gethash (cdr score-n-q-t-seq) seq-dev-hash)))
		       (score (car score-n-q-t-seq))
		       dofilter?)
		   (when (or (< score (gethash :score-threshold filter-hash))
			     (< num-reads (gethash :min-num-reads filter-hash))
			     (when var-freq
			       (or (eql var-freq :infinity)
				   (< var-freq (gethash :min-variant-freq filter-hash))))
			     (when strand-bias
			       (> strand-bias (gethash :max-strand-bias
						       filter-hash)))
			     (< strand-prob (gethash :strand-prob filter-hash)))
		     (setq dofilter? t))
		   (unless dofilter?
		     (push score-n-q-t-seq new-scores)
		     (push var-freq new-var-freqs)
		     (push strand-bias new-strand-bias-scores)
		     (push strand-count new-strand-counts)
		     (push strand-prob new-strand-probs)
		     (push rightmost-pos new-rightmost-positions)
		     )
		   dofilter?))
	     scores-n-q-t-seqs
	     (or variant-freqs-by-score (make-list (length scores-n-q-t-seqs)))
	     strand-bias-scores strand-counts strand-probs rightmost-positions) ;;)

      (setq scores-n-q-t-seqs (reverse new-scores))
      (setq variant-freqs-by-score (reverse new-var-freqs))
      (setq strand-bias-scores (reverse new-strand-bias-scores))
      (setq strand-counts (reverse new-strand-counts))
      (setq strand-probs (reverse new-strand-probs))
      (setq rightmost-positions (reverse new-rightmost-positions))

      (setq overall-filtered? (if scores-n-q-t-seqs nil t))    
      )))

(defun print-filter-settings (filter-hash stream)
  (let ((filters '(:min-mapq :score-threshold :min-num-reads :min-variant-freq 
		   :max-strand-bias :strand-prob))
	(descriptions '("Min. mapq for a read to be considered for variant calling."
			"Min score required for each variant found at a location."
			"Min number of reads for each variant."
			"Min read frequency for each variant."
			"Max strand bias for each variant."
			"Min binomial strand probability for each variant."))
	filters-n-values)
    (setq filters-n-values
	  (mapcar
	   #'(lambda (filter description)
	       (let ((filter-value (gethash filter filter-hash)))
		 (when filter-value
		   (format nil "##FILTER=<ID=~a,Description=\"~a=~a. ~a\">" filter filter filter-value description)
		   )))
	   filters descriptions))
    (setq filters-n-values (remove nil filters-n-values))
    (format stream "~{~a~^~%~}~%" filters-n-values)))

;; Memory cleanup after printing of variant is done
(defgeneric clean-uniq-seq-n-seq-hash (var-cand))
(defmethod clean-uniq-seq-n-seq-hash ((var-cand variant-candidate))
  (with-slots (seq-deviations seq-dev-hash scores-n-q-t-seqs rightmost-positions)
      var-cand
    (setq seq-dev-hash nil)
    (setq scores-n-q-t-seqs nil)
    (setq rightmost-positions nil)
    ))
    
(defun determine-map-qv-rms-n-0-counts (map-qvs top-n)
  (let ((mapq0-counts 0)
	relevant-map-qvs)
    (dolist (map-qvs-for-var (subseq map-qvs 0 top-n))
      (dolist (map-qv map-qvs-for-var)
	(when (eql 0 map-qv)
	  (incf mapq0-counts))
	(push map-qv relevant-map-qvs)))
    (values (root-mean-square relevant-map-qvs)
	    mapq0-counts)))
