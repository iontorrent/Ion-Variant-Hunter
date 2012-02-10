
;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

;; alignment-streamer.lisp
;; Takes bam records & alignments, and streams them through deviation and
;; variant candidate determination to end with output of the variant.

;; Functions that extract out relavent info from a bam record so that a flow space
;; alignment can be performed.

;; Extracts out part of flow intensities that was found to be
;; aligned to the genome
(defun get-number-flows-to-last-key (key-seq &optional (flow-order *flow-order*)
				     (num-flows 0))
  (if (> (length key-seq) 0)
      (let ((pos-key-base (position (char key-seq 0) flow-order
				    :start num-flows)))
	;; position is nil if not enough seq to cover key, i.e.
	;; (GET-NUMBER-FLOWS-TO-LAST-KEY "G" "GTCTGA" 6)
	;; Happens if called from calc flow numbers, where "flow-order" is
	;; adjusted during the alignment.
	;; i.e 15H179S8M in read COT5K:861:243
	(when pos-key-base
	  (get-number-flows-to-last-key (subseq key-seq 1) flow-order
					(1+ pos-key-base))))
      num-flows))

(defgeneric trim-n-orient-read-seq (bam))
(defmethod trim-n-orient-read-seq ((bam bam-record))
  (with-slots (cigar strand read-name seq)
      bam
    (unless (eql 0 (get-begin-soft-clipping cigar strand))
      (format *error-output* "Warning, soft clipping found at begining of read, and won't be aligned well, for read ~a.~%" read-name))
    (let ((end-trim
	   (get-begin-soft-clipping cigar strand t)))
      (subseq (if (eql strand :+)
		  seq
		  (reverse-complement seq)) 0 (- (length seq) end-trim)))))

(defgeneric get-relevant-flow-intensities (bam &optional key-seq flow-order))
(defmethod get-relevant-flow-intensities ((bam bam-record)
					  &optional
					  (key-seq "TCAG")
					  (flow-order *flow-order*))
  (with-slots (read-name flow-intensities strand seq cigar)
      bam
    ;;(format t "Working with read ~a~%" read-name)
    (let ((read-seq (trim-n-orient-read-seq bam))
	  (read-pos 0)
	  (cur-flow (get-number-flows-to-last-key key-seq flow-order))
	  relevant-intensities
	  (more-seq? t)
	  intensity
	  ;; take part of intensities that includes only last base of key
	  intensities)
      (setq intensities (subseq flow-intensities cur-flow))
      (decf cur-flow) ;; catch last base of key

      ;; first intensity is shared with key, so subtract that out.
      (when (< (car intensities) 50)
	(format *error-output* "Warning, missing flow intensity for last base of key in read ~a."
		read-name))
      (setf (car intensities) (- (car intensities) 100))
      (when (< (car intensities) 0) (setf (car intensities) 0))

      ;; loop seq and intensities
      (while more-seq?
	(setq intensity (pop intensities))
	(multiple-value-bind (num-bases residual)
	    (round-half-down intensity 100)
	  (let ((flow-order-nuc
		 (char flow-order (mod cur-flow
				       (length flow-order))))
		)
	    ;; Last flow could be present but number of bases (signal)
	    ;; could be trimmed.  So reset num-bases such that it does not
	    ;; go off the read sequence.
	    (when (>= (+ read-pos num-bases) (length read-seq))
	      (setq num-bases (- (length read-seq) read-pos 1))
	      (setq more-seq? nil)
	      ;;(when (> read-pos 60)
		;;(format t "Resetting num-bases to ~a.~%" num-bases))
	      )

	    #|;; debug statement
	    (when (or (> read-pos 60) (< read-pos 5))
	      (format t "flow-order-nuc:~a num-bases:~a read-pos:~a/~a cur-flow:~a/~a inten-size:~a~%"
		      flow-order-nuc
		      num-bases
		      read-pos (length read-seq)
		      cur-flow (length flow-intensities) (length intensities))) ;;|#
	     #| End of Debug statement |#

	    (dotimes (x num-bases)
	    #|  Debug statement
	      (when (> read-pos 60) ;;(- (length read-seq) 15))
		(format t "base-call:~a flow-order-nuc:~a inten:~a lgth:~a x:~a pos:~a pos-sum:~a~%"
			(when (<  (+ x read-pos) (length read-seq))
			  (char read-seq (+ x read-pos)))
			flow-order-nuc
			intensity
			(length read-seq) x read-pos (+ x read-pos))) ;; |#
	      #| End of Debug statement |#
	      (unless (eql (char read-seq (+ x read-pos))
			   flow-order-nuc)
		(error "Bam record ~a, Flow order nuc (~a) @ flow ~a w/ intensity ~a does not match base in seq (~a) @ position ~a" (read-name bam)
		       flow-order-nuc cur-flow intensity
		       (char read-seq (+ x read-pos)) (+ x read-pos))))
	    (incf read-pos num-bases)

	    #|  Debug statement 
	    (when (> read-pos 100)
	      (format t "~a < ~a ?~%" read-pos (length read-seq))
	      (format t "~a == 50 ?~%" residual)
	      (format t "~a == ~a ?~%" (when (< read-pos (length read-seq))
					 (char read-seq read-pos)) flow-order-nuc)) ;; |#
	    #| End of Debug statement |#

	    ;; check if x50 intensity base is actully there
	    (when (and (< read-pos (length read-seq))
		       (eql residual 50)
		       (eql (char read-seq read-pos)
			    flow-order-nuc))
	      ;; It actually does have extra base (as it was rounded up)
	      (incf read-pos)
	      (incf num-bases))

	    ;;loop increment & termination check
	    (incf cur-flow) ;; goto next flow in loop
	    (push intensity relevant-intensities)

	    (unless  (< read-pos (length read-seq))
	      (setq more-seq? nil))
	    (unless (< cur-flow (- (length flow-intensities) 1))
	      (setq more-seq? nil))

	    )
	  )) ;; end of while loop

      (reverse relevant-intensities))))

;; Determine which flow each position of the alignment is from.
;; Value of nil means an insertion to represent the reference was done.
(defgeneric calc-flow-numbers-lisp-aligner (aligner key-seq))
(defmethod calc-flow-numbers-lisp-aligner ((aligner flow-space-aligner) key-seq)
  ;;(setq *debug-aligner* aligner)
  (with-slots (flow-order symbols
			  flow-numbers)
      aligner
    ;; flow-order of aligner object is the flow bases in the full alignment (including insertions)
    (let ((cur-flow-num (get-number-flows-to-last-key key-seq flow-order)))
      (unless cur-flow-num  ;; See note in get-number-flows to see when this is nil.
	;; Causes flow-numbers to not be filled.  TODO, probably ok.
	(return-from calc-flow-numbers-lisp-aligner nil))
      (decf cur-flow-num)
      (setq flow-numbers
	    (mapcar #'(lambda (flow-base symbol)
			(unless (eql *align-del* symbol)
			  flow-base
			  #|  ;; debugging check of flow order
			  (unless (eql flow-base
				       (char *flow-order* (mod cur-flow-num (length *flow-order*))))
			    (format t "Warning ~a base in alignment doesn't match ~a base @ pos ~a in flow order.~%" flow-base (char *flow-order* (mod cur-flow-num (length *flow-order*))) (mod cur-flow-num (length *flow-order*))))  ;;|#
			  (incf cur-flow-num)))
		    flow-order symbols)))))

;; The following helper functions are for calc-flow-numbers-java-aligner
(defgeneric make-align-flow-array (aligner))
(defmethod make-align-flow-array ((aligner flow-space-aligner))
  (with-slots (flow-order symbols) ;; actually this is the align-flows
      aligner
    (let ((align-flow-array (make-array (length flow-order) :fill-pointer 0)))
      (mapc #'(lambda (flow-base align-symbol)
		(unless (eql align-symbol *align-del*)  ;;bases only in ref won't be in flow order
		  (vector-push flow-base align-flow-array)))
	    flow-order symbols)
      align-flow-array)))

(defun get-start-flow-number (align-flow-array flow-order-seq key-clip-offset)
  ;; Exact align flow-order from java aligner
  ;; "AGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTAC"
  ;; with flow order sequence
  ;; "TACGTACGTCTGAGCATCGATCGATGTACAGC"

  (let (;;(align-flows (make-array (length flow-order) :inital-contents flow-order))
	(flow-start-num (- key-clip-offset 1))  ;; 0 based & consider last base of key
	(req-num-matches (min (length flow-order-seq)
			      (length align-flow-array)))
	(num-matches 0))
    (while (< num-matches req-num-matches)
      (setq num-matches 0)  ;; start over again
      (let (;;(align-flow-num 0)
	    (i 0))
	;; This loop sees how many matches there are
	(while (and (< i (length align-flow-array))
		    (eql
		     (char flow-order-seq (mod (+ flow-start-num i) (length flow-order-seq)))
		     (aref align-flow-array i)))
	  (incf i)
	  (incf num-matches))
	(incf flow-start-num) ;; increase the flow start number
	(when (> flow-start-num (+ key-clip-offset (* 2 (length flow-order-seq))))
	  (return-from get-start-flow-number :cannot-align)
	  )))
    (decf flow-start-num)))

(defgeneric get-first-unclipped-flow (bam num-flows-key))
(defmethod get-first-unclipped-flow ((bam bam-record) num-flows-key)
  (with-slots (flow-intensities cigar strand)
      bam
    (let ((flow-array
	   (make-array (length flow-intensities) :initial-contents flow-intensities))
	  (clip-length (get-read-clipping cigar strand))
	  (cur-flow num-flows-key))
      (while (and (< cur-flow (length flow-intensities))
		  (> clip-length 0))
	(when (eql cur-flow num-flows-key)
	  (incf clip-length)) ;; add one for last base of key [which will be removed by decf right below
	(decf clip-length (round-half-down (/ (aref flow-array cur-flow) 100)))
	(incf cur-flow)
	;;(format t "~a ~a~%" clip-length cur-flow)
	)
      (when (< clip-length 0) ;; clipped in middle of hp
	(decf cur-flow))  ;;go back to flow that was clipped
      cur-flow)))

(defgeneric calc-flow-numbers-java-aligner (aligner key-seq flow-order-seq &optional num-flows-key))
(defmethod calc-flow-numbers-java-aligner ((aligner flow-space-aligner) key-seq flow-order-seq &optional num-flows-key)

  ;; purely for effeciency, need only calculate this once during execution
  (setq num-flows-key (or num-flows-key
			  (get-number-flows-to-last-key key-seq flow-order-seq)))

  (with-slots (flow-order symbols
			  flow-numbers)
      aligner
    ;; flow-order of aligner object is the flow bases in the full alignment (including insertions)
    ;; flow-order-seq is the flow order of the ion system that repeats over and over again

      (let ((start-flow-number (get-start-flow-number (make-align-flow-array aligner)
						      flow-order-seq
						      (get-first-unclipped-flow (bam aligner) num-flows-key))))
	(when (keywordp start-flow-number)
	  (format *error-output* "WARNING, ~a happened when attempting to register positions of alignment flow to flow order bases.  Alignment flows: ~{~a~^~}, Flow order: ~{~a~^~}~%" start-flow-number flow-order flow-order-seq))
	(decf start-flow-number)
	(setq flow-numbers
	      (mapcar #'(lambda (flow-base symbol)
			  ;; deleted bases won't have a position in the read
			  (unless (eql *align-del* symbol)
			    flow-base
			    (incf start-flow-number)))
		      flow-order symbols)))))

;; From bam record, make a flow space alignment record
;; Used only for LISP flow aligner method only
(defgeneric do-alignment (bam  &key query-flow-intensities
			       target-base-seq
			       flow-order
			       key-seq
			       phase-penalty
			       end-local?))
(defmethod do-alignment ((bam bam-record) 
			 &key query-flow-intensities
			 target-base-seq
			 (flow-order *flow-order*)
			 (key-seq "TCAG")
			 (phase-penalty 1)
			 (end-local? t))
  (with-slots (strand ref-seq reverse-complement-ref
		      flow-intensities is-mapped?)
      bam
    (unless is-mapped?
      (return-from do-alignment nil)) ;; can't realign if wasn't aligned to begin with.
    (let ((target-base-seq
	   (or target-base-seq
	       (if (eql strand :+)
		   ref-seq
		   reverse-complement-ref)))
	  (query-flow-seq
	   (or query-flow-intensities
	       (mapcar #'(lambda (x)
			   (/ x 100))
		       (get-relevant-flow-intensities bam)))))
      (setq query-flow-seq
	    (let ((counter -1))
	      (mapcar #'(lambda (counts)
			  (incf counter)
			  (cons (char flow-order (mod counter (length flow-order)))
				counts))
		      query-flow-seq)))
      (let ((aligner (make-flow-space-aligner query-flow-seq target-base-seq
					      flow-order key-seq bam)))
	(init-gap-penalties aligner phase-penalty)
	(calc-score-matrix aligner phase-penalty)
	(let ((best-score (calc-best-score aligner end-local?)))
	  (trace-path-back aligner best-score)
	  (calc-flow-numbers-lisp-aligner aligner key-seq)
	  (values aligner best-score))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Alignment streaming
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; alignment-streamer class
(defclass alignment-streamer ()
  ((bam-file :accessor bam-file :initarg :bam-file)
   (sam-file :accessor sam-file :initarg :sam-file)

   ;; stream from a file or a std out running process
   (cur-in-stream :accessor cur-in-stream)

   ;; :lisp is internal method, :java is Nil's method for doing flow space alignment
   (aligner-method :accessor aligner-method :initarg :aligner-method)

   ;; Which process (i.e. making bam or making deviation) the streamer is at currently
   (cur-process :accessor cur-process :initarg :initialization)

   ;; current bam and align objects
   (cur-bam :accessor cur-bam :initform nil)
   (cur-align :accessor cur-align)
   (cur-deviation :accessor cur-deviation :initform nil)
   (last-bam :accessor last-bam :initform nil)

   ;; Cache so that they can be properly resorted from
   ;; bam order to sequence deviation location order.
   (seq-deviation-cache :accessor seq-deviation-cache :initarg :seq-deviation-cache)
   (temp-unsorted-cache :accessor temp-unsorted-cache :initform nil)
   (do-variant-calling? :accessor do-variant-calling? :initform nil)

   ;;output streams for variants and individual deviants
   (deviant-strm :accessor deviant-strm)
   (variant-strm :accessor variant-strm)

   ;; settings (which could be auto set)
   (settings-hash :accessor settings-hash)
   ))

;; default values and settings
(defun set-value-or-default (settings-hash setting value default-value)
  (setf (gethash setting settings-hash)
	(or value default-value)))

(defun set-filters (settings-hash &key
		    min-mapq
		    score-threshold
		    min-num-reads
		    min-variant-freq
		    max-strand-bias
		    strand-prob
		    max-intensity-per-base
		    )
  (set-value-or-default settings-hash :min-mapq min-mapq min-mapq)
  (set-value-or-default settings-hash :score-threshold  score-threshold  score-threshold)
  (set-value-or-default settings-hash :min-num-reads    min-num-reads    min-num-reads)
  (set-value-or-default settings-hash :min-variant-freq min-variant-freq min-variant-freq)
  (set-value-or-default settings-hash :max-strand-bias  max-strand-bias max-strand-bias)
  (set-value-or-default settings-hash :strand-prob      strand-prob      strand-prob)
  (set-value-or-default settings-hash :max-intensity-per-base-hash
			(make-base-max-intensity-hash max-intensity-per-base) nil)
  )

(defun set-settings-hash-for-java (settings-hash &key
				   java-bin
				   java-bin-args
				   fs-align-jar
				   fs-align-jar-args
				   fs-align-num-threads
				   fs-align-range
				   bam-file
				   sam-file
				   )
  ;; java location
  (set-value-or-default settings-hash :java-bin
			java-bin
			"/usr/bin/java")

  (set-value-or-default settings-hash :java-bin-args
			(argument-underscore-to-dash java-bin-args) nil)

  ;; jar location (defaults to same location as this binary's)
  (set-value-or-default settings-hash :fs-align-jar
			fs-align-jar
			(full-path-at-exec-directory "SamToFlowSpace.jar"))

  ;; set the number of threads
  (set-value-or-default settings-hash :fs-align-num-threads
			fs-align-num-threads 4)
  (set-value-or-default settings-hash :fs-align-range
			fs-align-range nil)

  ;; make complete list of arguments for the java call
  (let ((non-settable-arg-list
	 (list
	  (format nil "REFERENCE=~a" (ref-file (gethash :ref-reader settings-hash)))
	  (format nil "INPUT=~a" (or bam-file sam-file))
	  "QUIET_STDERR=true"
	  ))
	settings-args
	;; default settings
	(num-threads (gethash :fs-align-num-threads settings-hash)))
    (when fs-align-jar-args
      (unless (listp fs-align-jar-args)
	(setq fs-align-jar-args (list fs-align-jar-args))))
    (setq settings-args
	  (cons (format nil "NUM_THREADS=~a" num-threads)
		fs-align-jar-args))
    (when fs-align-range
      (setq fs-align-range (format nil "RANGE=~a" fs-align-range))
      (setq settings-args
	    (cons fs-align-range settings-args)))

    (setf (gethash :fs-align-jar-args settings-hash)
	  (append non-settable-arg-list settings-args))
  )
  settings-hash)

(defun set-settings-hash-for-python (settings-hash &key
				     python-bin
				     python-py)
  ;; python binary location
  (set-value-or-default settings-hash :python-bin
			python-bin "/usr/bin/python")

  ;; python py script location (defaults to same location as this binary's)
  (set-value-or-default settings-hash :python-py
			python-py
			(full-path-at-exec-directory "samRegionOverlap.py")))

(defun decide-align-method (aligner-method sam-file bam-file)
  (case aligner-method
    (:auto (if bam-file :java (if sam-file :lisp
				  (error "Must have either sam or bam file.~%"))))
    (:lisp (unless sam-file
	     (error  "Cannot use lisp flow space aligner without a sam file.~%"))
	   :lisp)
    (:java :java)))

(defun make-alignment-streamer (&key
				sam-file
				bam-file
				(aligner-method :auto)
				reference-file
				align-file
				key-seq
				flow-order

				max-deviation-cache-size
				vh-num-threads
				;; filters
				min-mapq
				score-threshold min-num-reads min-variant-freq
				max-strand-bias strand-prob
				max-intensity-per-base

				java-bin
				java-bin-args
				fs-align-jar
				fs-align-jar-args
				fs-align-num-threads
				fs-align-range

				python-bin
				python-py

				validate-alignments
				)
  (let ((streamer-obj (make-instance 'alignment-streamer
				     :bam-file bam-file
				     :sam-file sam-file
				     :aligner-method aligner-method
				     :seq-deviation-cache nil))
	)
    (with-slots (settings-hash aligner-method)
	streamer-obj      
      (setq aligner-method (decide-align-method aligner-method sam-file bam-file))
      (format t "Using ~a flow space alignment method.~%" aligner-method)

      ;; Make and set settings hash
      (setq settings-hash (make-hash-table))
      (setf (gethash :ref-reader settings-hash)
	    (make-reference-reader reference-file))
      (setf (gethash :align-file settings-hash) align-file)
      (setf (gethash :key-seq settings-hash) key-seq)
      (setf (gethash :flow-order settings-hash) flow-order)
      (setf (gethash :num-flows-key settings-hash) (get-number-flows-to-last-key key-seq flow-order))
      (setf (gethash :validate-alignments settings-hash) validate-alignments)
      (setf (gethash :vh-num-threads settings-hash) vh-num-threads)
      (setf (gethash :max-deviation-cache-size settings-hash) max-deviation-cache-size) ;;150000
      (set-filters settings-hash
		   :min-mapq min-mapq
		   :score-threshold score-threshold 
		   :min-num-reads  min-num-reads
		   :min-variant-freq min-variant-freq
		   :max-strand-bias max-strand-bias
		   :strand-prob strand-prob
		   :max-intensity-per-base max-intensity-per-base
		   )

      (set-settings-hash-for-java settings-hash
				  :java-bin java-bin
				  :java-bin-args java-bin-args
				  :fs-align-jar fs-align-jar
				  :fs-align-jar-args fs-align-jar-args
				  :fs-align-num-threads fs-align-num-threads
				  :fs-align-range fs-align-range
				  :bam-file bam-file
				  :sam-file sam-file
				  )
      (set-settings-hash-for-python settings-hash
				    :python-bin python-bin
				    :python-py python-py
      ))
    streamer-obj))

(defgeneric close-streamer (streamer))
(defmethod close-streamer ((streamer alignment-streamer))
  (with-slots (settings-hash)
      streamer

    ;; Close the flow space alignment output stream
    (when (gethash :align-stream settings-hash)
      (close (gethash :align-stream settings-hash)))))

(defgeneric print-generic-header (streamer stream file-type))
(defmethod print-generic-header ((streamer alignment-streamer) stream file-type)
  (let (header-parts)
    (unless (eql :variant file-type)
      (push (format nil "file-type=~a" file-type) header-parts))
    (when (eql :variant file-type)
      (push (format nil "fileformat=VCFv4.1") header-parts))
    (push (format nil "fileDate=~a" (format-current-date :no-dividers)) header-parts)
    (push (format nil "fileTime=~a" (format-current-time :no-dividers)) header-parts)
    (push (format nil "source=~a" "IonVariantHunterV0.0.2") header-parts)

    (with-slots (bam-file sam-file settings-hash)
	streamer
      (when (gethash :ref-reader settings-hash)
	(let ((ref-file (ref-file (gethash :ref-reader settings-hash))))
	  (push (format nil "reference=file://~a" ref-file) header-parts)))
      (when bam-file
	(push (format nil "bam=file://~a" bam-file) header-parts))
      (when sam-file
	(push (format nil "sam-file://~a" sam-file) header-parts))
      )
    (when (boundp '*build-date*)
      (flet ((sym-val (str) (symbol-value (intern str))))
	(push (remove #\- (format nil "buildDate=~a" (sym-val "*BUILD-DATE*"))) header-parts)
	(push (remove #\: (format nil "buildTime=~a" (sym-val "*BUILD-TIME*"))) header-parts)))
    
    (dolist (header-part (reverse header-parts))
      (format stream "##~a~%" header-part))))

(defun print-vcf-header-column-names (stream)
  (let ((col-names '(CHROM POS ID REF ALT QUAL FILTER INFO)))
    (format stream (format nil "#~~{~~a~~^~a~~}~%" #\Tab)
	    col-names)))


;; These files aren't needed but could be outputted for debugging purposes
(defgeneric open-n-print-headers-intermediate-files (streamer))
(defmethod open-n-print-headers-intermediate-files ((streamer alignment-streamer))
  (with-slots (settings-hash)
      streamer
    (let ((align-file (gethash :align-file settings-hash)))
      (when align-file ;; open a stream for flow space alignment printouts
	(let ((align-stream
	       (open align-file :direction :output :if-exists :rename)))
	  (print-generic-header streamer align-stream :flow-space-alignments)
	  (setf (gethash :align-stream settings-hash) align-stream))))))

(defgeneric check-cur-bam-align (streamer))
(defgeneric output-to-intermediate-files (streamer))
(defgeneric make-deviation-from-cur-bam-align (streamer))

(defgeneric do-variant-calling (streamer))
(defgeneric process-variant-candidates (streamer variant-candidates))
(defgeneric reset-cache-after-variant-calling (streamer))

(defgeneric stream-bam (streamer method-on-align))
(defgeneric stream-ext-java (streamer method-on-align))
(defgeneric stream-alignments (streamer))

(defmethod stream-alignments ((streamer alignment-streamer))
  (flet ((method-on-align ()
	   ;; get the current bam record
	   (check-cur-bam-align streamer)

	   ;; checks for and removes mismatches based on very large intensity values
	   (let ((intensity-ht (gethash :max-intensity-per-base-hash (settings-hash streamer))))
	     (alter-alignment-for-large-intensities (cur-align streamer) intensity-ht))

	   ;; output flow space alignments
	   (output-to-intermediate-files streamer)

	   ;; Find sequence deviations on a per read basis
	   (make-deviation-from-cur-bam-align streamer)

	   ;; After seq-deviation-cache reaches a certain point
	   ;; can do variant calling
	   (when (do-variant-calling? streamer)
	     (do-variant-calling streamer))

	   ;; Set last BAM file
	   (when (cur-bam streamer)
	     (setf (last-bam streamer) (cur-bam streamer)))))
    (let (num-records)
      (with-slots (aligner-method seq-deviation-cache)
	  streamer
	(setq num-records
	      (case aligner-method
		(:java (stream-ext-java streamer #'method-on-align))
		(:lisp (stream-bam streamer #'method-on-align))
		(t (error "Method ~a is not available as a streaming process.~%"
			  aligner-method))))
	;; do variant calling on whatever is left in the cache
	(when (seq-deviation-cache streamer)
	  (setq seq-deviation-cache (sort-deviations seq-deviation-cache))
	  (do-variant-calling streamer)
	  )
	(format t "Processed ~a records.~%" num-records)
	(setq seq-deviation-cache (reverse seq-deviation-cache))))))

(defun progress-message (counter cur-bam start-time &optional force-print)
  (when (or force-print
	    (eql 0 (mod counter 2500)))
    (let* ((cur-time (get-universal-time))
	   (elapsed-time (- cur-time start-time)))
      (format t "Processed ~a reads in ~a sec. (~,2f reads/sec)."
	      counter elapsed-time
	      (unless (eql 0 elapsed-time)
		(/ counter elapsed-time))))
    (when cur-bam
      (format t " Last read: ~a" (get-basic-attrib-string cur-bam)))
    (format t "~%")))

(defgeneric abort-message (streamer &optional output-line in-read-counter?))
(defmethod abort-message ((streamer alignment-streamer) &optional output-line in-read-counter?)
  (with-slots (cur-process cur-bam)
      streamer
    ;; Caught some sort of error or interrupt
    (unless (eql output-line :eof)
      (format *error-output* "WARNING: An interrupt or an error occurred when performing ~a.~%"
	      cur-process)
      (when output-line
	(format *error-output* "Output line: ~a~%" output-line))
      (when (and cur-bam (not in-read-counter?))
	(format *error-output*
		"This was at or near the processing of the record with a read name of ~a.~%"
		(read-name cur-bam))))))

#|
(defvar *debug-bam-records* nil)
(defvar *debug-aligns* nil)
(setq *debug-bam-records* nil)
(setq *debug-aligns* nil)
;; |#

;; Used only for LISP flow aligner (not typical)
(defmethod stream-bam ((streamer alignment-streamer) method-on-align)
  (let ((num-header-lines (determine-number-header-lines 
			   (sam-file streamer) :header-char #\@))
	;; Counters
	(counter 0)
	(num-bam-records 0)
	(num-aligned-bam-records 0)
	(num-flow-space-aligns 0)
	;; Start time for progress meter
	(start-time (get-universal-time))
	)
    (with-slots (sam-file cur-bam cur-align cur-deviation cur-process settings-hash)
	streamer
      (with-delimited-file (samLine sam-file num-header-lines)
	(unwind-protect
	     (progn
	       (setq cur-process :start-line)
	       (incf counter)

	       ;; make bam record
	       (setq cur-process :parsing-bam)
	       (setq cur-bam (make-bam-record samLine (gethash :ref-reader settings-hash)))
	       (progress-message counter cur-bam start-time)
	       (when cur-bam
		 (incf num-bam-records)
		 (when (is-mapped? cur-bam)
		   (incf num-aligned-bam-records)))
	       ;;(push cur-bam *debug-bam-records*)

	       ;; make align record
	       (setq cur-process :doing-alignment)
	       (setq cur-align (do-alignment cur-bam
				 :flow-order (gethash :flow-order settings-hash)
				 :key-seq (gethash :key-seq settings-hash) ))
	       ;;(push cur-align *debug-aligns*)

	       (when cur-align
		 (setf (dp-matrix cur-align) nil)  ;; takes up a large amount of memory & no need to keep after alignment done
		 ;;(sb-ext:gc :full t)  ;; may slow things down too much
		 (incf num-flow-space-aligns))
	       (setq cur-process :alignment-done)

	       ;; Call method to make sequence deviation from alignment
	       (funcall method-on-align)

	       (setq cur-process :bam-line-done)
	       )
	  (abort-message streamer)
	  )) ;; end of loop

      ;; Final progress message and printout of number of records
      (progress-message counter cur-bam start-time t)
      (format t "Number of bam records: ~a~%" num-bam-records)
      (format t "Number of aligned bam records: ~a~%" num-aligned-bam-records)
      (format t "Number of flow space alignments found: ~a~%" num-flow-space-aligns)
      (unless (eql num-aligned-bam-records num-flow-space-aligns)
	(format *error-output* "WARNING, some aligned bam records did not produce a flow space alignment.~%"))

      counter ;; return
      )))

;; Sets individual parts coming from the external java process
(defun parse-line-to-get-intensities (output-line)
  (mapcar #'(lambda (x)
	      (if (find (char x 0) "-+")
		  0
		  (/ (parse-integer x) 100)))
	  (parse-string output-line #\,)))

(defgeneric set-bam-align-record (streamer output-line line-counter))
(defmethod set-bam-align-record ((streamer alignment-streamer) output-line line-counter)
  (with-slots (settings-hash cur-bam cur-align)
      streamer
    (let ((record-line-num (mod line-counter 6)))
      (case record-line-num
	;; sam record
	(0 (setq cur-bam (make-bam-record (parse-string output-line #\Tab)
					  (gethash :ref-reader settings-hash))))
	;; flow intensities
	(1 (setq cur-align (make-instance 'flow-space-aligner :bam cur-bam))
	   (setf (query-counts cur-align)
		 (parse-line-to-get-intensities output-line)))
	;; align symbols
	(2 (setf (symbols cur-align)
		 (mapcar #'(lambda (sym-string)
			     (let ((symbol
				    (char sym-string 0)))
			       (if (eql #\Space symbol)
				   *align-mismatch*
				   symbol)))
			 (parse-string output-line #\,)))
	   ;; Trimming upstream may cause some mismatches
	   ;; on the edge of the alignments, but they are not
	   ;; going to have evidence for 
	   (when (eql (car (symbols cur-align)) *align-mismatch*)
	     (setf (car (symbols cur-align)) *align-match*))
	   (when (eql (car (last (symbols cur-align))) *align-mismatch*)
	     (setf (car (last (symbols cur-align))) *align-match*))
	   )
	;; target "intensities"
	(3 (setf (target-counts cur-align)
		 (parse-line-to-get-intensities output-line)))
	;; flow bases for entire alignment (flow-order here isn't quite the right name).
	(4 (setf (flow-order cur-align)
		 (mapcar #'(lambda (x)
			     (char x 0))
			 (parse-string output-line #\,)))
	   (calc-flow-numbers-java-aligner (cur-align streamer)
					   (gethash :key-seq settings-hash)
					   (gethash :flow-order settings-hash)
					   (gethash :num-flows-key settings-hash)))
	(5 (let* ((parts (parse-string output-line #\,)))
	     ;; Note, original java start pos is 0-based and includes that position
	     ;;    and the original java end pos is 0-base but excludes that position.
	     ;; Here, both the start and end pos are 1-base, inclusive.
	     (with-slots (target-start-pos target-end-pos)
		 cur-align
	       ;; 1 based
	       (setq target-start-pos (parse-integer (car parts)))
	       ;; 1 based, inclusive
	       (setq target-end-pos (parse-integer (second parts)))

	       ;; shift if insertion at begining or end
	       (when (eql *align-ins* (car (symbols cur-align)))
		 (if (eql :+ (strand cur-bam))
		     (incf target-start-pos)
		     ))
	       (when (eql *align-ins* (car (last (symbols cur-align))))
		 (if (eql :- (strand cur-bam))
		     (incf target-start-pos)
		     ))
	       )
	     (setf (align-score cur-bam) (parse-integer (third parts))))
	   :finished)
	;; Other values aren't possible b/c of the modular done above.
	))))

(defun run-java-program (settings-hash &optional (error-output "/dev/null"))
  error-output
  (let (;;(java-command (format nil "'/usr/bin/java ~{~a~^ ~}~a'"
	;;		      java-args
	;;		      (if error-output
	;;			  "" ;;(format nil " 2> ~a" error-output)
	;;			  "")))
	(java-bin (gethash :java-bin settings-hash))
	(java-bin-args (gethash :java-bin-args settings-hash))
	(java-jar (gethash :fs-align-jar settings-hash))
	(jar-args (gethash :fs-align-jar-args settings-hash))
	java-args
	)
    (setq java-args (append (when java-bin-args (list java-bin-args)) (list "-jar" java-jar) jar-args))
    (format t "Executing java command:~%~a ~{~a~^ ~}~%"
	    java-bin java-args)
    
    (sb-ext:run-program ;;"/bin/bash"
			;;(list "-c" java-command)
			java-bin ;;"/usr/bin/java"
			java-args
			;;"/bin/cat"
			;;'("/data/antioch/projects/ion_variant_calling/mixed-NA12878-NA19099/dev-try3/test.out")
			;;'("/data/antioch/projects/ion_variant_calling/mixed-NA12878-NA19099/concreteVariantCaller-try2/ANG-189-r125148/ANG-189-r125148.concreteVariantCalls.sam-fsalign")
			:wait nil       ;; Don't wait for command to finish
			:output :stream ;; but rather stream it here.
			:error :stream
			)))

(defun check-process-run-status-n-print-stderr (process)
  (let ((process-stderr (process-error process))
	(process-exit (process-exit-code process))
	std-error-lines) ;;standard error output lines
    (labels ((add-error-line (error-line)
	       (push error-line std-error-lines)))
      (loop for error-line =
	   (read-line process-stderr nil :eof)
	 while (not (eql :eof error-line))
	 do
	   (when (string= "java.lang.Exception"
			  (subseq error-line 0 (min 19 (length error-line))))
	     (add-error-line "****************************~%")
	     (add-error-line "** ERROR by java process: **~%")
	     (add-error-line "****************************~%"))
	   (add-error-line error-line)
	   ))
    (setq std-error-lines (reverse std-error-lines))
    (when std-error-lines
      (format *error-output* "Note, there was stderr output from external process with PID of ~a:~%" (process-pid process))
      (format *error-output* "~{~a~^~%~}~%" std-error-lines)
      (format *error-output* "End of external process error output.~%"))
    (close process-stderr)
    (unless (eql 0 process-exit)
      (error "ERROR, external process with pid ~a exited abnormally with an exit code of ~a.~%"
	     (process-pid process)
	     (process-exit-code process))))
    )

(defmethod stream-ext-java ((streamer alignment-streamer) method-on-align)
  (with-slots (settings-hash cur-bam cur-align cur-process)
      streamer
    (let* ((process
	    (run-java-program settings-hash))
	   (java-stdout (sb-ext:process-output process))
	   (line-counter 0) 
	   (record-counter 0)
	   (start-time (get-universal-time))
	   cur-line
	   )
      (print-process-info process)
      (setq cur-process :starting-java-process)
      (unwind-protect
	   (loop for output-line = 
		(setq cur-process :reading-from-java-process
		      cur-line (read-line-of-process process java-stdout)) ;; this will set output-line too
	      while (and (setq cur-process :loop-done)
			 (not (eql :terminated output-line))
			 (not (eql :eof output-line)))
	      do
	      ;;(unwind-protect
		(let (status)
		  ;;(setq cur-line output-line)
		  (setq cur-process :processing-java-output)
		  (setq status 
			(set-bam-align-record streamer output-line line-counter))
		  
		  (incf line-counter)
		  (when (eql :finished status)
		    (progress-message record-counter cur-bam start-time)
		    (incf record-counter)
		    ;;Have bam and aligner objects filled, so now do the deviant detection
		    (setq cur-process :finding-deviations)
		    (funcall method-on-align)
		    )
		  (setq cur-process :bam-line-done)
		  ))
	;; Done with loop, check if it terminated early
	(sleep 1)
	(abort-message streamer cur-line)
	(unless (eql cur-process :bam-line-done)
	  (when (process-alive-p process)
	    (format *error-output* "Note, killing java process ~a.~%" process)
	    (sb-ext:process-kill process 15)))
	)
	   
      (print-process-info process)
      (progress-message record-counter cur-bam start-time t)
      (check-process-run-status-n-print-stderr process)
      (when (process-alive-p process)
	(format *error-output* "Note, java process still alive even after task completion. Terminating process ~a.~%" process)
	(sb-ext:process-kill process 15))
      record-counter
      )))

;;;; Calling python to get spanning reads
;;;; query bam files for counting spanning reads

(defun run-python-program (settings-hash bam-file)
  settings-hash
  (let ((python-bin (gethash :python-bin settings-hash)) ;;"/usr/bin/python")
	(python-py (gethash :python-py settings-hash)) ;;'("/data/antioch/projects/ion_variant_calling/mixed-NA12878-NA19099/samRegionOverlap/tests/samRegionOverlap.py"
	(python-args (list bam-file)) ;; code just takes bam as a single argument
	)
    ;;"/data/antioch/projects/ion_variant_calling/mixed-NA12878-NA19099/samRegionOverlap/tests/ANG-190-r125150.sorted.bam")))
    (format t "Executing python command:~%~a ~a ~{~a~^ ~}~%"
	    python-bin python-py python-args)
    (sb-ext:run-program python-bin
			(cons python-py python-args)
			:wait nil
			:output :stream
			:error :stream
			:input :stream)))

(defgeneric extract-reads-over-regions (streamer variant-candidates))
(defmethod extract-reads-over-regions ((streamer alignment-streamer) variant-candidates)
  (with-slots (settings-hash bam-file)
      streamer
  (let* ((process (run-python-program settings-hash bam-file))
	 (python-stdout (sb-ext:process-output process))
	 (python-stdin (sb-ext:process-input process))
	 spanning-read-names-sets
	 spanning-read-counts
	 spanning-plus-strand-read-counts
	 current-read-names
	 current-counter
	 current-plus-strand-counter
	 in-record?
	 cur-line
	 )
    ;; (print-process-info process)

    ;; give python all the positions
    (dolist (var-candidate variant-candidates)
      (with-slots (ref-name leftmost-pos max-rightmost)
	  var-candidate
	(let ((start-pos leftmost-pos)
	      end-pos ;; (or max-rightmost leftmost-pos)))
	      (length-contig (length-contig
			      (gethash ref-name (fai-entry-hash (gethash :ref-reader settings-hash))))))
	  (setq start-pos (max start-pos 1))             ;; Ensure start is 1 or higher
	  (setq end-pos (or max-rightmost start-pos))    ;; End position setting
	  (setq start-pos (min start-pos length-contig)) ;; Check if it goes off the chromosome
	  (setq end-pos (min end-pos length-contig))     ;; Same for end position

	  (format python-stdin "~a,~a,~a~%" ref-name start-pos end-pos))))

    (finish-output python-stdin)
    (close python-stdin)

    ;; Now read the output
      
      (unwind-protect

    (loop for output-line = (setq ;; cur-process :counting-reads-for-region
				  cur-line (read-line-of-process process python-stdout))
       while (and ;;(setq cur-process :loop-done)
		  (not (eql :terminated output-line))
		  (not (eql :eof output-line)))
       do
	 ;;(format t "OUT: ~a~%" output-line)
	 ;;(setq cur-process :counting-reads-for-region)
	 (cond
	   ;; start of the sam record list
	   ((string= "RECORDS-START"
		     (car (parse-string output-line #\:)))

	    (setq current-read-names nil)
	    (setq current-counter 0)
	    (setq current-plus-strand-counter 0)
	    (setq in-record? t))

	   ;; end of the sam record list
	   ((string= "RECORDS-END" output-line)
	    (when in-record?
	      (push current-read-names spanning-read-names-sets)
	      (push current-counter spanning-read-counts)
	      (push current-plus-strand-counter spanning-plus-strand-read-counts)
	      )
	    (unless in-record?
	      (format *error-output* "WARNING, detected RECORDS-END without a matching RECORDS-START.~%"))
	    (setq in-record? nil))

	   ;; inside record.  Just counting now.
	   (in-record?
	    (push (subseq output-line 2) current-read-names)
	    (when (eql (char output-line 0) #\+)
	      (incf current-plus-strand-counter))
	    (incf current-counter))

	   ;; ignored output
	   (t ;;ignore line
	    ))
	 ;;(setq cur-process :read-counter-line-done)
       ;;(clear-input python-stdout)
       ;;(setq process (run-python-program nil))
	 ) ;; end of reading python
    ;; check status

	(close python-stdout)
	;; Check to see if process was terminated early
	(sleep 1)
	(abort-message streamer cur-line t)
	(when (process-alive-p process)
	  (format *error-output* "Note, python process still alive even after task completion. Terminating process ~a.~%" process)
	  (sb-ext:process-kill process 15))
	) ;;end of unwind protect

    (setq spanning-read-names-sets
	  (reverse spanning-read-names-sets))
    (setq spanning-read-counts
	  (reverse spanning-read-counts))
    (setq spanning-plus-strand-read-counts
	  (reverse spanning-plus-strand-read-counts))
    (unless (eql (length spanning-read-counts)
		 (length variant-candidates))
      (format *error-output* "CRITICAL WARNING.  Number of records extracted from python script does not match the number of variants.  Total read counts are not to be trusted.")
      )
    ;; now annotate tags
    (when variant-candidates
      (while (and variant-candidates spanning-read-counts)
	(calc-spanning-reads-n-freqs (pop variant-candidates) (pop spanning-read-counts)
				     (pop  spanning-plus-strand-read-counts)
				     (pop spanning-read-names-sets)
				     ))
      )

    ;; check the status of the external python process and print out stderr
    ;; (print-process-info process)
    (check-process-run-status-n-print-stderr process)

    ;;won't be any return value
    )))
;;;; End of python extraction of reads

(defmethod check-cur-bam-align ((streamer alignment-streamer))
  (let ((good-alignment? t))
    (with-slots (cur-align cur-bam)
	streamer
      (unless (and cur-align cur-bam)
	(setq good-alignment? nil)
	;; print out reason reason why if there is some concern about those alignments
	(cond (cur-bam
	       (when (is-mapped? cur-bam) ;; unmapped read is normal, so no message needed.
		 ;; Something else caused no alignment to be found.
		 (format *error-output* "Warning, alignment failure for bam ~a.~%" cur-bam)))
	      (cur-align
	       (format *error-output* "Warning, unexpected situation with align ~a.~%" cur-align))
	      (t
	       (format *error-output* "Warning, unexpected situation of no alignment and no bam record.~%"))))
      good-alignment?)))

(defmethod output-to-intermediate-files ((streamer alignment-streamer))
  (with-slots (cur-align settings-hash cur-bam)
      streamer

    ;; Print out current alignment
    (when cur-align
      (let ((align-stream (gethash :align-stream settings-hash)))
	(when align-stream
	  (with-slots (read-name cigar seq strand ref-seq
				 reverse-complement-ref)
	      cur-bam
	    (format align-stream "READ:")
	    (pretty-print cur-bam align-stream)
	    (format align-stream "CIGAR:~a~%"  (cigar-str cigar))
	    (format align-stream "SEQ(~a): '~a'~%" strand (if (eql :+ strand)
						     seq
						     (reverse-complement seq)))
	    (format align-stream "TARGET: '~a'~%" (if (eql :+ strand)
					     ref-seq reverse-complement-ref))
	    (format align-stream "SEQ/TARGET lengths: ~a/~a~%" (length seq) (length ref-seq)))
	  (format align-stream "FLOWSPACE ALIGNMENT:~%")
	  (print-alignment cur-align t align-stream)
	  (format align-stream "FLOW INTENSITIES: ~{~a~^,~}~%" (flow-intensities cur-bam))
	  )))
    ))

(defmethod make-deviation-from-cur-bam-align ((streamer alignment-streamer))
  (with-slots (cur-bam cur-align cur-deviation last-bam
		       seq-deviation-cache temp-unsorted-cache
		       do-variant-calling? cur-process aligner-method)
      streamer
    (let (
	  last-ref
	  )
      (setq last-ref (when last-bam
		       (ref-name last-bam)))
      (setq do-variant-calling? nil)

      ;; New chromosome/contig so go ahead and call variants on the last one
      (when last-bam  ;; check it isn't the first contig
	(setq do-variant-calling?
	      (not (and last-ref
			(string= last-ref (ref-name cur-bam))))))
      (when (> (length seq-deviation-cache)
	       (gethash :max-deviation-cache-size (settings-hash streamer)))
	;; When cache size gets too big, force variant calling, even though
	;; it could be in the middle of a variant, b/c of memory constraints.
	;; There is now logic to handle being in the middle of the variant,
	;; by keeping in the temp-cache deviations that are beyond the read
	;; start of the current bam file.
	;; This doesn't seem to happen for the current defaults, but if cache
	;; size is still too small to handle very high coverage areas, it will
	;; force all sequence deviations to be used for variant calling. The
	;; end result would be the same indel being called twice.
	(setq do-variant-calling? :cache-exceeded))

      #|
      (when (equal (read-name cur-bam) "XOU43:473:67")
	(setq *debug-align* cur-align)
	(setq *debug-bam* cur-bam))
      ;; |#

      ;; Correction of the flow space alignment positions, due to trimming of hp at the end of the read
      #|
      (setq cur-process :correcting-flow-space-positions)
      (when (and cur-align cur-bam)
	(multiple-value-bind (found-issues? hp-bases-adjustment did-print?)
	    (print-n-check-target-seq-n-positions cur-align)
	  (when found-issues?
	    (alter-alignment-target-positions cur-align hp-bases-adjustment)
	    ;; check again if fixed worked
	    (when (print-n-check-target-seq-n-positions cur-align did-print?)
	      (format t "There are still issues with this alignment.~%")))))
      |#

      ;; Find the current deviation set for the alignment
      (setq cur-process :finding-deviations)
      (setq cur-deviation
	    (when (and cur-align cur-bam)
	      (unless (filtered-by-bam-attributes? cur-bam (settings-hash streamer))
		(find-seq-deviations cur-align cur-bam aligner-method))))

      ;; Do multiple processes on this deviation list, then
      ;; check (and then remove) any that have issues.
      (when cur-deviation
	;; Merge, trim, and make leftmost
	(setq cur-process :merge-n-trim-deviations)
	(setq cur-deviation (merge-deviation-list cur-deviation))
	(setq cur-deviation (mapcar #'trim-seq-deviation cur-deviation))
	(setq cur-deviation (remove-deviations-with-issues cur-deviation)))
      (when cur-deviation ;; Still deviations left w/o issues?
	(setq cur-process :find-leftmost-deviations)
	(setq cur-deviation (mapcar #'make-leftmost cur-deviation))
	;;(setq cur-deviation (mapcar #'set-rightmost-position cur-deviation))
	(setq cur-process :remove-deviations-with-issues)
	(setq cur-deviation (remove-deviations-with-issues cur-deviation)))
      (when cur-deviation ;; Still deviations left w/o issues?
	;; Put in in the cache
	(setq cur-process :placing-deviation-in-cache)

	(if (not do-variant-calling?)
	    ;; same chromosome, or it's the start
	    (dolist (cur-deviation-item cur-deviation)
	      (push cur-deviation-item seq-deviation-cache))
	    ))
      (when do-variant-calling?
	;; This needs to be done regardless if the cur-deviation had issues or not, so
	;; late 2.0 change was to bring this code out of the cur-deviation conditional.
	    ;; different chromosome, so call variants on completed chromosome
	    (progn
	      (when temp-unsorted-cache
		(error "Programming error, temp-unsorted-cache in make-deviation-from-cur-bam-align, should be empty here."))
	      (format t "ref comp: ~a ~a~%" last-ref (ref-name cur-bam))

	      ;; sort 
	      (when seq-deviation-cache
		(setq seq-deviation-cache (sort-deviations seq-deviation-cache))
		(when (eql do-variant-calling? :cache-exceeded)
		  (multiple-value-setq (seq-deviation-cache temp-unsorted-cache)
		    (extract-from-sorted-list (- (ref-pos cur-bam) 1) seq-deviation-cache
					      :from-end t :test #'> :key #'absolute-target-pos)))
		(when (< (* 4 (length seq-deviation-cache)) (length temp-unsorted-cache))
		  (format *error-output* "Warning, there may be repeated indel calls because not enough of the cache could be cleared.~%")
		  (setq seq-deviation-cache (append seq-deviation-cache temp-unsorted-cache))
		  (setq temp-unsorted-cache nil))
		)
	      (dolist (cur-deviation-item cur-deviation)
		(push cur-deviation-item temp-unsorted-cache))
	      (setq do-variant-calling? t)
	      )))
    (setq cur-process :finished-w-deviation)))

(defmethod process-variant-candidates ((streamer alignment-streamer) variant-candidates)
  (with-slots () ;;cur-process)
      streamer
    (dolist (var-cand variant-candidates)
      ;; Overall score
      ;;(setq cur-process :finding-scores)
      (calc-deviation-score var-cand)
      
      ;; Get the uniq set of sequence deviations
      ;;(setq cur-process :grouping-variant-sequences)
      (group-sequence-deviations var-cand)
      
      ;;(setq cur-process :processing-variant)
      (determine-strand-counts-n-bias var-cand)
      
      ;; Set rightmost
      (when (> (length (seq-deviations var-cand)) 1)
	;;(setq cur-process :making-rightmost)
	(set-rightmost-position var-cand))
      )
    ;; Apply all filters all execpt by frequency
    (dolist (var-cand variant-candidates)
      (filter-variants var-cand (settings-hash streamer)))
    
    ;; Delete variants that are filtered
    (let ((message (make-filtered-variant-count-statement variant-candidates nil)))
      (setq variant-candidates (delete t variant-candidates :key #'overall-filtered?))
      (format t "~a# of variants left after filters (all except var. freq. & max strand bias) applied: ~a~%"
	      message
	      (length variant-candidates)))

    ;;(setq cur-process :finding-spanning-reads-python)
    (when variant-candidates
      (extract-reads-over-regions streamer variant-candidates))
    ;;(setq cur-process :filtering-variants)
    
    
    (dolist (var-cand variant-candidates)
      (calc-spanning-reads-n-freqs var-cand)
      ;; now have strandness of spanning reads, so can calculate strand bias score  
      (determine-strand-counts-n-bias var-cand)
      ;; Final filter on freq. and strand bias score
      (filter-variants var-cand (settings-hash streamer)))
    (make-filtered-variant-count-statement variant-candidates t)
    ))

(defmethod reset-cache-after-variant-calling ((streamer alignment-streamer))
  (with-slots (seq-deviation-cache temp-unsorted-cache do-variant-calling?)
      streamer
    (setq seq-deviation-cache temp-unsorted-cache)
    (setq temp-unsorted-cache nil)
    (setq do-variant-calling? nil)
    ))

;; This will do it all, it will wash your car, it will clean up your house, it
;; will even do variant calling.  All you have to do is ask your torrent server.
(defmethod do-variant-calling ((streamer alignment-streamer))
  (with-slots (deviant-strm variant-strm cur-process)
      streamer
    ;; Retrieve the deviations and then sort them
    (with-slots (seq-deviation-cache)
	streamer
      ;;(setq *debug-seq-deviations* (copy-list seq-deviation-cache));; FOR DEBUG ONLY
      (unless seq-deviation-cache
	(format t "Note, no reads in alignment cache. Skipping.~%")
	(reset-cache-after-variant-calling streamer)
	(return-from do-variant-calling nil))
      ;;sorting done upstream now
      ;;(setq seq-deviation-cache (sort-deviations seq-deviation-cache))
      (format t "Performing variant calling from record ~a to ~a~%"
	      (get-basic-attrib-string (bam (car seq-deviation-cache)))
	      (get-basic-attrib-string (bam (car (last seq-deviation-cache)))))
      (format t "Total number of deviations found here: ~a~%"
	      (length seq-deviation-cache))

      (setq cur-process :making-variant-candidates)

;;      (setq *debug-seq-deviations*
;;	    (copy-list seq-deviation-cache))

      ;; Make variant candidate list and score each one
      (let ((variant-candidates
	     (make-variant-candidate-list-from-seq-deviations seq-deviation-cache)))
	;;(setq *debug-variant-candidates* (copy-list variant-candidates)) ;; FOR DEBUG ONLY
	(format t "A total of ~a variant candidates discovered.~%" (length variant-candidates))
	;;So for each variant candidate
	(setq cur-process :processing-variant-candidates)
	(multi-thread-mutating-function  (gethash :vh-num-threads
						  (settings-hash streamer))
					 #'process-variant-candidates
					 streamer variant-candidates)
	(dolist (var-cand variant-candidates)
	  ;; Print out candidate
	  (unless (overall-filtered? var-cand)
	    (setq cur-process :printing-variant)
	    ;; (push var-cand *chr1-73.642M-var-candidates*) ;; DEBUG!
	    (print-variant-candidate var-cand variant-strm)
	    (print-variant-candidate var-cand deviant-strm t))

	  ;; Memory clean up
	  ;;(unless (eql *debug-mode* :debug) ;; DEBUG
	    (clean-uniq-seq-n-seq-hash var-cand)
	  ;;)
	  )))

    ;; reset cache for the next contig
    (reset-cache-after-variant-calling streamer)))

(defvar *debug-streamer*)

;; main function that uses the streamer to make calls
(defun make-variant-calls-using-streamer (&key
					  sam-file
					  bam-file
					  merged-file
					  variant-file
					  (aligner-method :auto)
					  reference-file
					  align-file
					  key-seq
					  flow-order

					  max-deviation-cache-size 
					  vh-num-threads

					  ;; filters
					  min-mapq
					  score-threshold min-num-reads
					  min-variant-freq
					  max-strand-bias strand-prob
					  max-intensity-per-base

					  java-bin
					  java-bin-args
					  fs-align-jar
					  fs-align-jar-args
					  fs-align-num-threads
					  fs-align-range
					  python-bin
					  python-py

					  validate-alignments
					  )

  (let ((align-streamer (make-alignment-streamer
			 :sam-file sam-file
			 :bam-file bam-file
			 :aligner-method aligner-method
			 :flow-order flow-order
			 :reference-file reference-file
			 :align-file align-file
			 :key-seq key-seq

			 :max-deviation-cache-size max-deviation-cache-size 
			 :vh-num-threads vh-num-threads

			 ;; filters
			 :min-mapq min-mapq
			 :score-threshold score-threshold 
			 :min-num-reads  min-num-reads
			 :min-variant-freq min-variant-freq
			 :max-strand-bias max-strand-bias
			 :strand-prob strand-prob
			 :max-intensity-per-base max-intensity-per-base

			 :java-bin java-bin
			 :java-bin-args java-bin-args
			 :fs-align-jar fs-align-jar
			 :fs-align-jar-args fs-align-jar-args
			 :fs-align-num-threads fs-align-num-threads
			 :fs-align-range fs-align-range
			 :python-bin python-bin
			 :python-py python-py

			 :validate-alignments validate-alignments
			 ))
	)
    (setq *debug-streamer* align-streamer)
    (open-n-print-headers-intermediate-files align-streamer)

    ;; open up files for output of individual deviants, and variant list
    (with-open-file (merged merged-file :direction :output :if-exists :rename)
	(with-open-file (variant variant-file :direction :output :if-exists :rename)
	  ;; set the streams
	  (setf (deviant-strm align-streamer) merged)
	  (setf (variant-strm align-streamer) variant)
	  ;; Headers
	  (print-generic-header align-streamer merged :merged-deviation)
	  (print-generic-header align-streamer variant :variant)
	  (print-filter-settings (settings-hash align-streamer) merged)
	  (print-filter-settings (settings-hash align-streamer) variant)
	  (print-vcf-header-column-names variant)

	  ;; stream and call variants!
	  (stream-alignments align-streamer)))
    (close-streamer align-streamer)

    t))

;;
;; Old way:  Read in all bam records into memory, and then make alignments from them.
;; 
;; with list of bam objects, make a file with the sequence deviations
(defun write-alignments-to-file (bam-or-aligner-objects outfile &key return-align-objects? (flow-order *flow-order*) retain-dp-matrix? (find-deviations? t))
  (let ((counter 0)
	(start-time (get-universal-time))
	align-objects
	bam
	aligner
	all-dev-merged
	all-dev-unmerged)
    (with-open-file (out outfile :direction :output :if-exists :rename)
      (length 
       (mapc #'(lambda (bam-or-aligner)
		 ;;(setq  *test08-align-object-05* bam-or-aligner)
		 (incf counter)

		 ;; allows for input of bams or aligner objects
		 (if (eql 'bam-record (class-name (class-of bam-or-aligner)))
		     (setq bam bam-or-aligner
			   aligner nil)
		     (setq bam (bam bam-or-aligner)
			   aligner bam-or-aligner))

		 (with-slots (read-name cigar seq strand ref-seq
					reverse-complement-ref)
		     bam
		   (format out "READ:~a~aCIGAR:~a~%" read-name #\Tab (cigar-str cigar))
		   (format out "SEQ(~a): '~a'~%" strand (if (eql :+ strand)
							    seq
							    (reverse-complement seq)))
		   (format out "TARGET: '~a'~%" (if (eql :+ strand)
						    ref-seq reverse-complement-ref))
		   (format out "SEQ/TARGET lengths: ~a/~a~%" (length seq) (length ref-seq))
		   (let ((alignment (or aligner
					(do-alignment bam :flow-order flow-order ))))
		     (unless retain-dp-matrix?
		       ;; dp-matrix takes a lot of space, and typically not needed
		       ;; after alignment is done in do-alignment
		       (setf (dp-matrix alignment) nil))
		     (unless alignment
		       (format out "*** No original alignment ***~%"))
		     (when alignment
		       (print-alignment alignment t out)
		       (when return-align-objects?
			 (push alignment align-objects))
		       (let ((seq-dev (and find-deviations?
					   (find-seq-deviations alignment bam))))
			 (when seq-dev
			   (print-seq-deviations seq-dev out)
			   (format out "Deviations after merge:~%")
			   ;; In a read, merge adjacent deviations and make
			   ;; it leftmost.
			   (let ((merged-leftmost-dev
				  (mapcar #'make-leftmost
					  (merge-deviation-list seq-dev))))
			     ;; Set the rightmost position
			     ;;(mapcar #'set-rightmost-position merged-leftmost-dev)

			     (print-seq-deviations merged-leftmost-dev out t)
			     (setq all-dev-merged (append all-dev-merged merged-leftmost-dev))
			     (setq all-dev-unmerged (append all-dev-unmerged seq-dev))
			     )
			   ))))
		   (format out "~%")
		   (when (eql 0 (mod counter 250))
		     (with-slots (read-name ref-name ref-pos ref-end-pos)
			 bam
		       (let* ((cur-time (get-universal-time))
			      (elapsed-time (- cur-time start-time)))
			      
			 (format t "Processed ~a reads in ~a sec. (~,2f reads/sec). Last read: ~a, ~a:~a-~a~%"
			       counter elapsed-time
			       (unless (eql 0 elapsed-time)
				 (/ counter elapsed-time))
			       read-name ref-name ref-pos ref-end-pos))))
		   ))
	     bam-or-aligner-objects))
      ;;(format out "~%Complete sorted deviation list:~%")
      ;;(setq all-dev (sort all-dev #'< :key #'absolute-target-pos))
      ;;(print-seq-deviations all-dev out t)
      )
    (format t "Finished processing ~a bam records.~%" (length bam-or-aligner-objects))
    (values
     (when return-align-objects?
       (format t "~a aligned objects are being returned.~%" (length align-objects))
       (reverse align-objects))
     (when find-deviations?
       (sort-deviations all-dev-merged)
       (sort-deviations all-dev-unmerged))
     )))

;; OLD non-streaming method
;; From sequence deviations, print and then calculate variant list
(defun calc-n-print-deviations-n-variants (merged-deviations unmerged-deviations
					   merged-file unmerged-file variant-file)
  (let (
	variant-candidates
	)
    ;; write merged, leftmost adjusted deviation list
    ;;     (with-open-file (merged merged-file :direction :output :if-exists :rename)
    ;;	(print-seq-deviations merged-deviations merged t))

    ;; write unmerged devation
    (with-open-file (unmerged unmerged-file :direction :output :if-exists :rename)
      (print-seq-deviations unmerged-deviations unmerged))

    ;; Make variant candidate list and score each one
    (setq variant-candidates
	  (make-variant-candidate-list-from-seq-deviations merged-deviations))

    ;; open up files for output of individual deviants, and variant list
    (with-open-file (merged merged-file :direction :output :if-exists :rename)
	(with-open-file (variant variant-file :direction :output :if-exists :rename)

	  ;;So for each variant candidate
	  (dolist (var-cand variant-candidates)
	    ;; Overall score
	    (calc-deviation-score var-cand)
	    ;; Get the uniq set of sequence deviations
	    (group-sequence-deviations var-cand)
	    ;; Print out candidate
	    (print-variant-candidate var-cand variant)
	    (print-variant-candidate var-cand merged t)
	    )))
    (values
     variant-candidates)))

;; threading utility function
(defun wait-for-all-threads-to-finish (threads)
  (let ((running-threads threads)
	non-alive-threads
	last-non-alive-threads)
    (while (not (eql 0 (length running-threads)))
      (setq running-threads
	    (remove nil running-threads :key #'sb-thread:thread-alive-p))
      (setq non-alive-threads (set-difference threads running-threads))
      (unless (eql (length non-alive-threads)
		   (length last-non-alive-threads))
	(mapc #'(lambda (finished-thread)
		  (format t "~a thread has completed.~%" finished-thread))
	      (set-difference non-alive-threads last-non-alive-threads))
	(format t "~a thread~a still running.~%" (length running-threads)
		(if (eql 1 (length running-threads)) "" "s"))
	(setq last-non-alive-threads non-alive-threads))
      (sleep 3)))
  t)

;; multithread version
;; not up to date
(defun write-alignments-to-file-multitread-by-contig (bam-objects outfile contigs
						      &key (wait-for-all-threads? t)
						      (flow-order *flow-order*))
  (let ((threads 
	 (mapcar 
	  #'(lambda (contig)
	      (let ((contig-bam-objects
		     (extract-contig-bam-objects-from-sorted bam-objects contig)))
		(when contig-bam-objects
		  (sb-thread:make-thread
		   (lambda ()
		     (multiple-value-bind (base extension)
			 (split-filename-base-and-extension outfile)
		       (write-alignments-to-file
			contig-bam-objects
			(string-concat base "-" contig extension)
			:flow-order (copy-seq flow-order)))
		     )
		   :name contig))))
	  contigs)))
    (setq threads (remove nil threads))
    (format t "~a threads created from ~a contigs.~%" (length threads) (length contigs))
    (when wait-for-all-threads?
      (wait-for-all-threads-to-finish threads)
      )
    threads))
