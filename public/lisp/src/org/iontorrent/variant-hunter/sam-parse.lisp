;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

(in-package :cl-user)

;; requires util.lisp
;; requires parse-util.lisp
;; requires reference.lisp

(defvar *reference-padding-left* 5
  "Number of extra bases to include in the reference beyond what is indicated in the CIGAR.")
(defvar *reference-padding-right* 10)

;; BAM header class
(defclass bam-header ()
  ((headers :accessor headers :initarg :headers)
   (sequences :accessor sequences :initarg :sequences)))

(defun make-bam-header (sam-file)
  (let (headers
	sequences)
    (with-delimited-file (sam-line sam-file)
      (unless (eql #\@ (char (car sam-line) 0))
	(break-loop))
      (when (string= (car sam-line) "@HD")
	(push (cdr sam-line) headers))
      (when (string= (car sam-line) "@SQ")
	(mapc #'(lambda (column)
		  (when (string= "SN:" (subseq column 0 3))
		    (push (subseq column 3) sequences)))
	      (cdr sam-line))))
    (make-instance 'bam-header
		   :headers headers
		   :sequences sequences)))

;; CIGAR class
(defclass cigar ()
  ((cigar-str :accessor cigar-str :initarg :cigar-str)
   (elements :accessor elements :initarg :elements)
   (reversed :accessor reversed :initarg :reversed)))

(defun digits-to-number (digits)
  (if (cdr digits)
      (+ (car digits) (* 10 (digits-to-number (cdr digits))))
      (car digits)))

(defun split-cigar-into-reversed-elements ( cigar-string )
  (let (reversed-elements
	digit-list)
    (dotimes (pos (length cigar-string))
      (let ((cur-char (char cigar-string pos)))
	(cond ((<= 48 (char-code cur-char) 57) ;; it's a number
	       (push (- (char-code cur-char) 48) digit-list))
	      ((find cur-char "MIDNSHP")
	       (unless digit-list
		 (error "Cigar (~a) MUST start with a number." cigar-string))
	       (push (cons (digits-to-number digit-list)
			   cur-char) reversed-elements)
	       (setq digit-list nil)))))
    reversed-elements))

(defun make-cigar ( cigar-string )
  (let ((reversed-elements (split-cigar-into-reversed-elements cigar-string)))
    (make-instance 'cigar
		   :cigar-str cigar-string
		   :elements (reverse reversed-elements)
		   :reversed reversed-elements)))

(defgeneric get-read-clipping (cigar-obj strand))
(defmethod get-read-clipping ((cigar-obj cigar) strand)
  ;; see how many flows to ignore from the begining
  (let ((elements (if (eql strand :+)
		      (elements cigar-obj)
		      (reversed cigar-obj)))
	(clip-length 0))
    (while elements
      (if (find (cdr (car elements)) "HS")
	  (progn (incf clip-length (car (car elements)))
		 (pop elements))
	  (setq elements nil) ;; no more clipping found
	  ))
    clip-length))

(defgeneric get-total-clip-length (cigar-obj clip-type))
(defmethod get-total-clip-length ((cigar-obj cigar) clip-type)
  (let ((clip-length 0)
	(cigar-letters (case clip-type
			 (:both "HS")
			 (:soft "S")
			 (:hard "H"))))
    (dolist (element (elements cigar-obj))
      (when (find (cdr element) cigar-letters)
	(incf clip-length (car element))))
    clip-length))

(defgeneric get-begin-soft-clipping (cigar-obj strand &optional from-end?))
(defmethod get-begin-soft-clipping ((cigar-obj cigar) strand &optional from-end?)
  ;; see how much of the begining of the read is not aligned
  (when from-end?
    (setq strand (if (eql strand :+) :- :+)))
  (let ((elements (if (eql strand :+)
		      (elements cigar-obj)
		      (reversed cigar-obj)))
	(clip-length 0))
    (while elements
       (if (find (cdr (car elements)) "HS")
	   (progn
	     (when (eql (cdr (car elements)) #\S)
	       (incf clip-length (car (car elements))))
	     (pop elements))
	   (setq elements nil) ;; no more clipping found
	   ))
    clip-length))

(defgeneric get-reference-span-length (cigar-obj))
(defmethod get-reference-span-length ((cigar-obj cigar))
  ;; from CIGAR determine how many reference bases are included
  ;; in this CIGAR (adds M, D, N, and P)
  (let ((ref-span-length 0))
    (mapc #'(lambda (element)
	      (when (find (cdr element) "MDNP")
		(incf ref-span-length (car element)))
	      )
	  (elements cigar-obj))
    ref-span-length))

(defgeneric get-read-length (obj))
(defmethod get-read-length ((cigar-obj cigar))
  (let ((read-length 0))
    (mapc #'(lambda (element)
	      (when (find (cdr element) "MISH")
		(incf read-length (car element))))
	  (elements cigar-obj))
    read-length))

(defmethod get-read-length ((flow-intensities cons))
   (let ((read-length 0))
     (mapc #'(lambda (flow)
	       (incf read-length (round-half-up (/ flow 100))))
	   (cdr flow-intensities))
     (- read-length 4)))  ;; minus the key seq length


;; This just trims off the 0 intensity values at the end of the flow-intensity sequence
(defgeneric trim-flow-intensities (flow-intensities))
(defmethod trim-flow-intensities ((flow-intensities cons))
  (let ((trimmed-list (reverse (cdr flow-intensities))))
    (while (eql 0 (car trimmed-list))
      (pop trimmed-list))
    (reverse trimmed-list)))


(defmethod print-object ((cigar-obj cigar) strm)
  ;;(call-next-method)
  (format strm "#<CIGAR ~a>" (cigar-str cigar-obj)))

;; BAM record class
(defclass bam-record ()
  ((read-name :accessor read-name :initarg :read-name)
   (flag :accessor flag :initarg :flag)
   (strand :accessor strand :initarg :strand)
   (is-mapped? :accessor is-mapped? :initarg :is-mapped?)
   (ref-name :accessor ref-name :initarg :ref-name)
   (ref-pos :accessor ref-pos :initarg :ref-pos)
   (ref-end-pos :accessor ref-end-pos :initarg :ref-end-pos)
   (mapq :accessor mapq :initarg :mapq)
   (cigar :accessor cigar :initarg :cigar)
   (mate-ref-name :accessor mate-ref-name :initarg :mate-ref-name)
   (mate-ref-pos :accessor mate-ref-pos :initarg :mate-ref-pos)
   (insert-size :accessor insert-size :initarg :insert-size)
   (seq :accessor seq :initarg :seq)
   (qual :accessor qual :initarg :qual)
   (flow-intensities :accessor flow-intensities :initarg :flow-intensities)

   ;; These are not part of the original bam, but are useful additions
   ;; that are attributes the same read represented by the bam record

   ;; Score (sum of the residauls off from target)
   ;; Negative values, closer to 0 is closer to target
   ;; Only in the java aligner version
   (align-score :accessor align-score :initform nil)

   ;; Reference sequence around the bam record
   (ref-seq-pos :accessor ref-seq-pos :initarg :ref-seq-pos)
   (ref-seq :accessor ref-seq :initarg :ref-seq)
   (reverse-complement-ref :accessor reverse-complement-ref
			   :initarg :reverse-complement-ref)))

(defmethod print-object ((bam bam-record) strm)
  (call-next-method)
  (format strm "(~a[~a])" (sb-pcl::slot-value-or-default bam 'read-name)
	 (sb-pcl::slot-value-or-default bam 'strand ))
  )

(defgeneric get-basic-attrib-string (bam))
(defmethod get-basic-attrib-string ((bam bam-record))
  (with-slots (read-name strand ref-name ref-pos ref-end-pos)
      bam
    (format nil "~a(~a), ~a:~a-~a" read-name strand ref-name ref-pos ref-end-pos)))

(defgeneric pretty-print (obj &optional strm))
(defmethod pretty-print ((bam bam-record) &optional (strm t))
  (format strm "~a~%" (get-basic-attrib-string bam)))

(defun extract-bit (bit-wanted num)
  (mod (ash num (- bit-wanted)) 2))

(defun extract-strand (flag)
  (let ((strand-num (extract-bit 4 flag)))
    (if (eql strand-num 0)
	:+
	:-)))

(defun extract-mapped (flag)
  (let ((mapped-bit (extract-bit 2 flag)))
    (if (eql mapped-bit 0)
	:mapped
	:unmapped)))

(defun make-bam-record (sam-line &optional reference-reader)
  (let ((read-name (pop sam-line))
	(flag (parse-integer (pop sam-line)))
	strand
	mapped-status
	(ref-name (pop sam-line))
	(ref-pos (pop sam-line)) ;; position of the start of the read
	ref-end-pos
	(mapq (pop sam-line))
	(cigar (make-cigar (pop sam-line)))
	(mate-ref-name (pop sam-line))
	(mate-ref-pos (pop sam-line))
	(insert-size (pop sam-line))
	seq ;;(pop sam-line))  ;; To save memory, don't store
	qual ;;(pop sam-line)) ;; To save memory, don't store
	flow-intensities
	ref-seq
	ref-seq-pos  ;; position at the start of the ref sequence that's present
	)
    ;;(setq total-clip-length (get-total-clip-length cigar :both))

    ;; To save memory, don't store seq and qual!
    (pop sam-line)
    (pop sam-line)
    

    (dolist (flag sam-line)
      (cond ((string= (subseq flag 0 2) "FZ")
	     (setq flow-intensities (subseq flag 5)) ; Faster than (third (parse-string flag #\:)))
	     )))
    (setq strand (extract-strand flag))
    (setq mapped-status (extract-mapped flag))
    (setq flow-intensities
	  (read-from-string (format nil "(~a)" (substitute #\Space #\, flow-intensities))))
    (unless (string= ref-pos "0")
      ;;(format t "ref-name/pos = '~a'/'~a'~%" ref-name ref-pos)
      (setq ref-pos (parse-integer ref-pos))
      (setq ref-end-pos (+ ref-pos (get-reference-span-length cigar) -1))

      (let ((start-pad (if (eql :- strand)
			   *reference-padding-right*
			   *reference-padding-left*))
	    (end-pad (if (eql :+ strand)
			 *reference-padding-right*
			 *reference-padding-left*)))
	(when reference-reader
	  (multiple-value-bind (seq pos)
	      (extract-reference reference-reader
				 ref-name ref-pos
				 ref-end-pos
				 start-pad end-pad
				 )
	    (setq ref-seq (string-upcase seq))
	    (setq ref-seq-pos pos)))))
    (make-instance 'bam-record
		   :read-name read-name
		   :flag flag
		   :strand strand
		   :is-mapped? (eql mapped-status :mapped)
		   :ref-name ref-name
		   :ref-pos ref-pos
		   :ref-end-pos ref-end-pos
		   :mapq (parse-integer mapq)
		   :cigar cigar
		   :mate-ref-name mate-ref-name
		   :mate-ref-pos mate-ref-pos
		   :insert-size insert-size
		   :seq seq
		   :qual qual
		   :flow-intensities flow-intensities
		   :ref-seq ref-seq
		   :ref-seq-pos ref-seq-pos
		   :reverse-complement-ref (reverse-complement ref-seq)
		   )))

(defgeneric extract-reference-from-bam-obj (obj start-pos end-pos))
(defmethod extract-reference-from-bam-obj ((bam bam-record) start-pos end-pos)
  (with-slots (ref-seq ref-seq-pos)
      bam
    (let ((start (- start-pos ref-seq-pos))
	  (end (- end-pos ref-seq-pos))
	  extracted-ref)
      (when (or (< start 0)
		(< end 0)
		(> start (length ref-seq))
		(>= end (length ref-seq)))
	(format t "Note, bam record ~a which has ref starting at ~a with a length of ~a does not have sequence between ~a and ~a.~%" bam ref-seq-pos (length ref-seq) start-pos end-pos)
	(return-from extract-reference-from-bam-obj :no-reference))
      (when (< (- end start) -1)
	(format t "Note, bam record ~a, tried to extract a negative sized ref. sequence (start: ~a, end: ~a~%" bam start end)
	;; Note, end can be 1 less than start, since that's just an empty reference.
	(return-from extract-reference-from-bam-obj :no-reference))

      (setq extracted-ref (make-array (1+ (- end start)) :element-type 'standard-char
				      :initial-element #\?))
      (loop for rel-pos from start to end do
	   (let ((new-ref-pos (- rel-pos start)))
	     (setf (char extracted-ref new-ref-pos)
		   (char ref-seq rel-pos))))
      extracted-ref
      )))

(defvar *sam-file-debug* "/data/antioch/projects/ion_variant_calling/output_w_flowspace.sam")
(setq *sam-file-debug* "../data/CFTR.exon8-9.1570s.sam")

(defvar *reference-debug* "../data/cftr_24amp_v2.fasta")

(defun make-objects-from-sam-file (&key
				   (sam-file *sam-file-debug*)
				   (reference-file *reference-debug*)
				   (keep-non-aligned? t)
				   )
  (let ((num-header-lines (determine-number-header-lines sam-file :header-char #\@))
	(counter 0)
	bam-objects
	(ref-reader (make-reference-reader reference-file)))
    (format t "Found ~a header lines.~%" num-header-lines)
    (with-delimited-file (samLine sam-file num-header-lines)
      (incf counter)
      (let ((bam (make-bam-record samLine ref-reader)))
	(when (or keep-non-aligned?
		  (is-mapped? bam))
	  (push bam bam-objects))))
    (format t "Processed ~a sam records.~%" counter)
    (reverse
     bam-objects)))

(defun extract-contig-bam-objects-from-sorted (bam-objects contig)
  ;; Assuming that the bam objects are grouped by contig, find the start and end positions
  ;; of that particular contig, and return the subseq of bam objects for that contig
  ;; nil start position indicates contig not found
  ;; nil end position with a non-nil start position indicates contig goes to the end of the file
  (let* ((start-pos (position contig bam-objects
			      :key #'ref-name :test #'equal))
	 (end-pos (when start-pos
		    (position contig bam-objects
			      :start start-pos
			      :key #'ref-name :test #'(lambda (x y)
							(not (equal x y)))))))
    (values
     (when start-pos
       (subseq bam-objects start-pos end-pos))
     start-pos
     end-pos)))

;; Print attributes from a set of bam records
(defvar *bam-printing-attributes*)
(setq *bam-printing-attributes*
      '(read-name cigar ref-pos ref-end-pos))

(defgeneric print-attributes-of-bam-records (bam-objects stream &optional attributes))

(defmethod print-attributes-of-bam-records ((bam-objects cons)
					    stream
					    &optional (attributes
						       *bam-printing-attributes*
						       ))
  (let ((slots-n-values (get-object-values-for-slot-list bam-objects attributes)))
    (dolist (slot-n-values slots-n-values)
      (let ((slot (car slot-n-values))
	    (values (cdr slot-n-values)))
	(cond ((eql slot 'read-name)
	       (setq slot 'reads-names))
	      ((eql slot 'ref-pos)
	       (setq slot 'read-start-positions))
	      ((eql slot 'ref-end-pos)
	       (setq slot 'read-end-positions)))
	(format stream "~a="
		(string-downcase (format nil "~a" slot)
				 :start 1))
	(cond ((eql slot 'cigar)
	       (format stream "~{~a~^,~}" (mapcar #'cigar-str values)))
	      (t
	       (format stream "~{~a~^,~}" values))))
      (format stream ";"))))

;; Filters against certain BAM properties
(defgeneric filtered-by-bam-attributes? (bam filter-hash))
(defmethod filtered-by-bam-attributes? ((bam bam-record) (filter-hash hash-table))
  (let (do-filter?)
    (with-slots (mapq)
	bam
      (when (< mapq (gethash :min-mapq filter-hash))
	(setq do-filter? t))
      )
    do-filter?))

