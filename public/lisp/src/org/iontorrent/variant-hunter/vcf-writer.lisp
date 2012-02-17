;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

(in-package :cl-user)

;; requires util.lisp
;; requires parse-util.lisp
;; requires seq-deviations.lisp

(defclass vcf-tag ()
  ((vcf-tag-category :accessor vcf-category :initarg :vcf-tag-category)
   (vcf-tag-id :accessor vcf-tag-id :initarg :vcf-tag-id)
   (vcf-tag-number :accessor vcf-tag-number :initarg :vcf-tag-number)
   (vcf-tag-datatype :accessor vcf-tag-datatype :initarg :vcf-tag-datatype)
   (vcf-tag-desc :accessor vcf-tag-desc :initarg :Vcf-tag-desc)))

(defun make-vcf-tag-obj (vcf-tag-category vcf-tag-id vcf-tag-number vcf-tag-datatype vcf-tag-desc)
  (make-instance 'vcf-tag
		 :vcf-tag-category vcf-tag-category
		 :vcf-tag-id vcf-tag-id
		 :vcf-tag-number vcf-tag-number
		 :vcf-tag-datatype vcf-tag-datatype
		 :vcf-tag-desc vcf-tag-desc))

(defun make-vcf-tag-objs (list-of-vcf-tag-attribs)
  (mapcar #'(lambda (vcf-tag-item-attribs)
	      (apply #'make-vcf-tag-obj vcf-tag-item-attribs))
	  list-of-vcf-tag-attribs))

;; printing tags for the header
(defgeneric print-vcf-tag-for-vcf-header (vcf-tag stream)) 
(defmethod print-vcf-tag-for-vcf-header ((vcf-tag vcf-tag) stream)
  (with-slots (vcf-tag-category vcf-tag-id vcf-tag-datatype
				vcf-tag-number vcf-tag-desc)
      vcf-tag
    (let ((vcf-tag-datatype-string
	   (format nil "~a" vcf-tag-datatype)))
      (setq vcf-tag-datatype-string
	    (format nil "~a~a"
		    (string-upcase (char vcf-tag-datatype-string 0))
		    (string-downcase (subseq vcf-tag-datatype-string 1))))
      (format stream "##~a=<ID=~a,Number=~a,Type=~a,Description=\"~a\">~%"
	      vcf-tag-category vcf-tag-id (or vcf-tag-number ".")
	      vcf-tag-datatype-string
	      vcf-tag-desc))))

(defgeneric print-vcf-tags-for-vcf-header (vcf-tags stream))
(defmethod print-vcf-tags-for-vcf-header ((vcf-tags cons) stream)
  ;; Takes a list of vcf-tags
  (dolist (vcf-tag vcf-tags)
    (print-vcf-tag-for-vcf-header vcf-tag stream)))

(defun make-hash-of-vcf-header-tags ()
  (let ((vcf-header-tag-hash (make-hash-table)))
    (labels ((populate-tag-hash (tag-set tag-infos)
	       (setf (gethash tag-set vcf-header-tag-hash)
		     (make-vcf-tag-objs tag-infos))))
      ;; Dealing with Score
      (populate-tag-hash
       :SCORE 
       '(
	 (:INFO "Bayesian_Score" 1 :FLOAT "Bayesian score")
	 (:INFO "Score" 1 :FLOAT "Flow based score")
	 ))
      ;; Dealing with other parts of INFO
      (populate-tag-hash
       :INFO 
       (append
	'(
	  (:INFO "Variant-freq" 1 :FLOAT "Number of variant reads / number of spanning reads")
	  (:INFO "Num-spanning-reads" 1 :INTEGER "Number of reads that span both LM and RM positions")
	  (:INFO "Num-spanning-by-strand" nil :STRING "Above split by +/- strand")
	  (:INFO "Num-spanning-ref-reads" 1 :INTEGER "Number of spanning reads that are reference")
	  (:INFO "Num-variant-reads" nil :STRING "List of number of variant reads")
	  (:INFO "Zygosity" nil :STRING "Zygosity call")
	  (:INFO "LMPos" 1 :INTEGER "Leftmost position")
	  (:INFO "RMPos" 1 :INTEGER "Rightmost position")
	  (:INFO "Variants" nil :STRING "All the variants found, reported as variant/reference, sorted by Variants-scores. Only the highest scoring one is reported in the ALT column.")
	  (:INFO "Variants-scores" nil :STRING "Score for each variant found")
	  (:INFO "Variants-RMPos" nil :STRING "Rightmost position of those individual variants")
	  (:INFO "Variants-freqs" nil :STRING "Variant frequencies")
	  (:INFO "Num-reads" nil :STRING "Number of reads for each variant")
	  (:INFO "Plus-minus-strand-counts" nil :STRING "Counts on each strand, +/-")
	  )
	(list
	 (list :INFO "Strand-bias-scores" nil :STRING (format nil "Strand bias score for each variant. Strand bias = 2 * abs ( varRatio - overallRatio), where varRatio and overallRatio are the proportion of plus strand hits for the variants and for all the spanning reads, respectively.  A minimum of ~a reads for both the variant and spanning reads required, otherwise value will be NIL." *min-reads-for-sb-bias-calculation*)))
	'(
	  (:INFO "SB" nil :STRING "Strand bias score for the single variant reported in the ALT column.")
	  (:INFO "Strand-probabilities" nil :STRING "Bionomial probability of each variant seq.")
	  (:INFO "Strand-probability" 1 :FLOAT "Overall bionomial probability")
	  (:INFO "Read-names" nil :STRING "Read names that called the variant")
	  (:INFO "Cigars" nil :STRING "CIGAR strings of those reads")
	  (:INFO "MapQs" nil :STRING "MapQs of those reads")
	  (:INFO "MQ" 1 :FLOAT "RMS mapping quality")
	  (:INFO "MQ0" 1 :INTEGER "Number of MAPQ == 0 reads covering this record")
	  )))
      ;; Dealing with Genotype
      (populate-tag-hash
       :GENOTYPE
       '(
	 (:FORMAT :GT 1 :STRING "Genotype")
	 (:FORMAT :GQ 1 :FLOAT "Genotype Quality")
	 (:FORMAT :GL nil :STRING "Genotype Likelihood, number of values is (#ALT+1)*(#ALT+2)/2")
	 (:FORMAT :DP 1 :INTEGER "Read Depth")
	 (:FORMAT :FDP 1 :INTEGER "Filtered Read Depth")
	 (:FORMAT :AD nil :INTEGER "Allelic depths for the REF and ALT alleles in the order listed in the ALT field")
	 (:FORMAT :AST nil :INTEGER "Allelic unique start positions for the REF and ALT alleles in the order listed in the ALT field")
	 (:FORMAT :ABQV nil :INTEGER "Allelic average base qv for the REF and ALT alleles in the order listed  in the ALT field")
	 ))
      vcf-header-tag-hash
      )))

;; Helper function to make code for above function
(defun make-vcf-structured-list  (tag-set column-name ids is-numbs dtypes descriptions)
  (labels ((add-quotes-if-string (item)
	     (cond ((not item)
		    (format nil "nil"))
		   ((stringp item)
		    (format nil "\"~a\"" item))
		   ((or (keywordp item)
			(symbolp item))
		    (format nil ":~a" item))
		   (t
		    item))))
    (format t "(populate-tag-hash~% :~a '(~%" tag-set)
    (mapc #'(lambda (id is-numb dtype description)
	      (setq dtype (or dtype "STRING"))
	      (format t "   (~{~a~^ ~})~%"
		      (mapcar #'add-quotes-if-string
			      (list
			       column-name 
			       id
			       is-numb
			       (intern (string-upcase dtype) :keyword)
			       description))))
	  ids is-numbs dtypes descriptions)
    (format t "   ))~%")))

;; Populate settings-hash with tags
(defun populate-settings-hash-with-vcf-tag-hash (settings-hash)
  (setf (gethash :vcf-tag-hash settings-hash)
	(make-hash-of-vcf-header-tags)))

;; Set global settings-hash with this
;; vcf-tag-format is a comma seperated list of any or all of
;; score,info,genotype,indiv-read-info
(defun set-vcf-tag-format (settings-hash vcf-tag-format)
  (setf (gethash :vcf-tag-format settings-hash) (parse-string-to-intern-list vcf-tag-format))
  (populate-settings-hash-with-vcf-tag-hash settings-hash)
  settings-hash)

;; Actually prints the vcf header for tags
(defun print-vcf-header-tag-description (settings-hash stream)

  (let ((vcf-tag-hash (gethash :vcf-tag-hash settings-hash)))
    (dolist (tag-set '(:score :info :genotype))
      (when (find tag-set (gethash :vcf-tag-format settings-hash))
	(print-vcf-tags-for-vcf-header (gethash tag-set vcf-tag-hash) stream)))
    ))

(defun print-vcf-header-column-names (settings-hash stream)
  (let ((col-names '(CHROM POS ID REF ALT QUAL FILTER INFO)))
    (when (find :genotype (gethash :vcf-tag-format settings-hash))
      (setq col-names (append col-names '(FORMAT "Sample_1"))))
    (format stream (format nil "#~~{~~a~~^~a~~}~%" #\Tab)
	    col-names)))

(defvar *vcf-list-size-limit* 250)

(defgeneric print-variant-candidate (var-cand stream settings-hash &optional with-deviants?))
(defmethod print-variant-candidate ((var-cand variant-candidate) stream
				    settings-hash
				    &optional with-deviants?)
  (with-slots (ref-name leftmost-pos seq-deviations score
			scores-n-q-t-seqs seq-dev-hash
			max-rightmost rightmost-positions
			num-plus-spanning-reads
			num-spanning-ref-reads
			num-spanning-reads variant-freq
			num-variant-reads
			variant-freqs-by-score
			zygosity-call)
      var-cand

    ;;(format t "Candidate with ~a reads on ~a @ position ~a~%"
    ;;(length seq-deviations) ref-name leftmost-pos)
    
    ;; CHROM
    (format stream "~a~a" ref-name #\Tab)
    ;; POS [Leftmost]
    (format stream "~a~a" leftmost-pos #\Tab)
    ;; ID
    (format stream "~a~a" "." #\Tab)

    ;; ToDo, for ref/alt alleles, only top scoring candidate listed
    (let* ((variants (mapcar #'cdr scores-n-q-t-seqs))
	   (top-pair (car variants))
	   (top-var (car top-pair))
	   (top-ref (second top-pair)))
      ;; REF
      ;; TODO: For InDels or larger structural variants, the reference String must include the base before the event (which must be reflected in the POS field).
      (format stream "~a~a" top-ref  #\Tab)
      ;; ALT, TODO, need to list other alternate alleles.
      (format stream "~a~a" top-var #\Tab))

    ;; QUAL, TODO, will be some function of score
    (format stream "~a~a" "." #\Tab)
    ;; FILTER
    (format stream "~a~a" "." #\Tab)

    ;; INFO
    (when (find :score (gethash :vcf-tag-format settings-hash))
      (format stream "Score=~3$;" score))

    (when (find :info (gethash :vcf-tag-format settings-hash))

    (format stream "Variant-freq=~3$;" (or variant-freq :unknown))
    (format stream "Num-spanning-reads=~a;" num-spanning-reads)
    (format stream "Num-spanning-by-strand=~a/~a;" num-plus-spanning-reads
	    (- num-spanning-reads num-plus-spanning-reads))
    (format stream "Num-spanning-ref-reads=~a;" num-spanning-ref-reads)
    (format stream "Num-variant-reads=~{~a~^,~};" num-variant-reads)
    (let ((orig-calculation
	   (mapcar #'(lambda (q-t-seq)
		       (length (gethash q-t-seq seq-dev-hash)))
		   (mapcar #'cdr scores-n-q-t-seqs))))
      (unless (equal num-variant-reads
		     orig-calculation)
	(format *error-output* "WARNING, old 0.1.2 Num-variant-reads calculation does not match current one.~%")
	(format *error-output* "num-variant-reads: ~{~a~^,~}~%" num-variant-reads)
	(format *error-output* "orig calculation:  ~{~a~^,~}~%" orig-calculation)
	))
    (format stream "Zygosity=~a;" zygosity-call)

    (format stream "LMPos=~a;RMPos=~a;"  leftmost-pos max-rightmost)
    (format stream "Variants=~{~{~a~^/~}~^,~};" (mapcar #'cdr scores-n-q-t-seqs))
    (format stream "Variants-scores=~{~3$~^,~};" (mapcar #'car scores-n-q-t-seqs))
    (format stream "Variants-RMPos=~{~a~^,~};" rightmost-positions)
    (format stream "Variants-freqs=~{~3$~^,~};" variant-freqs-by-score) 

    ;; strand counts and bias
    (with-slots (strand-bias-scores strand-counts strand-probs overall-strand-prob)
	var-cand
	;;(determine-strand-counts-n-bias var-cand)
      (format stream "Plus-minus-strand-counts=~{~{~a~^/~}~^,~};"
	      strand-counts)
      (format stream "Strand-bias-scores=~{~a~^,~};"
	      (mapcar #'(lambda (prob)
			  (format-scientific prob 4)) strand-bias-scores))
      (format stream "Strand-probabilities=~{~a~^,~};"
	      (mapcar #'(lambda (prob)
			  (format-scientific prob 4)) strand-probs))
      (format stream "SB=~a;" (format-scientific (car strand-bias-scores) 4))
      (format stream "Strand-probability=~a;" (format-scientific overall-strand-prob 4)))

    ;; Print out some BAM attributes
    (flet ((get-tag-values (accessor)
	     (mapcar
	      #'(lambda (q-t-seq)
		  (mapcar accessor
			  (mapcar #'bam (gethash q-t-seq seq-dev-hash))))
	      (mapcar #'cdr scores-n-q-t-seqs))))
      (flet ((print-bam-tag (pretty-name accessor)
	       (format stream "~a=~{~{~a~^/~}~^,~};"
		       pretty-name
		       (limit-size-lists-in-list
			(get-tag-values accessor) *vcf-list-size-limit*)
		       )))
	(when (find :indiv-read-info (gethash :vcf-tag-format settings-hash))
	  (print-bam-tag "Read-names" #'read-name)
	  (print-bam-tag "Cigars" #'(lambda (item) (cigar-str (cigar item))))
	  ;; mapq and MQ and MQ0
	  (print-bam-tag "MapQs" #'mapq))

	(multiple-value-bind (rms mapq0-counts)
	    ;; TODO, for now, take top variant only
	    (when (get-tag-values #'mapq)
	      (determine-map-qv-rms-n-0-counts (get-tag-values #'mapq) 1))
	  (format stream "MQ=~,4f;" rms)
	  (format stream "MQ0=~a;" mapq0-counts))
	)
      )
    )
    (when (find :genotype (gethash :vcf-tag-format settings-hash))
      ;; Format field (9th column)
      (format stream "~aGT:GQ:GL:DP:FDP:AD:AST:ABQV~a" #\Tab #\Tab)
      ;; Genotype field (10th column)
      (format stream "~a:" (zygosity-to-genotype-call zygosity-call))  ;; GT
      (format stream "~a:" ".") ;; GQ
      (format stream "~a:" ".") ;; GL
      (format stream "~a:" (+  num-spanning-ref-reads (car num-variant-reads))) ;; DP
      (format stream "~a:" ".") ;; FPD
      (format stream "~a,~a:" num-spanning-ref-reads (car num-variant-reads))   ;; AD
      ;; TODO, DP and AD only include counts for only one of the variant sequences, whereas
      ;; variant-freq includes all variant sequences.  This is because VAR column
      ;; only includes only one of the variants.
      (format stream "~a:" ".") ;; AST
      (format stream "~a" ".")  ;; ABQV
      )
    (when with-deviants?
      (format stream "~%")
      ;;(format stream "Score=~a, candidate with ~a reads on ~a @ position ~a~%"
      ;;        score
      ;;        (length seq-deviations) ref-name leftmost-pos)
      (print-seq-deviations seq-deviations stream t))
    (format stream "~%")
    ))


#|
Tag name               Description
--------               -----------
Variant-freq           Number of variant reads / number of spanning reads  0.006
Num-spanning-reads     Number of reads that span both LM and RM positions  2206
Num-spanning-by-strand Above split by +/- strand                           1045/1161
Num-spanning-ref-reads Number of spanning reads that are reference         2196 [guessed]
LMPos                  Leftmost position                                   115256516
RMPos                  Rightmost position                                  115256517
Variants               All the variants found, variant/reference           G/A,T/A,AG/A
Variants-scores        Score for each variant found                        6.881,3.991,1.543
Variants-RMPos         Rightmost position of those individual variants     
Variants-freqs         Variant frequencies                                 0.004,0.001,0.001
Num-reads              Number of reads for each variant                    8,3,2
Plus-minus-strand-counts     Counts on each strand, +/-                    4/4,1/2,2/0
Strand-probabilities   Bionomial probability of each variant seq.          1.0000e+0,1.0000e+0,5.0000e-1
Strand-probability     Overall bionomial probability                       7.91E-01
                       Bionomial probabilities are of limited use because very sensitive to strand overall probability.
Read-names             Read names that called the variant
                                                     30O8V:1079:936/30O8V:576:772/30O8V:235:53/30O8V:483:549/30O8V:805:780/30O8V:594:98/30O8V:44:460/30O8V:1066:1189,30O8V:920:1002/30O8V:162:635/30O8V:412:986,30O8V:294:510/30O8V:203:505
Cigars                 CIGAR strings of those reads  53H3S50M/55H51M/55H65M/53H53M/65M55H/54M52H/54M52H/65M55H,55H51M/52H54M/52M56H,14M1I40M53H/14M1I40M54H
MapQs                  MapQs of those reads          36/39/84/35/84/80/80/77,39/36/72,62/62
MQ                     In VCF spec                   67.8997
MQ0                    In VCF spec                   0

For Read-names, Cigars, and MapQs, items are for each read, and seperated by variant with comma (,) and within a variant by slash (/).

|#
