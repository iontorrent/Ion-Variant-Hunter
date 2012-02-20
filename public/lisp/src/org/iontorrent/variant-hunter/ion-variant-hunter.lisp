;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

(in-package :cl-user)

;; requires util.lisp
;; requires stats.lisp
;; requires parse-util.lisp
;; requires reference.lisp
;; requires flow-space.lisp
;; requires sam-parse.lisp
;; requires seq-deviations.lisp
;; requires alignment-streamer.lisp

(defun make-individual-filenames-for-hunter (base-output-filename)
  (let ((exts '("align" "merged.dev" "unmerged.dev" "vcf")))
    (values-list
     (mapcar #'(lambda (ext)
		 (format nil "~a.~a" base-output-filename ext))
	     exts))))

;; From sam file, make bam objects, and then take the objects
(defun write-alignment-n-deviation-lists-from-sam (&key
						   sam-file ;; input sam file
						   bam-file ;; java can take in bam
						   reference-file
						   base-output-filename
						   (flow-order *flow-order*)
						   (key-seq "TCAG")
						   (aligner-method :auto)

						   (max-deviation-cache-size 300000)
						   (vh-num-threads 4)

						   ;; vcf format
						   (vcf-tag-format "score,info,genotype")

						   ;; filters
						   (min-mapq 0)
						   (score-threshold 4)
						   (min-num-reads 3)
						   (min-variant-freq 1/25)
						   (max-strand-bias 1/4)
						   (strand-prob 0)
						   (max-intensity-per-base
						    "A:600,T:600,G:600,C:600")

						   java-bin
						   java-bin-args
						   fs-align-jar
						   fs-align-jar-args
						   fs-align-range
						   (fs-align-num-threads 4)
						   python-bin
						   python-py
						   (samtools-bin "/usr/bin/samtools")
						   
						   return-align-objects?
						   retain-dp-matrix?
						   (find-deviations? t)
						   (streaming-method? t)
						   (print-flow-space-aligns? nil)
						   return-values?

						   validate-alignments)
  (unless base-output-filename
    (error "ERROR: Please specify --base-output-filename~%"))

  (let (;; Setting output files
	align-file
	merged-file
	unmerged-file
	variant-file
	bam-objects
	variant-candidates

	align-objects
	merged-deviations
	unmerged-deviations

	)
    (multiple-value-setq (align-file merged-file unmerged-file variant-file)
      (make-individual-filenames-for-hunter base-output-filename))

    (format *error-output* "Current working directory: ~a~%" *default-pathname-defaults*)

    (when streaming-method?
      (multiple-value-setq (merged-deviations variant-candidates)
	(make-variant-calls-using-streamer :sam-file sam-file
					   :bam-file bam-file
					   :merged-file merged-file
					   :variant-file variant-file
					   :aligner-method aligner-method
					   :flow-order flow-order
					   :reference-file reference-file
					   :align-file (when print-flow-space-aligns?
							 align-file)
					   :key-seq key-seq
					   :flow-order flow-order

					   :max-deviation-cache-size max-deviation-cache-size
					   :vh-num-threads vh-num-threads

					   ;; vcf format
					   :vcf-tag-format vcf-tag-format

					   ;;filters
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
					   :samtools-bin samtools-bin

					   :validate-alignments validate-alignments
					   )))
    (unless streaming-method?
      ;; getting bam objects
      (setq bam-objects (make-objects-from-sam-file :sam-file sam-file
						    :reference-file reference-file
						    :keep-non-aligned? nil))
      ;; Perform realingment and write output files
      (multiple-value-setq (align-objects merged-deviations unmerged-deviations)
	(write-alignments-to-file bam-objects align-file
				  :return-align-objects? return-align-objects?
				  :flow-order flow-order
				  :retain-dp-matrix? retain-dp-matrix?
				  :find-deviations? find-deviations?))
      (when find-deviations?
	(setq variant-candidates
	      (calc-n-print-deviations-n-variants merged-deviations
						  unmerged-deviations
						  merged-file
						  unmerged-file
						  variant-file))))

    ;;finished
    (format *error-output* "Completed.~%")

    ;; return values
    (if return-values?
	(values align-objects
		bam-objects
		merged-deviations
		unmerged-deviations)
	t)))

(defun write-alignment-n-deviation-lists-from-sam-usage (&optional (stream t) more-info?)
  (format stream "Usage:
./ion-variant-hunter-core  --function write-alignment-n-deviation-lists-from-sam --bam-file {bam-file} --reference-file {ref-file} --base-output-filename {output}~%")
  (when more-info?
    (format stream "For more info, use --detailed-help or --default-value-help~%~%")))

(defun write-alignment-n-deviation-lists-from-sam-doc (key-value &optional (stream t))  
  (write-alignment-n-deviation-lists-from-sam-usage stream)
  (when (stringp key-value)
    (setq key-value (intern (string-upcase key-value) :keyword)))

  (labels ((print-info (keyword description)
	     (when (or (eql key-value :all)
		       (eql keyword key-value))
	       (let ((arg-name (string-downcase (format nil "~a" keyword))))
		 (format stream "--~a " arg-name)
		 (dotimes (x (- 24 (length arg-name)))
		   (format stream " "))
		 (format stream "~a~%"
			 description)))))
    (when (eql key-value :all)
      (format stream "~%Detailed description of each command line argument:~%")
      (format stream "For default values, use --detailed-help.~%~%")
      (format stream "Required options:~%"))

    (print-info :bam-file "Input BAM file for variant hunter. REQUIRED.")
    (print-info :reference-file "Input reference file.  REQUIRED.  Associated .fai and .dict also are required to be present.")
    (print-info :base-output-filename "Base filename for output files (i.e. vcf/align).  REQUIRED")

    (when (eql key-value :all)
      (format stream "~%Memory options:~%"))
    (print-info :max-deviation-cache-size "The maximum number of reads to hold before variant calling is forced on.  Normally this is done only when a different contig is reached.")
    (print-info :dynamic-space-size "Dynamic memory space in megabytes.")
    (print-info :vh-num-threads "Variant hunter number of threads.")

    (when (eql key-value :all)
      (format stream "~%VCF output options:~%"))
    (print-info :vcf-tag-format "Comma separated list of desired classes of tags to be printed to the VCF file, i.e. 'score,info,genotype,indiv-read-info'")
    (format stream "    score            Displays scores~%")
    (format stream "    info             Displays non-score INFO tags~%")
    (format stream "    genotype         FORMAT and GENOTYPE columns~%")
    (format stream "    indiv-read-info  Individual read information (i.e. name/CIGAR)~%")

    (when (eql key-value :all)
      (format stream "~%Filter options:~%"))
    (print-info :score-threshold "Min score required for any variant to be reported.")
    (print-info :min-num-reads "Min number of reads required for each variant.")
    (print-info :min-variant-freq "Min variant frequency.")
    (print-info :max-strand-bias "Max allowed strand bias score which currently ranges from 0 to 2.")
    (print-info :strand-prob "Min. bionomical strand probability.")
    (print-info :max-intensity-per-base "Max intensity value in order to call a variant that is two less than the reference or longer.")
    (print-info :min-mapq "Min. MAPQ value required for a read for it to be considered.")

    (when (eql key-value :all)
      (format stream "~%External program options:~%"))
    (print-info :java-bin "Full path of java executable, can use `which java`")
    (print-info :java-bin-args "Arguments for java placed before the -jar argument.")
    (print-info :fs-align-jar "Full filename with path of SamToFlowSpace.jar (defaults to the same location as ion-variant-hunter-core)")
    (print-info :fs-align-jar-args "Arguments for the flow space aligner jar package.")
    (print-info :fs-align-range "RANGE option to specify genomic region, i.e. chr4:1000-2000.")
    (print-info :fs-align-num-threads "Number of threads to give the java flow space aligner.")

    (print-info :python-bin "Full path of python executable, can use `which python`.")
    (print-info :python-py "Full filename with path of samRegionOverlap.py (defaults to the same location as ion-variant-hunter-core)")
    (print-info :samtools-bin "Full path of samtools executable used to get contig information for the VCF header, can use `which samtools`.")

    (when (eql key-value :all)
      (format stream "~%Debug options:~%"))
    (print-info :print-flow-space-aligns? "Print out all flow space alignments to a .align file.  Default is off.")
    (print-info :aligner-method          "Options: ,:java, ,:lisp, ,:auto. Auto defaults to SamToFlowSpace.jar for bam and internal method for sam.")

    (when (eql key-value :all)
      (format stream "~%"))
    ))
