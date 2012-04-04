;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

(in-package :cl-user)

(defmacro while (test &rest body)
  `(do ()
    ((not ,test))
    ,@body))

(defun string-concat ( a &rest args)
  (apply #'concatenate (cons 'string
			     (cons a args))))

(defun determine-number-header-lines (filespec &key
				      (blank-as-headers? t)
				      (header-char #\#))
  (let ((counter 0))
    (with-open-file (stream filespec :direction :input)
      (loop for line = (read-line stream nil)
	    while (and line
		       (if (eql 0 (length line))
			   blank-as-headers?
			   (eql (char line 0) header-char)))
	    do (incf counter)))
    counter))

(defun parse-string-old (string delim &optional check-quotes? parsed-items)
  (let ((delim-pos (position delim string))
	parsed-string
	rest-string)
    (if delim-pos
	(progn
	  (setq parsed-string (subseq string 0 delim-pos))
	  (setq rest-string (subseq string (1+ delim-pos)))
	  (parse-string-old rest-string delim check-quotes?
			    (cons parsed-string parsed-items))
	  )
	(progn
	  (reverse (cons string parsed-items))))))

(defun parse-string (string delim &optional check-quotes?)
  ;; This has only one subseq per part instead of two in the above function
  check-quotes? ;; ignored, but not implemented in old function either
  (let (delim-positions
	(start-pos 0))		
    (dotimes (pos (length string))
      (when (eql (char string pos) delim)
	(push pos delim-positions)))
    (push nil delim-positions)
    (mapcar #'(lambda (delim-pos)
		(let ((this-part
		       (subseq string start-pos delim-pos)))
		  (when delim-pos
		    (setq start-pos (1+ delim-pos))) ; next start pos
		  this-part
		))
	    (reverse delim-positions))))

(defmacro with-delimited-file ((line filespec &optional
				     skip-lines
				     (delimiter #\Tab)
				     check-quotes?)
			       &body forms)
  (let ((stream (gensym))
	(in (gensym))
	(x (gensym))
	(break? (gensym)))
    `(with-open-file (,stream ,filespec :direction :input)
      (let (,break?)
	(labels ((get-next-line (,in)
		   (let ((line (read-line ,in nil :eof)))
		     (if (eql line :eof)
			 :eof
			 (parse-string line ,delimiter ,check-quotes?))))
		 (break-loop () (setq ,break? t)))
	  (break-loop) ;; just to not have function deleted
	  (setq ,break? nil)
	  (when ,skip-lines
	    (dotimes (,x ,skip-lines)
	      (get-next-line ,stream)))

	  (labels ((get-another-line ()
		     (get-next-line ,stream)))
	    (do ((,line (get-another-line)
			(get-another-line)))
		((or (eql ,line :eof) ,break?))
	      ,@forms)))))))

;; process stream
(defun print-process-info (process &optional (stream t))
  (let ((process-infos '(sb-ext:process-pid
			 ;;sb-ext:process-input
			 ;;sb-ext:process-output
			 sb-ext:process-error
			 sb-ext:process-alive-p
			 sb-ext:process-status
			 sb-ext:process-exit-code
			 sb-ext:process-core-dumped)))
    (mapc #'(lambda (info)
	      (format stream "(~a ~a) => ~a~%" info process
		      (funcall info process)))
	  process-infos)))

(defun read-line-of-process (process process-stdout)
  (let ((exit-code (process-exit-code process)))
    (if (or (not exit-code)     ;; still running
	    (eql exit-code 0))  ;; finished by there may be lines left in cache
	(read-line process-stdout nil :eof)
	:terminated)))

;; Used this before, but could be problematic b/c of leftover items in cache.
;;(if (process-alive-p process)
;;(read-line java-stdout nil :eof))
;;    :terminated))

;; argument parsing
(defun argument-underscore-to-dash (argument-string)
  ;; replaces begining _ (underscore) with - (dash)
  ;; if actual underscore is desired, then two conseq. underscores at the begining will be replaced by one
  ;; "Begining" is the first character of any word (parsed by spaces).
  (when argument-string
    (let ((arg-parts (parse-string argument-string #\Space))
	  new-args)
      (setq new-args
	    (mapcar #'(lambda (arg-part)
			(cond (;; __ (double underscore), means literal starting single underscore
			       (and (> (length arg-part) 1) 
				    (eql #\_ (char arg-part 0))
				    (eql #\_ (char arg-part 1)))
			       (subseq arg-part 1))
			      (;; else starting underscore means dash
			       (and (> (length arg-part) 1)
				    (eql #\_ (char arg-part 0)))
			       (format nil "-~a" (subseq arg-part 1)))
			      (;; otherwise, unmodified
			       arg-part)))
		    arg-parts))
      (format nil "~{~a~^ ~}" new-args))))

(defun parse-string-to-intern-list (argument-string &optional (separator #\,)
				    (intern-type :keyword))
  (mapcar #'(lambda (arg-item)
	      (intern
	       arg-item
	       intern-type))
	  (parse-string (string-upcase argument-string) separator)))

;; filename utils
(defun split-filename-base-and-extension (filename)
  (let ((last-period-pos
	 (position #\. filename :from-end t)))
    (values
     (subseq filename 0 last-period-pos)
     (subseq filename last-period-pos))))

(defun get-exec-pathname ()
  (pathname (car *posix-argv*)))

(defun get-exec-directory ()
  (directory-namestring (get-exec-pathname)))

(defun get-current-directory ()
  (directory-namestring *default-pathname-defaults*))

(defun ensure-full-directory (pathspec)
  (if (eql #\\ (char pathspec 0))
      pathspec 
      (string-concat (get-current-directory) pathspec)))

(defun get-default-paths (&optional (platform :unix))
  (case platform
    (:unix
     (parse-string (sb-unix::posix-getenv "PATH") #\:))))

(defun find-binary-in-path (binary)
  (unless (find #\/ binary)
    (dolist (path (get-default-paths))
      (let ((binary-w-path (format nil "~a/~a" path binary)))
	(when (directory binary-w-path)
	  (return-from find-binary-in-path binary-w-path))))
    nil))

(defun bam-n-bai-present? (bam-file)
  (and (directory bam-file)
       (directory (format nil "~a.bai" bam-file))))

(defun check-bam-n-bai (bam-file)
  (unless (stringp bam-file)
    (error (format nil "ERROR: Specify BAM file directly without any comma (,).")))
  (unless (bam-n-bai-present? bam-file)
    (error (format nil "ERROR: BAM and BAM index (.bai) files are both required.  BAM given: '~a'~%"
		   bam-file)))
  t)

(defun check-bam-n-bai-non-fatal (bam-file)
  (unless (bam-n-bai-present? bam-file)
    (format *error-output* "CRITICAL WARNING: BAM and/or BAM index (.bai) files are missing.  BAM given: '~a'~%"
	    bam-file))
  t)

(defvar *starting-current-directory*
  (unless (boundp '*build-date*)
    (get-current-directory)))

(defun full-path-at-exec-directory (filename)
  ;; get a file that's in the same directory as the built version
  (let ((path (if (boundp '*build-date*)
		  ;; built version
		  (get-exec-directory)
		   ;; dev version
		  *starting-current-directory*)))
    (string-concat path filename)))

(defun find-binary-in-path-or-exec-directory (filename)
  (or 
   (unless (eql 0 (length (get-exec-directory)))
     (full-path-at-exec-directory filename))
   (find-binary-in-path filename)
   (full-path-at-exec-directory filename)
   ))

;; General util functions
(defun get-random-item (list)
  (nth (random (length list)) list))

(defun get-random-hash-item (hash)
  (let (hash-items)
    (maphash #'(lambda (key val)
		 (push (cons key val) hash-items))
	     hash)
    (get-random-item hash-items)))


(defun round-half-up (number &optional (divisor 1))
  (multiple-value-bind (val rem)
      (round number divisor)
    (if (eql 1/2 (/ rem divisor))
	(values (1+ val) (* -1/2 divisor))
	(values val rem))))

(defun round-half-down (number &optional (divisor 1))
  (multiple-value-bind (val rem)
      (round number divisor)
    (if (eql -1/2 (/ rem divisor))
	(values (1- val) (* 1/2 divisor))
	(values val rem))))
