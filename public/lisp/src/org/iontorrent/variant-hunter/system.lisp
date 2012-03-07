;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

(in-package :cl-user)

;;;;;;;;;;;;;;;;
;; set and load source files
;;;;;;;;;;;;;;;;

(defvar *src-files*)
(setq *src-files* (list  "util.lisp" "stats.lisp" "parse-util.lisp" "reference.lisp" "sam-parse.lisp" "flow-space.lisp" "seq-deviations.lisp" "vcf-writer.lisp" "alignment-streamer.lisp" "ion-variant-hunter.lisp"))

(defvar *vh-version* "0.1.3.beta")
(dolist (src-file *src-files*)
  (load src-file))

(defvar *build-date* "")
(defvar *build-time* "")

(defvar *from-built-exec* nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Evaluation loop (for interactive sessions)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defun splash ()
  (format t "~%Ion Variant Hunter v~a~%" *vh-version*)
  (format t "Binary: ~a~%" (car *posix-argv*))
  (format t "Build date: ~a~%Build time: ~a~%" *build-date* *build-time*)
)

;; Evaluation loop
(defun repl-prompt-ion-variant-hunter (stream)
  (fresh-line stream)
  (write-string "Ion Variant Hunter >> " stream))

(defun indel-read-eval-loop ()
  (format t "~%Welcome to the Variant Hunter interactive session.~%Type (command-line-help) to see command line options.~%Type '(quit)' below to quit.~%")
  (setq sb-impl::*repl-prompt-fun* #'repl-prompt-ion-variant-hunter)
  (sb-impl::toplevel-init) ;; processes command line options and then enters interactive loop
  ;;(sb-impl::toplevel-repl nil)
  0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Shell command line option parsing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun read-command-line-value (value)
  (if (and (> (length value) 1)
	   (eql #\, (char value 0)))
      (with-input-from-string (stream (subseq value 1))
	(read stream))
      value))

(defun process-command-line-options ()
  (let ((options (rest *posix-argv*))
	(parsed-options nil)
	(disable-debugger t)
	(binary-fn (car *posix-argv*)))
    (loop while options do
	 (let ((option (first options)))
           (flet ((pop-option ()
                    (if options
                        (pop options)
                        (sb-impl::startup-error
                         "unexpected end of command line options"))))
	     (cond ((eql (char option 0) #\-)
		    (let ((option (pop-option))
			  (value (when options
				   (unless (eql #\- (char (car options) 0))
				     (pop-option)))))
		      (setq value (or (read-command-line-value value) :no-arg-value))
		      (push (list option value) parsed-options)))
		   (t
		    (format t "*********~%WARNING '~a' on command line is ignored.~%*********~%"
			    (pop-option)))
		   ))))
    (when disable-debugger
      ;;(sb-ext:disable-debugger)
      )
    (when parsed-options
      (setq parsed-options (nreverse parsed-options))
      (format t "~%Program executed with:~%")
      (dolist (parsed-option parsed-options)
	(format t "~a ~a~a~%"
		(car parsed-option)
		(if (stringp (second parsed-option)) "" ",")
		(second parsed-option))))
    (values
     parsed-options
     binary-fn)))

(defun print-argument-list-for-function (function &optional (stream t))
  (let ((arg-list (sb-kernel:%simple-fun-arglist function)) ;;(sb-introspect:function-lambda-list  function))
	arg-default-relist
	(max-arg-string-length 0))
    (when (eql '&key (car arg-list))
      (setq arg-list (cdr arg-list)))

    (dolist (arg arg-list)
      (let (cmd-line-arg
	    default-value)
	(if (consp arg)
	    (setq cmd-line-arg (car arg)
		  default-value (second arg))
	    (setq cmd-line-arg arg))
	(setq cmd-line-arg (string-downcase (format nil "~a" cmd-line-arg)))
	(when (> (length cmd-line-arg) max-arg-string-length)
	  (setq max-arg-string-length (length cmd-line-arg)))
	(push (list cmd-line-arg default-value) arg-default-relist)))

    (format stream "Avaliable command line options and default values:~%")
    (format stream "For descriptions of these options use --detailed-help.~%")
    (dolist (newarg (reverse arg-default-relist))
      (let ((cmd-line-arg (car newarg))
	    (default-value (second newarg)))
	(format stream "--~a" cmd-line-arg)
	(format stream " {value} ")
	(format stream "~{~a~^~}[" (make-list (- max-arg-string-length (length (car newarg))) :initial-element " "))
	(unless (stringp default-value)
	  (format stream ","))
	(when (keywordp default-value)
	  (format stream ":"))
	(format stream "~a]~%" default-value)))
    (format stream "~%")))

(defun info (func-string)
  (print-argument-list-for-function (symbol-function (intern (string-upcase func-string)))))

(defun command-line-help (&optional (binary (car *posix-argv*)))
  (format t "~%~a help:~%~%" binary)
  (format t "Options:~%")
  (format t "  --interactive              Enter into an interactive session.~%")
  (format t "  --dynamic-space-size <MiB> Size of reserved dynamic space in megabytes.~%")
  (format t "  --help                     SBCL help page.~%")
  (format t "  -h/-?                      This help page.~%")
  (format t "  --detailed-help            Description of each command line argument.~%")
  (format t "  --default-value-help       Shows default values for each argument.~%")
  (format t "  --function                 Function to run.~%")
  (format t "~%")
  ;;(print-argument-list-for-function #'write-alignment-n-deviation-lists-from-sam)
  (write-alignment-n-deviation-lists-from-sam-usage)
  (format t "~%"))

(defun find-option (options search-options)
  (let (option-found?)
    (dolist (option search-options)
      (setq option-found?
	    (or option-found? 
		(find (format nil "-~a" option) options :key #'car :test #'equal)
		(find (format nil "--~a" option) options :key #'car :test #'equal)) ))
    option-found?))

(defvar *exclude-argument-list* '("dynamic-space-size"))

(defun strip-beginning-dashes (argument)
  (while (eql #\- (char argument 0))
    (setq argument (subseq argument 1)))
  argument)

(defun get-function-n-arguments-from-options (posix-options)
  (let* ((function-name (second (find-option posix-options '("function"))))
	 (exclude-argument-list (push "function"
				      *exclude-argument-list*))
	 key-value-list
	 arg-list)
    (when function-name
      ;;(format t "'~a'~%~a" function-name posix-options)
      (dolist (option posix-options)
	(setf (car option) (strip-beginning-dashes (car option)))
	;;(format t "find ~a:~a~%" option exclude-argument-list)
	(unless (find (car option) exclude-argument-list :test #'equal)
	  (unless (string= "ARG" (string-upcase (car option)))
	    (push (intern (string-upcase (car option)) :keyword) key-value-list)
	    (push (second option) key-value-list))
	  (when (string= "ARG" (string-upcase (car option)))
	    (push (second option) arg-list))
	  )))

    (values
     (intern (string-upcase function-name))
     (reverse arg-list)
     (reverse key-value-list))))

(defun apply-function-from-posix-options (posix-options &optional state-list)
  (let (func-return-value)
    (multiple-value-bind (func arg-list key-value-list)
	(get-function-n-arguments-from-options posix-options)
      (when (car state-list)
	(setf (car state-list) :posix-options-parsed))
      ;;(format t "func = ~a~%" func)
      (setq func-return-value (apply func (append arg-list key-value-list)))
      (when (car state-list)
	(setq func-return-value
	      (setf (car state-list) :function-completed))))
    func-return-value
    ))

;; To enable profiling, uncomment:
;; (require :sb-sprof)
;; then wrap the apply function below with the this:
;; (sb-sprof:with-profiling (:mode :alloc :loop nil :show-progress t :max-samples 100 :report :flat)
;; This is memory profiling.  For CPU, just remove :mode :alloc

(defun act-on-command-line-options (posix-options binary)
  (let (state)
    (sb-ext:disable-debugger)
    ;;(format t "posix-options: ~{~a~^,~}~%find-option: ~a~%"  posix-options (find-option posix-options '("function")))

    (unwind-protect
	 (progn
	   (when posix-options
	     (when (find-option posix-options '("detailed-help"))
	       (write-alignment-n-deviation-lists-from-sam-doc :all)
	       (setq state :help)
	       (quit)
	       )
	     (when (find-option posix-options '("default-value-help"))
	       (print-argument-list-for-function #'write-alignment-n-deviation-lists-from-sam)
	       (setq state :help)
	       (quit)
	       )
	     (when (find-option posix-options '("?" "h" "help"))
	       (command-line-help binary)
	       (setq state :help)
	       (quit)
	       )
	     (when (find-option posix-options '("interactive"))
	       (setq state :interactive))
	     (when (find-option posix-options '("function"))
	       (let (return-value)
		 (multiple-value-bind (func arg-list key-value-list)
		     (get-function-n-arguments-from-options posix-options)
		   (setq state :posix-options-parsed)
		   (setq return-value (apply func (append arg-list key-value-list)))
		   (setq state :function-completed)
		   (format t "return value of ~a: ~a~%" func return-value))))
	     )
	   (setq state (or state :help)))
      (let (has-errors?)
	(when (eql state :posix-options-parsed)
	  (format t "Error in running function.  Terminating.~%")
	  (setq has-errors? t))
	(unless state
	  (format t "Error detected in parsing command line arguments.  Terminating.~%")
	  (setq has-errors? t)
	)
	(when has-errors?
	  (quit :unix-status 1))))
    ;;(format t "state detected: ~a~%" state)
    (case state
      (:help
       (command-line-help binary)
       (when *from-built-exec*
	 (quit :unix-status 0)))
      (:function-completed
       (format t "Function completed.  Ending.~%")
       (when *from-built-exec*
	 (quit :unix-status 0)))
      (:interactive
       (format t "Entering interactive mode.~%"))
      )
    (sb-ext:enable-debugger)))

;;;;;;;;;;;;;;;;
;;; Build call
;;;;;;;;;;;;;;;;
;; Default function called on startup as defined by the :toplevel key below.
(defun splash-and-eval-loop ()
  (splash)
  (multiple-value-bind (options binary)
      (process-command-line-options)
    (act-on-command-line-options options binary))
    
  (indel-read-eval-loop))

;; Actual build command
(defun build (&optional exec-fn)
  (when (eql (length exec-fn) 0)
    (setq exec-fn "ion-variant-hunter-core"))

  (setq *print-length* 25)  ;; limit debug output 

  (setq *build-date* (format-current-date))
  (setq *build-time* (format-current-time))

  (format t "~%")
  (format t "***************************************")
  (splash)
  (format t "Now building '~a'.~%" exec-fn)
  (format t "***************************************~%")
  (format t "~%")
  (sb-impl::flush-standard-output-streams)

  (setq *from-built-exec* t)
  (sb-ext:save-lisp-and-die exec-fn :executable t
			    :toplevel #'splash-and-eval-loop))
