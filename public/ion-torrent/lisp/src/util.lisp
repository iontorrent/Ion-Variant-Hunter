;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

;;;
;;; Date/Time
;;;
(defun format-current-date (&optional format)
  (multiple-value-bind (sec min hr day month year)
      (get-decoded-time)
    sec min hr
    (let ((format-string 
	   (cond ((eql format :no-dividers)
		  "~4,'0d~2,'0d~2,'0d")
		 (t
		  "~4,'0d-~2,'0d-~2,'0d"))))
      (format nil format-string year month day))))

(defun format-current-time (&optional format)
  (multiple-value-bind (sec min hr day month year)
      (get-decoded-time)
    day month year
    (let ((format-string 
	   (cond ((eql format :no-dividers)
		  "~2,'0d~2,'0d~2,'0d")
		 (t
		  "~2,'0d:~2,'0d:~2,'0d"))))
      (format nil format-string hr min sec))))

;; Number format
(defun format-scientific (number &optional (sig-figures 3))
  ;; Prints out up to 15 sig figs.
  (unless (numberp number)
    (return-from format-scientific "NIL"))
  (setq number (coerce number 'double-float))
  (let ((control-string (format nil "~~,~ae" sig-figures)))
    (substitute #\e #\d
		(format nil control-string number))))

(defun format-scientific-list (list  &optional (sig-figures 3))
  (mapcar #'(lambda (number)
	      (format-scientific number sig-figures))
	  list))

;; Lists
(defun limit-size-list (lst max-size)
  (subseq lst 0 (min max-size (length lst))))

(defun limit-size-lists-in-list (lists max-size)
  (mapcar #'(lambda (list)
	      (limit-size-list list max-size))
	  lists))
	      

(defun split-list-round-robbin (list num-parts)
  (setq num-parts (min (length list) num-parts))
  (let ((lists-array (make-array num-parts :initial-element nil))
	(counter 0))
    (dolist (item list)
      (push item
	    (aref lists-array (mod counter num-parts)))
      (incf counter))
    (dotimes (i num-parts)
      (setf (aref lists-array i)
	    (reverse (aref lists-array i))))
    (mapcar #'reverse lists-array)))

(defun split-list-contiguous (list num-parts)
  (setq num-parts (min (length list) num-parts))
  (let ((size-each-part (floor (/ (length list) num-parts)))
	segments)
    (dotimes (part-num num-parts)
      (push
       (subseq list (* part-num size-each-part)
	       (unless (eql part-num (- num-parts 1))
		 (* (1+ part-num) size-each-part)))
       segments))
    segments))

(defun extract-from-sorted-list (item sequence &key from-end key test)
  "What elements in the sequence where the item is greater (test #'>) or lesser (test #'<)"
  (let* ((position (position item sequence :from-end from-end :test test :key key)))
    (when sequence
      (unless (funcall test item (if key (funcall key (car sequence)) (car sequence)))
	(return-from extract-from-sorted-list
	  (values nil sequence))))
    (if (numberp position)
	(values
	 (subseq sequence 0 (1+ position))
	 (subseq sequence (1+ position)))
	sequence)))
;; Example:
;; (extract-from-sorted-list 24 '(1 1 1 2 4 5 7 9 23 24 24 24 24 25 25 26 27 28 55 59) :from-end t :test #'>)
;; (1 1 1 2 4 5 7 9 23)
;; (24 24 24 24 25 25 26 27 28 55 59)

;; sequences
(defun safe-subseq (sequence start &optional end)
  (let (was-filled?)
    (when (< start 0)
      (let ((start-diff (- start)))
	(setq was-filled? t)
	(format t "Warning filling in ~a '?' at the begining of ~a~%" start-diff sequence)
	(setq sequence (format nil "~{~a~^~}~a"
			       (make-list start-diff :initial-element #\?)
			       sequence))
	(setq start 0)
	(when end
	  (setq end (+ start-diff end)))))
    (when (and end
	       (> end (length sequence)))
      (let ((end-diff (- end (length sequence))))
	(setq was-filled? t)
	(format t "Warning filling in ~a '?' at the end of ~a~%" end-diff sequence)
	(setq sequence (format nil "~a~{~a~^~}"
			       sequence
			       (make-list end-diff :initial-element #\?)))
	(setq end nil) ;; just go to the end
	))
    (values
     (subseq sequence start end)
     was-filled?)))

(defun replace-end-?s-w-subst (seq-with-?s seq-without)
  (when (eql (length seq-with-?s) (length seq-without))
    (let ((first-pos-without-? (position #\? seq-with-?s :test-not #'eql))
	  (last-pos-without-? (position #\? seq-with-?s :test-not #'eql :from-end t)))
      (unless first-pos-without-?
	(return-from replace-end-?s-w-subst seq-without))
      ;;(format t "~a ~a~%" first-pos-without-? last-pos-without-?)
      (dotimes (x first-pos-without-?)
	(setf (char seq-with-?s x) (char seq-without x)))
      (dotimes (x (- (length seq-without)
		     last-pos-without-? 1))
	(let ((pos (- (length seq-without) x 1)))
	  (setf (char seq-with-?s pos) (char seq-without pos))))
      ))
  seq-with-?s)

(defun make-hash-table-from-list (key-items &key (test #'eql)
				  (value t))
  (let ((ht (make-hash-table :test test :size (length key-items))))
    (dolist (key-item key-items)
      (setf (gethash key-item ht) value))
    ht))

(defun set-difference-size (list1 list2 &key (test #'eql))
  ;;set difference size of list1 - list2
  (let ((ht (make-hash-table-from-list list1 :test test))
	(num-list2-items-in-list1 0))
    (dolist (l2-item list2)
      (when (gethash l2-item ht)
	(setf (gethash l2-item ht) nil) ;; found once already
	(incf num-list2-items-in-list1)))
    (- (hash-table-count ht) num-list2-items-in-list1)))
;;;
;;; CLOS extensions
;;;

;; For a given class, returns for each slot:
;; 1. the method to read that slot,
;; 2. the symbol of each slot, and
;; 3. the initarg of each slot.
(defgeneric get-class-slot-properties (class))

(defmethod get-class-slot-properties ((class-symbol symbol))
  (get-class-slot-properties (find-class class-symbol)))

(defmethod get-class-slot-properties ((class class))
  (values
;;   (mapcar #'sb-pcl::slot-definition-reader-function
;;	   (sb-pcl::class-slots class))
   (mapcar #'sb-pcl::slot-definition-name
	   (sb-pcl::class-slots class))
   (mapcar #'car
	   (mapcar #'sb-pcl::slot-definition-initargs
		   (sb-pcl::class-slots class)))))

;; For a list of objects, returns as a list the name followed by the
;; list of values for each slot.
;; If slots is nil, then just return all values
(defgeneric get-object-values-for-slot-list (objs &optional slot-list))
(defmethod get-object-values-for-slot-list ((objects cons)
					    &optional slot-list)
  (multiple-value-bind (names)
      (get-class-slot-properties (class-of (car objects)))
    (setq slot-list
	  (or slot-list names))
    (remove
     nil
     (mapcar #'(lambda (name)
		 (when (find name slot-list)
		   (cons
		    name
		    (mapcar #'(lambda (object)
				(sb-pcl::slot-value-or-default object name))
			       ;;(funcall reader object)
			    objects))))
	     names))))

;;;
;;; Multithreading
;;;
;; http://www.sbcl.org/manual/#Threading

(defun multi-thread-mutating-function (num-threads function first-argument second-list-argument
				       &rest rest-arguments)
  ;;The second argument must be in the form of a list of objects
  ;;that will be modified by the function in some way (process)

  (setq num-threads (min num-threads (length second-list-argument)))
  (when (eql num-threads 1)
    ;; 1 thread, so don't multithread!
    (return-from multi-thread-mutating-function
      (apply function first-argument second-list-argument rest-arguments)))

  (let (second-arg-segments
	threads)
      
  (setq second-arg-segments (split-list-contiguous second-list-argument num-threads))

  (dolist (second-arg-segment second-arg-segments)
    (push
     (sb-thread:make-thread (lambda ()
			      (apply function first-argument second-arg-segment
				     rest-arguments)))
     threads))
  (format t "~a threads running.~%" (length threads))
  (dolist (thread (reverse threads))
    ;; TODO, check to see if error checking needs to be done.
    (sb-thread:join-thread thread))
  (format t "All threads completed.~%")))