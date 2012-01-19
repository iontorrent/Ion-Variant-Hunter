;; Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

(in-package :cl-user)

(defun median (lst)
  (let ((num-items (length lst))
        (s-lst (sort lst #'>)))
    (if (evenp num-items)
        (/ (+ (nth (/ num-items 2) s-lst)
              (nth (- (/ num-items 2) 1) s-lst)) 2)
        (nth (floor (/ num-items 2)) s-lst))))

(defun average (lst)
  (let ((sum 0))
    (dolist (item lst)
      (incf sum item))
    (when lst
      (/ sum (length lst)))))

(defun euclidean-distance (lst1 &optional lst2)
  ;; if lst1 is already the diff, or just what euclidena magnitude, lst2 is optional
  (setq lst2 (or lst2 (make-list (length lst1) :initial-element 0)))
  (let ((sum-of-squares 0d0))
    (mapc #'(lambda (item1 item2)
	      (incf sum-of-squares
		    (expt (- item1 item2) 2))) ;; ensure double precision
	  lst1 lst2)
    (sqrt sum-of-squares)))

;; Math
(defun ln(n)
  (log n))

(defun root-mean-square (list)
  (when list
    (let ((sum-squares 0))
      (dolist (x list)
	(incf sum-squares (expt x 2)))
      (sqrt (+ 0d0 (/ sum-squares (length list)))))))

;; Combinatorics

;; Factorial/Gamma approximation functions
;; Lanczos approximation

;; Regular factorial
(defun factorial (n &optional (accum 1))
  ;;(declare (type fixnum n accum))
  (if (zerop n)
      accum
      (factorial (- n 1) (* accum n))))

;; Gosper approximation
(defun ln-factorial (n)
  (+ (* 1/2 (ln (* pi (+ (* 2 n) 1/3))))
     (* n (ln n))
     (- n)))

;; Auto select between the above functions
(defun auto-ln-factorial (n)
  ;; For low n, just use regular factorial function
  ;; For higher ones, use Gosper approximation
  (if (<= n 4)
      (ln (* 1d0 (factorial n)))
      (ln-factorial n)))

;; Choose function
(defun ln-n-choose-k (n k)
  (- (auto-ln-factorial n)
	  (auto-ln-factorial k)
	  (auto-ln-factorial (- n k))))

(defun n-choose-k (n k)
  (exp (ln-n-choose-k n k)))

;; Binomial distribution
(defun binomial-probability (k-successes n-trials probability)
  (when (or (< n-trials k-successes)
	    (> probability 1))
    (error "Cannot find binomial probability with k=~a, n=~a, and p=~a~%"
	   k-successes n-trials probability))
  (exp
   (+ (ln-n-choose-k n-trials k-successes)
      (* k-successes (ln probability))
      (* (- n-trials k-successes) (ln (- 1 probability))))))
;; One test case: (BINOMIAL-PROBABILITY 777 1569 1/2)
;; Here, without the ln representation of the last two terms, an arithmetic
;;  error FLOATING-POINT-OVERFLOW occurs.

(defun cumalative-binomial-probability (k-successes n-trials probability)
  (let ((cumalative-prob 0)
	(k-min (min k-successes (- n-trials k-successes))))
    (dotimes (k (+ k-min 1))
      (incf cumalative-prob
	    (+
	     (binomial-probability k n-trials probability)
	     (binomial-probability k n-trials (- 1 probability))))
      )
    (when (eql 1/2 (/ k-successes n-trials))
      (decf cumalative-prob (binomial-probability k-successes
						  n-trials probability)))
    (if (> cumalative-prob 1)
	1
	cumalative-prob)))
;; This is not the best way of doing this, as it is slow and errors
;; in the calculation can accumulate, i.e. if the >1 check wasn't there. . .
;; (cumalative-binomial-probability 750 1500 0.50d0) => 1.0004612624184774d0

;; Also this can be estimated using gaussian distribution.
