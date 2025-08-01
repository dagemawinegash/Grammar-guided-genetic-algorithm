; atom space creation
!(bind! &grammar (new-space))

;;;;;;;;;; grammar definition ;;;;;;;;
!(add-atom &grammar (start expr 0))  
!(add-atom &grammar (expr (bop expr expr) 0))
!(add-atom &grammar (expr (uop expr) 1))
!(add-atom &grammar (expr term 2))
!(add-atom &grammar (bop AND 0))
!(add-atom &grammar (bop OR 1))
!(add-atom &grammar (uop NOT 0))
!(add-atom &grammar (term A 0))
!(add-atom &grammar (term B 1))
!(add-atom &grammar (term C 2))


; general variable binds
!(bind! TERMINALS (A B C AND OR NOT))
!(bind! POP_SIZE 5)
!(bind! GENES 5)
!(bind! CROSSOVER_RATE 0.9)
!(bind! MUTATION_RATE 0.9)


;;;;;; truth table definition ;;;;;;
!(bind! &truth_table_data (
    ((False False False) False)
    ((False False True) True)
    ((False True False) False)
    ((False True True) True)
    ((True False False) False)
    ((True False True) True)
    ((True True False) False)
    ((True True True) True)
))

;;;;;; general helper functions ;;;;;;

; append function helper
(= (append $list $ele) 
    (union-atom $list ($ele)))

; gets the Nth element of an expression (0-indexed)
(= (nth-elem $list 0) 
    (car-atom $list))

(= (nth-elem $list $x) 
    (nth-elem (cdr-atom $list) (- $x 1)))

;;;;;; genotype creation ;;;;;;
(= (initialize-genotype) 
    (pop-helper 0 ()))

(= (make-random $i $acc)
    (if (< $i GENES) 
        (make-random (+ $i 1) 
            (append $acc (random-int &rng 0 10))) 
        $acc))

(= (pop-helper $i $acc)
    (if (< $i POP_SIZE) 
        (pop-helper (+ $i 1) 
            (append $acc (make-random 0 ()))) 
        $acc))  

; !(initialize-genotype)


;;;;;genotype to phenotype conversion ;;;;;;

(= (is-terminal $element)
    (if (== ()
            (collapse
                (let $item (superpose TERMINALS)
                    (if-equal $item $element $element (empty))
                )
            )
        )
        False
        True
    )
)
;!(is-terminal A)

(= (is-expression $item) 
    (let $type (get-metatype $item)
        (if (== $type Expression)
            True
            False
        )
    ))
; !(is-expression (bop))

(= (generate-phenotype $genotype $expand $max_depth $current_depth $current_idx)
    (if (>= $current_depth $max_depth)
        ()
        (let* (
            ($gene-idx (% $current_idx GENES))
            ($gene-val (nth-elem $genotype $gene-idx))
            ($all-options (collapse (match &grammar ($expand $nxt $_) $nxt)))
            ($size-options (size-atom $all-options))
            ($cur-choice (% $gene-val $size-options))
            ($next (match &grammar ($expand $nxt $cur-choice) $nxt))
            ($is-terminal (is-terminal $next))
            ($is-expression (is-expression $next))
        )
        (if $is-terminal
            $next
            (if $is-expression
                (ifexpr $next $genotype $max_depth (+ $current_depth 1) (+ $current_idx 1))
                (generate-phenotype $genotype $next $max_depth (+ $current_depth 1) (+ $current_idx 1))
            )
        )
        ))
)

(= (ifexpr $expr $genotype $max_depth $current_depth $current_idx) 
    (if-decons-expr $expr $head $tail 
        (let*(
            ($cur (generate-phenotype $genotype $head $max_depth $current_depth $current_idx))
            ($nxt (ifexpr $tail $genotype $max_depth $current_depth (+ $current_idx 1)))
        ) 
        (cons-atom $cur $nxt))
        ()
    )
)

!(generate-phenotype (2 0 0 1 2) start 10 0 0); expect [(AND (NOT A) C)]
!(generate-phenotype (0 2 0) start 10 0 0) ; expect A
!(generate-phenotype (0 1 0 2 2) start 10 0 0) ; expect (NOT C)


;;;;;; phenotype evaluation function ;;;;;;

(= (evaluate-expression $expr $A $B $C)
    (evaluate-helper $expr $A $B $C))

(= (evaluate-helper $expr $A $B $C)
    (if (== $expr A)
        $A
        (if (== $expr B)
            $B
            (if (== $expr C)
                $C
                (let* (
                    ($is_expression (is-expression $expr))
                )
                    (if $is_expression
                        (evaluate-expr $expr $A $B $C)
                        $expr
                    )
                )
            )
        )
    ))

(= (evaluate-expr $expr $A $B $C)
    (if-decons-expr $expr $op $args
        (evaluate-operator $op $args $A $B $C)
        $expr
    ))

(= (evaluate-operator $op $args $A $B $C)
    (if (== $op AND)
        (if-decons-expr $args $first $rest
            (if-decons-expr $rest $second $remaining
                (let* (
                    ($first-val (evaluate-helper $first $A $B $C))
                    ($second-val (evaluate-helper $second $A $B $C))
                )
                    (and $first-val $second-val)
                )
                (invalid)
            )
            (invalid)
        )
        (if (== $op OR)
            (if-decons-expr $args $first $rest
                (if-decons-expr $rest $second $remaining
                    (let* (
                        ($first-val (evaluate-helper $first $A $B $C))
                        ($second-val (evaluate-helper $second $A $B $C))
                    )
                        (or $first-val $second-val)
                    )
                    (invalid)
                )
                (invalid)
            )
            (if (== $op NOT)
                (if-decons-expr $args $operand $rest
                    (if (== $rest ())
                        (let* (
                            ($operand-val (evaluate-helper $operand $A $B $C))
                        )
                            (not $operand-val)
                        )
                        (invalid)
                    )
                    (invalid)
                )
                (invalid)
            )
        )
    )
)

; !(evaluate-expression (AND (NOT A) (OR B C)) True False True) ; expected: False

;;;;;; phenotype fitness evaluation ;;;;;;

(= (fitness-evaluator $phenotype)
    (if (== (evaluate-expression $phenotype True False True) (invalid))
        0
        (fitness-helper $phenotype &truth_table_data 0 0)
    )
)

(= (fitness-helper $phenotype $truth-data $index $correct-count)
    (if (>= $index 8)
        $correct-count
        (let* (
            ($row (nth-elem $truth-data $index))
            ($inputs (car-atom $row))
            ($expected-list (cdr-atom $row))
            ($expected (car-atom $expected-list))
            ($A (nth-elem $inputs 0))
            ($B (nth-elem $inputs 1))
            ($C (nth-elem $inputs 2))
            ($actual (evaluate-expression $phenotype $A $B $C))
            ($is-correct (== $actual $expected))
            ($new-count (if $is-correct (+ $correct-count 1) $correct-count))
        )
            (fitness-helper $phenotype $truth-data (+ $index 1) $new-count)
        )
    )
)

;!(fitness-evaluator (AND (NOT A) C)) ; expected: 6/8
;!(fitness-evaluator (OR A C)) ; expected: 6/8
;!(fitness-evaluator C) ; expected: 8/8
;!(fitness-evaluator (AND A C)) ; expected: 6/8
;!(fitness-evaluator (INVALID A C)) ; expected: 0
;!(fitness-evaluator (AND A)) ; expected: 0
;!(fitness-evaluator (A OR B)) ; expected: 0


;;;;;;;;;;;;;;;;; selection functions ;;;;;;;;;;;;;;

;;;;; roulette stochastic selection function ;;;;;;

(= (roulette-helper $max_fitness $population $fitnesses)
    (let* (
        ($random-num (random-float &rng 0 1))
        ($idx (random-int &rng 0 (size-atom $population)))
        ($random-fitness (nth-elem $fitnesses $idx)))   
        (if (< $random-num (/ $random-fitness $max_fitness))
            (nth-elem $population $idx)
            (roulette-helper $max_fitness $population $fitnesses)))
)

(= (roulette-stochastic $population $fitnesses)
    (let* (
        ($max_fitness (max-atom $fitnesses))
    )
        (if (== $max_fitness 0)
            (nth-elem $population (random-int &rng 0 (size-atom $population)))
            (roulette-helper $max_fitness $population $fitnesses)
        )
    )
)

; !(roulette-stochastic ((1 2 3) (0 0 1)) (7 2))

;;;;;; roulette classic selection function ;;;;;;

(= (sum-list $lst $acc $i)
    (if (< $i (size-atom $lst))
        (sum-list $lst (+ $acc (nth-elem $lst $i)) (+ $i 1))
        $acc
    )
)

(= (roulette-classic $population $fitnesses)
    (let* (
        ($total_fitness (sum-list $fitnesses 0 0))
    )
        (if (== $total_fitness 0)
            (nth-elem $population (random-int &rng 0 (size-atom $population)))
            (let* (
                ($pick (* (random-float &rng 0 1) $total_fitness))
            )
                (roulette-classic-helper $population $fitnesses $pick 0 0)
            )
        )
    )
)

(= (roulette-classic-helper $population $fitnesses $pick $current $i)
    (if (< $i (size-atom $population))
        (let* (
            ($fitness_score (nth-elem $fitnesses $i))
            ($new_current (+ $current $fitness_score))
        )
            (if (> $new_current $pick)
                (nth-elem $population $i)
                (roulette-classic-helper $population $fitnesses $pick $new_current (+ $i 1))
            )
        )
        (nth-elem $population (- (size-atom $population) 1)) ; for fall back
    )
)

; !(roulette-classic ((1 2 3) (0 0 1)) (7 2))

;;;;;; tournament selection function ;;;;;;

(= (tournament-selection $population $fitnesses $k)
    (let* (
        ($selected (tournament-select-helper $population $fitnesses $k ()))
        ($best (car-atom $selected))
        ($best-val (index-atom $best 1))
        ($winner (tournament-best-helper $selected $best $best-val))
    )
        (car-atom $winner)
    )
)

(= (tournament-select-helper $population $fitnesses $k $acc)
    (if (< (size-atom $acc) $k)
        (let* (
            ($idx (random-int &rng 0 (size-atom $population)))
            ($ind (nth-elem $population $idx))
            ($fit (nth-elem $fitnesses $idx))
        )
            (tournament-select-helper $population $fitnesses $k (append $acc ($ind $fit)))
        )
        $acc
    )
)

(= (tournament-best-helper $lst $best $best-val)
    (if (== $lst ())
        $best
        (let* (
            ($cur (car-atom $lst))
            ($cur-val (index-atom $cur 1))
        )
            (if (> $cur-val $best-val)
                (tournament-best-helper (cdr-atom $lst) $cur $cur-val)
                (tournament-best-helper (cdr-atom $lst) $best $best-val)
            )
        )
    )
)

;!(tournament-selection ((1 2 3) (0 0 1) (2 1 0)) (7 2 1) 2)


;;;;;; elite selection function ;;;;;;

;; map fitness to all individuals in population
(= (map-fitness $population $fitnesses $acc $idx)
    (if (< $idx (size-atom $population))
        (map-fitness $population $fitnesses (append $acc (nth-elem $fitnesses $idx)) (+ $idx 1))
        $acc
    )
)

(= (take-first-n $lst $n)
    (take-first-n-helper $lst $n ())
)

(= (take-first-n-helper $lst $n $acc)
    (if (or (== $n 0) (== $lst ()))
        $acc
        (let* (
            ($head (car-atom $lst))
            ($tail (cdr-atom $lst))
        )
            (take-first-n-helper $tail (- $n 1) (append $acc $head))
        )
    )
)

;!(take-first-n (1 2 3) 2)

;; sort population by fitness
(= (sort-by-fitness $population $fitnesses)
    (sort-by-fitness-helper $population $fitnesses ())
)

(= (sort-by-fitness-helper $pop $fits $acc)
    (if (== $pop ())
        $acc
        (let* (
            ($max-fit (max-atom $fits))
            ($idx (index-of $fits $max-fit))
            ($best-ind (nth-elem $pop $idx))
            ($new-pop (remove-nth $pop $idx))
            ($new-fits (remove-nth $fits $idx))
        )
            (sort-by-fitness-helper $new-pop $new-fits (append $acc $best-ind))
        )
    )
)

;; find index of a value in a list
(= (index-of $lst $val)
    (index-of-helper $lst $val 0)
)

(= (index-of-helper $lst $val $idx)
    (if (== $lst ())
        -1
        (let* (
            ($first (car-atom $lst))
        )
            (if (== $first $val)
                $idx
                (index-of-helper (cdr-atom $lst) $val (+ $idx 1))
            )
        )
    )
)

;; remove nth element from a list
(= (remove-nth $lst $n)
    (remove-nth-helper $lst $n 0 ())
)

(= (remove-nth-helper $lst $n $idx $acc)
    (if (== $lst ())
        $acc
        (if (== $idx $n)
            (remove-nth-helper (cdr-atom $lst) $n (+ $idx 1) $acc)
            (remove-nth-helper (cdr-atom $lst) $n (+ $idx 1) (append $acc (car-atom $lst)))
        )
    )
)

;; return top 2 individuals
(= (elite-selection $population $fitnesses)
    (let* (
        ($sorted-pop (sort-by-fitness $population $fitnesses))
        ($top-2 (take-first-n $sorted-pop 2))
    )
        $top-2)
)

;!(elite-selection ((1 2 3) (4 5 6) (7 8 9) (0 1 2)) (8 7 6 5))


;;;;;;; single-point crossover function ;;;;;;

;; get elements from position n to end
(= (drop-first $lst $n)
    (if (== $n 0)
        $lst
        (if (== $lst ())
            ()
            (drop-first (cdr-atom $lst) (- $n 1))
        )
    )
)

;!(drop-first (1 2 3) 2)

(= (single-point-crossover $parent1 $parent2)
    (let* (
        ($len1 (size-atom $parent1))
        ($len2 (size-atom $parent2))
        ($min-len (if (< $len1 $len2) $len1 $len2))
        ;(() (println! $min-len))
        ($crossover-point (random-int &rng 1 $min-len))
        ;(() (println! $crossover-point))
        ($p1-first (take-first-n $parent1 $crossover-point))
        ;(() (println! $p1-first))
        ($p1-second (drop-first $parent1 $crossover-point))
        ($p2-first (take-first-n $parent2 $crossover-point))
        ($p2-second (drop-first $parent2 $crossover-point))
        ($child1 (union-atom $p1-first $p2-second))
        ($child2 (union-atom $p2-first $p1-second))
    )
        ($child1 $child2))
)

(= (single-point-crossover-main $parent1 $parent2)
    (let* (
        ($rand (random-float &rng 0 1))
    )
        (if (> $rand CROSSOVER_RATE)
            ($parent1 $parent2) 
            (single-point-crossover $parent1 $parent2)
        )
    )
)

;!(single-point-crossover-main (1 2 3 4 5) (6 7 8 9 10))
;!(single-point-crossover-main (1 2 3) (4 5 6 7 8))
;!(single-point-crossover-main (1 2 3 4 5 6) (7 8 9))


;;;;;;; random resetting mutation function ;;;;;;

(= (random-reset-mutation $genotype)
    (let* (
        ($rand (random-float &rng 0 1))
    )
        (if (< $rand MUTATION_RATE)
            (single-gene-mutation $genotype)
            $genotype
        )
    )
)

(= (single-gene-mutation $genotype)
    (let* (
        ($genotype-size (size-atom $genotype))
        ($random-index (random-int &rng 0 $genotype-size))
        ($new-gene (random-int &rng 0 (max-atom $genotype)))
    )
        (replace-gene-at-index $genotype $random-index $new-gene)
    )
)

(= (replace-gene-at-index $genotype $index $new-gene)
    (replace-gene-helper $genotype $index $new-gene 0 ())
)

(= (replace-gene-helper $genotype $index $new-gene $current-index $acc)
    (if (== $genotype ())
        $acc
        (let* (
            ($head (car-atom $genotype))
            ($tail (cdr-atom $genotype))
        )
            (if (== $current-index $index)
                (replace-gene-helper $tail $index $new-gene (+ $current-index 1) (append $acc $new-gene))
                (replace-gene-helper $tail $index $new-gene (+ $current-index 1) (append $acc $head))
            )
        )
    )
)

;!(random-reset-mutation (1 2 3 4 5))
;!(random-reset-mutation (0 1 2 3))




