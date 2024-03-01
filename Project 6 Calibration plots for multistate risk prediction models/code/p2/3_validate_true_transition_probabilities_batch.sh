cd "/mnt/bmh01-rds/mrc-multi-outcome/Project_6/code/p2/"
module load apps/gcc/R/4.1.3-gcc830
nohup Rscript 3_validate_true_transition_probabilities.R 1 > 3_validate_true_transition_probabilities_set1.out &
nohup Rscript 3_validate_true_transition_probabilities.R 2 > 3_validate_true_transition_probabilities_set2.out &
nohup Rscript 3_validate_true_transition_probabilities.R 3 > 3_validate_true_transition_probabilities_set3.out &
nohup Rscript 3_validate_true_transition_probabilities.R 4 > 3_validate_true_transition_probabilities_set4.out &
nohup Rscript 3_validate_true_transition_probabilities.R 5 > 3_validate_true_transition_probabilities_set5.out &
nohup Rscript 3_validate_true_transition_probabilities.R 6 > 3_validate_true_transition_probabilities_set6.out &