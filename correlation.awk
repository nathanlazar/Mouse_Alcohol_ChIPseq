#!/usr/bin/awk -f

## Calculates Pearson Correlation Coefficient in one-pass
## pseudocode from http://en.wikipedia.org/wiki/Talk:Correlation_and_dependence
## Code from http://www.starklab.org/data/bardet_natprotoc_2011/correlation.awk

## Input: 2 cols from stdin

BEGIN{
    sum_sq_x = 0
    sum_sq_y = 0
    sum_coproduct = 0
    N = 0
}

(NF==2){
    N++
    if(NR==1){
	mean_x = $1
	mean_y = $2
    }else{
	sweep = (NR - 1.0) / NR
	delta_x = $1 - mean_x
	delta_y = $2 - mean_y
	sum_sq_x += delta_x * delta_x * sweep
	sum_sq_y += delta_y * delta_y * sweep
	sum_coproduct += delta_x * delta_y * sweep
	mean_x += delta_x / NR
	mean_y += delta_y / NR
    }
}

END{
    if(N > 0){
	pop_sd_x = sqrt( sum_sq_x/N )
	pop_sd_y = sqrt( sum_sq_y/N )
	cov_x_y = sum_coproduct/N
    }
    if(pop_sd_x * pop_sd_y != 0){
	correlation = cov_x_y / (pop_sd_x * pop_sd_y)
	print " "correlation
    }else{
	print " ""-"
    }
}
