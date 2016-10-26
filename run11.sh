#!/bin/bash
for i in {1..4}
	do
		./gerard11 sample num_warmup=5000 num_samples=5000 thin=10 random seed=12345 \
			id=$i data file=data.R init=init11.R output file=output11_$i.csv refresh=1000 &
done
