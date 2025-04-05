#!/bin/bash

# Clear previous output
> output_groebner.txt

# Loop 1000 times
for i in {1..1000}; do
    # Run the Python script with arguments and append the output to output.txt
    python3 groebner_solver.py -q=5 -n=5 --same_dim --minimal >> output_groebner.txt
done