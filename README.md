# CLL798
# SOC Earthquake Sandpile Project

## Self-Organized Criticality in Earthquake Fault Networks

A computational study using the Bak-Tang-Wiesenfeld (BTW) sandpile model to demonstrate self-organized criticality in earthquake fault dynamics.

### Building and Running

```bash
# Compile
g++ -O2 -std=c++17 -o sandpile src/sandpile.cpp

# Run simulation
./sandpile --L 64 --grains 500000 --transient 50000 --outdir .

# Analyze results
pip install numpy matplotlib scipy
python3 src/analyze.py .
```

### Project Structure

```
├── src/
│   ├── sandpile.cpp      # C++ simulation
│   └── analyze.py        # Python analysis & plotting
├── report/
│   ├── report.tex        # LaTeX manuscript
│   └── report.pdf        # Compiled PDF
├── README.md
```
