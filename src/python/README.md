HTML page generation happens in two stages.

1. First, plots are generated from reports

   ```python
    ./contentgen.py <report dir> <destination dir>
   ```

2. Then, HTML file is generated

   ```python
     ./htmlgen.py <output filename>
   ```

   (ATM, can be generated independently. That might not hold true in the future.)
