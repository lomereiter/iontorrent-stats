HTML page generation happens in two stages.

1. First, plots are generated from reports

   ```python
    ./contentgen.py <report dir> <destination dir>
   ```

2. (TBD) Then, directory with HTML and images is generated

   ```python
     ./htmlgen.py <content dir> <output dir>
   ```

   (here `<content dir>` is `<destination dir>` of the
   previous command)
