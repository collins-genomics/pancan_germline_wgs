# G2C helper libraries
### Source code for G2C helper libraries in R or Python

Copyright (c) 2023-Present, [Ryan L. Collins](mailto:Ryan_Collins@dfci.harvard.edu), the [Collins Genomics Lab](https://labs.dana-farber.org/collins-genomics/), and the [Van Allen Lab](https://labs.dana-farber.org/vanallenlab/) at the Dana-Farber Cancer Institute.  

We provide easy access to commonly used functions throughout the G2C project via two libraries:  
* `G2CR` : an R companion library with G2C-relevant functions; depends heavily on [`RLCtools`](https://github.com/RCollins13/RLCtools), and  
* `g2cpy` : a Python companion package with G2C-relevant functions.

---  

## G2CR  

### Installation & invocation  

In an R session:  

```
> install.packages("G2CR_0.2.0.tar.gz", type="source", repos=NULL)
> library(G2CR)
```

### Documentation  

All functions in `G2CR` are documented using R standards, meaning that full help text can be accessed for each function by prefixing each function with `?`.  

For example, in an R session:  
```
> library(G2CR)
> ?explode.by.cancer
```

The same holds true for all functions in [`RLCtools`](https://github.com/RCollins13/RLCtools), which contains most of the ultra-generic custom functions applied throughout this codebase.  

---  

## g2cpy  

### Installation & invocation  

The `g2cpy` Python module can be installed from the command line as follows:  

```
git clone https://github.com/collins-genomics/pancan_germline_wgs
cd pancan_germline_wgs/src/g2cpy
pip install -e .
```

The full library (or functions within) can be loaded with `import` as standard for all Python modules, and all functions imported by the `g2cpy` module can be listed with `dir(g2cpy)` as normal.  

### Documentation  

The `g2cpy` package is currently undocumented. This will be addressed eventually over time. Please contact [Ryan Collins](mailto:Ryan_Collins@dfci.harvard.edu) with questions.  
