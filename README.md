This github page contains the source code (as well as some descriptions) for the "mrdd" package. This package implements the two-dimensional RD estimation described in Cheng (2020), which is available on aldencheng.com.
The "mrdd" R package can also be installed in R if the user has the "devtools" package installed, by running the command:
```
devtools::install_github("alden1505/rd_multiple_running_var_code")
```
The two main functions are mrdd.fn (which estimates the two-dimensional RD), and create_plot (which generates plots based on the estimates).

**References**

Cheng, Alden, 2020. "Regression Discontinuity Designs with Multiple Running Variables." Working Paper.
