# Replication Package for: Status Externalities in Education and Low Birth Rates in Korea (AER, 2024)
## May 2025

This replication package replicates part of the results in Kim, Tertilt and Yum (2024), using *Julia* as the unique programming language.

We want to stress that even tough the outcome is the same in styling, some values might differ. We are not aware of all the causes of such divergences.

Code can run for a while, in particular the estimation of the OLG model can be very time consuming. Check instructions to change the maximum number of iterations of the various loops.
## Authors

- Luca Emanuele Marcianò - PhD, Collegio Carlo Alberto
- Francesco Sbrana - PhD, Collegio Carlo Alberto
## Instructions to replicators

- Download the GitHub Repo
- Edit `src/PkgTemplates.jl` to change the directory: set "main_dir" to "yourdirectory/PkgTemplates"
- Run `src/PkgTemplates.jl`  
- In the julia command line run:
    - `using .sbrana`to load the package
    - `sbrana.run()` to run the package   

### Additional Instructions
The OLG model is defined in `src/Table45.jl`, 
as convergence might require a lot of time we provide informations to change the number of iterations of each loop:
- Outer loop : line *130* 
- Inner loop 1 : line *141*   
- Inner loop 2 : line *266*

The Bootstrapping of the standard errors for the IV regression is set at 999 repetition, to change this:
- open `src/Table3.jl` 
- go to lines *432, 433, 435, 435*
- set the parameter `B` to the desired number  

To run Tests go to `test/test.jl` and set your working directory/src/Custom_fct.jl in `include`. 
## Dataset list

| Data file | Source | Notes    |Provided |
|-----------|--------|----------|---------|
| data_clean.dta | KTY (2024) | Serves as input for ... | Yes |
| cpi.dta | KTY (2024) | Serves as input for ... | Yes |
| children_school.dta| KTY (2024) | Serves as input for ... | Yes |
| data_educost_detail_wide.dta| KTY (2024) | Serves as input for ... | Yes |
| Private_Education_Participation_Rate_by_School_Level.csv| KTY (2024) | Serves as input for ... | Yes |

## List of tables and figures

The provided code reproduces *Figures 1a* and *1b, Tables 1, 2, 3, 4* and *5* and two in text numbers, namely *lifetime share of education spending* and *income elasticity of investment*. 
Figure 2 could not be reproduced, since it requires access to *not publicly* available data (PEES dataset). It would be possible to access such data upon using a South Korean phone number. Since this dataset cannot be shared, it was not available in the original replication package.


| Figure/Table #    | Program                  | Line Number | Output file                      |
|-------------------|--------------------------|-------------|----------------------------------|
| Figure 1a         | Figure1.jl              |     213        | Figure_1_A.pdf                    |                    
| Figure 1b         | Figure1.jl              |    258         | Figure_1_B.pdf                    |     
| Table 1           | Table1.jl                |  92           | Table1.tex                       |
| Table 2           | Table2.jl                |  325         | table2.tex                         |  
| Income elasticity of investment | Table2.jl | 285 | Income_elasticity.csv |
| Lifetime share of education spending | Table2.jl | 189| Lifetime_share.csv |
| Table 3           | Table3.jl                |  496       | table3.tex                          |   
| Table 4           | Table45.jl               |   799       |  table4.tex                                   |   
| Table 5           | Table45.jl               |   853       |      table5.tex                               |   
# Replication Package for: Status Externalities in Education and Low Birth Rates in Korea (AER, 2024)
## May 2025

This replication package replicates part of the results in Kim, Tertilt and Yum (2024), using *Julia* as the unique programming language.

We want to stress that even tough the outcome is the same in styling, some values might differ. We are not aware of all the causes of such divergences.

Code can run for a while, in particular the estimation of the OLG model can be very time consuming. Check instructions to change the maximum number of iterations of the various loops.
## Authors

- Luca Emanuele Marcianò - PhD, Collegio Carlo Alberto
- Francesco Sbrana - PhD, Collegio Carlo Alberto
## Data availability and provenance statements
### Statement about rights

The author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.

### Summary of availability

- All data available in this replication package *are* publicly available.

### Details on each data source

- This paper uses data *directly coming* from the original replication package of Kim, Tertilt and Yum (2024). Only publicly available data are hence held inside this package. Refer to the original replication package for details on original data sources.

## List of tables and figures

The provided code reproduces *Figures 1a* and *1b, Tables 1, 2, 3, 4* and *5* and two in text numbers, namely *lifetime share of education spending* and *income elasticity of investment*. 
Figure 2 could not be reproduced, since it requires access to *not publicly* available data (PEES dataset). It would be possible to access such data upon using a South Korean phone number. Since this dataset cannot be shared, it was not available in the original replication package.


| Figure/Table #    | Program                  | Line Number | Output file                      |
|-------------------|--------------------------|-------------|----------------------------------|
| Figure 1a         | Figure1.jl              |     213        | Figure_1_A.pdf                    |                    
| Figure 1b         | Figure1.jl              |    258         | Figure_1_B.pdf                    |     
| Table 1           | Table1.jl                |  92           | Table1.tex                       |
| Table 2           | Table2.jl                |  325         | table2.tex                         |  
| Income elasticity of investment | Table2.jl | 285 | Income_elasticity.csv |
| Lifetime share of education spending | Table2.jl | 189| Lifetime_share.csv |
| Table 3           | Table3.jl                |  496       | table3.tex                          |   
| Table 4           | Table45.jl               |   799       |  table4.tex                                   |   
| Table 5           | Table45.jl               |   853       |      table5.tex                               |   

## References

- Kim, Seongeun, Michèle Tertilt, and Minchul Yum. 2024. "Status Externalities in Education and Low Birth Rates in Korea." American Economic Review 114 (6): 1576–1611.
## Data availability and provenance statements
### Statement about rights

The author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.

### Summary of availability

- All data available in this replication package *are* publicly available.

### Details on each data source

- This paper uses data *directly coming* from the original replication package of Kim, Tertilt and Yum (2024). Only publicly available data are hence held inside this package. Refer to the original replication package for details on original data sources.

## License

[MIT](https://choosealicense.com/licenses/mit/)


## References

- Kim, Seongeun, Michèle Tertilt, and Minchul Yum. 2024. "Status Externalities in Education and Low Birth Rates in Korea." American Economic Review 114 (6): 1576–1611.