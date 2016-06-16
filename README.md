# intermediate

`intermediate` is R package for eQTl/pQTL mediation analysis.

## Installation

You can install `intermediate` directly from Github using `devtools` package:

```S
    # install.packages("devtools")
    devtools::install_github("simecek/intermediate")
```

## Example

```
  # DOQTL liver protein expresion dataset
  data(Tmem68)
  
  # Let us mediate Tmem68 to other proteins
  med <- mediation.scan(target = Tmem68$target,
                        mediator = Tmem68$mediator,
                        annotation = Tmem68$annotation,
                        covar = Tmem68$covar,
                        qtl.geno = Tmem68$qtl.geno)
                        
  # Plot mediation results and identify the mediator                      
  plot(med)                        
  identify(med)
  
  # Interactive plot
  kplot(med)
```


## References

Joel M. Chick,	Steven C. Munger,	Petr Simecek,	Edward L. Huttlin,	Kwangbom Choi,	Daniel M. Gatti,	Narayanan Raghupathy,	Karen L. Svenson,	Gary A. Churchill	& Steven P. Gygi: "**Defining the consequences of genetic variation on a proteome-wide scale**", Nature (2016) doi:10.1038/nature18270



