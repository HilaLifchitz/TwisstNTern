<img src="logo.png" height="300pt" align="bottom">

## A method for analysing topology weights in a ternary framework

sean to write

## Papers
Stankowski et al 2023 is where we first used the TwisstNTern method to study patterns of tree discrodnace in _Littorina_ 

Stankowski, S., Z. B. Zagrodzka, M. Garlovsky, A. Pal, D. Shipilina, D Garcia Castillo, T. Broquet, E. Leader, J. Reeve, K. Johannesson, A. M. Westram, R. K. Butlin. 2023. Selection on many loci drove the origin and spread of a key innovation. _bioRxiv_ doi: https://doi.org/10.1101/2023.02.13.528213

## Running TwisstNtern
There are two ways that you can run TwisstNTern on your data

**Option 1: From the command line**

1) Download twisstntern.py into your working directory; this should also include your input file of topology weights
2) Load python
3) Import the package
```bash
Import twisstntern
```
4) Execute the code
```bash
res=twisstntern.run_analysis(inputfile.csv, granularity)
```
Replace `bash inputfile.csv` is the name of your input file of topology weights. Only .csv format it supported. 
`granularity` determines the number of subtriangles used in the local symmetry analysis. We provided 3 preset options `coarse`, `fine`, and `superfine`. 

The most appropriate granuliaty for your will depend on the number of genomic regions analysed and the dispersion of the data aross the ternary plot. For most WGS studies where the genome has been section into block or arbitrarty size or SNP number, the `fine` option is a good place to start. Emerging tools like the ARG can produce very large numbers of marginal trees, so the `superfine` options may be suitable. Those with reduced representation datasets may have fewer loci than the average WGS study, so may opt for the coarse option. Similarly, in systems where there is fairly high gene tree-sepecies tree concordance, the dsitribution may be restricted to a fairly small section of the triangle so a higher granularity may provide better resolution. 

## Input file


