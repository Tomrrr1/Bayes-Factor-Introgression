# Bayesian test of introgression

This repository implements the Bayesian test of introgression developed by Ji et al. (2023). The program calculates the Bayes factor via the Savage-Dickey density ratio using an MCMC sample under the MSci or MSC-M to evaluate the evidence for a proposed gene flow event. Note that the MSci models gene flow as a discrete event that occurred at some fixed time, with its magnitude estimated through the introgression probability, $\varphi$. In contrast, the MSC-M assumes continuous gene flow at a particular rate every generation, as denoted by the migration rate, $M$. Both the MSci and MSC-M are implemented in the Bayesian programe BPP (Flouri et al. 2018, 2020).

The Bayes factor representing the evidence for $H_1$ (gene flow) over $H_0$ (no gene flow) is approximately

$$
B_{10,\epsilon} = \frac{1-\mathbb{P}(ø|x)}{\mathbb{P}(ø|x)} \biggl/ \frac{1-\mathbb{P}(ø)}{\mathbb{P}(ø)} \approx \frac{\mathbb{P}(ø)}{\mathbb{P}(ø|x)},
$$

where $\mathbb{P}(ø|x)$ represents the proportion of MCMC samples in which $\varphi$ or $M$ are less than $\epsilon$, and $\mathbb{P}(ø)$ refers to the probability that a random value taken from the prior distribution on $\varphi$ or M is less than $\epsilon$ (given by the cumulative distribution function). Intuitively, $B_{10} \ll 1$ is evidence for $H_0$ (no gene flow) and $B_{10} \gg 1$ is evidence for $H_1$ (gene flow). 

# Installation and Usage

The program is currently distributed as an R script. To install and run do the following:

- Clone the repository or download the source files.
- Navigate to the directory containing the script.
- Make sure you have R installed on your machine. If not, download and install R from https://www.r-project.org/.

The `stats` package (>= 4.2.2) must be installed. This can be done with the following command:

```r
install.packages("stats")
```

The program accepts as input a BPP MCMC sample file. 

From the command line, the script can be run as follows:

```text
Rscript BF-script.R [function] [alpha] [beta] [epsilon] [file] [column_indices]

function          Either BF_Gamma or BF_Beta.
alpha             Numeric value for the alpha parameter of the prior distribution on M or Varphi.
beta              Numeric value for the beta parameter of the prior distribution on M or Varphi.
epsilon           Numeric value for the epsilon parameter.
file              Path to the MCMC sample file.
column_indices    Column indices, separated by spaces. Ranges can be specified using ':'.

If the BF_Gamma option is selected then 'alpha' is the shape and 'beta' is the rate of the prior distribution.
If the BF_Beta option is selected then both 'alpha' and 'beta' are shape parameters of the prior distribution.  
```

For example:

```r
Rscript BF-script.R BF_Gamma 2 10 0.01 /path/to/sample-mcmc.txt 25:30
```

would use the BF_Gamma function with alpha = 2 (shape), beta = 10 (rate), epsilon = 0.01, an MCMC sample file located at /path/to/mcmc_file.txt, and columns 25:30. The MCMC file used in this example, along with the output file, can be found in the 'test' folder.

The output of the program is a text file containing the calculated Bayes factor for each proposed gene flow event. A Bayes factor threshold of 100 means strong support for $H_1$ and rejection of $H_0$ (which is similar to hypothesis testing at the 1\% significance level).

The usage instructions can be viewed by typing `Rscript BF-script.R --help`.

# References

Flouri, T., Jiao, X., Rannala, B., and Yang, Z. 2018. Species tree inference with BPP using genomic sequences and the multispecies coalescent. *Mol. Biol. Evol.*, 35(10): 2585–2593.

Flouri, T., Jiao, X., Rannala, B., and Yang, Z. 2020. A Bayesian implementation of the multispecies coalescent model with introgression for phylogenomic analysis. *Mol. Biol. Evol.*, 37(4): 1211–1223.

Ji, J., Jackson, D. J., Leache, A. D., and Yang, Z. 2023. Power of Bayesian and heuristic tests to detect cross-species introgression with reference to gene flow in the *Tamias quadrivittatus* group of North American chipmunks. *Syst. Biol.*, page 10.1093/sysbio/syac077.

