# Bayesian test of introgression

This repository implements the Bayesian test of introgression developed by Ji et al. (2022). The program calculates the Bayes factor via the Savage-Dickey density ratio using an MCMC sample under the MSci or MSC-M to evaluate the evidence for a proposed gene flow event. Note that the MSci models gene flow as a discrete event that occurred at some fixed time, with its magnitude estimated through the introgression probability, $\varphi$. In contrast, the MSC-M assumes continuous gene flow at a particular rate every generation, as denoted by the migration rate, $M$.

The Bayes factor is defined as

$$
B_{10,\epsilon} = \frac{1-\mathbb{P}(ø|x)}{\mathbb{P}(ø|x)} \biggl/ \frac{1-\mathbb{P}(ø)}{\mathbb{P}(ø)} \approx \frac{\mathbb{P}(ø)}{\mathbb{P}(ø|x)},
$$

where $\mathbb{P}(ø|x)$ represents the proportion of MCMC samples in which $\varphi$ or $M$ are less than $\epsilon$, and $\mathbb{P}(ø)$ refers to the probability that a random value taken from the prior distribution on $\varphi$ or M is less than $\epsilon$ (given by the cumulative distribution function). Intuitively, $B_{10} \ll 1$ is evidence for $H_0$ (null model: no gene flow) and $B_{10} \gg 1$ is evidence for $H_1$ (alternative model: gene flow).

# Installation and Usage

The program is currently distributed as an R script. To install and run do the following:

- Clone the repository or download the source files.
- Navigate to the directory containing the script.
- Make sure you have R installed on your machine. If not, download and install R from https://www.r-project.org/.

The `stats` package (>= 4.2.2) must be installed. This can be done with the following command:

```r
install.packages("stats")
```

From the command line, the script can be run as follows:

```text
Rscript BF-script.R [function] [alpha] [beta] [epsilon] [file] [column_indices]

function          Either BF_Gamma or BF_Beta.
alpha             Numeric value for the alpha parameter.
beta              Numeric value for the beta parameter.
epsilon           Numeric value for the epsilon parameter.
file              Path to the MCMC sample file.
column_indices    Column indices, separated by spaces. Ranges can be specified using ':'.
```

For example:

```r
Rscript BF-script.R BF_Gamma 2 10 0.01 /path/to/sample-mcmc.txt 25:30
```
would use the BF_Gamma function with alpha = 2, beta = 10, epsilon = 0.01, an MCMC sample file located at /path/to/mcmc_file.txt, and columns 25:30. The MCMC file used in this example, along with the output file, can be found in the 'test' folder. The usage instructions can be viewed by typing `Rscript BF-script.R --help`.

# References

Jiayi Ji, Donavan J Jackson, Adam D Leach&eacute;, and Ziheng Yang. Power of bayesian and heuristic tests to detect cross-species introgression with reference to gene flow in the tamias quadrivittatus group of north american chipmunks. Systematic Biology, 2022. ISSN 1063-5157
