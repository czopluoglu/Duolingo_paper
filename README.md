**Repository Overview**
This repository contains the code and resources associated with the following paper:

    Zopluoglu, C., Lockwood, J.R. (under review). A Comparative Study of Item Response Theory Models for Mixed Discrete-Continuous Responses. Journal of Intelligence.

For a tutorial-style introduction to the analyses conducted in the paper, please visit: https://czopluoglu.github.io/Duolingo_paper/

**Directory Structure and Contents**

1. **/data/**:

    - Contains simulated datasets in both long and wide formats for the Beta, Simplex, and SB IRT models.

2. **/script/beta_irt/**:
    - R script for simulating the dataset based on the Beta IRT model.

    - Scripts for fitting the model in Stan and evaluating model fit using posterior predictive checks.

    - Includes the Stan model syntax (beta.stan) and its compiled version for the cmdstanr package.

3. **/script/simplex_irt/:**

    - R script for simulating the dataset based on the Simplex IRT model.

    - Scripts for fitting the model in Stan and evaluating model fit using posterior predictive checks.

    - Includes the Stan model syntax (simplex.stan) and its compiled version for the cmdstanr package.

4. **/script/sb_irt/:**

    - R script for simulating the dataset based on the SB IRT model.

    - Scripts for fitting the model in Stan and evaluating model fit using posterior predictive checks.

    - Includes the Stan model syntax (sb.stan) and its compiled version for the cmdstanr package.

5. **/docs/:**

    - Contains the Rmarkdown file and the rendered HTML for the tutorial page.