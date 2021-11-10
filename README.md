# Generalized Unscented Transform

## Description
This package contains the code that can be used to reproduce the results in our paper titled `A Generalized Unscented Transformation for Probability Distributions`. 

The paper can be found [here](https://arxiv.org/abs/2104.01958).

The code structure, run instructions, and any special notes are given below.

### Code Structure
- **`Case_study_1a.m`** (Generates the results of the first part of case study 1 in our paper)
- **`Case_study_1b.m`** (Generates the results of the second part of case study 1 in our paper)
- **`Case_study_1c.m`** (Generates the results of the last part of case study 1 in our paper)
- **`Case_study_2.m`** (Generates the results of case study 2 in our paper)
- **`Case_study_3.m`** (Generates case study 3 in our paper)
- **`Evaluate_sample_statistics.m`** (Evaluates the sample statistics of sigma points)
- **`Examples_1_and_2.m`** (Generates the results Examples 1 and 2 in our paper)

- **/Unscented Transforms**
    - `GenUT_Ensemble.m` (Generates the sigma points using our generalized unscented transform)
    - `HOSP.m` (Generates the sigma points using the higher order sigma point unscented transform)
    - `unscentedEnsemble.m` (Generates the sigma points using the scaled unscented transform)

- **/Distribution Moments**
    - `Beta_moments.m` (Generates the moments for a Beta distribution)
    - `Binomial_moments.m` (Generates the moments for a Binomial distribution)
    - `Exponential_moments.m` (Generates the moments for an Exponential distribution)
    - `Gamma_moments.m` (Generates the moments for a Gamma distribution)
    - `Gaussian_moments.m` (Generates the moments for a Gaussian distribution)
    - `Geometric_moments.m` (Generates the moments for a Geometric distribution)
    - `Negative_Binomial_moments.m` (Generates the moments for a Negative Binomial distribution)
    - `Poisson_moments.m` (Generates the moments for a Poisson distribution)
    - `Rayleigh_moments.m` (Generates the moments for a Rayleigh distribution)
    - `Weibull_moments.m` (Generates the moments for a Weibull distribution)

---

### Run Instructions
Below are the run instructions to generate the results shown in our paper, as well as to generate resuts for other applications. We note that the only files capable of running are the `Examples_1_and_2.m`, `Case_study_1a.m`, `Case_study_1b.m`, `Case_study_1c.m`, `Case_study_2.m`,  and `Case_study_2.3.m` files. Other `m-files` only run by supplying the appropriate input arguments. Please see the `m-file` you want to run for its inut argument description.

#### Run examples 1 and 2 in the paper
Run the `Examples_1_and_2.m` file to see the performance of the constrained and unconstrained versions of our sigma points.

#### Run case study 1 in the paper
1. Run the `Case_study_1a.m` file to see the performance on the quadratic function of a random variable of different probability distributions. 
2. Run the `Case_study_1b.m` file to see the performance on the trigonometric function of a random variable of different probability distributions.
3. Run the `Case_study_1c.m` file to see the performance on the trigonometric function of a Poisson and Weibull random variable.

#### Run case study 2 in the paper
Run the `Case_study_2.m` file to see the performance on the vector of trigonometric functions.

#### Run case study 3 in the paper
Run the `Case_study_3.m` file to see the performance on the infectious disease model.

### Generate sigma points using our Generalized Unscented Transform (GenUT)
Run the `GenUT_Ensemble.m` file. The minimum input arguments that need to be supplied are the mean and covariance.

To constrain the sigma points, the mean must fall within the lower bound (`lb`) and upper bound (`ub`) arguments. 

### Analyze sample statistics of sigma points
Run the `Evaluate_sample_statistics.m` file by passing in two arguments; the sigma points and the weights.
