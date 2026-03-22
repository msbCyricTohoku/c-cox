# C-Cox: High-Performance Cox Proportional Hazards Regression

C-Cox is a lightweight, standalone C implementation of the Cox Proportional Hazards (PH) model. It is designed for high-performance survival analysis on massive datasets where memory efficiency and computational speed are critical.

Recently updated, C-Cox now fully supports advanced clinical trial and Electronic Health Record (EHR) methodologies, including **Time-Varying Covariates**, **Clustered Longitudinal Data**, and **Competing Risks**. These features are essential for modern biostatistics.

---

## Project Files

* `ccox.c`: Legacy source (does not support time-varying covariates).
* `main.c`, `ccox_math.c`, `ccox.h`: Current source supporting time-varying covariates, clustered longitudinal data, and competing risks.

---

## Mathematical Foundations & Advanced Features

### 1. The Base Cox Model (Time-Varying Covariates)
The model estimates the effects of covariates on the hazard function over time. By supporting time-varying datasets, covariates $X_i(t)$ can dynamically change as time progresses:

$$\lambda(t | X) = \lambda_0(t) \exp\left(\sum_{i=1}^{p} \beta_i X_i(t)\right)$$

### 2. Robust Sandwich Variance (Clustered Data)
When patients have multiple clinic visits, their records are statistically correlated. A standard Cox model assumes every row represents a completely independent person, which artificially shrinks error margins and leads to false discoveries. C-Cox solves this by applying the **Lin-Wei Robust Sandwich Estimator** grouped by patient clusters:

$$V_{robust} = I^{-1} B I^{-1}$$

Where $I$ is the Information Matrix (Inverse Hessian) and $B$ is the variance of the cluster-grouped score residuals. This produces accurate standard errors (`Rob_SE`).

### 3. Competing Risks (Aalen-Johansen Cumulative Incidence, CIF)
Standard survival curves (like Kaplan-Meier) fail when competing risks exist (e.g., a patient dies of an unrelated cause and can no longer experience the target event). C-Cox isolates your specific target event and calculates the true absolute probability over time using the **Aalen-Johansen Estimator**:

$$\hat{F}_1(t) = \sum_{t_j \le t} \hat{S}(t_{j-1}) \frac{d_{1j}}{n_j}$$

Where $\hat{S}$ is the overall survival from all causes, ensuring the risk of the target event is not artificially inflated by impossible scenarios.

---

## Supported Dataset Style: Counting Process

C-Cox relies on the counting process (Start and Stop) format, which is the standard for longitudinal data in modern biostatistics. It allows users to update a patient's health data (like aging or fluctuating BMI) dynamically across multiple intervals by simply creating a new row per visit.

**Example CSV Format:**

| patient_id | start | stop | status | treatment | age | bmi |
|---|---|---|---|---|---|---|
| 1 | 0.00 | 4.64 | 0 | 1 | 65.20 | 24.92 |
| 1 | 4.64 | 5.27 | 0 | 1 | 69.84 | 25.10 |
| 1 | 5.27 | 7.40 | 1 | 1 | 70.47 | 24.01 |
| 2 | 0.00 | 2.50 | 2 | 0 | 55.00 | 30.00 |

* **patient_id**: Groups rows belonging to the same subject (clustering).
* **start / stop**: The exact time window the patient was observed.
* **status**: 
  * 0 = Censored (Routine visit, no event)
  * 1 = Target Event (e.g., Heart Attack)
  * 2 = Competing Risk (e.g., Death from other causes)
* **Covariates**: Variables like age and BMI update dynamically at every visit.

---

## Features

* Massive Parallelism.
* Binary Search.
* Minimal Dependencies.
* Tie Handling.
* Memory Efficient.

---

### Installation

To compile the library, use any standard C compiler (e.g., GCC or Clang) with OpenMP support.

### Install GSL dependencies (Ubuntu/Debian example)
`sudo apt update && sudo apt install -y libgsl-dev`

### Compile with OpenMP and GSL flags
`gcc -O3 -fopenmp main.c ccox_math.c -o ccox -lgsl -lgslcblas -lm`

### Alternatively, if you are using a Makefile:
`make`

### Run the model
`./ccox config.dat`

---

#### Configuration File (config.dat)

C-Cox is controlled entirely by a simple text configuration file. No hardcoding is required in the C source. This allows users to toggle advanced statistical features on or off instantly.

#### Parameters
`MAX_ITER=100`                 # Max Newton-Raphson iterations (default: 100)

`TOLERANCE=1E-7`               # Convergence tolerance (default: 1E-7)

#### Data Settings
`file=dummy_data.csv`  # Path to your CSV dataset

`n=100000`                     # Max number of rows to read from the CSV

`covno=3`                      # Number of predictive covariates to analyze

`covariates=treatment,age,bmi` # Exact column names of covariates separated by commas

#### Column Mappers
`start_col=start`             # Column name for interval start time (leave blank for standard Cox)

`stop_col=stop`                # Column name for interval end time

`status_col=status`            # Column name for the event status

`cluster_col=patient_id`       # Column name for patient grouping (leave blank for independent rows)

#### Advanced
`event_code=1`                 # The specific status integer that represents your target event

`rob_se=1`                     # 1 = Turn ON Robust Sandwich Estimator (Rob_SE), 0 = Standard SE

`cif_times=1.0,5.0,10.0`       # Comma separated times to output Aalen-Johansen CIF absolute risk %

---

`Developed by msb -- 2026`
