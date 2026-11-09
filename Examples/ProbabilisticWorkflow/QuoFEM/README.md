# Using ASREp with quoFEM

This example demonstrates how to use **ASREp** with **quoFEM** for uncertainty quantification and sensitivity analysis.

## Setup

1. Open the `main.py` script and modify the `ASREp_path` variable on **line 7** so that it points to your local ASREp installation directory.

   ```python
   ASREp_path = "/path/to/your/ASREp"
   ```

2. The `main.py` script serves as the interface between **quoFEM** and **ASREp**.

3. The `params.py` script defines templates for the uncertain variables used as inputs for Monte Carlo simulations.

## Running the Analysis

1. Launch **quoFEM**.

2. Load one of the provided input files using:

   **File → Open → Input File**

3. Select the appropriate JSON input file:

   * `quoFEM_input_forward_propagation.json`
     Example of **forward uncertainty quantification**.

   * `quoFEM_input_sensitivity.json`
     Example of **sensitivity analysis**.

4. Click **Run** to start the analysis.

## Results

After the analysis is complete, a summary of the simulation results can be found in:

```text
response.csv
```

The file is located in the **quoFEM working directory**.

## quoFEM Documentation

For instructions on downloading, installing, and using quoFEM, see the official documentation:

[quoFEM — NHERI SimCenter](https://simcenter.designsafe-ci.org/research-tools/quofem-application/)
