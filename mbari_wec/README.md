Example plotting for the MBARI WEC 2021 dataset. For more information about the MBARI WEC, please see https://doi.org/10.1007/s40722-021-00197-9.

## Getting started

1. Install [Anaconda](https://anaconda.org) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

2. Use the `environment.yml` file in the root directory of this repository to create a dedicated [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) and install the packages used in this project.

	```bash
	conda env create --file environment.yml
	```

	Note that the environment can be updated to match any changes in the `environment.yml` file by calling

	```bash
	conda env update -f environment.yml
	```

3. Activate the new environment.

	```bash
	conda activate mbari_wec_2021
	```

4. Download the data from MHK-DR: 

5. Call the plotting example script.

	```bash
	python mbari_wec_2021_example.py
	```
