# NEIVAv1.0

This repository consists of a series of emission factor (EF) databases and processed datasets, including the recommended EF dataset that has EFs for over 1500 constituents averaged over laboratory and field studies organized by 14 globally relevant fuel and fire types. NEIVA is documented in [NEIVAv1.0: Next-generation Emissions InVentory expansion of Akagi et al. version 1.0](https://egusphere.copernicus.org/preprints/2024/egusphere-2024-633/), currently under review.

### Citation
* NEIVA is documented in: (link). We ask that if you use this inventory, you cite this paper.
* If you would like, we encourage you to add any papers citing (link) to our NEIVA papers list.
 [NEIVA papers list](https://docs.google.com/spreadsheets/d/1uXLA59hYS1TJNgUj3USroiDX7IaCfrBNx_SZjSJkd6Q/edit#gid=0)

## Organization of the repository

 * [data](data): The database contents are available as '.sql' files. Additionally, files in CSV format are located within the data/csv/ directory.
 * [python_scripts](python_scripts): Functions used to create the datasets, add new data, query data.
 * [jupyter_notebooks](jupyter_notebooks): Jupyter notebooks demonstrating querying data, and adding new data and generating new datasets using the Python package, neivapy. Jupyter notebooks demonstrating querying data using MySQL.

## Contact
* Technical questions: Samiha Binte Shahid (sbint003@ucr.edu)
* General questions: Kelley Barsanti (barsanti@ucar.edu)