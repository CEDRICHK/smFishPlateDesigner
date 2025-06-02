# smFishPlateDesigner

The `smFishPlateDesigner` package offers a suite of tools for processing and analyzing high-throughput single-molecule fluorescence in situ hybridization (HT-smFish) data. It provides functionalities for creating PCR plate layouts, analyzing TIV dosage data, and more.

## Prerequisites

Before installing `smFishPlateDesigner`, ensure you have installed the following R packages: `tidyverse`, `magrittr`, `xlsx`, `readxl` and `targets`. You can install these packages using the following commands:

```r
install.packages(c("tidyverse", "magrittr", "xlsx", "readxl", "targets"))
```

## Installation from GitHub and the Tarball

You can install the latest version of `smFishPlateDesigner` from GitHub using the source tarball. First, download the tarball (`.tar.gz` file) from the GitHub releases. Then install the package using the `install.packages` function in R:

```r
# Replace "path/to" with the actual path to the tarball
install.packages("path/to/smFishPlateDesigner_X.Y.Z.tar.gz", repos = NULL, type = "source")
```

Make sure to replace `"path/to/smFishPlateDesigner_X.Y.Z.tar.gz"` with the actual path to the downloaded tarball.

## Getting Started with `_targets.R` Workflow

To streamline the analysis process, `smFishPlateDesigner` is designed to work with the `targets` package for creating reproducible workflows. Here's how to set it up:

1. **Ensure Prerequisite Packages Are Installed**: Follow the instructions above to install required packages.

2. **Create an Output Directory**: Ensure you have an `output/` directory in your project, as the package will save some results to files in this directory. 

3. **Use the `_targets.R` Script**: Instead of using source files directly, incorporate your analysis workflow within a `_targets.R` script. An example of this script can be found in the package's GitHub repository.

4. **Configure `tar_option()`**: In your `_targets.R` script, ensure to include the `smFishPlateDesigner` package within the `tar_option()` function to make it available during the execution of targets.

5. **Set Up Excel File Options**: If using a specific Excel file for your data, configure it as an option within your `_targets.R` script:

   ```r
   options(smFishPlateDesigner_excel_path = "path/to/your/excel_file.xlsx")
   ```
   This allows `smFishPlateDesigner` to access your Excel file during the analysis.
   
6. **Execute the Workflow**: With your `_targets.R` script configured, run the workflow using the following command in your R console:

   ```r
   library(targets)
   tar_make()
   ```
This setup ensures that your HT-smFish experiments' analysis is reproducible and efficiently managed.

### Further Information

For a detailed guide on using `smFishPlateDesigner`, including explanations of each target in the workflow, refer to the package vignette provided in the `man` folder or check the [targets documentation](https://docs.ropensci.org/targets/).

## Documentation

The `smFishPlateDesigner` package comes with detailed documentation located in the `doc` folder. These documents offer an in-depth understanding of how to use the package and its features.

You can find these documents in the form of HTML files in the `doc` folder of the package's [GitHub repository](https://github.com/cedrichk/smFishPlateDesigner/tree/main/doc).

## Issues and Contributions

Feel free to report any issues or feature requests on the [GitHub repository](https://github.com/cedrichk/smFishPlateDesigner). Contributions to the package are also welcome.


