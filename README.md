# smFishPlateDesigner
The `smFishPlateDesigner` package offers a suite of tools for processing and analyzing high-throughput single-molecule fluorescence in situ hybridization (HT-smFish) data. It provides functionalities for creating PCR plate layouts, analyzing TIV dosage data, and more.

#### Installation from GitHub and the Tarball
You can install the latest version of `smFishPlateDesigner` from GitHub using the source tarball. First, download the tarball (`.tar.gz` file) from the GitHub releases or clone/download the repository. Then install the package using the `install.packages` function in R:

```r
# Replace "path/to" with the actual path to the tarball
install.packages("path/to/smFishPlateDesigner_X.Y.Z.tar.gz", repos = NULL, type = "source")
```

Make sure to replace `"path/to/smFishPlateDesigner_X.Y.Z.tar.gz"` with the actual path to the downloaded tarball.

## Getting Started

### Basic Usage

The `smFishPlateDesigner` package is designed to work with the `targets` package to create reproducible workflows for HT-smFish experiments. To get started, you'll need to set up a few things:

1. **Create an Output Directory**: Ensure you have an `output/` directory in your project, as the package will save some results to files in this directory.

2. **Set Excel File Path**: If you have a specific Excel file for your data, set its path before running the workflow:

   ```r
   options(smFishPlateDesigner_excel_path = "path/to/your/excel_file.xlsx")
   ```

   If you do not have a specific file and wish to use the default embedded file, no action is required for the file path.

3. **Run the Workflow**: Execute the workflow using the following command in your R console:

   ```r
   library(targets)
   tar_make()
   ```

### Further Information

For a detailed guide on using `smFishPlateDesigner`, including explanations of each target in the workflow, refer to the package vignette provided in the `man` folder or check the [targets documentation](https://docs.ropensci.org/targets/).

## Documentation

The `smFishPlateDesigner` package comes with detailed documentation located in the `doc` folder. These documents offer an in-depth understanding of how to use the package and its features.

You can find these documents in the form of HTML files in the `doc` folder of the package's [GitHub repository](https://github.com/cedrichk/smFishPlateDesigner/tree/main/doc).

## Issues and Contributions

Feel free to report any issues or feature requests on the [GitHub repository](https://github.com/cedrichk/smFishPlateDesigner). Contributions to the package are also welcome.


