Install all dependencies with [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

    conda env create -f htqpcr_env.yml

Activate environment

    conda activate htqpcr_env

Start Jupyter Notebook

    jupyter-notebook

Open and Run the notebooks
- htqpcr_validation_figures.ipynb
- bm_rawdata_annotation.ipynb

Usage examples for the biomarkdataparser.py script can be found in the notebooks,
you can use the biomarkdataparser.py script in the terminal or in the notebooks.


__Command line options__

|Command line option [Input]|Abbreviation|Description|Default|
|---------------------------|----|-------------------|-------|
|--help|-h|Display help message|-|
|--inputfile [str]|-i|Path to csv file (BM qPCR data export)|-|
|--outdir [str]|-o|Path/Name of the output directory|current directory|
|--validation|-v|Create labels for the annotation of raw Cq data image (examples in bm_rawdata_annotation.ipynb)|-|
|--replicates [assay, sample, None]|-r|Technical replicates were made with assays (primer) or samples|assay|
|--number [int]|-n|Number of technical replicates measured|-|
|--typefilter [str ...]|-t|Filter sample types, e.g. NTC, Unknown, Standard|-|
|--samplefilter [str ...]|-s|Filter sample names\*|-|
|--removesamples [str ...]|-x|Remove samples that have the input string(s) in the sample name|-|
|--standardcurve|-c|Draw plot with standard curves|-|
|--gformat [eps, png, jpg]|-f|Choose an output format for the figures|png|
|--assay_species [str]|-a|Path to csv file (no header) with first column: assay name\* and the second column: the annotation (species) to display in the results (example: labelfiles/assay_species.csv)|-|
|--labelfile [str]|-l|Path to csv file (with header) with figure labels as columns (x, y) example: labelsfiles/cheese_labels.csv|-|
|--title [str]||Optional plot title|-|
|--transpose||Flip x and y axes in heatmap (does not affect labelfiles)|-|
|--italic [x, y ...]||Select axis with species names (italic label), (example with partial italic labels: labelfiles/inoculated_cheese_labels.csv|-|
|--figsize [str ...]||Figure size in inches, format: width height|-|
|--legendposition [str ...]||Adjust legend position, format: x y or x y width height|-|
|--deltact [str ...]||Delta Cq of samples measured with and without pre-amplification, give sample name strings that identifies the two groups, format: preamp id no preamp id (e.g. "nopreamp" "preamp" )|-|
|--preamp [str ...]||Include qualitative pre-amplification data in results. First element preamp id, second element no preamp id|-|
|--datatoplot [cq, copy, None]||Select data to display in the plots|-|
|--color||Use color in plots|-|
|--rawdata||Create plots with raw data (show replicates separately)|-|
|--reverse||Sort sample names in reverse order|-|
|--outfile [str]||filename for the output figure|-|
|--axeslabels||Display axes labels in heatmaps|-|

\*Make sure you use the exact same sample names and assay names as in the input csv file ([Fluidigm](https://www.fluidigm.com/) Real-Time PCR Analysis software export file)
