## Proteomics Pipeline
Code for proteomic analysis of prefrontal and hippocampal samples after chronic unpredictable stress regime to examine the effect of microglial HMGB1 and sex on the brain's response to chronic stress. The all_analyses.py file contains all code necessary for the identification and visualization of proteins significantly upregulated or downregulated between different experimental conditions. The all_analyses.py file is divided into 10 different functions.
## Description of Functions Contained in all_analyses.py
# The following 4 functions are strictly for the visualization of proteins in the form of different figures:

-create_volcano(volcano_frame, comparison, exceptions=None): Creates a volcano plot from a provided dataframe (volcano_frame) of untransformed p-values and log2(fold changes). Labels the top 15 most significant upregulated and downregulated proteins. Excludes labelling proteins provided in an iterable exceptions argument.

-create_cluster(cluster_array, group_1, group_2, comparison, names): Creates a clustermap from a provided array of log10(raw intensity) values (cluster_array). The slices of the raw data file must also be provided as the group_1 and group_2 arguments. The comparison argument is a string in the form "Group 1 v. Group 2". Lastly, names is the names of the proteins making up the rows of cluster_array.

-create_venn(file): Creates Venn diagrams from lists of proteins provided as the columns of the Excel spreadsheet file. This file should be in the same format as Venn Diagram Data Set.xlsx.

-volcano_log(a): Function that computes the base 10 logarithm of an integer a. In create_volcano.py, this function is applied to the pandas.Series containing the untransformed p-values of volcano_frame.

# The following 6 functions read raw data from raw data spreadsheets organized in the format of Franklin Male Analysis Results.xlsx and Franklin Female Analysis Results.xlsx. These sheets contain rows representing proteins and a hierarchically organized header dividing the columns representing each animal into different experimental groups. Each cell provides raw log10(intensity) values of the protein of the corresponding row for the animal of the corresponding column. Descriptions of each function can be found below:

-genetic_comparison(male_file, female_file): Identifies proteins significantly upregulated or downregulated as a result of the knockout of microglial HMGB1 for each sex and brain region. This function creates the required data structures and passes them to create_volcano() and create_cluster() for the creation of figures. The arguments male_file and female_file should be the absolute paths of Franklin Male Analysis Results.xlsx and Franklin Female Analysis Results.xlsx respectively.

-sex_comparison(male_file, female_file): Identifies proteins significantly upregulated or downregulated by stress between sexes for each brain region and genetic condition. This function creates the required data structures and passes them to create_volcano() and create_cluster() for the creation of figures. The arguments male_file and female_file should be the absolute paths of Franklin Male Analysis Results.xlsx and Franklin Female Analysis Results.xlsx respectively.

-stress_comparison(file, sex): Identifies proteins significantly upregulated or downregulated between stress and control groups for each sex, brain region, and genetic condition. This function creates the required data structures and passes them to create_volcano() and create_cluster() for the creation of figures. The arguments male_file and female_file should be the absolute paths of Franklin Male Analysis Results.xlsx and Franklin Female Analysis Results.xlsx respectively.

-resilience_comparison(male_file, female_file): Identifies proteins significantly upregulated or downregulated between WT control groups and KO stress groups for each brain region and sex to determine the effective resilience induced by the genetic knockout. This function creates the required data structures and passes them to create_volcano() and create_cluster() for the creation of figures. The arguments male_file and female_file should be the absolute paths of Franklin Male Analysis Results.xlsx and Franklin Female Analysis Results.xlsx respectively.

-get_group_vals(group, protein): Called by the four functions above to get the values corresponding to a certain protein for the provided experimental group.

-compare_groups(group_1, group_2): Also called by the same four functions to compare two groups whose values were obtained by get_group_vals(). Outputs the p-value of an independent t-test, the log2(fold change), the direction of regulation, and the values of the two groups.
# Instructions for Use
Simply press run and follow the prompts that show up in the Python console. You should also have a folder named "Figures" in the same directory location where all_analyses.py is located.
