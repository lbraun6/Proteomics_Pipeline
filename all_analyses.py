import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib_venn import *
import scipy.stats as stats
from adjustText import adjust_text
import seaborn as sns
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter

# Update font for all figures
font = {'family': 'georgia'}
plt.rc('font', **font)


def genetic_comparison(male_file, female_file):
    # MALE -------------------------------------------------------------------------------------------------------------
    # Read raw data Excel sheet into a dataframe with a multilevel header. Index column is the protein name
    male_raw_intensities = pd.read_excel(male_file, "Raw log10 intensities", header=[0, 1, 2, 3], index_col=3)
    male_raw_intensities.replace("Missing Value", np.nan, inplace=True)

    # Assign each control and stress group subframe to its own variable (includes all animals within the group as
    # columns)
    male_pfc_ko_ctl_protein_expression = male_raw_intensities["PFC", "HMGB1 KO", "Control"]
    male_pfc_wt_ctl_protein_expression = male_raw_intensities["PFC", "WT", "Control"]
    male_hip_ko_ctl_protein_expression = male_raw_intensities["HIP", "HMGB1 KO", "Control"]
    male_hip_wt_ctl_protein_expression = male_raw_intensities["HIP", "WT", "Control"]

    # FEMALE -----------------------------------------------------------------------------------------------------------
    # Read raw data Excel sheet into a dataframe with a multilevel header. Index column is the protein name
    female_raw_intensities = pd.read_excel(female_file, "Raw log10 intensities", header=[0, 1, 2, 3], index_col=3)
    female_raw_intensities.replace("Missing Value", np.nan, inplace=True)

    # Assign each control and stress group subframe to its own variable (includes all animals within the group as
    # columns)
    female_pfc_ko_ctl_protein_expression = female_raw_intensities["PFC", "HMGB1 KO", "Control"]
    female_pfc_wt_ctl_protein_expression = female_raw_intensities["PFC", "WT", "Control"]
    female_hip_ko_ctl_protein_expression = female_raw_intensities["HIP", "HMGB1 KO", "Control"]
    female_hip_wt_ctl_protein_expression = female_raw_intensities["HIP", "WT", "Control"]

    # Initialize lists to hold proteins identified as significant for at least one sex as well as the direction of their
    # regulation and their fold change
    id_male_pfc = []
    id_female_pfc = []
    id_male_hip = []
    id_female_hip = []

    reg_male_pfc = []
    reg_female_pfc = []
    reg_male_hip = []
    reg_female_hip = []

    fc_male_pfc = []
    fc_female_pfc = []
    fc_male_hip = []
    fc_female_hip = []

    p_male_pfc = []
    p_female_pfc = []
    p_male_hip = []
    p_female_hip = []

    male_pfc_cluster = []
    female_pfc_cluster = []
    male_hip_cluster = []
    female_hip_cluster = []

    male_pfc_names = []
    female_pfc_names = []
    male_hip_names = []
    female_hip_names = []

    for protein in male_raw_intensities.index:
        male_pfc = compare_groups(get_group_vals(male_pfc_ko_ctl_protein_expression, protein),
                                  get_group_vals(male_pfc_wt_ctl_protein_expression, protein))
        if isinstance(male_pfc, list):
            id_male_pfc += [protein]
            p_male_pfc += [male_pfc[0]]
            reg_male_pfc += [male_pfc[1]]
            fc_male_pfc += [male_pfc[2]]
            if male_pfc[0] <= 0.05:
                list_1 = male_pfc[3]
                list_2 = male_pfc[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                male_pfc_cluster += [list_1 + list_2]
                male_pfc_names += [protein]

        male_hip = compare_groups(get_group_vals(male_hip_ko_ctl_protein_expression, protein),
                                  get_group_vals(male_hip_wt_ctl_protein_expression, protein))
        if isinstance(male_hip, list):
            id_male_hip += [protein]
            p_male_hip += [male_hip[0]]
            reg_male_hip += [male_hip[1]]
            fc_male_hip += [male_hip[2]]
            if male_hip[0] <= 0.05:
                list_1 = male_hip[3]
                list_2 = male_hip[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                male_hip_cluster += [list_1 + list_2]
                male_hip_names += [protein]

    for protein in female_raw_intensities.index:
        female_pfc = compare_groups(get_group_vals(female_pfc_ko_ctl_protein_expression, protein),
                                    get_group_vals(female_pfc_wt_ctl_protein_expression, protein))
        if isinstance(female_pfc, list):
            id_female_pfc += [protein]
            p_female_pfc += [female_pfc[0]]
            reg_female_pfc += [female_pfc[1]]
            fc_female_pfc += [female_pfc[2]]
            if female_pfc[0] <= 0.05:
                list_1 = female_pfc[3]
                list_2 = female_pfc[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                female_pfc_cluster += [list_1 + list_2]
                female_pfc_names += [protein]

        female_hip = compare_groups(get_group_vals(female_hip_ko_ctl_protein_expression, protein),
                                    get_group_vals(female_hip_wt_ctl_protein_expression, protein))
        if isinstance(female_hip, list):
            id_female_hip += [protein]
            p_female_hip += [female_hip[0]]
            reg_female_hip += [female_hip[1]]
            fc_female_hip += [female_hip[2]]
            if female_hip[0] <= 0.05:
                list_1 = female_hip[3]
                list_2 = female_hip[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                female_hip_cluster += [list_1 + list_2]
                female_hip_names += [protein]

    male_pfc_volcano = pd.DataFrame({"Protein": id_male_pfc, "Fold Change": fc_male_pfc, "p-Value": p_male_pfc})
    female_pfc_volcano = pd.DataFrame({"Protein": id_female_pfc, "Fold Change": fc_female_pfc, "p-Value": p_female_pfc})
    male_hip_volcano = pd.DataFrame({"Protein": id_male_hip, "Fold Change": fc_male_hip, "p-Value": p_male_hip})
    female_hip_volcano = pd.DataFrame({"Protein": id_female_hip, "Fold Change": fc_female_hip, "p-Value": p_female_hip})

    create_volcano(male_pfc_volcano, "Male PFC KO Control v. Male PFC WT Control")
    create_cluster(male_pfc_cluster, male_pfc_ko_ctl_protein_expression, male_pfc_wt_ctl_protein_expression,
                   "Male PFC KO Control v. Male PFC WT Control", male_pfc_names)

    create_volcano(female_pfc_volcano, "Female PFC KO Control v. Female PFC WT Control")
    create_cluster(female_pfc_cluster, female_pfc_ko_ctl_protein_expression, female_pfc_wt_ctl_protein_expression,
                   "Female PFC KO Control v. Female PFC WT Control", female_pfc_names)

    create_volcano(male_hip_volcano, "Male HIP KO Control v. Male HIP WT Control")
    create_cluster(male_hip_cluster, male_hip_ko_ctl_protein_expression, male_hip_wt_ctl_protein_expression,
                   "Male HIP KO Control v. Male HIP WT Control", male_hip_names)

    create_volcano(female_hip_volcano, "Female HIP KO Control v. Male HIP WT Control")
    create_cluster(female_hip_cluster, female_hip_ko_ctl_protein_expression, female_hip_wt_ctl_protein_expression,
                   "Female HIP KO Control v. Female HIP WT Control", female_hip_names)


def sex_comparison(male_file, female_file, venn_file):
    # MALE -------------------------------------------------------------------------------------------------------------
    # Read raw data Excel sheet into a dataframe with a multilevel header. Index column is the protein name
    male_raw_intensities = pd.read_excel(male_file, "Raw log10 intensities", header=[0, 1, 2, 3], index_col=3)
    male_raw_intensities.replace("Missing Value", np.nan, inplace=True)

    # Assign each control and stress group subframe to its own variable (includes all animals within the group as
    # columns)
    male_pfc_ko_ctl_protein_expression = male_raw_intensities["PFC", "HMGB1 KO", "Control"]
    male_pfc_wt_ctl_protein_expression = male_raw_intensities["PFC", "WT", "Control"]
    male_hip_ko_ctl_protein_expression = male_raw_intensities["HIP", "HMGB1 KO", "Control"]
    male_hip_wt_ctl_protein_expression = male_raw_intensities["HIP", "WT", "Control"]

    male_pfc_ko_stress_protein_expression = male_raw_intensities["PFC", "HMGB1 KO", "Stress"]
    male_pfc_wt_stress_protein_expression = male_raw_intensities["PFC", "WT", "Stress"]
    male_hip_ko_stress_protein_expression = male_raw_intensities["HIP", "HMGB1 KO", "Stress"]
    male_hip_wt_stress_protein_expression = male_raw_intensities["HIP", "WT", "Stress"]

    # FEMALE -----------------------------------------------------------------------------------------------------------
    # Read raw data Excel sheet into a dataframe with a multilevel header. Index column is the protein name
    female_raw_intensities = pd.read_excel(female_file, "Raw log10 intensities", header=[0, 1, 2, 3], index_col=3)
    female_raw_intensities.replace("Missing Value", np.nan, inplace=True)

    # Assign each control and stress group subframe to its own variable (includes all animals within the group as
    # columns)
    female_pfc_ko_ctl_protein_expression = female_raw_intensities["PFC", "HMGB1 KO", "Control"]
    female_pfc_wt_ctl_protein_expression = female_raw_intensities["PFC", "WT", "Control"]
    female_hip_ko_ctl_protein_expression = female_raw_intensities["HIP", "HMGB1 KO", "Control"]
    female_hip_wt_ctl_protein_expression = female_raw_intensities["HIP", "WT", "Control"]

    female_pfc_ko_stress_protein_expression = female_raw_intensities["PFC", "HMGB1 KO", "Stress"]
    female_pfc_wt_stress_protein_expression = female_raw_intensities["PFC", "WT", "Stress"]
    female_hip_ko_stress_protein_expression = female_raw_intensities["HIP", "HMGB1 KO", "Stress"]
    female_hip_wt_stress_protein_expression = female_raw_intensities["HIP", "WT", "Stress"]

    # Initialize lists to hold proteins identified as significant for at least one sex as well as the direction of their
    # regulation and their fold change
    id_pfc_ko = []
    id_pfc_wt = []
    id_hip_ko = []
    id_hip_wt = []

    reg_pfc_ko = []
    reg_pfc_wt = []
    reg_hip_ko = []
    reg_hip_wt = []

    fc_pfc_ko = []
    fc_pfc_wt = []
    fc_hip_ko = []
    fc_hip_wt = []

    p_pfc_ko = []
    p_pfc_wt = []
    p_hip_ko = []
    p_hip_wt = []

    pfc_ko_cluster = []
    pfc_wt_cluster = []
    hip_ko_cluster = []
    hip_wt_cluster = []

    pfc_ko_names = []
    pfc_wt_names = []
    hip_ko_names = []
    hip_wt_names = []

    # Iterate through every protein in the raw data file, select the values within matching control and stress groups
    # for each sex, compare these groups to determine if there are significant differences in protein expression for at
    # least one sex. If there are, perform a comparison between the stress conditions of the sexes and store the protein
    # name, direction of regulation, and fold change
    for protein in male_raw_intensities.index:
        if protein in female_raw_intensities.index:
            # PFC KO
            male_pfc_ko = compare_groups(get_group_vals(male_pfc_ko_stress_protein_expression, protein),
                                         get_group_vals(male_pfc_ko_ctl_protein_expression, protein))
            female_pfc_ko = compare_groups(get_group_vals(female_pfc_ko_stress_protein_expression, protein),
                                           get_group_vals(female_pfc_ko_ctl_protein_expression, protein))

            if isinstance(male_pfc_ko, list) and isinstance(female_pfc_ko, list):
                if male_pfc_ko[0] <= 0.05 or female_pfc_ko[0] <= 0.05:
                    pfc_ko = compare_groups(get_group_vals(male_pfc_ko_stress_protein_expression, protein),
                                            get_group_vals(female_pfc_ko_stress_protein_expression, protein))
                    if isinstance(pfc_ko, list):
                        id_pfc_ko += [protein]
                        p_pfc_ko += [pfc_ko[0]]
                        reg_pfc_ko += [pfc_ko[1]]
                        fc_pfc_ko += [pfc_ko[2]]
                        if pfc_ko[0] <= 0.05:
                            list_1 = pfc_ko[3]
                            list_2 = pfc_ko[4]
                            pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                            pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                            pfc_ko_cluster += [list_1 + list_2]
                            pfc_ko_names += [protein]

            # PFC WT
            male_pfc_wt = compare_groups(get_group_vals(male_pfc_wt_stress_protein_expression, protein),
                                         get_group_vals(male_pfc_wt_ctl_protein_expression, protein))
            female_pfc_wt = compare_groups(get_group_vals(female_pfc_wt_stress_protein_expression, protein),
                                           get_group_vals(female_pfc_wt_ctl_protein_expression, protein))

            if isinstance(male_pfc_wt, list) and isinstance(female_pfc_wt, list):
                if male_pfc_wt[0] <= 0.05 or female_pfc_wt[0] <= 0.05:
                    pfc_wt = compare_groups(get_group_vals(male_pfc_wt_stress_protein_expression, protein),
                                            get_group_vals(female_pfc_wt_stress_protein_expression, protein))
                    if isinstance(pfc_wt, list):
                        id_pfc_wt += [protein]
                        p_pfc_wt += [pfc_wt[0]]
                        reg_pfc_wt += [pfc_wt[1]]
                        fc_pfc_wt += [pfc_wt[2]]
                        if pfc_wt[0] <= 0.05:
                            list_1 = pfc_wt[3]
                            list_2 = pfc_wt[4]
                            pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                            pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                            pfc_wt_cluster += [list_1 + list_2]
                            pfc_wt_names += [protein]

            # HIP KO
            male_hip_ko = compare_groups(get_group_vals(male_hip_ko_stress_protein_expression, protein),
                                         get_group_vals(male_hip_ko_ctl_protein_expression, protein))
            female_hip_ko = compare_groups(get_group_vals(female_hip_ko_stress_protein_expression, protein),
                                           get_group_vals(female_hip_ko_ctl_protein_expression, protein))

            if isinstance(male_hip_ko, list) and isinstance(female_hip_ko, list):
                if male_hip_ko[0] <= 0.05 or female_hip_ko[0] <= 0.05:
                    hip_ko = compare_groups(get_group_vals(male_hip_ko_stress_protein_expression, protein),
                                            get_group_vals(female_hip_ko_stress_protein_expression, protein))
                    if isinstance(hip_ko, list):
                        id_hip_ko += [protein]
                        p_hip_ko += [hip_ko[0]]
                        reg_hip_ko += [hip_ko[1]]
                        fc_hip_ko += [hip_ko[2]]
                        if hip_ko[0] <= 0.05:
                            list_1 = hip_ko[3]
                            list_2 = hip_ko[4]
                            pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                            pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                            hip_ko_cluster += [list_1 + list_2]
                            hip_ko_names += [protein]

            # HIP WT
            male_hip_wt = compare_groups(get_group_vals(male_hip_wt_stress_protein_expression, protein),
                                         get_group_vals(male_hip_wt_ctl_protein_expression, protein))
            female_hip_wt = compare_groups(get_group_vals(female_hip_wt_stress_protein_expression, protein),
                                           get_group_vals(female_hip_wt_ctl_protein_expression, protein))

            if isinstance(male_hip_wt, list) and isinstance(female_hip_wt, list):
                if male_hip_wt[0] <= 0.05 or female_hip_wt[0] <= 0.05:
                    hip_wt = compare_groups(get_group_vals(male_hip_wt_stress_protein_expression, protein),
                                            get_group_vals(female_hip_wt_stress_protein_expression, protein))
                    if isinstance(hip_wt, list):
                        id_hip_wt += [protein]
                        p_hip_wt += [hip_wt[0]]
                        reg_hip_wt += [hip_wt[1]]
                        fc_hip_wt += [hip_wt[2]]
                        if hip_wt[0] <= 0.05:
                            list_1 = hip_wt[3]
                            list_2 = hip_wt[4]
                            pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                            pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                            hip_wt_cluster += [list_1 + list_2]
                            hip_wt_names += [protein]

    pfc_wt_volcano = pd.DataFrame({"Protein": id_pfc_wt, "Fold Change": fc_pfc_wt,
                                   "p-Value": p_pfc_wt})
    pfc_ko_volcano = pd.DataFrame({"Protein": id_pfc_ko, "Fold Change": fc_pfc_ko,
                                   "p-Value": p_pfc_ko})
    hip_wt_volcano = pd.DataFrame({"Protein": id_hip_wt, "Fold Change": fc_hip_wt,
                                   "p-Value": p_hip_wt})
    hip_ko_volcano = pd.DataFrame({"Protein": id_hip_ko, "Fold Change": fc_hip_ko,
                                   "p-Value": p_hip_ko})

    create_volcano(pfc_wt_volcano, "Male PFC WT Stress v. Female PFC WT Stress")
    create_cluster(pfc_wt_cluster, male_pfc_wt_stress_protein_expression, female_pfc_wt_stress_protein_expression,
                   "Male PFC WT Stress v. Female PFC WT Stress", pfc_wt_names)
    create_venn(venn_file, "Male PFC WT v. Female PFC WT")

    create_volcano(pfc_ko_volcano, "Male PFC KO Stress v. Female PFC KO Stress")
    create_cluster(pfc_ko_cluster, male_pfc_ko_stress_protein_expression, female_pfc_ko_stress_protein_expression,
                   "Male PFC KO Stress v. Female PFC KO Stress", pfc_ko_names)
    create_venn(venn_file, "Male PFC KO v. Female PFC KO")

    create_volcano(hip_wt_volcano, "Male HIP WT Stress v. Female HIP WT Stress")
    create_cluster(hip_wt_cluster, male_hip_wt_stress_protein_expression, female_hip_wt_stress_protein_expression,
                   "Male HIP WT Stress v. Female HIP WT Stress", hip_wt_names)
    create_venn(venn_file, "Male HIP WT v. Female HIP WT")

    create_volcano(hip_ko_volcano, "Male HIP KO Stress v. Female HIP KO Stress")
    create_cluster(hip_ko_cluster, male_hip_ko_stress_protein_expression, female_hip_ko_stress_protein_expression,
                   "Male HIP KO Stress v. Female HIP KO Stress", hip_ko_names)
    create_venn(venn_file, "Male HIP KO v. Female HIP KO")

    # with pd.ExcelWriter(male_file[:male_file.rfind("\\")] + "\\" + "Male v. Female Volcano Plot Dataset.xlsx") as writer:
    #     pfc_wt_volcano.to_excel(writer, sheet_name="PFC WT", index=False)
    #     pfc_ko_volcano.to_excel(writer, sheet_name="PFC KO", index=False)
    #     hip_wt_volcano.to_excel(writer, sheet_name="HIP WT", index=False)
    #     hip_ko_volcano.to_excel(writer, sheet_name="HIP KO", index=False)


def stress_comparison(file, sex, venn_file):
    # Read raw data Excel sheet into a dataframe with a multilevel header. Index column is the protein name
    raw_intensities = pd.read_excel(file, "Raw log10 intensities", header=[0, 1, 2, 3], index_col=3)
    raw_intensities.replace("Missing Value", np.nan, inplace=True)

    # Assign each control and stress group subframe to its own variable (includes all animals within the group as
    # columns)
    pfc_ko_ctl_protein_expression = raw_intensities["PFC", "HMGB1 KO", "Control"]
    pfc_ko_stress_protein_expression = raw_intensities["PFC", "HMGB1 KO", "Stress"]

    pfc_wt_ctl_protein_expression = raw_intensities["PFC", "WT", "Control"]
    pfc_wt_stress_protein_expression = raw_intensities["PFC", "WT", "Stress"]

    hip_ko_ctl_protein_expression = raw_intensities["HIP", "HMGB1 KO", "Control"]
    hip_ko_stress_protein_expression = raw_intensities["HIP", "HMGB1 KO", "Stress"]

    hip_wt_ctl_protein_expression = raw_intensities["HIP", "WT", "Control"]
    hip_wt_stress_protein_expression = raw_intensities["HIP", "WT", "Stress"]

    # Initialize lists to hold the analyzed proteins as well as the direction of their regulation and their fold change
    id_pfc_ko = []
    id_pfc_wt = []
    id_hip_ko = []
    id_hip_wt = []

    reg_pfc_ko = []
    reg_pfc_wt = []
    reg_hip_ko = []
    reg_hip_wt = []

    fc_pfc_ko = []
    fc_pfc_wt = []
    fc_hip_ko = []
    fc_hip_wt = []

    p_pfc_ko = []
    p_pfc_wt = []
    p_hip_ko = []
    p_hip_wt = []

    pfc_ko_cluster = []
    pfc_wt_cluster = []
    hip_ko_cluster = []
    hip_wt_cluster = []

    pfc_ko_names = []
    pfc_wt_names = []
    hip_ko_names = []
    hip_wt_names = []

    # Initialize lists for scatter plots
    scatter_id_pfc = []
    scatter_id_hip = []

    scatter_pfc_wt = []
    scatter_pfc_ko = []
    scatter_hip_wt = []
    scatter_hip_ko = []

    # Iterate through every protein in the raw data file, select the values within matching control and stress groups,
    # compare these groups to calculate p-values and store the protein name, direction of regulation, and fold change
    for protein in raw_intensities.index:
        pfc_ko = compare_groups(get_group_vals(pfc_ko_stress_protein_expression, protein),
                                get_group_vals(pfc_ko_ctl_protein_expression, protein))
        if isinstance(pfc_ko, list):
            id_pfc_ko += [protein]
            p_pfc_ko += [pfc_ko[0]]
            reg_pfc_ko += [pfc_ko[1]]
            fc_pfc_ko += [pfc_ko[2]]
            if pfc_ko[0] <= 0.05:
                list_1 = pfc_ko[3]
                list_2 = pfc_ko[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                pfc_ko_cluster += [list_1 + list_2]
                pfc_ko_names += [protein]

        pfc_wt = compare_groups(get_group_vals(pfc_wt_stress_protein_expression, protein),
                                get_group_vals(pfc_wt_ctl_protein_expression, protein))
        if isinstance(pfc_wt, list):
            id_pfc_wt += [protein]
            p_pfc_wt += [pfc_wt[0]]
            reg_pfc_wt += [pfc_wt[1]]
            fc_pfc_wt += [pfc_wt[2]]
            if pfc_wt[0] <= 0.05:
                list_1 = pfc_wt[3]
                list_2 = pfc_wt[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                pfc_wt_cluster += [list_1 + list_2]
                pfc_wt_names += [protein]

        hip_ko = compare_groups(get_group_vals(hip_ko_stress_protein_expression, protein),
                                get_group_vals(hip_ko_ctl_protein_expression, protein))
        if isinstance(hip_ko, list):
            id_hip_ko += [protein]
            p_hip_ko += [hip_ko[0]]
            reg_hip_ko += [hip_ko[1]]
            fc_hip_ko += [hip_ko[2]]
            if hip_ko[0] <= 0.05:
                list_1 = hip_ko[3]
                list_2 = hip_ko[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                hip_ko_cluster += [list_1 + list_2]
                hip_ko_names += [protein]

        hip_wt = compare_groups(get_group_vals(hip_wt_stress_protein_expression, protein),
                                get_group_vals(hip_wt_ctl_protein_expression, protein))
        if isinstance(hip_wt, list):
            id_hip_wt += [protein]
            p_hip_wt += [hip_wt[0]]
            reg_hip_wt += [hip_wt[1]]
            fc_hip_wt += [hip_wt[2]]
            if hip_wt[0] <= 0.05:
                list_1 = hip_wt[3]
                list_2 = hip_wt[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                hip_wt_cluster += [list_1 + list_2]
                hip_wt_names += [protein]

        # For Scatter Plots
        if isinstance(pfc_wt, list) and isinstance(pfc_ko, list):
            scatter_id_pfc += [protein]
            scatter_pfc_wt += [pfc_wt[2]]
            scatter_pfc_ko += [pfc_ko[2]]

        if isinstance(hip_wt, list) and isinstance(hip_ko, list):
            scatter_id_hip += [protein]
            scatter_hip_wt += [hip_wt[2]]
            scatter_hip_ko += [hip_ko[2]]

    pfc_wt_volcano = pd.DataFrame({"Protein": id_pfc_wt, "Fold Change": fc_pfc_wt, "p-Value": p_pfc_wt})
    pfc_ko_volcano = pd.DataFrame({"Protein": id_pfc_ko, "Fold Change": fc_pfc_ko, "p-Value": p_pfc_ko})
    hip_wt_volcano = pd.DataFrame({"Protein": id_hip_wt, "Fold Change": fc_hip_wt, "p-Value": p_hip_wt})
    hip_ko_volcano = pd.DataFrame({"Protein": id_hip_ko, "Fold Change": fc_hip_ko, "p-Value": p_hip_ko})

    pfc_scatter = pd.DataFrame({"Protein": scatter_id_pfc, "WT Fold Change": scatter_pfc_wt,
                                "KO Fold Change": scatter_pfc_ko})
    hip_scatter = pd.DataFrame({"Protein": scatter_id_hip, "WT Fold Change": scatter_hip_wt,
                                "KO Fold Change": scatter_hip_ko})

    if sex == "Male":
        create_scatter(pfc_scatter, "Male PFC WT v. Male PFC KO")
        create_scatter(hip_scatter, "Male HIP WT v. Male HIP KO")

        create_venn(venn_file, "Male PFC WT v. Male PFC KO")
        create_venn(venn_file, "Male HIP WT v. Male HIP KO")

        create_volcano(pfc_wt_volcano, "Male PFC WT Stress v. Male PFC WT Control", exceptions=["Camkk2"])
        create_cluster(pfc_wt_cluster, pfc_wt_stress_protein_expression, pfc_wt_ctl_protein_expression,
                       "Male PFC WT Stress v. Male PFC WT Control", pfc_wt_names)

        create_volcano(pfc_ko_volcano, "Male PFC KO Stress v. Male PFC KO Control")
        create_cluster(pfc_ko_cluster, pfc_ko_stress_protein_expression, pfc_ko_ctl_protein_expression,
                       "Male PFC KO Stress v. Male PFC KO Control", pfc_ko_names)

        create_volcano(hip_wt_volcano, "Male HIP WT Stress v. Male HIP WT Control")
        create_cluster(hip_wt_cluster, hip_wt_stress_protein_expression, hip_wt_ctl_protein_expression,
                       "Male HIP WT Stress v. Male HIP WT Control", hip_wt_names)

        create_volcano(hip_ko_volcano, "Male HIP KO Stress v. Male HIP KO Control")
        create_cluster(hip_ko_cluster, hip_ko_stress_protein_expression, hip_ko_ctl_protein_expression,
                       "Male HIP KO Stress v. Male HIP KO Control", hip_ko_names)

    elif sex == "Female":
        create_scatter(pfc_scatter, "Female PFC WT v. Female PFC KO")
        create_scatter(hip_scatter, "Female HIP WT v. Female HIP KO")

        create_venn(venn_file, "Female PFC WT v. Female PFC KO")
        create_venn(venn_file, "Female HIP WT v. Female HIP KO")

        create_volcano(pfc_wt_volcano, "Female PFC WT Stress v. Female PFC WT Control")
        create_cluster(pfc_wt_cluster, pfc_wt_stress_protein_expression, pfc_wt_ctl_protein_expression,
                       "Female PFC WT Stress v. Female PFC WT Control", pfc_wt_names)

        create_volcano(pfc_ko_volcano, "Female PFC KO Stress v. Female PFC KO Control")
        create_cluster(pfc_ko_cluster, pfc_ko_stress_protein_expression, pfc_ko_ctl_protein_expression,
                       "Female PFC KO Stress v. Female PFC KO Control", pfc_ko_names)

        create_volcano(hip_wt_volcano, "Female HIP WT Stress v. Female HIP WT Control")
        create_cluster(hip_wt_cluster, hip_wt_stress_protein_expression, hip_wt_ctl_protein_expression,
                       "Female HIP WT Stress v. Female HIP WT Control", hip_wt_names)

        create_volcano(hip_ko_volcano, "Female HIP KO Stress v. Female HIP KO Control")
        create_cluster(hip_ko_cluster, hip_ko_stress_protein_expression, hip_ko_ctl_protein_expression,
                       "Female HIP KO Stress v. Female HIP KO Control", hip_ko_names)

    # with pd.ExcelWriter(file[:file.rfind("\\")] + "\\" + sex + " Volcano Plot Dataset.xlsx") as writer:
    #     pfc_wt_volcano.to_excel(writer, sheet_name="PFC WT", index=False)
    #     pfc_ko_volcano.to_excel(writer, sheet_name="PFC KO", index=False)
    #     hip_wt_volcano.to_excel(writer, sheet_name="HIP WT", index=False)
    #     hip_ko_volcano.to_excel(writer, sheet_name="HIP KO", index=False)


def resilience_comparison(male_file, female_file):
    # MALE -------------------------------------------------------------------------------------------------------------
    # Read raw data Excel sheet into a dataframe with a multilevel header. Index column is the protein name
    male_raw_intensities = pd.read_excel(male_file, "Raw log10 intensities", header=[0, 1, 2, 3], index_col=3)
    male_raw_intensities.replace("Missing Value", np.nan, inplace=True)

    # Assign each control and stress group subframe to its own variable (includes all animals within the group as
    # columns)
    male_pfc_ko_stress_protein_expression = male_raw_intensities["PFC", "HMGB1 KO", "Stress"]
    male_pfc_wt_ctl_protein_expression = male_raw_intensities["PFC", "WT", "Control"]
    male_hip_ko_stress_protein_expression = male_raw_intensities["HIP", "HMGB1 KO", "Stress"]
    male_hip_wt_ctl_protein_expression = male_raw_intensities["HIP", "WT", "Control"]

    # FEMALE -----------------------------------------------------------------------------------------------------------
    # Read raw data Excel sheet into a dataframe with a multilevel header. Index column is the protein name
    female_raw_intensities = pd.read_excel(female_file, "Raw log10 intensities", header=[0, 1, 2, 3], index_col=3)
    female_raw_intensities.replace("Missing Value", np.nan, inplace=True)

    # Assign each control and stress group subframe to its own variable (includes all animals within the group as
    # columns)
    female_pfc_ko_stress_protein_expression = female_raw_intensities["PFC", "HMGB1 KO", "Stress"]
    female_pfc_wt_ctl_protein_expression = female_raw_intensities["PFC", "WT", "Control"]
    female_hip_ko_stress_protein_expression = female_raw_intensities["HIP", "HMGB1 KO", "Stress"]
    female_hip_wt_ctl_protein_expression = female_raw_intensities["HIP", "WT", "Control"]

    # Initialize lists to hold proteins identified as significant for at least one sex as well as the direction of their
    # regulation and their fold change
    id_male_pfc = []
    id_female_pfc = []
    id_male_hip = []
    id_female_hip = []

    reg_male_pfc = []
    reg_female_pfc = []
    reg_male_hip = []
    reg_female_hip = []

    fc_male_pfc = []
    fc_female_pfc = []
    fc_male_hip = []
    fc_female_hip = []

    p_male_pfc = []
    p_female_pfc = []
    p_male_hip = []
    p_female_hip = []

    male_pfc_cluster = []
    female_pfc_cluster = []
    male_hip_cluster = []
    female_hip_cluster = []

    male_pfc_names = []
    female_pfc_names = []
    male_hip_names = []
    female_hip_names = []

    for protein in male_raw_intensities.index:
        male_pfc = compare_groups(get_group_vals(male_pfc_ko_stress_protein_expression, protein),
                                  get_group_vals(male_pfc_wt_ctl_protein_expression, protein))
        if isinstance(male_pfc, list):
            id_male_pfc += [protein]
            p_male_pfc += [male_pfc[0]]
            reg_male_pfc += [male_pfc[1]]
            fc_male_pfc += [male_pfc[2]]
            if male_pfc[0] <= 0.05:
                list_1 = male_pfc[3]
                list_2 = male_pfc[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                male_pfc_cluster += [list_1 + list_2]
                male_pfc_names += [protein]

        male_hip = compare_groups(get_group_vals(male_hip_ko_stress_protein_expression, protein),
                                  get_group_vals(male_hip_wt_ctl_protein_expression, protein))
        if isinstance(male_hip, list):
            id_male_hip += [protein]
            p_male_hip += [male_hip[0]]
            reg_male_hip += [male_hip[1]]
            fc_male_hip += [male_hip[2]]
            if male_hip[0] <= 0.05:
                list_1 = male_hip[3]
                list_2 = male_hip[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                male_hip_cluster += [list_1 + list_2]
                male_hip_names += [protein]

    for protein in female_raw_intensities.index:
        female_pfc = compare_groups(get_group_vals(female_pfc_ko_stress_protein_expression, protein),
                                    get_group_vals(female_pfc_wt_ctl_protein_expression, protein))
        if isinstance(female_pfc, list):
            id_female_pfc += [protein]
            p_female_pfc += [female_pfc[0]]
            reg_female_pfc += [female_pfc[1]]
            fc_female_pfc += [female_pfc[2]]
            if female_pfc[0] <= 0.05:
                list_1 = female_pfc[3]
                list_2 = female_pfc[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                female_pfc_cluster += [list_1 + list_2]
                female_pfc_names += [protein]

        female_hip = compare_groups(get_group_vals(female_hip_ko_stress_protein_expression, protein),
                                    get_group_vals(female_hip_wt_ctl_protein_expression, protein))
        if isinstance(female_hip, list):
            id_female_hip += [protein]
            p_female_hip += [female_hip[0]]
            reg_female_hip += [female_hip[1]]
            fc_female_hip += [female_hip[2]]
            if female_hip[0] <= 0.05:
                list_1 = female_hip[3]
                list_2 = female_hip[4]
                pd.Series(list_1).fillna(np.nanmean(list_1)).tolist()
                pd.Series(list_2).fillna(np.nanmean(list_2)).tolist()
                female_hip_cluster += [list_1 + list_2]
                female_hip_names += [protein]

    male_pfc_volcano = pd.DataFrame({"Protein": id_male_pfc, "Fold Change": fc_male_pfc, "p-Value": p_male_pfc})
    female_pfc_volcano = pd.DataFrame({"Protein": id_female_pfc, "Fold Change": fc_female_pfc, "p-Value": p_female_pfc})
    male_hip_volcano = pd.DataFrame({"Protein": id_male_hip, "Fold Change": fc_male_hip, "p-Value": p_male_hip})
    female_hip_volcano = pd.DataFrame({"Protein": id_female_hip, "Fold Change": fc_female_hip, "p-Value": p_female_hip})

    create_volcano(male_pfc_volcano, "Male PFC KO Stress v. Male PFC WT Control")
    create_cluster(male_pfc_cluster, male_pfc_ko_stress_protein_expression, male_pfc_wt_ctl_protein_expression,
                   "Male PFC KO Stress v. Male PFC WT Control", male_pfc_names)

    create_volcano(female_pfc_volcano, "Female PFC KO Stress v. Female PFC WT Control")
    create_cluster(female_pfc_cluster, female_pfc_ko_stress_protein_expression, female_pfc_wt_ctl_protein_expression,
                   "Female PFC KO Stress v. Female PFC WT Control", female_pfc_names)

    create_volcano(male_hip_volcano, "Male HIP KO Stress v. Male HIP WT Control")
    create_cluster(male_hip_cluster, male_hip_ko_stress_protein_expression, male_hip_wt_ctl_protein_expression,
                   "Male HIP KO Stress v. Male HIP WT Control", male_hip_names)

    create_volcano(female_hip_volcano, "Female HIP KO Stress v. Female HIP WT Control")
    create_cluster(female_hip_cluster, female_hip_ko_stress_protein_expression, female_hip_wt_ctl_protein_expression,
                   "Female HIP KO Stress v. Female HIP WT Control", female_hip_names)


def get_group_vals(group, protein):
    vals = []
    for x in group.columns:
        vals += [group[x][protein]]
    return vals


def compare_groups(group_1, group_2):
    stat, p = stats.ttest_ind(group_1, group_2)
    group_1_raw = []
    group_2_raw = []

    for x in group_1:
        group_1_raw += [10 ** x]
    for x in group_2:
        group_2_raw += [10 ** x]

    if np.mean(group_1) > np.mean(group_2):
        fc = np.log2(stats.mstats.gmean(group_1_raw) / stats.mstats.gmean(group_2_raw))
        return [p, "up", fc, group_1, group_2]

    elif np.mean(group_1) < np.mean(group_2):
        fc = np.log2(stats.mstats.gmean(group_1_raw) / stats.mstats.gmean(group_2_raw))
        return [p, "down", fc, group_1, group_2]

    else:
        return None


def create_volcano(volcano_frame, comparison, exceptions=None):
    # Font Information: family=georgia, size=16 (annotations=14), color=black
    volcano_frame["p-Value"] = volcano_frame["p-Value"].apply(volcano_log)
    up_reg = {"x": [], "y": [], "id": []}
    down_reg = {"x": [], "y": [], "id": []}
    un_reg = {"x": [], "y": [], "id": []}
    fig, ax = plt.subplots()
    fig.set_size_inches(7, 7)
    plotted_up = {}
    plotted_down = {}

    for x in range(len(volcano_frame.index)):

        # Adjust protein ID so it includes just protein name
        protein_id = str(volcano_frame["Protein"][x])
        if protein_id.find("GN=") >= 0:
            protein_id = protein_id.split()
            protein_id = protein_id[-3]
            protein_id = protein_id[3:]
        elif protein_id.find("|") >= 0:
            protein_id = protein_id.split("|")[2]
            protein_id = protein_id.split("_")[0]

        # Record location of HMGB1 separate from other proteins
        if protein_id == "Hmgb1":
            hmgb1_x = volcano_frame["Fold Change"][x]
            hmgb1_y = volcano_frame["p-Value"][x]

        # Record locations of all significantly up-regulated proteins
        if volcano_frame["p-Value"][x] >= (-1 * np.log10(0.05)) and volcano_frame["Fold Change"][x] > 0:
            up_reg["x"] += [volcano_frame["Fold Change"][x]]
            up_reg["y"] += [volcano_frame["p-Value"][x]]
            up_reg["id"] += [protein_id]

        # Record locations of all significantly down-regulated proteins
        elif volcano_frame["p-Value"][x] >= (-1 * np.log10(0.05)) and volcano_frame["Fold Change"][x] < 0:
            down_reg["x"] += [volcano_frame["Fold Change"][x]]
            down_reg["y"] += [volcano_frame["p-Value"][x]]
            down_reg["id"] += [protein_id]

        # Record locations of all non-significantly regulated proteins
        elif volcano_frame["p-Value"][x] < (-1 * np.log10(0.05)):
            un_reg["x"] += [volcano_frame["Fold Change"][x]]
            un_reg["y"] += [volcano_frame["p-Value"][x]]

    # Eliminate isoforms
    up_reg = pd.DataFrame(up_reg)
    down_reg = pd.DataFrame(down_reg)
    up_reg.drop_duplicates(subset=["id"], inplace=True)
    down_reg.drop_duplicates(subset=["id"], inplace=True)

    # Plot significantly up-regulated proteins in blue with annotations
    sizes = [50] * len(up_reg["x"])
    ax.scatter(up_reg["x"], up_reg["y"], sizes, c="b", marker=".")
    up_len = len(up_reg["x"])
    if len(up_reg["x"]) > 10:
        up_reg = up_reg.nlargest(10, "y", "all")
    up_updated = False
    up_reg = up_reg.nsmallest(10, "x", "all")
    for x in up_reg.index:
        if not up_reg["id"][x] in plotted_up.keys():
            plotted_up[up_reg["id"][x]] = ax.annotate(up_reg["id"][x], ha="center",
                                                      xy=(up_reg["x"][x], up_reg["y"][x]),
                                                      xytext=(0, 5),
                                                      textcoords="offset points",
                                                      arrowprops={"arrowstyle": "-", "color": "black"},
                                                      fontsize=14)
            up_updated = False
        else:
            if not up_updated:
                up_dupes = up_reg["id"].to_list().count(up_reg["id"][x])
                plotted_up[up_reg["id"][x]].set_text(up_reg["id"][x] + " (" + str(up_dupes) + ")")
                up_updated = True

    # Plot significantly down-regulated proteins in red with annotations
    sizes = [50] * len(down_reg["x"])
    ax.scatter(down_reg["x"], down_reg["y"], sizes, c="r", marker=".")
    down_len = len(down_reg["x"])
    if len(down_reg["x"]) > 10:
        down_reg = down_reg.nlargest(10, "y", "all")
    down_updated = False
    down_reg = down_reg.nlargest(10, "x", "all")
    for x in down_reg.index:
        if not down_reg["id"][x] in plotted_down.keys():
            plotted_down[down_reg["id"][x]] = ax.annotate(down_reg["id"][x], ha="center",
                                                          xy=(down_reg["x"][x], down_reg["y"][x]),
                                                          xytext=(0, 5),
                                                          textcoords="offset points",
                                                          arrowprops={"arrowstyle": "-", "color": "black"},
                                                          fontsize=14)
            down_updated = False
        else:
            if not down_updated:
                down_dupes = down_reg["id"].to_list().count(down_reg["id"][x])
                plotted_down[down_reg["id"][x]].set_text(down_reg["id"][x] + " (" + str(down_dupes) + ")")
                down_updated = True

    # Adjust the labels of any exceptions that were passed to the function
    if exceptions is not None:
        for x in exceptions:
            if x in plotted_up.keys():
                plotted_up[x].set_text("")
            elif x in plotted_down.keys():
                plotted_down[x].set_text("")

    # Plot non-significantly regulated proteins in gray
    sizes = [30] * len(un_reg["x"])
    ax.scatter(un_reg["x"], un_reg["y"], sizes, c="gray", marker=".")

    # Highlight HMGB1 in yellow
    if "hmgb1_x" in locals():
        ax.scatter(x=hmgb1_x, y=hmgb1_y, c="yellow", alpha=0.5, marker="o")
        ax.annotate("Hmgb1", xy=(hmgb1_x, hmgb1_y), xytext=(0, 5), textcoords="offset points",
                    arrowprops={"arrowstyle": "-", "color": "black"}, ha="center", fontsize=14)

    # Change axis limits to avoid text overlap/cut-off
    y_bot, y_top = ax.get_ylim()
    x_left, x_right = ax.get_xlim()

    ax.set_ylim(0, y_top + 0.3)
    ax.set_xlim(x_left - 0.2, x_right + 0.2)

    # Label the number of up- and down-regulated proteins at the top of each section
    obj1 = ax.text(x=0.01, y=0.99, s="Down-regulated: " + str(down_len), ha="left", va="top", fontsize=16,
                   transform=ax.transAxes)
    obj2 = ax.text(x=0.99, y=0.99, s="Up-regulated: " + str(up_len), ha="right", va="top", fontsize=16,
                   transform=ax.transAxes)
    obj3 = ax.text(x=0.01, y=0.01, s="Total: " + str(len(volcano_frame)), ha="left", va="bottom", fontsize=16,
                   transform=ax.transAxes)

    # Figure title and axis titles
    ax.set_xlabel(r"$log_{2}(Fold Change)$", fontsize=16)
    ax.set_ylabel(r"$-log_{10}(p-Value)$", fontsize=16)
    plt.suptitle(comparison, fontsize=16)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    # Formatting
    obj3 = ax.axhline(y=(-1 * np.log10(0.05)), linestyle="--", c="black")
    ax.axvline(x=0, linestyle="-", c="black")

    ax.spines["bottom"].set_color("black")
    ax.spines["top"].set_color("black")
    ax.spines["left"].set_color("black")
    ax.spines["right"].set_color("black")

    ax.spines["bottom"].set_linewidth(2)
    ax.spines["top"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)
    ax.spines["right"].set_linewidth(2)

    fig.tight_layout()
    plt.savefig("Figures" + "\\" + "Volcano " + comparison + ".png")
    plt.show()


def create_cluster(cluster_array, group_1, group_2, comparison, names):
    # Font Information: family=georgia, size=16 (14 for y-ticks), color=black
    col_names = []
    row_names = []
    name_1 = comparison.split(" v. ")[0]
    name_2 = comparison.split(" v. ")[1]

    remove_term = []
    for descriptor in name_1.split():
        if descriptor in name_2.split():
            remove_term += [descriptor]

    split_name_1 = name_1.split()
    split_name_2 = name_2.split()
    for term in remove_term:
        split_name_1.remove(term)
        split_name_2.remove(term)
    name_1 = " ".join(split_name_1)
    name_2 = " ".join(split_name_2)

    for x in names:
        protein_id = str(x)
        if protein_id.find("GN=") >= 0:
            protein_id = protein_id.split()
            protein_id = protein_id[-3]
            protein_id = protein_id[3:]
        elif protein_id.find("|") >= 0:
            protein_id = protein_id.split("|")[2]
            protein_id = protein_id.split("_")[0]

        row_names += [protein_id]

    print(row_names)
    print(len(row_names))
    for x in range(1, len(group_1.columns) + 1):
        col_names += [name_1 + " " + str(x)]

    for x in range(1, len(group_2.columns) + 1):
        col_names += [name_2 + " " + str(x)]

    # Eliminate isoforms
    cluster_frame = pd.DataFrame(cluster_array, columns=col_names, index=row_names)
    cluster_frame = cluster_frame[~cluster_frame.index.duplicated(keep='first')]

    # Make Cluster Map
    c = sns.clustermap(cluster_frame, col_cluster=False, z_score=0, cmap=sns.color_palette("vlag_r", as_cmap=True),
                       cbar_kws=dict(orientation='horizontal'), dendrogram_ratio=0.15, xticklabels=False,
                       figsize=(6, 8), norm=mpl.colors.Normalize(vmin=-2, vmax=2))
    plt.setp(c.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=16)
    plt.subplots_adjust(bottom=0.1)
    plt.subplots_adjust(right=0.80)
    plt.suptitle(comparison, x=0.5, fontsize=18)

    # Color Bar Formatting
    c.ax_col_dendrogram.set_visible(False)
    dendro_box = c.ax_col_dendrogram.get_position()
    dendro_box.y1 = dendro_box.y0 + 0.06
    dendro_box.y0 += 0.03
    dendro_box.x0 += 0.2
    dendro_box.x1 -= 0.2
    c.cax.set_position(dendro_box)
    c.cax.xaxis.set_ticks_position("bottom")
    c.ax_cbar.set_title("z-score of $log_{10}(intensity)$", fontsize=16)
    for spine in c.ax_cbar.spines:
        c.ax_cbar.spines[spine].set_color('black')
        c.ax_cbar.spines[spine].set_linewidth(2)

    # Figure Outline and Dotted Delimiter Line
    ax = c.ax_heatmap
    ax.add_patch(Rectangle((0, 0), len(cluster_frame.columns), len(cluster_frame.index), fill=False, edgecolor='black',
                           lw=5))
    ax.add_patch(Rectangle((0, 0), len(group_1.columns), len(cluster_frame.index),
                           fill=False, edgecolor='black', lw=2, linestyle="--"))

    # Group Bars and Group Labels
    if len(cluster_frame.index) >= 130:
        # Left Group
        ax.add_patch(Rectangle((0.2, len(cluster_frame.index) + 3), len(group_1.columns) - 0.2, 0, clip_on=False,
                               fill=False, edgecolor='black', lw=3, linestyle="-"))
        ax.text(len(group_1.columns) / 2, len(cluster_frame.index) + 10, name_1, ha="center", fontsize=20)

        # Right Group
        ax.add_patch(Rectangle((len(group_1.columns) + 0.2, len(cluster_frame.index) + 3), len(group_2.columns) - 0.2,
                               0, clip_on=False, fill=False, edgecolor='black', lw=3, linestyle="-"))
        ax.text(len(group_1.columns) + len(group_2.columns) / 2, len(cluster_frame.index) + 10, name_2, ha="center", fontsize=20)

    elif 130 > len(cluster_frame.index) >= 50:
        # Left Group
        ax.add_patch(Rectangle((0.2, len(cluster_frame.index) + 1.4), len(group_1.columns) - 0.2, 0, clip_on=False,
                               fill=False, edgecolor='black', lw=3, linestyle="-"))
        ax.text(len(group_1.columns) / 2, len(cluster_frame.index) + 4.5, name_1, ha="center", fontsize=20)

        # Right Group
        ax.add_patch(Rectangle((len(group_1.columns) + 0.2, len(cluster_frame.index) + 1.4), len(group_2.columns) - 0.2,
                               0, clip_on=False, fill=False, edgecolor='black', lw=3, linestyle="-"))
        ax.text(len(group_1.columns) + len(group_2.columns) / 2, len(cluster_frame.index) + 4.5, name_2, ha="center", fontsize=20)

    elif 50 > len(cluster_frame.index) > 25:
        # Left Group
        ax.add_patch(Rectangle((0.2, len(cluster_frame.index) + 0.8), len(group_1.columns) - 0.2, 0, clip_on=False,
                               fill=False, edgecolor='black', lw=4, linestyle="-"))
        ax.text(len(group_1.columns) / 2, len(cluster_frame.index) + 2.5, name_1, ha="center", fontsize=20)

        # Right Group
        ax.add_patch(Rectangle((len(group_1.columns) + 0.2, len(cluster_frame.index) + 0.8), len(group_2.columns) - 0.2,
                               0, clip_on=False, fill=False, edgecolor='black', lw=4, linestyle="-"))
        ax.text(len(group_1.columns) + len(group_2.columns) / 2, len(cluster_frame.index) + 2.5, name_2, ha="center", fontsize=20)

    elif 25 >= len(cluster_frame.index) >= 14:
        # Left Group
        ax.add_patch(Rectangle((0.2, len(cluster_frame.index) + 0.4), len(group_1.columns) - 0.2, 0, clip_on=False,
                               fill=False, edgecolor='black', lw=4, linestyle="-"))
        ax.text(len(group_1.columns) / 2, len(cluster_frame.index) + 1.2, name_1, ha="center", fontsize=20)

        # Right Group
        ax.add_patch(Rectangle((len(group_1.columns) + 0.2, len(cluster_frame.index) + 0.4), len(group_2.columns) - 0.2,
                               0, clip_on=False, fill=False, edgecolor='black', lw=4, linestyle="-"))
        ax.text(len(group_1.columns) + len(group_2.columns) / 2, len(cluster_frame.index) + 1.2, name_2, ha="center", fontsize=20)

    else:
        # Left Group
        ax.add_patch(Rectangle((0.2, len(cluster_frame.index) + 0.2), len(group_1.columns) - 0.2, 0, clip_on=False,
                               fill=False, edgecolor='black', lw=4, linestyle="-"))
        ax.text(len(group_1.columns) / 2, len(cluster_frame.index) + 0.6, name_1, ha="center", fontsize=20)

        # Right Group
        ax.add_patch(Rectangle((len(group_1.columns) + 0.2, len(cluster_frame.index) + 0.2), len(group_2.columns) - 0.2,
                               0, clip_on=False, fill=False, edgecolor='black', lw=4, linestyle="-"))
        ax.text(len(group_1.columns) + len(group_2.columns) / 2, len(cluster_frame.index) + 0.6, name_2, ha="center", fontsize=20)

    # Save the significant proteins
    # cluster_frame.to_csv("Figures" + "\\" + "Protein Lists" + "\\" + "Significant Proteins " + comparison + ".csv")

    # Save and Show Figure
    plt.savefig("Figures" + "\\" + "Cluster " + comparison + ".png")
    plt.show()


def create_venn(file, comparison):
    sig_proteins = pd.read_excel(file, "Sheet1")
    name_1 = comparison.split(" v. ")[0]
    name_2 = comparison.split(" v. ")[1]
    list_1 = list(sig_proteins[name_1])
    list_2 = list(sig_proteins[name_2])

    # Find NaN count in list_1
    nan_count = 0

    for x in range(len(list_1)):
        list_1[x] = str(list_1[x])
        if list_1[x] == "nan":
            nan_count += 1

    for x in range(nan_count):
        list_1.remove("nan")

    # Find NaN count in list_2
    nan_count = 0

    for x in range(len(list_2)):
        list_2[x] = str(list_2[x])
        if list_2[x] == "nan":
            nan_count += 1

    for x in range(nan_count):
        list_2.remove("nan")

    # Set colors and label text based on comparison
    colors = ("gray", "black")
    labels = ("Unique to...", "Unique to...")

    if "Male" in name_1.split() and "Female" in name_2.split():
        colors = ("teal", "pink")
        labels = ("Unique to Male", "Unique to Female")
    elif "WT" in name_1.split() and "KO" in name_2.split():
        if "Male" in name_1.split():
            colors = ("teal", "green")
            labels = ("Unique to WT", "Unique to KO")
        elif "Female" in name_1.split():
            colors = ("pink", "magenta")
            labels = ("Unique to WT", "Unique to KO")

    # Create Venn Diagram
    fig, ax = plt.subplots()
    fig.set_size_inches(4.5, 4.5)
    v1 = venn2([set(list_1), set(list_2)], labels, ax=ax, set_colors=colors)
    c1 = venn2_circles([set(list_1), set(list_2)])

    # Venn diagram customization for all Venn diagrams
    c1[0].set_edgecolor("white")
    c1[1].set_edgecolor("white")
    plt.suptitle(comparison, y=0.9, fontsize=16)
    for text in v1.set_labels:
        text.set_fontsize(16)

    # Find proteins in shared circle
    shared = list(set(list_1) & set(list_2))
    nan_count = 0
    share_count = 0
    for x in range(len(shared)):
        shared[x] = str(shared[x])
        if shared[x] == "nan":
            nan_count += 1
        else:
            share_count += 1

    for x in range(nan_count):
        shared.remove("nan")

    # Customize annotation and other properties based on whether diagram has a shared region
    if len(shared) == 0:
        # Annotate
        plt.annotate("No shared proteins", xy=(0.5, 0.5), xytext=(0, -125),
                     ha='center', xycoords='figure fraction', textcoords='offset points', fontsize=16)
        v1.subset_labels[0].set_fontsize(16)
        v1.subset_labels[1].set_fontsize(16)

        # Move labels to the top of circles
        l1 = v1.get_label_by_id("A")
        x, y = l1.get_position()
        l1.set_position((x + 0.3, -y + 0.1))
        l2 = v1.get_label_by_id("B")
        x, y = l2.get_position()
        l2.set_position((x - 0.3, -y + 0.1))

    else:
        # Adjust annotation based on length
        shared = ", ".join(shared)

        if len(shared) < 3:
            shared = "None"
        if len(shared) > 25:
            shared = shared[:shared.find(",", 25) + 1] + "\n" + shared[shared.find(",", 25) + 1:]
        if len(shared) > 50:
            shared = shared[:shared.find(",", 50) + 1] + "\n" + shared[shared.find(",", 50) + 1:]
        if len(shared) > 75:
            shared = shared[:shared.find(",", 75) + 1] + "\n" + shared[shared.find(",", 75) + 1:]

        # Annotate
        plt.annotate(shared, xy=v1.get_label_by_id('11').get_position() - np.array([0, 0.05]), xytext=(0, -125),
                     ha='center', textcoords='offset points', arrowprops=dict(arrowstyle='->', color='black'), fontsize=16)

        v1.get_label_by_id("11").set_text(str(share_count))
        for text in v1.subset_labels:
            text.set_fontsize(16)

        # Move labels to the top of circles
        l1 = v1.get_label_by_id("A")
        x, y = l1.get_position()
        l1.set_position((x + 0.2, -y + 0.1))
        l2 = v1.get_label_by_id("B")
        x, y = l2.get_position()
        l2.set_position((x - 0.2, -y + 0.1))

    # Show and save Venn diagram
    plt.savefig("Figures" + "\\" + "Venn Diagram " + comparison + ".png")
    plt.show()


def create_scatter(scatter_frame, comparison, exceptions=None):
    q1 = {"x": [], "y": [], "id": []}
    q2 = {"x": [], "y": [], "id": []}
    q3 = {"x": [], "y": [], "id": []}
    q4 = {"x": [], "y": [], "id": []}
    center = {"x": [], "y": [], "id": []}
    fig, ax = plt.subplots()
    fig.set_size_inches(7, 7)
    plotted = {}

    for x in range(len(scatter_frame.index)):

        # Adjust protein ID so it includes just protein name
        protein_id = str(scatter_frame["Protein"][x])
        if protein_id.find("GN=") >= 0:
            protein_id = protein_id.split()
            protein_id = protein_id[-3]
            protein_id = protein_id[3:]
        elif protein_id.find("|") >= 0:
            protein_id = protein_id.split("|")[2]
            protein_id = protein_id.split("_")[0]

        # Record location of HMGB1 separate from other proteins
        if protein_id == "Hmgb1":
            hmgb1_x = scatter_frame["WT Fold Change"][x]
            hmgb1_y = scatter_frame["KO Fold Change"][x]

        # Record locations of all proteins with positive fold change for both WT and KO animals
        if (scatter_frame["WT Fold Change"][x] >= 1 and scatter_frame["KO Fold Change"][x] >= 0) or \
                (scatter_frame["WT Fold Change"][x] >= 0 and scatter_frame["KO Fold Change"][x] >= 1):
            q1["x"] += [scatter_frame["WT Fold Change"][x]]
            q1["y"] += [scatter_frame["KO Fold Change"][x]]
            q1["id"] += [protein_id]

        # Record locations of all proteins with negative fold change for WT and positive fold change for KO animals
        elif (scatter_frame["WT Fold Change"][x] <= -1 and scatter_frame["KO Fold Change"][x] >= 0) or \
                (scatter_frame["WT Fold Change"][x] <= 0 and scatter_frame["KO Fold Change"][x] >= 1):
            q2["x"] += [scatter_frame["WT Fold Change"][x]]
            q2["y"] += [scatter_frame["KO Fold Change"][x]]
            q2["id"] += [protein_id]

        # Record locations of all proteins with negative fold change for both WT and KO animals
        elif (scatter_frame["WT Fold Change"][x] <= -1 and scatter_frame["KO Fold Change"][x] <= 0) or \
                (scatter_frame["WT Fold Change"][x] <= 0 and scatter_frame["KO Fold Change"][x] <= -1):
            q3["x"] += [scatter_frame["WT Fold Change"][x]]
            q3["y"] += [scatter_frame["KO Fold Change"][x]]
            q3["id"] += [protein_id]

        # Record locations of all proteins with positive fold change for WT and negative fold change for KO animals
        elif (scatter_frame["WT Fold Change"][x] >= 1 and scatter_frame["KO Fold Change"][x] <= 0) or \
                (scatter_frame["WT Fold Change"][x] >= 0 and scatter_frame["KO Fold Change"][x] <= -1):
            q4["x"] += [scatter_frame["WT Fold Change"][x]]
            q4["y"] += [scatter_frame["KO Fold Change"][x]]
            q4["id"] += [protein_id]

        # Record locations of all proteins with fold changes below threshold
        elif -1 < scatter_frame["WT Fold Change"][x] < 1 and -1 < scatter_frame["KO Fold Change"][x] < 1:
            center["x"] += [scatter_frame["WT Fold Change"][x]]
            center["y"] += [scatter_frame["KO Fold Change"][x]]

    # Plot all proteins with colors/annotations according to their quadrant
    for quadrant, color, marker, markersize in zip([center, q1, q2, q3, q4],
                                                   ["gray", "blue", "orange", "red", "orange"],
                                                   [".", ".", ".", ".", "."],
                                                   [30, 50, 50, 50, 50]):
        sizes = [markersize] * len(quadrant["x"])
        ax.scatter(quadrant["x"], quadrant["y"], sizes, c=color, marker=marker)
        # if not quadrant == center:
        #     for x in range(len(quadrant["id"])):
        #         plotted[quadrant["id"][x]] = ax.annotate(quadrant["id"][x], ha="center",
        #                                                  xy=(quadrant["x"][x], quadrant["y"][x]), xytext=(0, 5),
        #                                                  arrowprops={"arrowstyle": "-", "color": "black"},
        #                                                  textcoords="offset points")

    # Adjust the labels of any exceptions that were passed to the function
    if exceptions is not None:
        for x in exceptions:
            if x in plotted.keys():
                plotted[x].set_text("")

    # Highlight HMGB1 in yellow
    if "hmgb1_x" in locals():
        ax.scatter(x=hmgb1_x, y=hmgb1_y, c="yellow", alpha=0.5, marker="o")
        ax.annotate("Hmgb1", xy=(hmgb1_x, hmgb1_y), xytext=(0, 5), textcoords="offset points",
                    arrowprops={"arrowstyle": "-", "color": "black"}, ha="center", fontsize=14)

    # Figure title and axis titles
    ax.set_xlabel(r"$log_{2}(WT-Fold Change)$", fontsize=16)
    ax.set_ylabel(r"$log_{2}(KO-Fold Change)$", fontsize=16)
    plt.suptitle(comparison, fontsize=18)

    # Formatting
    obj3 = ax.plot([-1, 1, 1, -1, -1], [1, 1, -1, -1, 1], linestyle=(0, (5, 10)), c="gray")

    # Formatting axis ticks
    x_min = int(scatter_frame["WT Fold Change"].min() // 1)
    x_max = int((scatter_frame["WT Fold Change"].max() + 1) // 1)

    y_min = int(scatter_frame["KO Fold Change"].min() // 1)
    y_max = int((scatter_frame["KO Fold Change"].max() + 1) // 1)

    x_ticks = range(x_min, x_max + 1)
    y_ticks = range(y_min, y_max + 1)

    plt.xticks(x_ticks, fontsize=16)
    plt.yticks(y_ticks, fontsize=16)

    # Add x=0 and y=0 lines
    ax.axvline(x=0, linestyle="-", c="black")
    ax.axhline(y=0, linestyle="-", c="black")

    # Add arrows with annotations above each
    ax.text(x=0.225, y=-0.205, s="Stress-induced\ndownregulation in WT", ha="center", va="center", fontsize=16,
            transform=ax.transAxes)
    ax.annotate('', xy=(0, -0.15), xycoords='axes fraction', xytext=(0.45, -0.15), arrowprops=dict(arrowstyle="->", color='k'))
    ax.text(x=0.775, y=-0.205, s="Stress-induced\nupregulation in WT", ha="center", va="center", fontsize=16,
            transform=ax.transAxes)
    ax.annotate('', xy=(0.55, -0.15), xycoords='axes fraction', xytext=(1, -0.15), arrowprops=dict(arrowstyle="<-", color='k'))

    ax.text(x=-0.205, y=0.225, s="Stress-induced\ndownregulation in KO", ha="center", va="center", fontsize=16,
            transform=ax.transAxes, rotation="vertical")
    ax.annotate('', xy=(-0.15, 0), xycoords='axes fraction', xytext=(-0.15, 0.45), arrowprops=dict(arrowstyle="->", color='k'))
    ax.text(x=-0.205, y=0.775, s="Stress-induced\nupregulation in KO", ha="center", va="center", fontsize=16,
            transform=ax.transAxes, rotation="vertical")
    ax.annotate('', xy=(-0.15, 0.55), xycoords='axes fraction', xytext=(-0.15, 1), arrowprops=dict(arrowstyle="<-", color='k'))

    # Add % annotation in quadrant 2
    percent = str(round(((len(q2["x"]) + len(q4["x"])) / (len(q1["x"]) + len(q2["x"]) + len(q3["x"]) + len(q4["x"]))) * 100, 1))
    ax.text(x=0.01, y=0.99, s=percent + "%\nHMGB1 Dependent", ha="left", va="top", fontsize=16, transform=ax.transAxes)

    ax.spines["bottom"].set_color("black")
    ax.spines["top"].set_color("black")
    ax.spines["left"].set_color("black")
    ax.spines["right"].set_color("black")

    ax.spines["bottom"].set_linewidth(2)
    ax.spines["top"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)
    ax.spines["right"].set_linewidth(2)

    with pd.ExcelWriter("Figures" + "\\" + "Protein Lists" + "\\" + "Q2+Q4 of Scatter Plot Dataset " + comparison + ".xlsx") as writer:
        pd.DataFrame(q2).to_excel(writer, sheet_name="Q2", index=False)
        pd.DataFrame(q4).to_excel(writer, sheet_name="Q4", index=False)

    fig.tight_layout()
    plt.subplots_adjust(bottom=0.2, left=0.2)

    plt.savefig("Figures" + "\\" + "Scatter " + comparison + ".png")
    plt.show()


def volcano_log(a):
    return -1 * np.log10(a)


analysis = input("Examine the effect of genetic manipulation, stress condition, sex, experimental resilience profile, "
                 "or all? Type 'Genetic', 'Stress', 'Sex', 'Resilience', or All: ")
if analysis == "Stress":
    venn_file_input = input(r"Input path of Venn diagram data file: ")
    file_input = input(r"Input path of raw data file: ")
    sex_input = input(r"Input sex of mice: ")
    stress_comparison(file_input, sex_input, venn_file_input)
elif analysis == "Sex":
    venn_file_input = input(r"Input path of Venn diagram data file: ")
    male_file_input = input(r"Input path of male raw data file: ")
    female_file_input = input(r"Input path of female raw data file: ")
    sex_comparison(male_file_input, female_file_input, venn_file_input)
elif analysis == "Genetic":
    male_file_input = input(r"Input path of male raw data file: ")
    female_file_input = input(r"Input path of female raw data file: ")
    genetic_comparison(male_file_input, female_file_input)
elif analysis == "Resilience":
    male_file_input = input(r"Input path of male raw data file: ")
    female_file_input = input(r"Input path of female raw data file: ")
    resilience_comparison(male_file_input, female_file_input)
elif analysis == "All":
    male_file_input = r"C:\Users\Luke\Desktop\College\Research\Dr. Franklin\Proteomics Pipeline\Franklin Male Analysis Results.xlsx"
    female_file_input = r"C:\Users\Luke\Desktop\College\Research\Dr. Franklin\Proteomics Pipeline\Franklin Female Analysis Results.xlsx"
    venn_file_input = r"C:\Users\Luke\Desktop\College\Research\Dr. Franklin\Proteomics Pipeline\Venn Diagram Data Set.xlsx"
    stress_comparison(male_file_input, "Male", venn_file_input)
    stress_comparison(female_file_input, "Female", venn_file_input)
    sex_comparison(male_file_input, female_file_input, venn_file_input)
    genetic_comparison(male_file_input, female_file_input)
    resilience_comparison(male_file_input, female_file_input)

# C:\Users\Luke\Desktop\College\Research\Dr. Franklin\Proteomics Pipeline\Franklin Male Analysis Results.xlsx
# C:\Users\Luke\Desktop\College\Research\Dr. Franklin\Proteomics Pipeline\Franklin Female Analysis Results.xlsx
# C:\Users\Luke\Desktop\College\Research\Dr. Franklin\Proteomics Pipeline\Venn Diagram Data Set.xlsx
