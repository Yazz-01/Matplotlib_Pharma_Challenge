#Matplotlib_Pharma_Challenge

By using Matplotlib and Pandas I helped to Pharmaceuticals Inc. to look for potential treatments for squamous cell carcinoma (SCC) done by an Analysis of Drug Regimens on Mice. This contains a script analysing pharmaceutical data from mice and tumor volumes. The goal was to analyze the mice and the drug regimens in the study. The main drug analysis is done on Capumulin.

1. Resources Required
Resources files containing the csv files.
#Study data files
mouse_metadata_path = "Resources/Mouse_metadata.csv"
study_results_path = "Resources/Study_results.csv"

2. #Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
from scipy.stats import sem
import numpy as np
import csv

3. #Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

#Combine the data into a single dataset
merge_df= pd.DataFrame.merge(mouse_metadata,study_results,
                    on="Mouse ID", how="outer")
#Display the data table for preview
merge_df.head(5)


#4.Cleaning the data

#Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicated_mouse=merge_df.loc[merge_df.duplicated(subset=["Mouse ID","Timepoint"]),"Mouse ID"].unique()
print(f"The duplicated mouse in the study is {duplicated_mouse}")


#Get all the data for the duplicate mouse ID. 
clean_study_df= merge_df.drop_duplicates(subset=["Mouse ID","Timepoint"],keep='first', inplace=False)
clean_study_df.loc[clean_study_df["Mouse ID"]=="g989"]

#Checking for null values 
na=clean_study_df.isnull().values.any()

#Checking for the sum of null values 
null= clean_study_df.isnull().sum()

#Checking the number of mice in the clean DataFrame.
num_mice =len(clean_study_df["Mouse ID"])
print(f"The final number of mice in the study with no duplicated is {num_mice}")


5. Statistics Summary

#Generate a summary statistics table of mean, median, variance, standard deviation, and SEM 
#of the tumor volume for each regimen

mean_tumor = clean_study_df.groupby(["Drug Regimen"])["Tumor Volume (mm3)"].mean()
median_tumor = clean_study_df.groupby(["Drug Regimen"])["Tumor Volume (mm3)"].median()
var_tumor = clean_study_df.groupby(["Drug Regimen"])["Tumor Volume (mm3)"].var()
std_tumor = clean_study_df.groupby(["Drug Regimen"])["Tumor Volume (mm3)"].std()
sem_tumor = clean_study_df.groupby(["Drug Regimen"])["Tumor Volume (mm3)"].sem()

statistcs_df= pd.DataFrame({"Mean":mean_tumor, "Median":median_tumor,"Variance":var_tumor, "Stan_Dev":std_tumor, "SEM":sem_tumor})


Double checking the results by doing it again from a different method
#Drug regime per Tumor volume statistics: in a single line mean, median, variance, standard deviation and sem in one line
stats_tumor = clean_study_df.groupby(["Drug Regimen"])["Tumor Volume (mm3)"].agg([np.mean,np.median,np.var,np.std,st.sem])

6. Line chart
#Generate a bar plot showing the total number of measurements taken per drug regimen using pandas.
total_measures = clean_study_df.groupby(["Drug Regimen"]).count()["Tumor Volume (mm3)"]

#Plotting total measurements
plot_total_measures = total_measures.plot(kind='bar',rot=90, figsize=(10,5))
plt.title("Total Number of Measurements per Drug Regime")
plt.ylabel("Total Measurements per Drug Regime")
plt.show()

#Plotting total measurements
plot_total_measures = total_measures.plot(kind="bar",rot=90, figsize=(10,5))
plt.title("Total Number of Measurements per Drug Regime")
plt.ylabel("Total Measurements per Drug Regime")

7. Pie chart

#Find the distribution of female versus male mice using pandas
total_mice_grouped= clean_study_df.groupby(["Sex","Mouse ID"]).size()
total_mice= len(total_mice_grouped)
print(f"The total number of mice is {total_mice}")
#--------Male mice
male_mice = total_mice_grouped["Male"].count()
percent_male_mice = male_mice/total_mice*100
print(f"The number of Male mice is {male_mice} and the percentage is {percent_male_mice}")

#--------Female mice
female_mice = total_mice_grouped["Female"].count()
percent_female_mice = female_mice/total_mice*100
print(f"The number of Female mice is {female_mice} and the percentage is {percent_female_mice}")

#----- Summary of mice in the study
male_female_df = pd.DataFrame({"Distribution of mice":[male_mice,female_mice],"Percentage":[percent_male_mice,percent_female_mice]},
                               index=["Male", "Female"])
#Generate a pie plot showing the distribution of female versus male mice using pandas
explode = (0.1, 0)
colors = ['blue', 'orange']
male_female_df_plot = male_female_df.plot.pie(y='Distribution of mice',figsize=(5,5), colors = colors, startangle=140, explode = explode, shadow = True, autopct="%1.1f%%")


8. Final TUmor Volume with Capomulin
#Calculate the final tumor volume for the Top four Drug Regimenes "Capomulin", "Ramicane", "Infubinol", "Ceftamin"

max_tumor_vol= clean_study_df.groupby("Mouse ID").max()
max_tumor_vol_df= max_tumor_vol.reset_index() 

#Filtering the data by Drug regimen vs Timepoint
cap_df=clean_study_df.loc[(clean_study_df["Drug Regimen"]=="Capomulin")&(clean_study_df["Timepoint"]==45)]

9. Boxplot

#Merging information of the four top drug regimen and the max tumor volume vs timepoint 
top_tumor_df= max_tumor_vol_df[["Mouse ID","Timepoint"]].merge(clean_study_df,on= ["Mouse ID","Timepoint"],how= "left")

#Identifying the Tumor Volume per the top four Drug Regimen
capomulin_vol= top_tumor_df.loc[top_tumor_df["Drug Regimen"]=="Capomulin"]["Tumor Volume (mm3)"]

ramicane_vol= top_tumor_df.loc[top_tumor_df["Drug Regimen"]=="Ramicane"]["Tumor Volume (mm3)"]

infubinol_vol= top_tumor_df.loc[top_tumor_df["Drug Regimen"]=="Infubinol"]["Tumor Volume (mm3)"]

ceftamin_vol= top_tumor_df.loc[top_tumor_df["Drug Regimen"]=="Ceftamin"]["Tumor Volume (mm3)"]

tumor_volumes = [capomulin_vol,ramicane_vol,infubinol_vol,ceftamin_vol]

#Listing top drug remines to label dataframe
best_treat = ['Capomulin', 'Ramicane', 'Infubinol','Ceftamin']

#Generate a box plot of the final tumor volume of each mouse across four regimens of interest
plt.boxplot(tumor_volumes, labels=best_treat)
plt.ylim(10, 80)
plt.show()

10. Line and Scatter Plot
#Generate a line plot of time point versus tumor volume for a mouse treated with Capomulin
time_tumor = clean_study_df[clean_study_df["Mouse ID"].isin(["j119"])]

time_tumor_df = [["Mouse ID", "Timepoint", "Tumor Volume (mm3)"]]

line_plot = time_tumor.reset_index()

line_plot_graph = line_plot[["Mouse ID", "Timepoint", "Tumor Volume (mm3)"]]

lines = line_plot_graph.plot.line()


#Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen

cap_scatter = cap_df[["Mouse ID","Weight (g)", "Tumor Volume (mm3)"]]

cap_sorted = cap_scatter.sort_values(["Weight (g)"], ascending=True)

cap_scat_plot = cap_scatter.reset_index()

cap_group_weight = cap_scat_plot.groupby("Weight (g)")["Tumor Volume (mm3)"].mean()

cap_group_plot = pd.DataFrame(cap_group_weight).reset_index()

cap_scatter = cap_group_plot.plot(kind='scatter', x='Weight (g)', y='Tumor Volume (mm3)', grid = True, figsize= (8,8))


11. Correlation and Regression

#Calculate the correlation coefficient and linear regression model for mouse weight
x_values = cap_group_plot["Weight (g)"]
y_values = cap_group_plot["Tumor Volume (mm3)"]
(slope, intercept, rvalue, pvalue, stderr) = st.linregress(x_values, y_values)
regress_values = x_values * slope + intercept
line_eq = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))
plt.scatter(x_values, y_values)
plt.plot(x_values,regress_values,"r-")
plt.annotate(line_eq,(6,10),fontsize=10,color="red")
plt.xlabel("Weight")
plt.ylabel("Tumor Volume")
plt.title("Weight vs Tumor Vol")
plt.show()

