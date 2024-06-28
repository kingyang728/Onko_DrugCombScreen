# About Onko_DrugCombScreen

## Purpose
Onko_DrugCombScreen is a powerful tool designed to assist researchers in identifying potential drug combination treatments for targeted cancers compared to comparison cohorts. It integrates various data analysis and visualization tools to analyze the potential efficacy and synergy of drug combinations.

## Features
- **Data Input**: Users can upload custom target and comparison cohorts' variants data, as well as cell line data in CSV format, for drug recommendation analysis.
- **Analysis Tools**: The app offers several analysis tools, including volcano plots, heatmaps, circle plots, and alluvial diagrams, to visualize drug combination recommendations.
- **Interactive Visualizations**: Users can interact with the plots to explore the data more deeply, such as hovering over points for more details or selecting specific drugs to highlight.
- **Data Export**: Results from the analyses can be downloaded for further investigation or reporting.
- **Test Data Analysis**: Click the `Test_Data` button to run a sample analysis with test data.

## How to Use
1. **Upload Data**: Begin by uploading your target and comparison data files using the file input option on the left panel, along with any relevant cell line data.
2. **Select Analysis Parameters**:  Choose your target cancer disease, set the percentage threshold, and select other parameters such as test type (one-sided or two-sided) and response type (positive, negative, etc.) as needed.
3. **Run Analysis**: Click "Submit" button to process the data. Explore the generated plots and tables to examine potential drug combinations.
4. **Download Results**: Use the download buttons to save your results for downstream analysis or publication.

## Analysis Reuslts Descriptions
- **Volcano Plot**: Visualizes significant recommended drug combinations, where each node represents a potential drug combination. Nodes vary in size according to occurrence percentage and shape to indicate evidence level.
- **Heatmap**: Shows the proportion of recommended drug combinations for the target cancer group, with customizable percentage thresholds and colors.
- **Circle Plot**: Displays recommended drug proportions and co-recommendation connections.
- **Alluvial Diagram**: Represents connections from drugs to genes and variants to Sample IDs within the target cancer group, allowing for drug selection customization.
- **UpsetPlot**: Illustrates proportions and intersections of recommended drugs, providing insights into overlap and distinct recommendations.
- **BarPlot**: Depicts the percentage of patient counts for a single recommended drug in the target cancer group.
- **Primary_Table**: Offers a drug recommendation table for the target cancer group.
- **Comparison_Table**: Provides a drug recommendation table for the comparison group.
- **Cellline_Table**: Shows a drug recommendation table for the cell line group.
- **DrugComb_Analysis_Table**: Features an analysis of candidate drug combinations, incorporating Fisher's exact tests and cell line data.

## Imprint
Georg-August-Universität Göttingen  
Public Law Foundation  
Universitätsmedizin Göttingen  
Robert-Koch-Str. 40  
37075 Göttingen  
Postal address: 37099 Göttingen  
represented by the Management Board  
For more details, visit our [Imprint Page](https://bioinformatics.umg.eu/imprint/).

## Contact
For support or more information, please reach out to us at [jingyu.yang@bioinf.med.uni-goettingen.de].

Thank you for using Onko_DrugCombScreen.