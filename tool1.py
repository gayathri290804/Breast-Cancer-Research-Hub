import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import requests
from io import BytesIO
import webbrowser
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors



# Heatmapper
def heatmapper():
    st.title("Heatmapper:")
    st.subheader("Heatmapper is an advanced tool widely used in breast cancer research for visualizing and analyzing data, such as gene expression levels, protein interactions, or other biological datasets. Heatmaps generated by Heatmapper can highlight patterns and correlations in complex datasets, making it a valuable resource for identifying biomarkers, pathways, and molecular signatures associated with breast cancer.")
    
    # File uploader
    uploaded_file = st.file_uploader("Upload Dataset for Heatmap", type=["xlsx", "csv", "txt", "tsv"])
    
    if uploaded_file:
        # Read the uploaded file
        if uploaded_file.name.endswith(".csv"):
            df = pd.read_csv(uploaded_file)
        elif uploaded_file.name.endswith(".tsv") or uploaded_file.name.endswith(".txt"):
            df = pd.read_csv(uploaded_file, delimiter='\t')
        else:
            df = pd.read_excel(uploaded_file)
        
        st.write("### Dataset Preview")
        st.dataframe(df)

        # Column selection
        st.write("### Select columns for Heatmap:")
        columns = st.multiselect("Columns", df.columns)

        if len(columns) > 1:
            # Compute correlation matrix
            corr_matrix = df[columns].corr()

            # Plot heatmap
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", ax=ax)
            plt.title("Correlation Heatmap")

            st.pyplot(fig)

            # Add download button
            st.write("### Download Heatmap as PNG")
            buf = BytesIO()
            fig.savefig(buf, format="png")
            buf.seek(0)
            st.download_button(
                label="Download Heatmap",
                data=buf,
                file_name="heatmap.png",
                mime="image/png"
            )






