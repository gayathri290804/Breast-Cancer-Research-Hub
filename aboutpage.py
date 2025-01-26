import streamlit as st

def about_page():
    st.write("""
    ## About Breast Cancer Research Hub:

    The **Breast Cancer Research Hub** is a comprehensive platform designed to integrate and provide valuable insights into breast cancer research. This initiative includes a curated database containing information on genes, proteins, molecular structures, pathways, post-translational modifications, subcellular localizations, and other essential data relevant to breast cancer.

    ### Key Features of the Breast Cancer Research Hub:

    - **Gene and Protein Information**: Access detailed information about genes and proteins associated with breast cancer.
    - **Structural and Functional Insights**: Explore molecular structures, pathways, and functional annotations.
    - **Integrated Tools**: The platform incorporates six powerful tools to support researchers and clinicians:
        - **Heatmapper**: For visualizing gene or protein expression data.
        - **SNP Tool**: To analyze single nucleotide polymorphisms (SNPs) relevant to breast cancer.
        - **SWISS ADME**: For evaluating the pharmacokinetics and drug-likeness of potential therapeutic compounds.
        - **PubChem**: A database for accessing chemical information and potential drug candidates.
        - **Breast Cancer Risk Tool**: For assessing individual risk factors for breast cancer.
        - **Clustal Omega**: For performing multiple sequence alignments of gene or protein sequences.

    ### Mission and Vision:

    Our mission is to create a unified hub that empowers researchers with easy access to critical breast cancer data and tools. By integrating diverse data sources and providing user-friendly tools, we aim to:
    - Accelerate the pace of breast cancer research.
    - Enhance understanding of genetic and molecular mechanisms.
    - Support the development of new and effective therapeutic strategies.

    ### Advancing Breast Cancer Research:

    With its holistic approach, the Breast Cancer Research Hub is a valuable resource for researchers, clinicians, and students aiming to deepen their understanding of breast cancer and contribute to the fight against this disease.

    Together, we can push the boundaries of research and bring hope to millions affected by breast cancer.
    """)

# Example of how the about_page function would be used in a Streamlit app
if __name__ == "__main__":
    st.title("Breast Cancer Research Hub")
    about_page()
