import streamlit as st
import os
from tool1 import heatmapper
from tool2 import snp_ncbi
from tool4 import pubchem
from tool3 import swissadme_tool
from tool5 import clustalo
from tool6 import risk_assessment
from Homepage import home_page
from excel import upload_excel
from aboutpage import about_page
from feedback import feed_back

# Inject custom CSS for styling
def apply_custom_css():
    st.markdown(
        """
        <style>
        /* Center and spread tabs across the screen */
        .stTabs [data-baseweb="tab"] {
            flex: 1; /* Make tabs evenly spread */
            text-align: center; /* Center-align text */
            background-color: transparent !important; /* Remove background color */
        }
        /* Increase font size and make tab labels bold */
        .stTabs [data-baseweb="tab"] > div {
            font-size: 18px; /* Adjust font size */
            font-weight: bold !important; /* Make it bold */
            color: black; /* Set text color for clarity */
            background-color: transparent !important; /* Ensure no background color for labels */
        }
        </style>
        """,
        unsafe_allow_html=True,
    )

def main():
    # Apply custom CSS
    apply_custom_css()

    # Display title and image
    col1, col2 = st.columns([7, 1])

    with col1:
        st.title("Breast Cancer Research Hub")

    with col2:
        # Check if the image exists using a relative path
        image_path = os.path.join("extra", "breast-cancer.png")
        if os.path.exists(image_path):
            st.image(image_path, width=500)
        else:
            st.error("Image not found!")

    # Horizontal Navigation Bar
    tabs = st.tabs(["Home", "Database Viewer", "Tools", "About", "Feedback"])

    with tabs[0]:
        st.subheader("Welcome to the Breast Cancer Research Hub")
        home_page()

    with tabs[1]:
        st.subheader("Database Viewer")
        upload_excel()

    with tabs[2]:
        tools_menu()

    with tabs[3]:
        about_page()

    with tabs[4]:
        feed_back()

# Tools Submenu Logic
def tools_menu():
    # Horizontal Tabs for Tools
    tool_tabs = st.tabs([  
        "Heatmapper", "SNP-NCBI", "SwissADME", "PubChem", 
        "Clustal Omega", "Breast Cancer Risk Tool"
    ])

    with tool_tabs[0]:
        heatmapper()
    with tool_tabs[1]:
        snp_ncbi()
    with tool_tabs[2]:
        swissadme_tool()
    with tool_tabs[3]:
        pubchem()
    with tool_tabs[4]:
        clustalo()
    with tool_tabs[5]:
        risk_assessment()

if __name__ == "__main__":
    main()
