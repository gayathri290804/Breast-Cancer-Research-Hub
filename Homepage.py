import streamlit as st
import pandas as pd
from openpyxl import load_workbook
import os

def home_page():
    # Title of the home page
    st.title("Welcome to the Breast Cancer Research Hub")

    # Subheader and Text about Breast Cancer
    st.subheader("What is Breast Cancer?")
    st.write("""
    Breast cancer is a type of cancer that begins in the cells of the breast. It can occur in both men and women, though it is far more common in women. 
    The cells in the breast tissue divide and grow uncontrollably, forming a lump or tumor. Early detection is crucial for improving survival rates.
    """)

    # Add key facts about Breast Cancer
    st.subheader("Key Facts about Breast Cancer")
    st.write("""
    - Breast cancer is the most common cancer among women globally.
    - Early detection through regular screenings can significantly improve outcomes.
    - Common risk factors include age, genetic predisposition, and lifestyle factors.
    - Treatments include surgery, radiation therapy, chemotherapy, hormone therapy, and targeted therapy.
    """)

    # Add an image related to Breast Cancer
    st.subheader("Breast Cancer Awareness")

    # Use relative path and add fallback handling
    image_path = "extra/Breast_cancer.png"
    if os.path.exists(image_path):
        st.image(image_path, caption="Breast Cancer Awareness")
    else:
        st.warning("Image not found. Please check the file path or ensure the image is included in the deployment.")

    # Add a call-to-action
    st.subheader("Join Us in Advancing Breast Cancer Research")
    st.write("""
    Whether you are a researcher, clinician, or student, the Breast Cancer Research Hub offers valuable resources to aid your work.
    Together, we can make a difference in the fight against breast cancer.
    """)

    # Footer
    st.markdown("""
    ---
    © 2025 Breast Cancer Research Hub
    """)

# Run the app
if __name__ == "__main__":
    home_page()
