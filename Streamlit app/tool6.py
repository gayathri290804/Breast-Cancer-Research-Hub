import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import requests
from io import BytesIO
import webbrowser
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


# Breast Cancer Risk Assessment
def risk_assessment():
    st.title("Breast Cancer Risk Assessment:")
    st.subheader("Breast Cancer Risk Assessment involves evaluating various factors, such as genetics, family history, lifestyle, and environmental influences, to estimate an individual's likelihood of developing breast cancer. This assessment is crucial for improving patient outcomes by enabling early intervention and tailored management plans.")

    # Collect input data
    age = st.number_input("Enter your age:", min_value=18, max_value=120)
    family_history = st.radio("Do you have a family history of breast cancer?", ("Yes", "No"))
    genetic_testing = st.radio("Have you undergone genetic testing for breast cancer susceptibility?", ("Yes", "No"))
    hormone_therapy = st.radio("Have you ever taken hormone replacement therapy?", ("Yes", "No"))
    alcohol_consumption = st.radio("Do you consume alcohol regularly?", ("Yes", "No"))
    smoking = st.radio("Do you smoke?", ("Yes", "No"))
    physical_activity = st.radio("Do you engage in regular physical activity?", ("Yes", "No"))
    bmi = st.number_input("Enter your Body Mass Index (BMI):", min_value=10.0, max_value=50.0)
    childbirth = st.radio("Have you had children?", ("Yes", "No"))
    breastfeeding = st.radio("Did you breastfeed your children?", ("Yes", "No"))
    diet = st.radio("Do you follow a healthy diet?", ("Yes", "No"))

    # Implement risk calculation logic (simplified example)
    risk_score = 0

    # Age factor
    if age > 50:
        risk_score += 2
    else:
        risk_score += 1

    # Family history factor
    if family_history == "Yes":
        risk_score += 3

    # Genetic testing factor
    if genetic_testing == "Yes":
        risk_score += 4

    # Hormone therapy factor
    if hormone_therapy == "Yes":
        risk_score += 2

    # Alcohol consumption factor
    if alcohol_consumption == "Yes":
        risk_score += 1

    # Smoking factor
    if smoking == "Yes":
        risk_score += 2

    # Physical activity factor
    if physical_activity == "No":
        risk_score += 1

    # BMI factor
    if bmi > 25:
        risk_score += 2

    # Childbirth factor
    if childbirth == "No":
        risk_score += 2

    # Breastfeeding factor
    if breastfeeding == "No":
        risk_score += 1

    # Diet factor
    if diet == "No":
        risk_score += 1

    # Display result
    if risk_score < 5:
        st.success("Your risk of developing breast cancer is low.")
        st.info("Note: A low risk score indicates that while your likelihood of developing breast cancer is low, regular screenings and a healthy lifestyle are still recommended.")
    elif 5 <= risk_score < 10:
        st.warning("Your risk of developing breast cancer is moderate.")
        st.info("Note: A moderate risk score suggests that you should consult with a healthcare professional for further evaluation and consider lifestyle modifications to lower your risk.")
    else:
        st.error("Your risk of developing breast cancer is high.")
        st.info("Note: A high risk score indicates a significant likelihood of developing breast cancer. Immediate consultation with a healthcare provider is strongly recommended, along with genetic testing or regular monitoring.")

    st.info(f"Total risk score: {risk_score}")

# Run the tool in the Streamlit app
if __name__ == '__main__':
    risk_assessment()
