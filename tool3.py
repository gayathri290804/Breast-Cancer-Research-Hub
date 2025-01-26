import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import requests
from io import BytesIO
import webbrowser
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


def calculate_admet(smiles):
    """
    Calculate basic ADMET properties from a SMILES string using RDKit.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"Error": "Invalid SMILES string"}

        # Calculate molecular descriptors
        properties = {
            "Molecular Weight": round(Descriptors.MolWt(mol), 2),
            "LogP (lipophilicity)": round(Crippen.MolLogP(mol), 2),
            "TPSA (Topological Polar Surface Area)": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "Number of H-Bond Donors": Descriptors.NumHDonors(mol),
            "Number of H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
            "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
            "Number of Rings": Descriptors.RingCount(mol),
            "Absorption": round(100 - rdMolDescriptors.CalcTPSA(mol) * 0.1, 2),  # Simplified calculation
            "Distribution": "Moderate" if Crippen.MolLogP(mol) < 5 else "Low",
            "Metabolism": "Likely Hepatic" if Descriptors.NumHDonors(mol) + Descriptors.NumHAcceptors(mol) > 5 else "Minimal",
            "Excretion": "Renal" if Descriptors.MolWt(mol) < 500 else "Biliary",
            "Toxicity": "Low" if rdMolDescriptors.CalcTPSA(mol) < 140 else "High",
        }
        return properties
    except Exception as e:
        return {"Error": str(e)}

def swissadme_tool():
    """
    Streamlit interface for calculating ADMET properties from SMILES.
    """
    st.title("ADMET Property Predictor:")

    st.subheader("SwissADME and Lipinski's Rule play a pivotal role in breast cancer drug discovery by evaluating the pharmacokinetics, drug-likeness, and medicinal chemistry properties of potential compounds. SwissADME predicts absorption, distribution, metabolism, excretion, and toxicity (ADMET) parameters. Lipinski's Rule of Five ensures that drug candidates meet criteria such as molecular weight (<500 Dalton), hydrophobicity (LogP <5), and limits on hydrogen bond donors and acceptors, enhancing their suitability as oral drugs.")

    st.subheader("Input SMILES String")
    smiles_string = st.text_input("Enter SMILES String:(e.g., CC(=O)OC1=CC=CC=C1C(=O)O)")

    if st.button("Calculate ADMET Properties"):
        if smiles_string:
            # Calculate properties
            admet_properties = calculate_admet(smiles_string)

            # Display results
            if "Error" in admet_properties:
                st.error(admet_properties["Error"])
            else:
                st.subheader("Predicted ADMET Properties")
                for key, value in admet_properties.items():
                    st.write(f"**{key}:** {value}")
        else:
            st.warning("Please enter a valid SMILES string.")

if __name__ == "__main__":
    swissadme_tool()
