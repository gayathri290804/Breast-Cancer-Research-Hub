import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import requests
from io import BytesIO
import webbrowser
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors



#clustral omega
def clustalo():
    st.title("Multiple Sequence Alignment:")
    st.subheader("Clustal Omega is a tool used in breast cancer research to align protein sequences, helping identify conserved regions and mutations in cancer-related proteins like HER2 or estrogen receptors. It aids in understanding structural differences that influence cancer progression and drug resistance, supporting the development of targeted therapies.")


    # User inputs
    st.write("Enter your details and sequence information below.")

    # Input email
    email = st.text_input("Enter your email address", placeholder="example@example.com")

    # Select sequence type (Protein/DNA)
    sequence_type = st.selectbox("Select sequence type", ["protein", "dna"])

    # Upload file or input sequences manually
    st.write("Upload a file containing sequences in FASTA format or input the sequences directly.")
    uploaded_file = st.file_uploader("Upload FASTA File", type=["txt", "fasta"])
    
    fasta_sequences = ""
    if uploaded_file is not None:
        fasta_sequences = uploaded_file.read().decode("utf-8")
        st.text_area("FASTA Sequences", value=fasta_sequences, height=200, disabled=True)
    else:
        fasta_sequences = st.text_area("Input FASTA Sequences", placeholder="Enter sequences in FASTA format", height=200)

    # Submit for alignment
    if st.button("Run Clustal Omega Alignment"):
        if not email.strip():
            st.error("Please provide a valid email address.")
        elif not fasta_sequences.strip():
            st.error("Please provide FASTA sequences for alignment.")
        else:
            st.info("Submitting sequences to Clustal Omega API...")
            try:
                # Clustal Omega API Endpoint
                url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run/"

                params = {
                    "email": email,
                    "sequence": fasta_sequences,
                    "stype": sequence_type,
                }
                response = requests.post(url, data=params)

                if response.status_code == 200:
                    job_id = response.text.strip()
                    st.success(f"Job submitted successfully. Job ID: {job_id}")

                    # Create clickable link for the results page
                    result_url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/out"
                    st.markdown(f"[Click here to view the results](https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{job_id}/out)")

                    # Optionally fetch results automatically
                    st.info("Fetching alignment results...")
                    result_response = requests.get(result_url)

                    if result_response.status_code == 200:
                        alignment_result = result_response.text
                        st.subheader("Alignment Results")
                        st.text_area("Aligned Sequences", alignment_result, height=300)
                    else:
                        st.warning("The results may not be ready yet. Use the link above to check later.")
                else:
                    st.error("Failed to submit job to Clustal Omega API. Please try again.")
            except Exception as e:
                st.error(f"An error occurred: {e}")

if __name__ == "__main__":
    clustalo()

