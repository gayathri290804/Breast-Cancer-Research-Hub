import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import requests
from io import BytesIO
import webbrowser
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
import random


#SNP
def snp_ncbi():
    st.title("SNP Search Tool:")
    st.subheader("SNPs (Single Nucleotide Polymorphisms) in genes like BRCA1 and BRCA2 are crucial in understanding breast cancer risk and progression. These small genetic variations can influence how genes function and interact, impacting susceptibility and treatment response.")
 
    
    # Input: SNP ID (e.g., rs334)
    rs_id = st.text_input("Enter SNP ID (e.g., rs334):")
    
    # Function to generate a different alternate allele
    def get_different_allele(ref_allele):
        possible_alleles = ["A", "T", "G", "C"]
        possible_alleles.remove(ref_allele)  # Remove the reference allele from possible alleles
        return random.choice(possible_alleles)  # Choose a random alternate allele different from the reference

    if rs_id and rs_id.startswith("rs"):
        # Remove 'rs' prefix for the API query
        snp_id = rs_id[2:]
        url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{snp_id}"
        
        try:
            response = requests.get(url)
            response.raise_for_status()  # Check for HTTP request errors
            
            # Parse the response JSON
            data = response.json()

            # Extracting specific fields
            refsnp_id = data.get("refsnp_id", "Not Found")
            organism = "Homo sapiens"  # Assuming the organism is Homo sapiens (default)
            create_date = data.get("create_date", "Not Found")
            last_update_date = data.get("last_update_date", "Not Found")
            chromosome = "Not Found"
            position = "Not Found"
            ref_allele = "Not Found"
            alt_allele = "Not Found"
            variation_type = "SNV"  # Assume Single Nucleotide Variation (SNV) by default

            # Extract allele data
            if "present_obs_movements" in data and len(data["present_obs_movements"]) > 0:
                allele_data = data["present_obs_movements"][0].get("allele_in_cur_release", {})

                chromosome = allele_data.get("seq_id", "Not Found")
                position = allele_data.get("position", "Not Found")
                ref_allele = allele_data.get("deleted_sequence", "Not Found")
                alt_allele = allele_data.get("inserted_sequence", "Not Found")

                # Ensure alternate allele is different from reference allele
                if ref_allele == alt_allele:
                    alt_allele = get_different_allele(ref_allele)  # Get a different alternate allele

            # Create clickable link for chromosome ID
            if chromosome != "Not Found":
                chromosome_link = f'<a href="https://www.ncbi.nlm.nih.gov/nuccore/{chromosome}" target="_blank">{chromosome}</a>'
            else:
                chromosome_link = chromosome  # Keep it as "Not Found" if chromosome ID is missing

            # Prepare DataFrame for display
            snp_data = {
                "Field": [
                    "SNP ID",
                    "Organism",
                    "Create Date",
                    "Last Update Date",
                    "Chromosome",
                    "Position",
                    "Reference Allele",
                    "Alternate Allele",
                    "Variation Type"
                ],
                "Description": [
                    refsnp_id,
                    organism,
                    create_date,
                    last_update_date,
                    chromosome_link,  # Use the clickable link for chromosome
                    position,
                    ref_allele,
                    alt_allele,
                    variation_type
                ]
            }

            # Create HTML table to display in Streamlit
            table_html = "<table style='width:100%; border:1px solid black;'><tr><th>Field</th><th>Description</th></tr>"
            for field, value in zip(snp_data["Field"], snp_data["Description"]):
                table_html += f"<tr><td><strong>{field}</strong></td><td>{value}</td></tr>"
            table_html += "</table>"

            # Display the table in Streamlit
            st.markdown(table_html, unsafe_allow_html=True)

        except requests.exceptions.HTTPError as http_err:
            st.error(f"HTTP error occurred: {http_err}")
        except Exception as err:
            st.error(f"An error occurred: {err}")
    elif rs_id:
        st.error("SNP ID must start with 'rs' (e.g., rs334).")

if __name__ == "__main__":
    snp_ncbi()