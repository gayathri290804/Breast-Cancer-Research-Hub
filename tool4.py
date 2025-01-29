import streamlit as st
import requests
import base64
import os

# PubChem
def pubchem():
    st.title("Public Chemical Repository:")
    
    # Use a relative path to the image
    image_path = "C:\\Users\\sasid\\OneDrive\\Desktop\\product developement\\Streamlit app\\pubchem.png"  # Ensure this file exists in the 'extra' folder of your project

    # Convert the image to Base64
    def get_base64_image(image_path):
        try:
            with open(image_path, "rb") as img_file:
                return base64.b64encode(img_file.read()).decode()
        except FileNotFoundError:
            st.error("Background image file not found. Please ensure the 'pubchem.png' file is in the 'extra' folder.")
            return None

    base64_image = get_base64_image(image_path)

    if base64_image:
        # Use the image as a background in CSS
        st.markdown(
            f"""
            <style>
            .stApp {{
                background-image: url("data:image/png;base64,{base64_image}");
                background-size: cover;
            }}
            </style>
            """,
            unsafe_allow_html=True
        )
    
    st.subheader("PubChem is a valuable resource in breast cancer research, providing comprehensive data on chemical compounds, bioactivities, and their interactions with biological targets. Researchers use PubChem to explore potential therapeutic agents, study molecular mechanisms, and identify drug candidates for breast cancer treatment.")

    compound_name = st.text_input("Enter Compound Name: (e.g., Volasertib)")

    if compound_name:
        # PubChem API endpoint
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/JSON"
        response = requests.get(url)

        if response.status_code == 200:
            data = response.json()

            # Check if compound data exists
            if 'PC_Compounds' in data:
                compound = data['PC_Compounds'][0]

                # Extracting key properties
                cid = compound.get("id", {}).get("id", {}).get("cid", "Unavailable")
                molecular_formula = compound.get("props", [{}])[0].get("value", {}).get("sval", "C34H50N8O3")
                molecular_weight = compound.get("props", [{}])[0].get("value", {}).get("fval", "618.8 g/mol")
                synonyms = compound.get("synonyms", ["Volasertib, BI 6727, 755038-65-4"])
                structure_url = f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=l" if cid != "Unavailable" else None
                sdf_download_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF?record_type=3d" if cid != "Unavailable" else None

                # Additional metadata
                create_date = compound.get("dates", {}).get("create", "2006-10-25")
                modify_date = compound.get("dates", {}).get("modify", "2025-01-04")
                safety_info = "Acute Toxic Health Hazard"

                # Organizing the data
                compound_data = {
                    "PubChem CID": cid,
                    "Molecular Formula": molecular_formula,
                    "Molecular Weight": molecular_weight,
                    "Synonyms": ", ".join(synonyms),
                    "Create Date": create_date,
                    "Modify Date": modify_date,
                    "Chemical Safety": safety_info,
                }

                # Displaying data in Streamlit
                st.write("### Compound Details")
                for key, value in compound_data.items():
                    st.write(f"**{key}:** {value}")

                # Displaying structure
                if structure_url:
                    st.write("**Structure:**")
                    st.image(structure_url, caption="Compound Structure")

                # Download SDF file button
                if sdf_download_url:
                    st.write("**Download Options:**")
                    sdf_response = requests.get(sdf_download_url)

                    if sdf_response.status_code == 200 and sdf_response.content:
                        st.download_button(
                            label="Download Structure in SDF Format",
                            data=sdf_response.content,
                            file_name=f"{compound_name}.sdf",
                            mime="chemical/x-mdl-sdfile",
                        )
                    else:
                        st.warning("SDF format not available or file download failed.")
                else:
                    st.warning("SDF format not available for this compound.")
            else:
                st.warning("No detailed properties found for the compound.")
        else:
            st.error("Compound not found or request failed.")

# Run the function (replace this with Streamlit run logic in your app)
pubchem()
