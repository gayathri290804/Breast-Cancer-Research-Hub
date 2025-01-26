import streamlit as st
import pandas as pd
from openpyxl import load_workbook
import streamlit.components.v1 as components






def feed_back():
    
    st.write("""
    ## We Value Your Feedback! ✨
    
    Your feedback helps us improve the Breast Cancer Research Hub. Please take a moment to share your thoughts, suggestions, or any issues you encountered while using the app.

    ### How would you rate your experience?
    """)

    # Feedback Rating (Scale of 1 to 5)
    rating = st.slider("Rate the app", min_value=1, max_value=5, step=1)

    # Feedback message input
    feedback_message = st.text_area("Your feedback", "Type your feedback here...")

    if st.button("Submit Feedback"):
        # You can store this feedback in a file or database, but for now, we'll just display it.
        st.write(f"**Thank you for your feedback!**")
        st.write(f"**Your rating**: {rating}/5")
        st.write(f"**Your message**: {feedback_message}")

        # You can add logic to save this feedback to a file or database
        # For example: save_feedback_to_file(rating, feedback_message)





