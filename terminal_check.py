import pandas as pd
import matplotlib.pyplot as plt
import re
from fpdf import FPDF

# Load Excel file using pandas
df = pd.read_excel('/Users/asameerpradhan/Downloads/ter_check.xlsx')

# Filter rows where the "New_residue" column contains "TER", "terminal", "ter", or "Ter"
filtered_df = df[df['New_residue'].astype(str).str.contains(r'TER|terminal|ter', case=False, na=False)]

# Count the number of occurrences of "Terminal" in the filtered dataframe
terminal_count = filtered_df['New_residue'].astype(str).str.count(r'TER|terminal|ter', flags=re.IGNORECASE).sum()

# Create PDF
pdf = FPDF()
pdf.add_page()
pdf.set_font("Arial", size=12)
pdf.cell(200, 10, txt=f"{terminal_count} Terminal residues found", ln=True)
pdf.output("terminal_count.pdf")

print(f'Number of times "Terminal" has appeared: {terminal_count}')
