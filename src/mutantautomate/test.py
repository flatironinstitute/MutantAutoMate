import unittest
from combined_single import search_residue  # Import the function to be tested
import requests
import re

# Define the extract_pdb_ids function
def extract_pdb_ids(text):
    # Define a regular expression pattern to match "database":"PDB","id":"<>"
    pattern = r'"database":"PDB","id":"([^"]+)"'

    # Use regex to find all matches of the pattern
    matches = re.findall(pattern, text)

    return matches

class TestPDBOutput(unittest.TestCase):

    def test_pdb_output(self):
        # Input for testing (SHANK3 D26)
        gene_name = "SHANK3"
        residue_info = "D26"
        
        # Call the search_residue function
        matching_isoforms = search_residue(residue_info[0], int(residue_info[1:]), gene_name)

        # Check if matching_isoforms is not empty
        self.assertTrue(matching_isoforms, "Matching isoforms should not be empty")

        # Ensure that each value in matching_isoforms is 6 characters long by cutting them to 6 characters
        trimmed_isoforms = [isoform[:6] for isoform in matching_isoforms]

        # Assert that "Q9BYB0" is in the list of matching isoforms
        self.assertIn("Q9BYB0", trimmed_isoforms, "Expected 'Q9BYB0' in the list of matching isoforms")

        # Extract PDB IDs from the UniProt API for the first isoform
        if matching_isoforms:
            isoform_id = "Q9BYB0"
            print(isoform_id)
            url = f"https://rest.uniprot.org/uniprotkb/{isoform_id}"
            response = requests.get(url)
            if response.status_code == 200:
                response_text = response.text
                pdb_ids = extract_pdb_ids(response_text)
                self.assertTrue(pdb_ids, "Extracted PDB IDs should not be empty")
                print("PDB IDs extracted from UniProt:")
                for pdb_id in pdb_ids:
                    print(pdb_id)

                # Check if "6CPK" is in the list of PDB IDs
                self.assertIn("6CPK", pdb_ids, "Expected '6CPK' in the list of PDB IDs")
            else:
                self.fail(f"Failed to retrieve UniProt data. Status code: {response.status_code}")

if __name__ == '__main__':
    unittest.main()
