"""DEPECRATED

This script scrapes the LvEDGE website at https://lvedge.bu.edu/.
A specific SPU ID is used to search for every corresponding gene.
The pfam domain of each gene is then printed to the console, preceded by the gene name.
This code also comes with an object that can be imported to perform the same functionality elsewhere.

Input: SPU ID (ex. SPU_003582)
Output: 'Gene name / SPU ID: pfam domain' for each gene, outputed to console

Copyright (C) <RELEASE_YEAR_HERE> Bradham Lab

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

Correspondence: Correspondence: 
    Cynthia A. Bradham - cbradham@bu.edu - *
    Anthony B. Garza   - abgarza@bu.edu  - **
    Stephanie P. Hao   - sphao@bu.edu    - **
    Yeting Li          - yetingli@bu.edu - **
    Nofal Ouardaoui    - naouarda@bu.edu - **

    *  - Principle Investigator
    ** - Software Developers
"""


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Imports
from selenium import webdriver
from selenium.webdriver.common.by import By
import requests
from io import StringIO
from Bio import SeqIO
import time


class LvEdgeScraper:
    """
    Object with functionality for scraping LvEDgE website at https://lvedge.bu.edu/cgi-bin/lvedge/main.py

    Attributes
    ----------
        url: str
            URL of LvEDGE website.
        is_headless: bool
            Run chrome in headless mode? Default is True.
        try_count: int
            Number of times to try to retrieve an spu ID from the website before giving up. Default is 10.

    Methods
    -------
        open_search: None
            Opens LvEDGE website to the search tab
        search: List[List[str]]
            Search for a specific SPU ID and return the gene table
        __del__: None
            Close the chrome window when the object is deleted
        close: None
            Easier name to use for user to close the chrome window than __del__ :)
    """

    url = 'https://lvedge.bu.edu/cgi-bin/lvedge/main.py'
    def __init__(self, is_headless = True, try_count = 10):
        """Constructor
        
        When built, the search tab is opened in the LvEDGE website automatically.
        """

        # Options for the chrome window
        self.chrome_options = webdriver.ChromeOptions()

        # If is_headless is True, then run chrome in headless mode, i.e., no window mode.
        if is_headless:
            self.chrome_options.add_argument('--headless')

        # Number of times to try to retrieve an spu ID
        self.try_count = try_count
        
        # Open the search tab
        self.open_search()
        

    def open_search(self):
        """
        Opens LvEDGE website to the search tab
        """
        
        # Open chrome and add options
        self.driver = webdriver.Chrome(options=self.chrome_options)

        # Go to the LvEDGE website home page
        self.driver.get(self.url)

        # Click on search tab
        search_tab = self.driver.find_element(By.LINK_TEXT, 'Search')
        search_tab.click()



    def search(self, spu_id):
        """
        Search for a specific SPU ID and return the gene table.

        Given that the search tab is open, the SPU ID is typed in and the search button is clicked.
        A list of genes would typically appear in the browser.
        This list of scraped and returned as a table (list of lists).

        Parameters
        ---------
            spu_id: str
                SPU ID to search for.
        
        Returns
        -------
            List[List[str]]
                Table of genes returned by search.
                First row is the header row.
        """

        r = [] # List to return
  
        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Search for SPU ID and get results
        try:
            # Try to get the gene table multiple times (query may fail)
            count = 0
            table_count = 0
            while count < self.try_count and table_count < 3:

                        
                # Grab the form to fill out the SPU ID
                spu_id_form = self.driver.find_element(By.NAME, 'SPU_id')

                # Type in SPU ID
                spu_id_form.send_keys(spu_id)

                # Button to submit SPU ID
                submit_button = self.driver.find_element(By.ID, 'search_genes')
                submit_button.click()

                # Get gene table, i.e., table with all of the information per gene
                gene_table_list = self.driver.find_elements(By.TAG_NAME, 'table')
                table_count = len(gene_table_list)

                count += 1

                if table_count < 3:
                    time.sleep(10)  # Sleep for 5 seconds before trying again

            if table_count < 3:
                return r
            
            # Get the gene table
            gene_table = gene_table_list[2]

            #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
            # For each gene returned by search, return the gene's row and include the fasta sequence instead of the link


            # Get all of the rows in the gene table
            gene_rows = gene_table.find_elements(By.TAG_NAME, 'tr')

            # Add to output variable 'r' each row
            for row in gene_rows:

                # Get all of the columns in the current row
                columns = row.find_elements(By.TAG_NAME, 'td')

                if len(columns) == 0:
                    continue

                fasta_link = columns[3].find_element(By.TAG_NAME, 'a').get_attribute('href')
                fasta_text = requests.get(fasta_link).text
                fasta = list(SeqIO.parse(StringIO(fasta_text), 'fasta'))

                # Transform each value inside the columns into a string
                text_columns = [column.text for column in columns]

                text_columns[3] = fasta

                # Add the row to the output variable 'r'
                r.append(text_columns)
        
        except Exception as e:
            print(e)
            input(spu_id)

        return r
    


    def __del__(self):
        """
        Close the chrome window when the object is deleted.
        """
        self.driver.quit()


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# STAND ALONE SCRIPT
if __name__ == "__main__":

    import sys

    if len(sys.argv) != 2 or sys.argv[1] in ['-h', '--help']:
        print('Usage: python lvedge_scraper.py <SPU ID>')
        sys.exit()

    spu_id = sys.argv[1]

    # Initialize LvEDGEScraper object
    scraper = LvEdgeScraper()

    # Search for SPU ID
    gene_rows = scraper.search(spu_id)

    # Get the PFAM Domain of each gene
    for row in gene_rows:

        gene_name = row[0]
        pfam_domain = row[9]
        pfam_domain = ', '.join([x for x in pfam_domain.split(',') if x != 'NA' and x != ''])

        print(f'{gene_name} / {spu_id}: {pfam_domain}')

