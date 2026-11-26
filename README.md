# ZoopData_Expl
Import and explore zooplankton data gathered by the Environment Agency since 2022 around the English coast.

# To-do list
* Generate report investigating drivers of change in inshore zooplankton data
* What is the bigger picture?
  + are some regions different to others?
* consider regional and temporal differences in taxon diversity and density
* how do these parameters link to environmental pressures? (nutrients, temperature, salinity)
* Do relative proportions of zooplankton lifeforms show regional differences? (e.g., Large:Small copepods)
* **Update 'carbon content' variable names to specify units.  e.g., currently have names like 'md_carbTot_m3'.  Should be 'md_ug_carbTot_m3'**
* **look at variability of 'clustered' sites (Solent, Kent, (Northumberland?))**

# Approach
* Place returned analysis files in the folder location indicated in the `hldfol` object identified in the `R/set_meta.R` script (not uploaded to git)
* Run script: `00_imp_raw_lab_data.R` to extract returned lab analysis data and consolidate into a single `results_[date].csv` file within the `hldfol` folder
* Copy the extracted data from `results_[date].csv` into the data file at: `datfol/processedData/[Date]_MBA_Returns_Amalgamated_USE.xlsx`
  * Add the extracted data to the sheet `outR04` and ensure that lookups in columns from Col S onwards are filled
  * Drag the lookups in rows in sheet `outR04_LF` to include the values just pasted.  This will add lifeform, taxonomic, and carbon content information to the new data.
  * Search the `outR04_LF` data for #N/A values indicating failed lookups. Missing lookup values will need to be added to the relevant sheets:
    * `taxWide` - high resolution taxonomic information taken from the World Register of Marine Species
    * `TaxonomicRaw` - taxonomic information presenting names to display (deals with inconsistencies/typos on lab returns)
    * `LifeFormsUSE` - PLET Master lifeform lookup table, with additional lifeforms added for taxa not on master list
