import streamlit as st
import sqlite3
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

### SQL Database Scheme query

csv_path = r'cell-count.csv'
db_path = r'Loblaw_Bio_cell_count.db'
if os.path.exists(db_path):
    os.remove(db_path)
    
st.set_page_config(page_title='Loblaw Bio: Immune Cell Population Analysis')
st.title('Loblaw Bio: Clinical Trial Dashboard')

### Part 1 ================================================================================================================================


def initialize_database(csv_path, db_path):
	conn = sqlite3.connect(db_path)
	cursor = conn.cursor()
	cursor.execute('PRAGMA foreign_keys = ON;')

	df_raw = pd.read_csv(csv_path)
	df_raw.to_sql('cell-count', conn, if_exists='replace', index=False)

	sql_script = """
	/* =========================================================================================================
	DATA CHECK
	Preliminary checks to ensure data robustness. Commented out in final product to keep execution clean.
	=========================================================================================================
	*/

	/*
	-- Ensure Primary Key has no duplicates
	SELECT sample, COUNT(*) 
	FROM "cell-count" 
	GROUP BY sample
	HAVING COUNT(*) > 1;

	-- Ensure no NULLs in data for critical columns
	SELECT COUNT(*) AS missing_ids
	FROM "cell-count"
	WHERE subject IS NULL OR subject = '';
	*/

	/* 
	=========================================================================================================
	PRODUCTION SCHEMA 
		Table cleaning, separation, and normalization
	=========================================================================================================
	*/

	-- Subject Table Processing
	-- Define new table
	DROP TABLE IF EXISTS subject_table;
	CREATE TABLE subject_table (
		subject TEXT PRIMARY KEY,
		condition TEXT,
		age INTEGER,
		sex TEXT,
		treatment TEXT,
		response TEXT
	);

	-- Data insertion, select only distinct values, and normalize subject primary key via padding
	INSERT INTO subject_table (subject, condition, age, sex, treatment, response)
	SELECT DISTINCT 'sbj' || PRINTF('%04d', CAST(SUBSTR(subject, 4) AS INTEGER)), condition, age, sex, treatment, response
	FROM "cell-count";

	/*
	-- sanity check for duplicates
	SELECT 
		COUNT(*) AS total_rows, 
		COUNT(DISTINCT subject) AS unique_subjects
	FROM subject_table;
	*/


	-- Project Table Processing

	-- Check to see if project corresponds to one label each
	/*
	SELECT 
		project, 
		COUNT(DISTINCT project) AS label_count,
		GROUP_CONCAT(DISTINCT sample_type) AS labels_found
	FROM "cell-count"
	GROUP BY project;
	*/

	DROP TABLE IF EXISTS project_table;
	CREATE TABLE project_table (
		project TEXT PRIMARY KEY,
		sample_type TEXT
	);

	INSERT INTO project_table (project, sample_type)
	SELECT DISTINCT project, sample_type
	FROM "cell-count";


	-- Sample Table Processing

	-- Create and define sample table
	DROP TABLE IF EXISTS sample_table;
	CREATE TABLE sample_table (
		sample TEXT PRIMARY KEY,
		subject TEXT,
		project TEXT,
		time_from_treatment_start INTEGER,
		b_cell REAL,
		cd8_t_cell REAL,
		cd4_t_cell REAL,
		nk_cell REAL,
		monocyte REAL,
		FOREIGN KEY (subject) REFERENCES subject_table(subject),
		FOREIGN KEY (project) REFERENCES project_table(project)
	);

	INSERT INTO sample_table (sample, subject, project, time_from_treatment_start, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte)
	SELECT sample, 'sbj' || PRINTF('%04d', CAST(SUBSTR(subject, 4) AS INTEGER)), project, time_from_treatment_start, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte
	FROM "cell-count";
	"""

	cursor.executescript(sql_script)
	conn.commit()
	print(f'Database scheme successful')
	return conn

conn = initialize_database(csv_path, db_path)

### Part 2 ================================================================================================================================
#Retrieve table from database
df_samples = pd.read_sql_query('SELECT * FROM sample_table', conn)
df_samples_total = df_samples.copy()
df_samples_total['total_count'] = df_samples_total[['b_cell', 'cd8_t_cell', 'cd4_t_cell', 'nk_cell', 'monocyte']].sum(axis = 1)

#Melt relevant columns and some minor edits for intuitiveness
melt_cols = ['b_cell', 'cd8_t_cell', 'cd4_t_cell', 'nk_cell', 'monocyte']
df_total_melt = pd.melt(df_samples_total, 
                        id_vars = [c for c in df_samples_total.columns if c not in melt_cols],
                        value_vars = melt_cols,
                        var_name='population',
                        value_name='count'
                        )

df_total_melt = df_total_melt.sort_values(by='sample')
df_total_melt = df_total_melt.reset_index(drop=True)
df_total_melt['percentage'] = df_total_melt['count'] / df_total_melt['total_count']
summary_table = df_total_melt[['sample', 'total_count', 'population', 'count', 'percentage']]


### Part 3 ================================================================================================================================
#Retrieve table from database
with sqlite3.connect('Loblaw_Bio_cell_count.db') as conn:
    df_subjects = pd.read_sql_query('SELECT * FROM subject_table', conn)
    df_project = pd.read_sql_query('SELECT * FROM project_table', conn)

df_master = pd.merge(df_total_melt, df_subjects, on='subject', how='inner')
df_master = pd.merge(df_master, df_project, on='project', how='inner')

# Filtering according to requirements
filtered_df = df_master[
    (df_master['treatment'] == 'miraclib') & 
    (df_master['condition'] == 'melanoma') &
    (df_master['sample_type'] == 'PBMC')
]

# Analysis wrapped with Streamlit for dashboard. 
tab1, tab2, tab3, tab4 = st.tabs([
    'Project Documentation', 
    'Initial Analysis', 
    'Statistical Analysis', 
    'Data Subset Analysis'
])

with tab2:
	st.header('Part 2: Initial Analysis - Data Overview')
	st.dataframe(summary_table, width='stretch')
     
	st.markdown("""
	### Notes
             Similar to the logic behind keeping 'sbj' included in the numbers, 'sample' is still kept in the values to allow for a more flexible scale-up, 
             as there is no guarantee that the prefix used will always be 'sample', especially when introducing data/samples from external vendors/labs when scaling up. 
             Keeping 'sample' in the value will help distinguish between different prefix labels in the future while still allowing the values to be sorted in numerical order.
	""")
    
with tab3:
    st.header('Part 3: Statistical Analysis')
    
    # We remove the columns here so the plot takes the full width
    # Row 1: The Plot
    st.subheader('Population Distribution Boxplot')
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.boxplot(
        data=filtered_df, 
        x='population', 
        y='percentage', 
        hue='response',
        palette={'yes': '#1ff802', 'no': '#ff1900'},
        ax=ax
    )
    plt.xlabel("Immune Cell Population", fontsize=12)
    plt.ylabel("Frequency (Relative %)", fontsize=12)
    st.pyplot(fig)

    st.divider()

    # Row 2: Summary Table and Explanation
    # We create new columns here to place the table and text side-by-side below the plot
    col_table, col_text = st.columns([1, 1])

    with col_table:
        st.subheader('Statistical Significance Summary')

        populations = filtered_df['population'].unique()
        stats_results = []

        for cell in populations:
            response_yes = filtered_df[(filtered_df['population'] == cell) & (filtered_df['response'] == 'yes')]
            response_no = filtered_df[(filtered_df['population'] == cell) & (filtered_df['response'] == 'no')]

            t_stat, p_val = stats.ttest_ind(response_yes['percentage'], response_no['percentage'], equal_var=False)

            verdict = 'Significant' if p_val < 0.05 else 'NOT Significant'

            stats_results.append({
                'Population': cell,
                'p_value': p_val,
                'verdict': verdict
            })

        stats_df = pd.DataFrame(stats_results)

        st.dataframe(
            stats_df,
            hide_index=True,
            width='stretch',
            column_config={
                "p_value": st.column_config.NumberColumn("p-value", format="%.4f"),
                "verdict": st.column_config.TextColumn("Verdict")
            }
        )

    with col_text:
        st.markdown("### Explanation of Results")
        st.write("""
        The boxplot above visualizes the distribution of cell frequencies for responders (green)
        versus non-responders (red). The table summarizes Welch's t-tests for each population.

        Based on the results, the cd4_t_cell population is statistically significant, while
        b_cell, nk_cell, monocyte, and cd8_t_cell are not.
        This suggests cd4_t_cell may be a promising biomarker for further investigation.
        """)

### Part 4 ================================================================================================================================

with tab4:
    st.header('Part 4: Baseline Data Subset Analysis (t=0)')
    # Filtering based on problem statement
    time_zero_df = filtered_df[(filtered_df['time_from_treatment_start'] == 0)]
    unique_subjects_df = time_zero_df.drop_duplicates(subset=['subject'])

    c1, c2, c3 = st.columns(3)

    with c1:
        st.write('Samples per Project')
        unique_samples_df = time_zero_df.drop_duplicates(subset=['sample'])
        st.write(unique_samples_df['project'].value_counts())

    with c2:
        st.write('Responders vs Non-Responders')
        st.write(unique_subjects_df['response'].value_counts())

    with c3:
        st.write('Male vs Female (Baseline)')
        st.write(unique_subjects_df['sex'].value_counts())

    st.divider()

    st.subheader('Extra Problem Result')
    

    ### Extra Problem
    extra_problem_df = time_zero_df[(time_zero_df['population'] == 'b_cell')]
    extra_problem_df = extra_problem_df[(extra_problem_df['sex'] == 'M')]
    extra_problem_df = extra_problem_df[(extra_problem_df['response'] == 'yes')]
    mean_b_cell = extra_problem_df['count'].mean()
    st.info(f'Average B-cells (2 sig figs): {mean_b_cell:.2g}')
    

with tab1:
    st.header('Project Documentation & Rationale')
    st.markdown("""
	# Teiko_Loblaw_Bio

	## How to run
	Dependencies are listed in requirements.txt. They can be installed by running ```pip install -r requirements.txt``` in the terminal.
	Please run the code through ```streamlit run Analysis.py```

	Dashboard link: https://teikoloblawbio-wbwtpe8un8du2xwf7usoy9.streamlit.app/

	## Relational Database Schema and Rationale
	Given the cell-count.csv file, the raw data is separated into three tables through SQL: 
	- subject_table: Stores unique patient metadata (subject [Primary Key], Condition, Age, Sex, Treatment, Response).
	- project_table: Stores unique project identifiers and associated sample types (project [Primary Key], sample type).
	- sample_table: Stores individual sample observations, including timepoints and cell counts. It uses Foreign Keys to link back to the parent tables.

	This scheme has normalized tables that fulfill the requirements of the third normal form (3NF) tables. As a result, the scheme avoids redundancy, ensuring that any data is stored only once and is connected between tables through foreign keys, linking every sample record to its respective subject and project.

	This efficient layout scales well, as 3NF structures help minimize data storage requirements by eliminating the repetition of metadata, leaving more space for high-volume sample observations. Furthermore, this allows any data correction to happen at one point, instead of multiple, such as a correction on subject data would only require an edit in the subject_table once, instead of multiple times in the original table. Lastly, should the company introduce additional products/services, the new data can simply be tied to the subject data to keep consistent data between services. 

	### Notes on data management rationale:
	- 'sbj' was purposefully left behind and not cleaned (sbj0001 -> 0001) with scaling in mind. To account for more data entries in the future, where 'sbj' is not guaranteed and clinicians can assign a different prefix to the subjects, 'sbj' was not removed to avoid ambiguity and allow such changes when scaling up.
	- Subject values are padded (sbj001 -> sbj0001) to ensure alphanumeric sorting and create a more intuitive table organization for clinicians/clients to review. Depending on the scale of the data, the amount of padding could be easily or dynamically adjusted to accommodate the larger number of subjects.
	- Sanity and data integrity checks are provided but commented out in the query to ensure efficient execution.
	- Although the current analysis retrieves data from all tables for comprehensive analysis, when the number of tables and samples increases in the future, this scheme is still optimal to enable selective retrieval and reduces memory overhead compared to processing the raw table. 

	Given the flexibility of the modified database schema, there are a few interesting analyses that might be possible as the data sama sample size scales up, such as:
	- Changes in cell count over time as new treatment is introduced at various points in time.
	- Cross-comparison between old project data and newer projects covering different treatments or services.
	- Training a predictive model to predict response based on cell counts and other newly gathered features.

	## Code Structure

	The code features a hybrid approach, where the relational database is built through SQL prompts, while the table retrieval and data analysis are done through Python (pandas), and the interface design was done through Streamlit. 

	- Data Management (SQL): The cell-count.csv raw data was processed by a SQL prompt, separating it into three tables (mentioned above) and storing them in a .db file to help enforce Relational Integrity through Foreign Key constraints. SQL provides an efficient storage and retrieval process, which will be crucial as the database size scales up. 
	- Analytical Layer (Python/Pandas): Afterward, for each part of the problem given, I retrieved the necessary and relevant data from the database and converted each to a pandas dataframe. Further processing is done through pandas in Python to tailor to the needs of every question. The code is dynamically coded to ensure that further additions or edits would be a simple change, such as by having my analysis logic (melting, grouping, and calculating frequencies) driven by the column headers and unique values present in the database. The boxplot for part 3 is created through seaborn, and the significance test is done through Welch's T-test that iterates through all available cell populations, setting the significance threshold at <0.05, a commonly used value. 
	- Presentation (Streamlit): Lastly, all the analysis was wrapped into a presentable format for Yah D’yada through Streamlit. This helps create an interactive dashboard to organize and show all the analysis results, highlighting only the important results without all the analysis clutter. 
	""")
