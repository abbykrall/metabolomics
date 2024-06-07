# %%
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import os
from pandastable import Table, TableModel
import sys
import xlsxwriter

# Dynamically construct the class_ion_data_path using the current working directory
base_dir = os.path.dirname(os.path.abspath(__file__))  # Get the directory of the current script
class_ion_data_path = os.path.join(base_dir, 'mets_w_classes_hek_model.xlsx')


# Screen size variable width x height in px
app_size = "1600x1100"


# Scoring contribution of variables
peak_area_impact = 0.45  # Positive impact
iqr_diff_impact = 0.1  # Negative impact
variability_impact = 0.2  # Negative impact
outlier_impact = 0.3  # Negative impact


def save_figure_as_image(fig):
    file_path = filedialog.asksaveasfilename(
        defaultextension='.svg', 
        filetypes=[
            ('SVG files', '*.svg'), 
            ('PNG files', '*.png'), 
            ('JPEG files', '*.jpg'),
            ('PDF files', '*.pdf')
        ]
    )
    if file_path:
        file_ext = os.path.splitext(file_path)[1].lower()
        if file_ext in ['.png', '.jpg', '.jpeg']:
            fig.savefig(file_path, format=file_ext[1:], dpi=300)  # High resolution for raster formats
        elif file_ext == '.pdf':
            fig.savefig(file_path, format='pdf')
        else:
            fig.savefig(file_path, format='svg')
        messagebox.showinfo("Save Image", f"Figure saved as {file_ext.upper()} file: {file_path}")


class MetaboliteAnalysisApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Metabolite Standard Peak Area Analysis")
        self.geometry(app_size)
        
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill='both', expand=True)

        # Load stored data at initialization
        self.stored_data = self.load_stored_data()

        # Load class and ion data
        self.class_ion_data = self.load_class_ion_data(class_ion_data_path)
        self.create_class_ion_mappings()

        self.setup_met_std_peak_area_check_ui()

        self.protocol("WM_DELETE_WINDOW", self.on_closing)


    def on_closing(self):
        self.quit()
        self.destroy()
        sys.exit(0)
    

    def apply_custom_theme_for_table(self, table):
        options = {
            'cellbackgr': '#e6e6e6',
            'cellforegr': '#000000',
            'rowselectedcolor': '#d9d9d9',
            'colsselectedcolor': '#d9d9d9',
            'grid_color': '#a6a6a6',
            'text_color': '#000000',
            'background_color': '#ffffff',
            'font': ('Arial', 12),
            'rowheight': 25,
            'colheadercolor': '#f2f2f2',
            'floatprecision': 2,
            'align': 'center'
        }
        for key, value in options.items():
            setattr(table, key, value)
        table.redraw()


    def set_column_widths(self, table, df):
        """Set column widths based on the content of the dataframe."""
        for col in df.columns:
            max_len = max(df[col].astype(str).apply(len).max(), len(col)) + 2  # Add padding
            table.columnwidths[col] = max_len * 8  # Adjust the multiplier as needed
        table.redraw()

    def setup_met_std_peak_area_check_ui(self):
        self.tab_peak_area_check = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_peak_area_check, text='Std Peak Area Check')

        self.file_path_frame_peak_area = tk.Frame(self.tab_peak_area_check)
        self.file_path_frame_peak_area.grid(row=0, column=0, sticky='ew', pady=15)

        self.file_path_entry_peak_area = tk.Entry(self.file_path_frame_peak_area)
        self.file_path_entry_peak_area.grid(row=0, column=0, sticky='ew')
        self.upload_button_peak_area = tk.Button(self.file_path_frame_peak_area, text="Upload Met Std Peak Area Data", command=self.upload_std_peak_area_file)
        self.upload_button_peak_area.grid(row=0, column=1, padx=5)

        self.save_button_peak_area = tk.Button(self.file_path_frame_peak_area, text="Save Uploaded Data", command=self.save_uploaded_data)
        self.save_button_peak_area.grid(row=0, column=2, padx=5)

        self.std_met_peak_area_stats = tk.Frame(self.tab_peak_area_check)
        self.std_met_peak_area_stats.grid(row=1, column=0, sticky='nsew', pady=15)

        self.table_frame = tk.Frame(self.std_met_peak_area_stats)
        self.table_frame.grid(row=0, column=0, sticky='nsew')

        # Create scrollbars
        self.vsb = tk.Scrollbar(self.table_frame, orient="vertical", command=self._on_vertical_scroll)
        self.hsb = tk.Scrollbar(self.table_frame, orient="horizontal", command=self._on_horizontal_scroll)

        # Initialize pandastable with scrollbars
        self.table = Table(self.table_frame, dataframe=pd.DataFrame(), showtoolbar=False, showstatusbar=False, editable=False)

        # Attach scrollbars to the table
        self.table.configure(yscrollcommand=self.vsb.set, xscrollcommand=self.hsb.set)
        self.vsb.grid(row=0, column=1, sticky='ns')
        self.hsb.grid(row=1, column=0, sticky='ew')

        # Apply custom theme
        self.apply_custom_theme_for_table(self.table)

        self.table.show()

        # Bind the double-click event
        self.table.bind("<Double-Button-1>", self.on_double_click)

        self.class_groups = {}
        self.ion_groups = {}

        # Configure grid weights
        self.file_path_frame_peak_area.columnconfigure(0, weight=1)  # Make the first column expandable
        self.tab_peak_area_check.grid_rowconfigure(1, weight=1)
        self.tab_peak_area_check.grid_columnconfigure(0, weight=1)
        self.std_met_peak_area_stats.grid_rowconfigure(0, weight=1)
        self.std_met_peak_area_stats.grid_columnconfigure(0, weight=1)

        # Summary frame
        self.summary_frame = tk.Frame(self.tab_peak_area_check)
        self.summary_frame.grid(row=2, column=0, sticky='ew', pady=10)

    
    def save_uploaded_data(self):
        if not hasattr(self, 'df'):
            messagebox.showerror("Error", "No data has been uploaded.")
            return

        popup = tk.Toplevel()
        popup.title("Save Uploaded Data")
        popup.geometry("300x200")

        tk.Label(popup, text="Enter Date (YYYYMMDD):").pack(pady=5)
        date_entry = tk.Entry(popup)
        date_entry.pack(pady=5)

        tk.Label(popup, text="Enter Column ID (##), only the number:").pack(pady=5)
        col_id_entry = tk.Entry(popup)
        col_id_entry.pack(pady=5)

        def save():
            date = date_entry.get()
            col_id = col_id_entry.get()
            if not date or not col_id:
                messagebox.showerror("Error", "Both fields are required.")
                return

            filename = f"{date}_col_id_{col_id}_hek_area_edited.xlsx"
            filepath = os.path.join(base_dir, 'hek_stored_runs', filename)

            # Save the DataFrame to Excel with the specified sheet name
            with pd.ExcelWriter(filepath, engine='xlsxwriter') as writer:
                self.df.to_excel(writer, sheet_name='PoolAfterDF', index=True)

            messagebox.showinfo("Success", f"File saved as {filename}")
            popup.destroy()

        tk.Button(popup, text="Save", command=save).pack(pady=20)


    def _on_vertical_scroll(self, *args):
        self.table.yview(*args)

    def _on_horizontal_scroll(self, *args):
        self.table.xview(*args)

    def load_class_ion_data(self, filepath):
        class_ion_data = pd.read_excel(filepath)
        class_ion_data.columns = class_ion_data.columns.str.strip()
        class_ion_data['metabolite'] = class_ion_data['metabolite'].str.replace('|', ';').str.replace('"', '').str.strip()
        class_ion_data['class'] = class_ion_data['class'].str.replace('"', '').str.strip()
        class_ion_data['ion'] = class_ion_data['ion'].str.replace('"', '').str.strip()
        class_ion_data.set_index('metabolite', inplace=True)
        return class_ion_data

    def create_class_ion_mappings(self):
        self.class_mapping = self.class_ion_data['class'].to_dict()
        self.ion_mapping = self.class_ion_data['ion'].to_dict()

    def clean_dataframe_index(self, df):
        df.index = df.index.str.replace('|', ';').str.replace('"', '').str.strip()
        return df


    def get_color_tag(self, score):
        try:
            score = float(score)
        except ValueError:
            return None
        if score > 0.5:
            return 'lightgreen'
        elif score < 0:
            return 'lightcoral'
        else:
            return 'yellow'


    def upload_std_peak_area_file(self):
        analysis_fpath = filedialog.askopenfilename()
        if not analysis_fpath:
            return

        self.file_path_entry_peak_area.delete(0, tk.END)
        self.file_path_entry_peak_area.insert(0, analysis_fpath)

        if self.is_excel_file(analysis_fpath):
            xls = pd.ExcelFile(analysis_fpath)
            if 'PoolAfterDF' in xls.sheet_names:
                self.df = pd.read_excel(analysis_fpath, sheet_name='PoolAfterDF', index_col='Compound')
                self.df = self.clean_dataframe_index(self.df)  # Clean the DataFrame index
                self.df_normalized = self.df.div(self.df.loc['trifluoromethanesulfonate'])
                self.update_table()
            else:
                messagebox.showerror("Error", "'PoolAfterDF' sheet not present in the Excel file.")
        else:
            messagebox.showerror("Error", "Selected file is not a valid Excel file.")

    def is_excel_file(self, fpath):
        return fpath.endswith(('.xls', '.xlsx'))
    
    def load_stored_data(self):
        directory = os.path.join(base_dir, 'hek_stored_runs')
        stored_data = {}
        if not os.path.exists(directory):
            print(f"Directory not found: {directory}")
            return {}  # Return an empty dictionary if the directory does not exist
        for file in os.listdir(directory):
            if file.endswith('.xlsx'):
                path = os.path.join(directory, file)
                try:
                    df_original = pd.read_excel(path, index_col='Compound')
                    if 'trifluoromethanesulfonate' in df_original.index:
                        df_normalized = df_original.div(df_original.loc['trifluoromethanesulfonate'])
                        date_label, col_id = file.split('_')[0], file.split('_')[3]  # Extract date and col_id
                        stored_data[date_label] = (df_original, df_normalized, col_id)
                except Exception as e:
                    print(f"Failed to load {file}: {e}")
        return stored_data
    
    def calculate_rsd(self, data):
        """Calculate the Relative Standard Deviation (RSD), handling cases where the mean is zero."""
        mean = data.mean()
        if mean == 0:
            return 0  # Return 0 or some other appropriate value instead of NaN
        return (data.std() / mean) * 100
        
    def calculate_iqr_and_range(self, data):
        """Calculates the interquartile range and min-max of the given data."""
        q75, q25 = np.percentile(data, [75 ,25])
        iqr = q75 - q25
        data_min = np.min(data)
        data_max = np.max(data)
        min_max_range = f"{data_min:.2e}-{data_max:.2e}"
        return iqr, min_max_range
    
    def calculate_scores(self):
        # Calculate scores for original and normalized data
        original_scores = self.calculate_data_scores(self.df, normalized=False)
        normalized_scores = self.calculate_data_scores(self.df_normalized, normalized=True)

        # Normalize the scores
        # original_scores = self.normalize_scores(original_scores)
        # normalized_scores = self.normalize_scores(normalized_scores)

        # Combine original and normalized scores into a single dictionary
        scores = {compound: (original_scores.get(compound, 'ND'), normalized_scores.get(compound, 'ND'))
                for compound in self.df.index}
        return scores

    def calculate_data_scores(self, df, normalized):
        stored_means = {}
        stored_iqrs = {}
        internal_standard = 'trifluoromethanesulfonate'

        for date_label, (df_orig, df_norm, col_id) in self.stored_data.items():
            df_to_use = df_norm if normalized else df_orig
            for compound in df_to_use.index:
                if compound not in stored_means:
                    stored_means[compound] = []
                    stored_iqrs[compound] = []
                values = df_to_use.loc[compound].values.flatten()
                stored_means[compound].extend(values)
                q75, q25 = np.percentile(values, [75, 25])
                iqr = q75 - q25
                stored_iqrs[compound].append(iqr)

        # Calculate mean and IQR of the stored data
        stored_means = {compound: np.mean(values) for compound, values in stored_means.items()}
        stored_iqrs = {compound: np.mean(values) for compound, values in stored_iqrs.items()}

        scores = {}
        for compound in df.index:
            if compound in stored_means:
                values = df.loc[compound].values.flatten()
                mean_peak_area = np.mean(values)
                iqr_diff = abs(np.percentile(values, 75) - np.percentile(values, 25) - stored_iqrs[compound])

                # Calculate relative measures
                mean_peak_area_score = mean_peak_area / stored_means[compound] if stored_means[compound] != 0 else mean_peak_area
                iqr_diff_score = iqr_diff / stored_iqrs[compound] if stored_iqrs[compound] != 0 else iqr_diff

                # Internal standard adjustment for normalized data
                if normalized and internal_standard in df.index:
                    internal_standard_ratio = df.loc[internal_standard].mean() / self.stored_data_mean(internal_standard, normalized)
                    mean_peak_area_score *= internal_standard_ratio

                # Calculate variability impact
                variability_impact = self.calculate_variability_impact(values)

                # Detect and penalize outliers
                outlier_impact = self.detect_outliers(values)

                # Combine the scores with weights: mean positively, IQR difference negatively, variability negatively, outliers negatively
                peak_area_contrib = mean_peak_area_score * peak_area_impact  # Positive impact
                iqr_diff_contrib = -iqr_diff_score * iqr_diff_impact  # Negative impact
                variability_contrib = -variability_impact * variability_impact  # Negative impact
                outlier_contrib = -outlier_impact * outlier_impact  # Negative impact

                scores[compound] = peak_area_contrib + iqr_diff_contrib + variability_contrib + outlier_contrib
        return scores

    def stored_data_mean(self, compound, normalized):
        means = []
        for date_label, (df_orig, df_norm, col_id) in self.stored_data.items():
            df_to_use = df_norm if normalized else df_orig
            if compound in df_to_use.index:
                means.append(df_to_use.loc[compound].mean())
        return np.mean(means) if means else 1  # Return 1 if no means found to avoid division by zero


    def calculate_variability_impact(self, data):
        # Calculate the standard deviation as a measure of variability
        return data.std() / data.mean() if data.mean() != 0 else data.std()


    def detect_outliers(self, data):
        # Calculate the number of outliers based on IQR
        q75, q25 = np.percentile(data, [75, 25])
        iqr = q75 - q25
        lower_bound = q25 - 1.5 * iqr
        upper_bound = q75 + 1.5 * iqr
        outliers = [x for x in data if x < lower_bound or x > upper_bound]
        outlier_impact = len(outliers) / len(data)  # Proportion of outliers
        return outlier_impact


    def normalize_scores(self, scores):
        min_score = min(scores.values())
        max_score = max(scores.values())
        range_score = max_score - min_score
        if range_score == 0:
            return {compound: 0.5 for compound in scores}
        return {compound: (score - min_score) / range_score for compound, score in scores.items()}


    def update_table(self, sort_by=None):
        scores = self.calculate_scores()  # Calculate scores if scoring is applied

        data = []
        columns = [
            'Compound', 'Class', 'Ion', 'Mean', 'min-max', 'RSD', 'IQR',
            'N Mean', 'N min-max', 'N RSD', 'N IQR', 'Score', 'N Score'
        ]

        if sort_by == 'class':
            sorted_compounds = sorted(self.df.index, key=lambda x: self.class_mapping.get(x.replace('|', ';').replace('"', '').strip(), ''))
        elif sort_by == 'ion':
            sorted_compounds = sorted(self.df.index, key=lambda x: self.ion_mapping.get(x.replace('|', ';').replace('"', '').strip(), ''))
        else:
            sorted_compounds = self.df.index

        for compound in sorted_compounds:
            original_data = self.df.loc[compound]
            normalized_data = self.df_normalized.loc[compound]

            original_rsd = self.calculate_rsd(original_data)
            normalized_rsd = self.calculate_rsd(normalized_data)

            original_iqr, original_min_max = self.calculate_iqr_and_range(original_data)
            normalized_iqr, normalized_min_max = self.calculate_iqr_and_range(normalized_data)

            score = scores[compound][0] if scores and compound in scores else "N/A"
            n_score = scores[compound][1] if scores and compound in scores else "N/A"

            try:
                formatted_score = f"{float(score):.3f}" if score != "N/A" else "N/A"
            except ValueError:
                formatted_score = "N/A"

            try:
                formatted_n_score = f"{float(n_score):.3f}" if n_score != "N/A" else "N/A"
            except ValueError:
                formatted_n_score = "N/A"

            values = [
                compound.replace('"', ''),
                self.class_mapping.get(compound, 'N/A'),
                self.ion_mapping.get(compound, 'N/A'),
                f"{original_data.mean():.3e}",
                original_min_max,
                f"{original_rsd:.2f}",
                f"{original_iqr:.3e}",
                f"{normalized_data.mean():.3e}",
                normalized_min_max,
                f"{normalized_rsd:.2f}",
                f"{normalized_iqr:.3e}",
                formatted_score,
                formatted_n_score
            ]
            data.append(values)

        df_display = pd.DataFrame(data, columns=columns)

        if self.table is None:
            self.table = Table(self.table_frame, dataframe=df_display, showtoolbar=True, showstatusbar=True)
            self.table.show()
        else:
            self.table.updateModel(TableModel(df_display))
            self.table.redraw()

        # Set column widths
        self.set_column_widths(self.table, df_display)

        # Apply color to the 'Score' and 'N Score' columns
        self.apply_color_to_scores()

        # Update summary
        self.update_summary(df_display)


    def apply_color_to_scores(self):
        # Define the color function for scores
        def color_score(value):
            if value == "N/A":
                return None
            value = float(value)
            if value > 0.5:
                return "lightgreen"
            elif value < 0:
                return "lightcoral"
            else:
                return "yellow"

        # Create a mask for the dataframe
        mask_score = pd.Series(index=self.table.model.df.index, dtype=object)
        mask_n_score = pd.Series(index=self.table.model.df.index, dtype=object)

        for row in range(len(self.table.model.df)):
            score_value = self.table.model.df.iloc[row]['Score']
            n_score_value = self.table.model.df.iloc[row]['N Score']
            
            mask_score.iloc[row] = color_score(score_value)
            mask_n_score.iloc[row] = color_score(n_score_value)

        # Apply the color mask to the table
        self.table.setColorByMask(col='Score', mask=mask_score.notnull(), clr=mask_score)
        self.table.setColorByMask(col='N Score', mask=mask_n_score.notnull(), clr=mask_n_score)
        self.table.redraw()


    def update_summary(self, df_display):
        # Clear existing widgets in summary frame
        for widget in self.summary_frame.winfo_children():
            widget.destroy()

        summary_text = tk.Text(self.summary_frame, height=10, wrap='word', state='normal')
        summary_text.pack(fill='both', expand=True)

        class_summary = self.calculate_summary_by_color(df_display, 'Class')
        ion_summary = self.calculate_summary_by_color(df_display, 'Ion')

        summary_text.insert(tk.END, "Class Summary:\n")
        summary_text.insert(tk.END, class_summary.to_string())
        summary_text.insert(tk.END, "\n\nIon Summary:\n")
        summary_text.insert(tk.END, ion_summary.to_string())

        summary_text.config(state='disabled')


    def calculate_summary_by_color(self, df_display, group_by_column):
        def safe_apply_color_tag(value):
            try:
                return self.get_color_tag(float(value))
            except ValueError:
                return None

        summary = df_display.groupby(group_by_column)[['Score', 'N Score']].apply(
            lambda x: pd.Series({
                'Green Original': (x['Score'].apply(safe_apply_color_tag) == 'lightgreen').sum(),
                'Yellow Original': (x['Score'].apply(safe_apply_color_tag) == 'yellow').sum(),
                'Red Original': (x['Score'].apply(safe_apply_color_tag) == 'lightcoral').sum(),
                'Green Normalized': (x['N Score'].apply(safe_apply_color_tag) == 'lightgreen').sum(),
                'Yellow Normalized': (x['N Score'].apply(safe_apply_color_tag) == 'yellow').sum(),
                'Red Normalized': (x['N Score'].apply(safe_apply_color_tag) == 'lightcoral').sum()
            })
        )
        return summary


    def on_double_click(self, event):
        row_clicked = self.table.get_row_clicked(event)
        col_clicked = self.table.get_col_clicked(event)
        column_name = self.table.model.df.columns[col_clicked]
        
        compound = self.table.model.df.iloc[row_clicked]['Compound']
        
        if column_name in ["Mean", "min-max"]:  
            data = self.df.loc[compound]
            title_suffix = "Original Replicates"
            self.plot_data_points_scatter(data, compound, title_suffix)

        elif column_name in ["N Mean", "N min-max"]:
            data = self.df_normalized.loc[compound]
            title_suffix = "Normalized Replicates"
            self.plot_data_points_scatter(data, compound, title_suffix)

        elif column_name in ["RSD", "IQR"]:
            data = self.df.loc[compound]
            title_suffix = "Original"
            self.plot_data_points(data, compound, title_suffix)

        elif column_name in ["N RSD", "N IQR"]:
            data = self.df_normalized.loc[compound]
            title_suffix = "Normalized"
            self.plot_data_points(data, compound, title_suffix)

        elif column_name in ["Score", "N Score"]:
            if column_name == "Score":
                title_suffix = "Original with Stored Data"
                normalized = False
            else:
                title_suffix = "Normalized with Stored Data"
                normalized = True
            self.plot_box_w_stored(compound, title_suffix, normalized)


    def plot_box_w_stored(self, compound_name, title_suffix, normalized):
        if not self.stored_data:
            messagebox.showerror("Error", "Stored data is not available.")
            return

        data_lists = []
        labels = []
        col_id_list = []

        # Prepare data from stored files
        stored_data_entries = []
        for date_label, (df_orig, df_norm, col_id) in self.stored_data.items():
            df_to_use = df_norm if normalized else df_orig
            if compound_name in df_to_use.index:
                data = df_to_use.loc[compound_name].values.flatten()
                stored_data_entries.append((date_label, data, col_id))

        # Sort data by date label to ensure consistent order
        stored_data_entries.sort(key=lambda x: x[0])

        # Append sorted data
        for date_label, data, col_id in stored_data_entries:
            data_lists.append(data)
            labels.append(f"{date_label} ({col_id})")
            col_id_list.append(col_id)

        # Add the new data last
        if compound_name in (self.df_normalized if normalized else self.df).index:
            current_data = (self.df_normalized if normalized else self.df).loc[compound_name].values.flatten()
            data_lists.append(current_data)
            labels.append('New')
            col_id_list.append('New')  # This label will appear last

        # Ensure all entries in data_lists are 1D arrays
        data_lists = [np.array(data).flatten() for data in data_lists]

        # Generate plot if there is any data to show
        if data_lists:
            popup = tk.Toplevel()
            popup.title(f"{title_suffix} - {compound_name}")
            popup.geometry("1300x1000")

            fig, ax = plt.subplots()

            # Generate a color palette
            unique_col_ids = list(set(col_id_list))
            color_palette = plt.get_cmap("tab20")(np.linspace(0, 1, len(unique_col_ids)))
            col_id_to_color = {col_id: color_palette[i] for i, col_id in enumerate(unique_col_ids)}

            boxprops = dict(linestyle='-', linewidth=2)
            medianprops = dict(linestyle='-', linewidth=1, color='firebrick')

            for i, (data, col_id) in enumerate(zip(data_lists, col_id_list)):
                box = ax.boxplot(data, positions=[i], widths=0.6, patch_artist=True, 
                                boxprops=boxprops, medianprops=medianprops)
                for patch in box['boxes']:
                    patch.set_facecolor(col_id_to_color[col_id])

            ax.set_xticklabels(labels, rotation=45)
            ax.set_title(f"{title_suffix} - {compound_name}")
            ax.set_ylabel('Normalized Peak Area Values' if normalized else 'Original Peak Area Values')
            ax.set_xlabel('Date (Column ID)')

            plt.tight_layout()
            canvas = FigureCanvasTkAgg(fig, master=popup)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

            save_button = tk.Button(popup, text="Save Graph", command=lambda: save_figure_as_image(fig))
            save_button.pack(side=tk.LEFT)

            tk.Button(popup, text="Close", command=popup.destroy).pack(side=tk.RIGHT)
        else:
            messagebox.showinfo("Data Unavailable", f"No data available for {compound_name}.")


    def plot_data_points_scatter(self, data, compound_name, title_suffix):
        # Create a popup window for the scatter plot
        popup = tk.Toplevel()
        popup.title(f"{compound_name} - {title_suffix}")
        popup.geometry("1000x1000")

        fig, ax = plt.subplots()
        # Plot data points
        ax.scatter(range(len(data)), data, color='blue', alpha=0.7, label=f'{compound_name} data')

        # Customizing the plot
        ax.set_title(f"{compound_name} - {title_suffix}")
        ax.set_ylabel("Values")
        ax.set_xlabel("Sample Index")
        plt.xticks(rotation=45)

        # Creating a canvas as a matplotlib backend
        canvas = FigureCanvasTkAgg(fig, master=popup)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        save_button = tk.Button(popup, text="Save as SVG", command=lambda: save_figure_as_image(fig))
        save_button.pack(side=tk.LEFT)

        # Add a close button to the popup
        tk.Button(popup, text="Close", command=popup.destroy).pack(side=tk.RIGHT)


    def plot_data_points(self, data, compound_name, title_suffix):
        # Prepare data by grouping by the initial part of the column name
        groups = {}
        for col in data.dropna().index:
            date = col.split('-')[0]  # Assuming date is the first part before '-HEK-std'
            if date not in groups:
                groups[date] = []
            groups[date].append(data[col])

        # Setup the popup window
        popup = tk.Toplevel()
        popup.title(f"{compound_name} - {title_suffix}")
        popup.geometry("1000x1000")

        fig, ax = plt.subplots()
        # Create boxplot for each group of data
        box_data = [groups[date] for date in sorted(groups)]
        bp = ax.boxplot(box_data, tick_labels=sorted(groups.keys()), notch=True, vert=True, patch_artist=True, showfliers=True)

        # Customize the boxplot appearance
        for box in bp['boxes']:
            # Set edge color and fill with a more transparent color
            box.set(color='#1f77b4', linewidth=2)
            box.set(facecolor='#1f77b4', alpha=0.5)  # Set transparency

        for whisker in bp['whiskers']:
            whisker.set(color='#1f77b4', linewidth=2)

        for cap in bp['caps']:
            cap.set(color='#1f77b4', linewidth=2)

        for median in bp['medians']:
            median.set(color='yellow', linewidth=2)  # Set medians to yellow for visibility

        for flier in bp['fliers']:
            flier.set(marker='o', color='#e7298a', alpha=0.9)  # Outliers visible as pink dots

        # Add individual data points on the plot for clarity
        for i, line in enumerate(groups):
            y_data = groups[line]
            x_data = np.random.normal(1 + i, 0.02, size=len(y_data))  # Add some jitter to the x-axis
            ax.plot(x_data, y_data, 'r.', alpha=0.5)  # Points are plotted as red dots with transparency

        ax.set_title(f"{compound_name} - {title_suffix}")
        ax.set_ylabel("Peak Area")
        ax.set_xlabel("Date Run")
        plt.xticks(rotation=45)
        plt.grid(True)

        # Create a canvas as a matplotlib backend
        canvas = FigureCanvasTkAgg(fig, master=popup)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        save_button = tk.Button(popup, text="Save as SVG", command=lambda: save_figure_as_image(fig))
        save_button.pack(side=tk.LEFT)

        # Add a close button to the popup
        tk.Button(popup, text="Close", command=popup.destroy).pack(side=tk.RIGHT)


    

if __name__ == "__main__":
    app = MetaboliteAnalysisApp()
    app.mainloop()


# %%



