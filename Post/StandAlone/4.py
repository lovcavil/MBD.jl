import os
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import Tk, Label, Button, Listbox, MULTIPLE, END, filedialog, Checkbutton, IntVar, Entry
from tkinter import ttk
from multiprocessing import Pool, Manager, cpu_count, Value, Lock

# Default colors and line styles
default_colors = [
    "#606c38",  # Green
    "#2a9d8f",
    "#e9c46a",
    "#d62828",
    "#219ebc",
    '#104680',  # Dark blue
    '#317CB7',  # Medium blue
    '#6DADD1',  # Light blue
    '#B6D7E8',  # Pale blue
    '#F6B293',  # Pale orange
    '#DC6D57',  # Medium red
    '#B72230'   # Dark red
]
default_line_styles = ['-', '--', '-.', ':', (0, (1, 10)), (0, (5, 10)), (0, (3, 5, 1, 5)),
                       (0, (5, 1)), (0, (3, 10, 1, 10)), (0, (3, 1, 1, 1))]

# Create a global lock
lock = Lock()

def plot_and_save(params):
    csv_files, column, output_folder, colors, line_styles, progress_counter, total_figures, dpi = params
    plt.figure(figsize=(12, 8))
    
    for i, file in enumerate(csv_files):
        df = pd.read_csv(file)
        if column in df.columns:
            plt.plot(df['Time'], df[column], label=os.path.basename(file),
                     color=colors[i % len(colors)], linestyle=line_styles[i % len(line_styles)])

    plt.title(f'{column} vs Time')
    plt.xlabel('Time')
    plt.ylabel(column)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'{output_folder}/{column}_vs_Time.png', dpi=dpi)
    plt.close()

    with lock:
        progress_counter.value += 1
        print(f"Progress: {progress_counter.value}/{total_figures}")

class CSVPlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("CSV Plotter")

        self.notebook = ttk.Notebook(root)
        self.tab1 = ttk.Frame(self.notebook)
        self.tab2 = ttk.Frame(self.notebook)
        self.tab3 = ttk.Frame(self.notebook)
        self.notebook.add(self.tab1, text='Main')
        self.notebook.add(self.tab2, text='Settings')
        self.notebook.add(self.tab3, text='Resolution')
        self.notebook.pack(expand=True, fill='both')

        # Tab 1: Main
        self.label = Label(self.tab1, text="Select a folder to list CSV files:")
        self.label.pack()

        self.folder_button = Button(self.tab1, text="Browse", command=self.browse_folder)
        self.folder_button.pack()

        self.refresh_button = Button(self.tab1, text="Refresh", command=self.refresh_folder)
        self.refresh_button.pack()

        self.csv_listbox = Listbox(self.tab1, selectmode=MULTIPLE)
        self.csv_listbox.pack(fill='both', expand=True)

        self.plot_button = Button(self.tab1, text="Plot Selected CSVs", command=self.plot_selected_csvs)
        self.plot_button.pack()

        self.plot_button2 = Button(self.tab1, text="Plot Selected CSVs2", command=self.plot_selected_csvs2)
        self.plot_button2.pack()

        self.subfolder_label = Label(self.tab1, text="Subfolder Name:")
        self.subfolder_label.pack()

        self.subfolder_entry = Entry(self.tab1)
        self.subfolder_entry.pack()

        self.column_file_label = Label(self.tab1, text="CSV File with Column Names (in selected folder):")
        self.column_file_label.pack()

        self.column_file_entry = Entry(self.tab1)
        self.column_file_entry.insert(0, "col.csv")  # Default value for the CSV column file
        self.column_file_entry.pack()

        self.progress = ttk.Progressbar(self.tab1, orient="horizontal", length=300, mode="determinate")
        self.progress.pack()

        self.csv_files = []
        self.current_folder = None

        self.output_folder = './plots'
        os.makedirs(self.output_folder, exist_ok=True)

        # Tab 2: Settings
        self.color_vars = [IntVar(value=1) for _ in default_colors]
        self.style_vars = [IntVar(value=1) for _ in default_line_styles]

        Label(self.tab2, text="Select Colors:").pack()
        for i, color in enumerate(default_colors):
            Checkbutton(self.tab2, text=color, variable=self.color_vars[i]).pack(anchor='w')

        Label(self.tab2, text="Select Line Styles:").pack()
        for i, style in enumerate(default_line_styles):
            Checkbutton(self.tab2, text=str(style), variable=self.style_vars[i]).pack(anchor='w')

        # Tab 3: Resolution
        Label(self.tab3, text="Set Resolution (DPI):").pack()
        self.dpi_entry = Entry(self.tab3)
        self.dpi_entry.insert(0, "300")  # Default DPI value
        self.dpi_entry.pack()

    def browse_folder(self):
        folder_selected = filedialog.askdirectory()
        if folder_selected:
            self.current_folder = folder_selected
            self.update_file_list()

    def refresh_folder(self):
        if self.current_folder:
            self.update_file_list()

    def update_file_list(self):
        if self.current_folder:
            self.csv_files = [os.path.join(self.current_folder, f) for f in os.listdir(self.current_folder) if f.endswith('.csv')]
            self.csv_listbox.delete(0, END)
            for csv_file in self.csv_files:
                self.csv_listbox.insert(END, os.path.basename(csv_file))

    def plot_selected_csvs(self):
        selected_indices = self.csv_listbox.curselection()
        selected_files = [self.csv_files[i] for i in selected_indices]

        if not selected_files:
            print("No files selected!")
            return

        subfolder = self.subfolder_entry.get()
        output_folder = os.path.join(self.output_folder, subfolder)
        os.makedirs(output_folder, exist_ok=True)

        df_example = pd.read_csv(selected_files[0])
        columns_to_plot = [col for col in df_example.columns if col != 'Time']

        selected_colors = [default_colors[i] for i, var in enumerate(self.color_vars) if var.get() == 1]
        selected_styles = [default_line_styles[i] for i, var in enumerate(self.style_vars) if var.get() == 1]

        total_figures = len(columns_to_plot)
        self.progress["maximum"] = total_figures

        manager = Manager()
        self.progress_counter = manager.Value('i', 0)

        try:
            dpi = int(self.dpi_entry.get())
        except ValueError:
            dpi = 300  # Default DPI if invalid input

        params_list = [(selected_files, column, output_folder, selected_colors, selected_styles, self.progress_counter, total_figures, dpi) for column in columns_to_plot]

        self.pool = Pool(processes=cpu_count())
        self.pool.map_async(plot_and_save, params_list, callback=self.on_done)

        self.check_progress()

    def plot_selected_csvs2(self):
        selected_indices = self.csv_listbox.curselection()
        selected_files = [self.csv_files[i] for i in selected_indices]

        if not selected_files:
            print("No files selected!")
            return

        subfolder = self.subfolder_entry.get()
        output_folder = os.path.join(self.output_folder, subfolder)
        os.makedirs(output_folder, exist_ok=True)

        column_file_name = self.column_file_entry.get()
        column_file_path = os.path.join(self.current_folder, column_file_name)
        if not os.path.isfile(column_file_path):
            print("Invalid column file path!")
            return

        df_columns = pd.read_csv(column_file_path)
        columns_to_plot = df_columns.iloc[:, 0].tolist()

        selected_colors = [default_colors[i] for i, var in enumerate(self.color_vars) if var.get() == 1]
        selected_styles = [default_line_styles[i] for i, var in enumerate(self.style_vars) if var.get() == 1]

        total_figures = len(columns_to_plot)
        self.progress["maximum"] = total_figures

        manager = Manager()
        self.progress_counter = manager.Value('i', 0)

        try:
            dpi = int(self.dpi_entry.get())
        except ValueError:
            dpi = 300  # Default DPI if invalid input

        params_list = [(selected_files, column, output_folder, selected_colors, selected_styles, self.progress_counter, total_figures, dpi) for column in columns_to_plot]

        self.pool = Pool(processes=cpu_count())
        self.pool.map_async(plot_and_save, params_list, callback=self.on_done)

        self.check_progress()

    def check_progress(self):
        self.progress["value"] = self.progress_counter.value
        if self.progress_counter.value < self.progress["maximum"]:
            self.root.after(100, self.check_progress)
        else:
            print(f"Plots saved in the folder: {self.output_folder}")

    def on_done(self, _):
        self.pool.close()
        self.pool.join()

if __name__ == "__main__":
    root = Tk()
    csv_plotter = CSVPlotter(root)
    root.mainloop()
