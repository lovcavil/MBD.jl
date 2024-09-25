import os
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import Tk, Label, Button, Listbox, MULTIPLE, END, filedialog, Checkbutton, IntVar, Entry
from tkinter import ttk
from multiprocessing import Pool, Manager, cpu_count, Value, Lock

# Default colors and line styles
default_colors = [
    (16/256, 70/256, 128/256),   
    (49/256, 124/256, 183/256),    
    (109/256, 173/256, 209/256),     
    (182/256, 215/256, 232/256),       
    (246/256, 178/256, 147/256),
    (220/256, 109/256, 87/256),
    (183/256, 34/256, 48/256)
]
default_line_styles = ['-', '--', '-.', ':', (0, (1, 10)), (0, (5, 10)), (0, (3, 5, 1, 5)), 
                       (0, (5, 1)), (0, (3, 10, 1, 10)), (0, (3, 1, 1, 1))]

# Create a global lock
lock = Lock()

def plot_and_save(params):
    csv_files, column, output_folder, colors, line_styles, progress_counter, total_figures = params
    plt.figure()
    
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
    plt.savefig(f'{output_folder}/{column}_vs_Time.png')
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
        self.notebook.add(self.tab1, text='Main')
        self.notebook.add(self.tab2, text='Settings')
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

        self.subfolder_label = Label(self.tab1, text="Subfolder Name:")
        self.subfolder_label.pack()

        self.subfolder_entry = Entry(self.tab1)
        self.subfolder_entry.pack()

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
            Checkbutton(self.tab2, text=str(color), variable=self.color_vars[i]).pack(anchor='w')

        Label(self.tab2, text="Select Line Styles:").pack()
        for i, style in enumerate(default_line_styles):
            Checkbutton(self.tab2, text=str(style), variable=self.style_vars[i]).pack(anchor='w')

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
        
        params_list = [(selected_files, column, output_folder, selected_colors, selected_styles, self.progress_counter, total_figures) for column in columns_to_plot]

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
