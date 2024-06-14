import os
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import Tk, Label, Button, Listbox, MULTIPLE, END, filedialog
from multiprocessing import Pool

# Function to plot and save figures
def plot_and_save(params):
    csv_files, column, output_folder = params
    plt.figure()
    
    for file in csv_files:
        df = pd.read_csv(file)
        if column in df.columns:
            plt.plot(df['Time'], df[column], label=os.path.basename(file))
    
    plt.title(f'{column} vs Time')
    plt.xlabel('Time')
    plt.ylabel(column)
    plt.legend()
    plt.grid(True)
    plt.savefig(f'{output_folder}/{column}_vs_Time.png')
    plt.close()

class CSVPlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("CSV Plotter")

        self.label = Label(root, text="Select a folder to list CSV files:")
        self.label.pack()

        self.folder_button = Button(root, text="Browse", command=self.browse_folder)
        self.folder_button.pack()

        self.refresh_button = Button(root, text="Refresh", command=self.refresh_folder)
        self.refresh_button.pack()

        self.csv_listbox = Listbox(root, selectmode=MULTIPLE)
        self.csv_listbox.pack(fill='both', expand=True)

        self.plot_button = Button(root, text="Plot Selected CSVs", command=self.plot_selected_csvs)
        self.plot_button.pack()

        self.csv_files = []
        self.current_folder = None

        self.output_folder = './plots'
        os.makedirs(self.output_folder, exist_ok=True)

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

        df_example = pd.read_csv(selected_files[0])
        columns_to_plot = [col for col in df_example.columns if col != 'Time']

        params_list = [(selected_files, column, self.output_folder) for column in columns_to_plot]

        with Pool(processes=os.cpu_count()) as pool:
            pool.map(plot_and_save, params_list)

        print(f"Plots saved in the folder: {self.output_folder}")

if __name__ == "__main__":
    root = Tk()
    csv_plotter = CSVPlotter(root)
    root.mainloop()
