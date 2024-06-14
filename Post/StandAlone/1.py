import os
import pandas as pd
import matplotlib.pyplot as plt
from tkinter import Tk, Label, Button, Listbox, MULTIPLE, END, filedialog
from plot1 import *
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
            self.csv_files = [f for f in os.listdir(self.current_folder) if f.endswith('.csv')]
            self.csv_listbox.delete(0, END)
            for csv_file in self.csv_files:
                self.csv_listbox.insert(END, csv_file)

    def plot_selected_csvs(self):
        selected_indices = self.csv_listbox.curselection()
        selected_files = [self.csv_files[i] for i in selected_indices]

        for file in selected_files:
            file_path = os.path.join(self.current_folder, file)
            #print(file_path)
            df = pd.read_csv(file_path)
            # df.plot()
            # plt.title(file)
            # plt.show()

if __name__ == "__main__":
    root = Tk()
    csv_plotter = CSVPlotter(root)
    root.mainloop()
