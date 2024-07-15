import os
import pandas as pd
from tkinter import Tk, Label, Button, Listbox, MULTIPLE, END, filedialog, Entry
from tkinter import ttk

class CSVColumnSelector:
    def __init__(self, root):
        self.root = root
        self.root.title("CSV Column Selector")

        self.label = Label(root, text="Select a CSV file:")
        self.label.pack()

        self.file_button = Button(root, text="Browse", command=self.browse_file)
        self.file_button.pack()

        self.refresh_button = Button(root, text="Refresh", command=self.refresh_file)
        self.refresh_button.pack()

        self.column_listbox = Listbox(root, selectmode=MULTIPLE)
        self.column_listbox.pack(fill='both', expand=True)

        self.save_selection_button = Button(root, text="Save Selection", command=self.save_selection)
        self.save_selection_button.pack()

        self.load_selection_button = Button(root, text="Load Selection", command=self.load_selection)
        self.load_selection_button.pack()

        self.progress = ttk.Progressbar(root, orient="horizontal", length=300, mode="determinate")
        self.progress.pack()

        self.selection_label = Label(root, text="Filename to Save/Load Selection:")
        self.selection_label.pack()

        self.selection_entry = Entry(root)
        self.selection_entry.pack()

        self.csv_file = None
        self.df = None

    def browse_file(self):
        file_selected = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if file_selected:
            self.csv_file = file_selected
            self.update_column_list()

    def refresh_file(self):
        if self.csv_file:
            self.update_column_list()

    def update_column_list(self):
        if self.csv_file:
            self.df = pd.read_csv(self.csv_file)
            self.column_listbox.delete(0, END)
            for column in self.df.columns:
                self.column_listbox.insert(END, column)

    def save_selection(self):
        selected_indices = self.column_listbox.curselection()
        selected_columns = [self.df.columns[i] for i in selected_indices]

        if not selected_columns:
            print("No columns selected!")
            return

        selection_file = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if selection_file:
            with open(selection_file, 'w') as f:
                for column in selected_columns:
                    f.write(f"{column}\n")
            print(f"Selection saved to {selection_file}")

    def load_selection(self):
        selection_file = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if not selection_file or not os.path.exists(selection_file):
            print("Invalid filename or file does not exist!")
            return

        with open(selection_file, 'r') as f:
            selected_columns = [line.strip() for line in f.readlines()]

        self.column_listbox.selection_clear(0, END)
        for i, column in enumerate(self.df.columns):
            if column in selected_columns:
                self.column_listbox.selection_set(i)
        print(f"Selection loaded from {selection_file}")

if __name__ == "__main__":
    root = Tk()
    csv_column_selector = CSVColumnSelector(root)
    root.mainloop()
