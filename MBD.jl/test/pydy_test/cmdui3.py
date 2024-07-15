import tkinter as tk
from tkinter import messagebox
import subprocess
import os
import tempfile

class CommandApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Command Line Executor")

        self.commands = []
        self.text_boxes = []
        absolute_script_path = os.path.abspath(__file__)
        script_directory = os.path.dirname(absolute_script_path)
        self.default_file_path = os.path.join(script_directory, "commands.txt")
        self.file_path = self.default_file_path

        self.update_button = tk.Button(root, text="Update Commands", command=self.update_commands)
        self.update_button.pack()

        self.save_button = tk.Button(root, text="Save Commands", command=self.save_commands)
        self.save_button.pack()

        self.frame = tk.Frame(root)
        self.frame.pack(fill=tk.BOTH, expand=True)

        # Load default commands
        self.load_commands()

    def load_commands(self):
        if os.path.isfile(self.file_path):
            with open(self.file_path, 'r') as file:
                self.commands = file.readlines()

            self.text_boxes = []
            for widget in self.frame.winfo_children():
                widget.destroy()

            for command in self.commands:
                self.create_command_entry(command.strip())

    def update_commands(self):
        self.load_commands()

    def create_command_entry(self, command):
        command_frame = tk.Frame(self.frame)
        command_frame.pack(pady=2, fill=tk.X)

        text_box = tk.Text(command_frame, height=5, width=100)
        text_box.pack(side=tk.LEFT, padx=5)
        text_box.insert(tk.END, command)
        self.text_boxes.append(text_box)

        run_button = tk.Button(command_frame, text="Run", command=lambda tb=text_box: self.run_command(tb))
        run_button.pack(side=tk.RIGHT, padx=5)

    def run_command(self, text_box):
        command = text_box.get("1.0", tk.END).strip()
        try:
            # Create a temporary batch file to run the full command sequence
            with tempfile.NamedTemporaryFile(delete=False, suffix='.bat', mode='w') as batch_file:
                batch_file.write(command)
                batch_file_path = batch_file.name

            # Run the batch file in a new Windows Terminal tab
            subprocess.Popen(f'wt -w 0 nt -- cmd /c "{batch_file_path}"', shell=True)
        except Exception as e:
            messagebox.showerror("Command Error", f"Command: {command}\n\nError:\n{str(e)}")

    def save_commands(self):
        if self.file_path:
            with open(self.file_path, 'w') as file:
                for text_box in self.text_boxes:
                    file.write(text_box.get("1.0", tk.END).strip() + '\n')
            messagebox.showinfo("Success", "Commands saved successfully!")
        else:
            messagebox.showerror("Error", "No file loaded to save commands!")

if __name__ == "__main__":
    root = tk.Tk()
    app = CommandApp(root)
    root.mainloop()
