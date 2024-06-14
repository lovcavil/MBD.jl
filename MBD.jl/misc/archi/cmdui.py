import tkinter as tk
from tkinter import filedialog, messagebox
import subprocess
import os

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
        # new_file_path = filedialog.askopenfilename(filetypes=[("Text files", "*.txt")])
        # if new_file_path:
        #     self.file_path = new_file_path
        # else:
        #     self.file_path = self.default_file_path
        
        self.load_commands()

    def create_command_entry(self, command):
        command_frame = tk.Frame(self.frame)
        command_frame.pack(pady=2, fill=tk.X)

        text_box = tk.Entry(command_frame, width=80)
        text_box.pack(side=tk.LEFT, padx=5)
        text_box.insert(0, command)
        self.text_boxes.append(text_box)

        run_button = tk.Button(command_frame, text="Run", command=lambda tb=text_box: self.run_command(tb))
        run_button.pack(side=tk.RIGHT, padx=5)

    def run_command(self, text_box):
        command = text_box.get()
        try:
            # Open a new command window and run the command without waiting for it to complete
            subprocess.Popen(f'start cmd /c "{command}"', shell=True)
            messagebox.showinfo("Command Started", f"Command: {command}\n\nThe command has been started in a new window.")
        except Exception as e:
            messagebox.showerror("Command Error", f"Command: {command}\n\nError:\n{str(e)}")

    def save_commands(self):
        if self.file_path:
            with open(self.file_path, 'w') as file:
                for text_box in self.text_boxes:
                    file.write(text_box.get() + '\n')
            messagebox.showinfo("Success", "Commands saved successfully!")
        else:
            messagebox.showerror("Error", "No file loaded to save commands!")

if __name__ == "__main__":
    root = tk.Tk()
    app = CommandApp(root)
    root.mainloop()
