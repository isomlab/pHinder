import tkinter as tk
from tkinter import filedialog, ttk
import os

class FilePathWidget:
    def __init__(self, root, process_file=None, options_label="Options:", extra_options={}):
        self.root = root
        self.file_path = tk.StringVar()
        self.save_path = tk.StringVar()
        self.process_file = process_file if process_file else self.default_process_file
        self.options_label = options_label
        self.extra_options = extra_options

        self.frame = ttk.Frame(self.root)
        self.frame.pack(pady=10, padx=10, fill="x")

        self.hidden_row_frame = None
        self.option_vars = {}  # Changed to a dictionary
        self.create_widgets()

    def create_widgets(self):
        # File Path Row
        ttk.Label(self.frame, text="File Path:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        file_path_entry = ttk.Entry(self.frame, textvariable=self.file_path, width=50)
        file_path_entry.grid(row=0, column=1, padx=5, pady=5)

        browse_file_button = ttk.Button(self.frame, text="Browse", command=self.browse_and_process_file)
        browse_file_button.grid(row=0, column=2, padx=5, pady=5)

        # Save Path Row
        ttk.Label(self.frame, text="Save Path:").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        save_path_entry = ttk.Entry(self.frame, textvariable=self.save_path, width=50)
        save_path_entry.grid(row=1, column=1, padx=5, pady=5)

        browse_save_button = ttk.Button(self.frame, text="Browse", command=self.browse_save_path)
        browse_save_button.grid(row=1, column=2, padx=5, pady=5)

        # Hidden Row for File Specific Options
        self.hidden_row_frame = ttk.Frame(self.frame)
        self.hidden_row_frame.grid(row=2, column=0, columnspan=3, sticky="w", padx=5, pady=5)
        self.hidden_row_frame.grid_remove()

    def browse_and_process_file(self):
        file_path = filedialog.askopenfilename(title="Select a File")
        if file_path:
            self.file_path.set(file_path)
            self.save_path.set(os.path.dirname(file_path))  # Set default save path to the directory of the file
            options = self.process_file(file_path)
            if options:
                self.show_file_specific_options(self.options_label, options)

    def browse_save_path(self):
        initial_dir = self.save_path.get() if self.save_path.get() else os.getcwd()
        save_path = filedialog.askdirectory(title="Select Save Directory", initialdir=initial_dir)
        if save_path:
            self.save_path.set(save_path)

    def show_file_specific_options(self, label, options):
        # Clear previous options
        for widget in self.hidden_row_frame.winfo_children():
            widget.destroy()

        self.option_vars = {}  # Clear dictionary

        # Add label to the left of options
        ttk.Label(self.hidden_row_frame, text=label).grid(row=0, column=0, padx=5, pady=5, sticky="w")

        # Populate new options with checkboxes (selected by default), max 10 per row
        row = 0
        col = 1
        for option in options:
            var = tk.IntVar(value=1)
            self.option_vars[option] = var  # Store in dictionary
            ttk.Checkbutton(self.hidden_row_frame, text=option, variable=var).grid(row=row, column=col, padx=5, pady=2, sticky="w")
            col += 1
            if col >= 11:  # Account for label in column 0
                col = 1
                row += 1

        if self.extra_options:
            for option in self.extra_options:
                var = tk.IntVar(value=self.extra_options[option])
                self.option_vars[option] = var  # Store in dictionary
                ttk.Checkbutton(self.hidden_row_frame, text=option, variable=var).grid(row=row, column=col, padx=5, pady=2, sticky="w")
                col += 1
                if col >= 11:  # Account for label in column 0
                    col = 1
                    row += 1

        # Make the hidden row visible
        self.hidden_row_frame.grid()

    def get_values(self):
        if not self.option_vars:
            return None
        return {option: var.get() for option, var in self.option_vars.items()}

    def get_file_path(self):
        return self.file_path.get()

    def get_save_path(self):
        return self.save_path.get()

    def default_process_file(self, file_path):
        # Default placeholder for processing file
        print(f"Default processing file: {file_path}")
        # Example: Return file-specific options
        return ["Option 1: Example", "Option 2: Placeholder"]

    def report_values(self):
        values = self.get_values()
        if values:
            print("Selected Options:", values)
        else:
            print("No options available.")

if __name__ == "__main__":
    def custom_process_function(file_path):
        print(f"Custom processing file: {file_path}")
        # Example: Return file-specific options based on file
        return ["Option 1: Custom", "Option 2: File Details", "Option 3: More Info"]

    root = tk.Tk()
    root.title("File Path Widget")

    app = FilePathWidget(root, custom_process_function)

    ttk.Button(root, text="Report Values", command=app.report_values).pack(pady=10)

    root.mainloop()
