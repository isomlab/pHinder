import tkinter as tk
from tkinter import ttk

class IntegerInputWidget:
    def __init__(self, root):
        self.root = root

        self.frame = ttk.Frame(self.root)
        self.frame.pack(pady=5, padx=10, fill="x")

        self.input_var = tk.StringVar()
        self.has_typed = False

        self.create_widgets()

    def create_widgets(self):
        # Input Label and Field in the Same Row
        ttk.Label(self.frame, text="Enter an Integer:").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.entry = ttk.Entry(self.frame, textvariable=self.input_var, width=30)
        self.entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        self.entry.bind("<KeyRelease>", self.validate_input)
        self.entry.bind("<Leave>", self.on_mouse_leave)

        # Error Label
        self.error_label = ttk.Label(self.frame, text="", foreground="red", font=("Arial", 10))
        self.error_label.grid(row=1, column=1, sticky="w", padx=5, pady=0)

    def on_mouse_leave(self, event):
        """Validate when mouse leaves the field."""
        if self.input_var.get().strip():  # Only validate if input exists
            self.validate_input()

    def validate_input(self, event=None):
        """Validate the integer input."""
        value = self.input_var.get()
        try:
            int_value = int(value)
            self.error_label.config(text="")  # Clear error message
            print(f"Valid integer: {int_value}")
        except ValueError:
            self.error_label.config(text="Invalid input. Please enter a valid integer.")

    def get_values(self):
        """Retrieve the validated integer value or indicate invalid state."""
        value = self.input_var.get()
        try:
            return int(value)
        except ValueError:
            return "Invalid input"

class FloatInputWidget:
    def __init__(self, root):
        self.root = root

        self.frame = ttk.Frame(self.root)
        self.frame.pack(pady=5, padx=10, fill="x")

        self.input_var = tk.StringVar()

        self.create_widgets()

    def create_widgets(self):
        # Input Label and Field in the Same Row
        ttk.Label(self.frame, text="Enter a Float:").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.entry = ttk.Entry(self.frame, textvariable=self.input_var, width=30)
        self.entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        self.entry.bind("<KeyRelease>", self.validate_input)
        self.entry.bind("<Leave>", self.on_mouse_leave)

        # Error Label
        self.error_label = ttk.Label(self.frame, text="", foreground="red", font=("Arial", 10))
        self.error_label.grid(row=1, column=1, sticky="w", padx=5, pady=0)

    def on_mouse_leave(self, event):
        """Validate when mouse leaves the field."""
        if self.input_var.get().strip():  # Only validate if input exists
            self.validate_input()

    def validate_input(self, event=None):
        """Validate the float input."""
        value = self.input_var.get()
        try:
            float_value = float(value)
            self.error_label.config(text="")  # Clear error message
            print(f"Valid float: {float_value}")
        except ValueError:
            self.error_label.config(text="Invalid input. Please enter a valid float.")

    def get_values(self):
        """Retrieve the validated float value or indicate invalid state."""
        value = self.input_var.get()
        try:
            return float(value)
        except ValueError:
            return "Invalid input"

class StringInputWidget:
    def __init__(self, root):
        self.root = root

        self.frame = ttk.Frame(self.root)
        self.frame.pack(pady=5, padx=10, fill="x")

        self.input_var = tk.StringVar()

        self.create_widgets()

    def create_widgets(self):
        # Input Label and Field in the Same Row
        ttk.Label(self.frame, text="Enter a String:").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.entry = ttk.Entry(self.frame, textvariable=self.input_var, width=30)
        self.entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        self.entry.bind("<KeyRelease>", self.validate_input)
        self.entry.bind("<Leave>", self.on_mouse_leave)

        # Error Label
        self.error_label = ttk.Label(self.frame, text="", foreground="red", font=("Arial", 10))
        self.error_label.grid(row=1, column=1, sticky="w", padx=5, pady=0)

    def on_mouse_leave(self, event):
        """Validate when mouse leaves the field."""
        if self.input_var.get().strip():  # Only validate if input exists
            self.validate_input()

    def validate_input(self, event=None):
        """Validate the string input."""
        value = self.input_var.get()
        if value.strip():
            self.error_label.config(text="")  # Clear error message
            print(f"Valid string: {value.strip()}")
        else:
            self.error_label.config(text="Invalid input. Please enter a valid string.")

    def get_values(self):
        """Retrieve the validated string value or indicate invalid state."""
        value = self.input_var.get()
        if value.strip():
            return value.strip()
        return "Invalid input"

class CheckButtonWidget:
    def __init__(self, root, label="Check Option"):
        self.root = root

        self.frame = ttk.Frame(self.root)
        self.frame.pack(pady=5, padx=10, fill="x")

        self.input_var = tk.BooleanVar()

        self.create_widgets(label)

    def create_widgets(self, label):
        # Checkbutton
        ttk.Label(self.frame, text=label).grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.checkbutton = ttk.Checkbutton(self.frame, variable=self.input_var)
        self.checkbutton.grid(row=0, column=1, sticky="w", padx=5, pady=2)

    def get_values(self):
        """Retrieve the value of the checkbutton."""
        return self.input_var.get()

class ComboBoxWidget:
    def __init__(self, root, label="Select an Option", options=None):
        self.root = root
        self.options = options if options else []

        self.frame = ttk.Frame(self.root)
        self.frame.pack(pady=5, padx=10, fill="x")

        self.input_var = tk.StringVar()

        self.create_widgets(label)

    def create_widgets(self, label):
        # ComboBox Label and Field
        ttk.Label(self.frame, text=label).grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.combobox = ttk.Combobox(self.frame, textvariable=self.input_var, values=self.options, state="readonly")
        self.combobox.grid(row=0, column=1, sticky="w", padx=5, pady=2)

    def get_values(self):
        """Retrieve the selected value from the combobox."""
        return self.input_var.get()

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Input Widgets")

    int_widget = IntegerInputWidget(root)
    float_widget = FloatInputWidget(root)
    string_widget = StringInputWidget(root)
    check_widget = CheckButtonWidget(root, label="Enable Feature")
    combo_widget = ComboBoxWidget(root, label="Choose an Option", options=["Option 1", "Option 2", "Option 3"])

    root.mainloop()
