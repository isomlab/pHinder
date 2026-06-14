import tkinter as tk
from tkinter import LabelFrame, Button, Frame
from tkinter import ttk
import json

class DynamicOptionWidget:
    def __init__(self, root, title, options, save_file="options_state.json"):
        self.root = root
        self.title = title
        self.options = options
        self.save_file = save_file
        self.vars = {option: tk.IntVar(value=1) for option in options}
        self.visible = False

        self.frame = Frame(self.root)
        self.frame.pack(pady=10, padx=10, fill='x')

        self.options_frame = LabelFrame(self.frame, text=self.title, padx=5, pady=5)

        Button(self.frame, text=f"Toggle {self.title}", command=self.toggle_visibility).pack(anchor='w')

        self.load_state()
        self.create_options()

    def create_options(self):
        for widget in self.options_frame.winfo_children():
            widget.destroy()

        for i, option in enumerate(self.options):
            ttk.Checkbutton(
                self.options_frame,
                text=option,
                variable=self.vars[option]
            ).grid(row=i // 4, column=i % 4, sticky='w', padx=5, pady=2)

    def toggle_visibility(self):
        if self.visible:
            self.options_frame.pack_forget()
        else:
            self.options_frame.pack(pady=5, fill='x')
        self.visible = not self.visible

    def save_state(self):
        try:
            state = {option: var.get() for option, var in self.vars.items()}
            with open(self.save_file, "w") as f:
                json.dump(state, f)
        except Exception as e:
            print(f"Error saving state: {e}")

    def load_state(self):
        try:
            with open(self.save_file, "r") as f:
                state = json.load(f)
                for option, value in state.items():
                    if option in self.vars:
                        self.vars[option].set(value)
        except FileNotFoundError:
            print("No saved state found. Using default values.")
        except Exception as e:
            print(f"Error loading state: {e}")

    def get_values(self):
        """Return the current state of options as a dictionary."""
        return {option: var.get() for option, var in self.vars.items()}

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Collapsible Options Example")

    # Apply a modern theme
    style = ttk.Style()
    style.theme_use("clam")

    options_list = ["Option 1", "Option 2", "Option 3", "Option 4"]
    widget = DynamicOptionWidget(root, "Example Options", options_list)

    Button(root, text="Save State", command=widget.save_state).pack(pady=10)
    Button(root, text="Print Values", command=lambda: print(widget.get_values())).pack(pady=10)
    Button(root, text="Exit", command=root.quit).pack(pady=10)

    root.mainloop()
