import tkinter as tk
from tkinter import ttk

class CollapsiblePane(ttk.Frame):
    """
    A collapsible pane widget that can be embedded into a GUI layout.
    """
    def __init__(self, parent, title="title", *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        # Toggle button
        self.title = title
        self.toggle_button = ttk.Button(self, text=self.title, command=self.toggle)
        self.toggle_button.pack(fill="x")

        # Container for the collapsible content
        self.container = ttk.Frame(self)
        self.container.pack(fill="x", expand=True)

        # Start collapsed
        self.is_collapsed = True
        self.container.pack_forget()

    def toggle(self):
        """Toggle the visibility of the container."""
        if self.is_collapsed:
            self.container.pack(fill="x", expand=True)
            self.toggle_button.config(text="Hide")
        else:
            self.container.pack_forget()
            self.toggle_button.config(text=self.title)
        self.is_collapsed = not self.is_collapsed

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Collapsible Pane Embedded Example")

    # Example usage
    collapsible_pane = CollapsiblePane(root, title="Example Collapsible Pane")
    collapsible_pane.pack(fill="x", padx=10, pady=10)

    # Adding content to the collapsible pane
    ttk.Label(collapsible_pane.container, text="This is a collapsible pane.").pack(pady=10)
    ttk.Button(collapsible_pane.container, text="Button inside pane").pack(pady=10)

    root.mainloop()
