import tkinter as tk

class TerminalOutputWidget(tk.Frame):
    def __init__(self, parent, title="Terminal Output", *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        # Get the root window reference
        root = self.winfo_toplevel()

        # # Calculate default height as 80% of the screen height
        # screen_height = root.winfo_screenheight()
        # widget_height = int(screen_height * 0.8)

        # # Calculate default width as 50% of the screen width
        # screen_width = root.winfo_screenwidth()
        # widget_width = int(screen_width * 0.6)

        # # Set default size of the root window
        # root.geometry(f"{widget_width}x{widget_height}")

        # Create a label for the widget
        self.label = tk.Label(self, text=title, anchor="center")
        self.label.pack(pady=5)

        # Create the text widget
        self.text_widget = tk.Text(self, wrap="word", state="disabled", bg="black", fg="white")
        self.text_widget.pack(fill="both", expand=True, padx=5, pady=5)

        # Configure the scrollbar
        self.scrollbar = tk.Scrollbar(self, command=self.text_widget.yview)
        self.text_widget.configure(yscrollcommand=self.scrollbar.set)
        self.scrollbar.pack(side="right", fill="y")

    def write(self, message):
        """Write a message to the text widget."""
        self.text_widget.configure(state="normal")
        self.text_widget.insert("end", message)
        self.text_widget.see("end")
        self.text_widget.configure(state="disabled")

    def clear(self):
        """Clear the contents of the text widget."""
        self.text_widget.configure(state="normal")
        self.text_widget.delete(1.0, "end")
        self.text_widget.configure(state="disabled")

    def flush(self):
        """Flush method for compatibility with sys.stdout."""
        pass

# Example usage
if __name__ == "__main__":
    root = tk.Tk()
    root.title("Terminal Output Widget Example")

    # Create the terminal output widget
    terminal_widget = TerminalOutputWidget(root, title="Custom Terminal Output")
    terminal_widget.pack(fill="both", expand=True, padx=10, pady=10)

    # Redirect stdout to the widget
    import sys
    sys.stdout = terminal_widget

    # Example output
    print("This is an example of redirecting terminal output.")
    print("The text appears in the widget.")

    root.mainloop()
