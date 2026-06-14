import tkinter as tk
from tkinter import LabelFrame, Button, Frame
from tkinter import ttk

class AminoAcidSelectionWidget:
    def __init__(self, root, title, default_selections=None):
        """
        Initialize the Amino Acid Selection Widget.

        Parameters:
        - root: The root Tkinter widget.
        - title: Title of the widget.
        - default_selections: Optional list of amino acids to be selected by default.
        """
        self.amino_acids = [
            "Alanine (A)", "Valine (V)", "Leucine (L)", "Isoleucine (I)", "Methionine (M)", "Phenylalanine (F)", "Tryptophan (W)", "Proline (P)",
            "Serine (S)", "Threonine (T)", "Asparagine (N)", "Glutamine (Q)", "Tyrosine (Y)", "Cysteine (C)",
            "Aspartic Acid (D)", "Glutamic Acid (E)", "Lysine (K)", "Arginine (R)", "Histidine (H)"
        ]

        self.root = root
        self.title = title
        self.amino_vars = {amino: tk.IntVar(value=0) for amino in self.amino_acids}
        self.group_vars = {}
        self.visible = False

        # Apply default selections if provided
        if default_selections:
            for amino in default_selections:
                if amino in self.amino_vars:
                    self.amino_vars[amino].set(1)

        self.frame = Frame(self.root)

        self.create_amino_acid_groups()

    def create_amino_acid_groups(self):
        grouped_amino_acids = {
            "Hydrophobic": ["Alanine (A)", "Valine (V)", "Leucine (L)", "Isoleucine (I)", "Methionine (M)", "Phenylalanine (F)", "Tryptophan (W)", "Proline (P)"],
            "Polar": ["Serine (S)", "Threonine (T)", "Asparagine (N)", "Glutamine (Q)", "Tyrosine (Y)", "Cysteine (C)"],
            "Acidic": ["Aspartic Acid (D)", "Glutamic Acid (E)"],
            "Basic": ["Lysine (K)", "Arginine (R)", "Histidine (H)"]
        }

        for group, acids in grouped_amino_acids.items():
            group_frame = LabelFrame(self.frame, text=group, padx=10, pady=10)
            group_frame.pack(fill="x", padx=10, pady=5)

            # Determine if all acids in the group are selected
            group_var_value = int(all(self.amino_vars[acid].get() == 1 for acid in acids))
            group_var = tk.IntVar(value=group_var_value)
            self.group_vars[group] = group_var

            group_checkbox = ttk.Checkbutton(group_frame, text="Select All", variable=group_var, command=lambda g=group, a=acids: self.toggle_amino_group(g, a))
            group_checkbox.grid(row=0, column=0, columnspan=4, sticky="w", padx=10, pady=5)

            for i, acid in enumerate(acids):
                checkbox = ttk.Checkbutton(group_frame, text=acid, variable=self.amino_vars[acid], command=lambda g=group, a=acids: self.update_group_state(g, a))
                checkbox.grid(row=(i // 4) + 1, column=i % 4, sticky="w", padx=10, pady=5)

    def toggle_amino_group(self, group, acids):
        state = self.group_vars[group].get()
        for acid in acids:
            self.amino_vars[acid].set(state)

    def update_group_state(self, group, acids):
        # Update the group checkbox state based on individual selections
        all_selected = all(self.amino_vars[acid].get() == 1 for acid in acids)
        self.group_vars[group].set(int(all_selected))

    def get_values(self):
        amino_code_map = {
            "Alanine (A)": "Ala", "Valine (V)": "Val", "Leucine (L)": "Leu", "Isoleucine (I)": "Ile", "Methionine (M)": "Met",
            "Phenylalanine (F)": "Phe", "Tryptophan (W)": "Trp", "Proline (P)": "Pro", "Serine (S)": "Ser", "Threonine (T)": "Thr",
            "Asparagine (N)": "Asn", "Glutamine (Q)": "Gln", "Tyrosine (Y)": "Tyr", "Cysteine (C)": "Cys", "Aspartic Acid (D)": "Asp",
            "Glutamic Acid (E)": "Glu", "Lysine (K)": "Lys", "Arginine (R)": "Arg", "Histidine (H)": "His"
        }
        return {amino_code_map[amino]: var.get() for amino, var in self.amino_vars.items()}

if __name__ == "__main__":
    def toggle_visibility():
        if app.visible:
            app.frame.pack_forget()
        else:
            app.frame.pack(pady=10, padx=10)
        app.visible = not app.visible

    root = tk.Tk()
    root.title("Amino Acid Widget")

    default_selections = ["Alanine (A)", "Leucine (L)", "Glutamic Acid (E)"]

    app = AminoAcidSelectionWidget(root, "Amino Acids", default_selections=default_selections)

    toggle_button = ttk.Button(root, text="Select Amino Acids", command=toggle_visibility)
    toggle_button.pack(pady=5)

    # Initially collapse the widget
    app.frame.pack_forget()

    root.mainloop()
