"""
Common type aliases for types in  bioseqkit.
"""

seq_tuple = tuple[str, str]
seq_dict = dict[str, seq_tuple]
bounds = int | float | tuple[int | float, int | float]