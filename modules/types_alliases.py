"""
Common type aliases for types in  bioseqkit.
"""

SeqTuple = tuple[str, str]
SeqDict = dict[str, SeqTuple]
Num = int | float
Bounds = int | float | tuple[int | float, int | float]
