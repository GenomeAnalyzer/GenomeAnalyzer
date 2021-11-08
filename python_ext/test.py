import FASTAOP
import os
import json

b1 = int("0100101010010010", 2)
b2 = int("1101101000010110", 2)

print(f"Hamming difference between {b1} and {b2}")

print(FASTAOP.hamming(b1, b2))