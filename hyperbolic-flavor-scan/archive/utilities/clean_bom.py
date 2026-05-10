with open("CORE_MASTER_v5.tex", "rb") as f:
    data = f.read()

bad = [
    b'\xef\xbb\xbf',        # BOM
    b'\xe2\x94\x82',        # box-drawing |
    b'\xe2\x94\x80',        # box-drawing -
    b'\xe2\x95\x9e',        # box-drawing misc
    b'\xc3\x83\xc2\xa2',   # double-encoded
    b'\xc2\x9d',            # U+009D
    b'\x9d',
]

original = data
for b in bad:
    data = data.replace(b, b'')

# Fix em-dashes
data = data.replace(b'\xe2\x80\x94', b'---')
data = data.replace(b'\xe2\x80\x93', b'--')

print(f"Removed {len(original)-len(data)} bad bytes")

with open("CORE_MASTER_v5.tex", "wb") as f:
    f.write(data)
print("Saved.")
