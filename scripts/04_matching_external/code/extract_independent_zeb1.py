#!/usr/bin/env python3
"""Stream one ZEB1 column from a wide dependency matrix."""
import csv
import sys
import urllib.request

url, dest = sys.argv[1], sys.argv[2]
req = urllib.request.Request(url, headers={"User-Agent": "research/1.0"})
with urllib.request.urlopen(req, timeout=300) as resp, open(dest, "w", newline="", encoding="utf-8") as out:
    text = (line.decode("utf-8", "replace") for line in resp)
    reader = csv.reader(text)
    header = next(reader)
    # first column is the row id
    idx = None
    for i, name in enumerate(header):
        gene = name.split(" (")[0].strip().strip('"')
        if gene == "ZEB1":
            idx = i
            break
    if idx is None:
        print("ZEB1_NOT_IN_HEADER", len(header), file=sys.stderr)
        found = None
        for row in reader:
            if row and row[0].split(" (")[0].strip() == "ZEB1":
                found = row
                break
        if found is None:
            sys.exit(2)
        writer = csv.writer(out)
        writer.writerow(["id"] + header[1:])
        writer.writerow(found)
    else:
        writer = csv.writer(out)
        writer.writerow(["id", "ZEB1"])
        for row in reader:
            if len(row) > idx:
                writer.writerow([row[0], row[idx]])
print("wrote", dest, file=sys.stderr)
