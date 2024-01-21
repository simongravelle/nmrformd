#!/usr/bin/env python
# coding: utf-8

zotero_file = open("from-zotero.bib", "r")
cleaned_file = open("../source/references.bib", "w")

for line in zotero_file:

    if "journaltitle =" in line:
        line = line.replace("journaltitle", "journal")

    if ("date =" in line) & (("urldate =" not in line)):
        if "-" in line:
            year = line.split("-")[0].split("{")[1]
        else:
            year = line.split("}")[0].split("{")[1]
        assert len(year) == 4
        line = "  year = {"+year+"},\n"
    
    if ("shortjournal =" not in line) & \
        ("keywords =" not in line) & \
        ("abstract =" not in line) & \
        ("urldate =" not in line) & \
        ("langid =" not in line) & \
        ("file =" not in line):

        cleaned_file.write(line)

zotero_file.close()
cleaned_file.close()